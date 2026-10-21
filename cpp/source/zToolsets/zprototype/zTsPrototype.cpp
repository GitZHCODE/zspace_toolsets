#include <headers/zToolsets/zprototype/zTsPrototype.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace zSpace {
// Everything in this anonymous namespace is a plain helper function, not
// a class member - internal linkage only, invisible outside this file,
// not declared anywhere in zTsPrototype.h. They take the zTsPrototype
// instance (whose data members are public) as an explicit parameter
// instead of relying on implicit member access.
namespace {
// Plain (positions, polyCounts, polyConnects) buffer used by the
// crease-subdivision / flood-fill pipeline before a real zObjMesh is built.
struct PrototypeRawMesh {
  zPointArray positions;
  zIntArray pCounts;
  zIntArray pConnects;
};

int64_t edgeKey(int a, int b) {
  int lo = std::min(a, b), hi = std::max(a, b);
  return (int64_t(lo) << 32) | uint32_t(hi);
}

// One distinct colour per id via the golden-ratio conjugate hue
// step, so any number of ids stays spread out and never repeats.
// Shared by computeFaceBlocks and computeBlock so both colour the
// same id identically.
zColor colorForId(int id) {
  const float goldenRatioConjugate = 0.6180339887f;
  float hue = std::fmod(id * goldenRatioConjugate, 1.0f);
  float sat = 0.65f;
  float val = 0.90f;

  float hh = hue * 6.0f;
  int sector = int(hh) % 6;
  float f = hh - std::floor(hh);
  float p = val * (1.0f - sat);
  float q = val * (1.0f - sat * f);
  float t = val * (1.0f - sat * (1.0f - f));

  float r, g, b;
  switch (sector) {
  case 0:
    r = val;
    g = t;
    b = p;
    break;
  case 1:
    r = q;
    g = val;
    b = p;
    break;
  case 2:
    r = p;
    g = val;
    b = t;
    break;
  case 3:
    r = p;
    g = q;
    b = val;
    break;
  case 4:
    r = t;
    g = p;
    b = val;
    break;
  default:
    r = val;
    g = p;
    b = q;
    break;
  }
  return zColor(r, g, b, 1.0f);
}

// Checks every scalar setting that doesn't depend on the input
// graph having been parsed yet (see computePlane/computePiers) -
// pier-spacing / in-bounds checks against the DERIVED piers happen
// separately, in validatePiers, once computePiers has actually run.
bool validateSettings(zTsPrototype &self) {
  PrototypeSettings &i_Settings = self.i_Settings;

  zFnGraph fnGraph(self.i_Graph);
  if (fnGraph.numVertices() < 2) {
    self._Status = "zTsPrototype::i_Graph must contain at least the two "
                   "top-chord endpoints.";
    return false;
  }
  if (fnGraph.numEdges() < 1) {
    self._Status = "zTsPrototype::i_Graph must contain at least one edge.";
    return false;
  }
  if (i_Settings.ColumnHeight <= 0.0f) {
    self._Status = "PrototypeSettings.ColumnHeight must be positive.";
    return false;
  }
  if (i_Settings.ArchHeightCoefficient <= 0.0f) {
    self._Status = "PrototypeSettings.ArchHeightCoefficient must be positive.";
    return false;
  }
  if (i_Settings.ColumnWidth <= 0.0f) {
    self._Status = "PrototypeSettings.ColumnWidth must be positive.";
    return false;
  }
  if (i_Settings.BlockSize <= 0.0f) {
    self._Status = "PrototypeSettings.BlockSize must be positive.";
    return false;
  }
  if (i_Settings.SubdivisionCount < 0) {
    self._Status = "PrototypeSettings.SubdivisionCount must not be negative.";
    return false;
  }

  return true;
}

// Fits the plane i_Graph's own vertex positions lie in (via Newell's
// method - robust for a possibly-non-convex, near-planar point set,
// without needing to know which 3 points are safely non-collinear),
// then places the local origin at that plane's own bottom-left
// corner - the graph's own axis-aligned (in local x/z) bounding box
// minimum - so every local coordinate this produces (and everything
// computePiers through computeFaceId then builds from them) is
// >= 0, the same way the old fixed-world-XZ-plane convention always
// had x = 0 at the left edge and z = 0 at the ground. Every graph
// vertex is then projected into these local 2D coordinates
// (self._LocalGraphPositions, index-parallel to i_Graph's own vertex
// list; x along the plane's xAxis, z along its yAxis, y always 0 here
// - the thickness axis, the plane's normal, only comes in once
// computeOffsetMeshes' offset math is re-enabled). remapToGlobal
// (called right before the output meshes are built) is the only
// place local and global ever meet again.
//
// The plane's own local "vertical" (yAxis) is world Z
// projected into the plane, so a graph drawn in (or close to) the
// familiar world XZ elevation keeps behaving exactly as before;
// this falls back to world Y if the plane itself is (near)
// horizontal, where that projection would vanish.
void computePlane(zTsPrototype &self) {
  zFnGraph fnGraph(self.i_Graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  int n = int(positions.size());

  zPoint normalSum(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < n; i++) {
    const zPoint &a = positions[i];
    const zPoint &b = positions[(i + 1) % n];
    normalSum.x += (a.y - b.y) * (a.z + b.z);
    normalSum.y += (a.z - b.z) * (a.x + b.x);
    normalSum.z += (a.x - b.x) * (a.y + b.y);
  }

  float normalLen =
      std::sqrt(normalSum.x * normalSum.x + normalSum.y * normalSum.y +
                normalSum.z * normalSum.z);
  zPoint normal = (normalLen > 1e-8f)
                      ? zPoint(normalSum.x / normalLen, normalSum.y / normalLen,
                               normalSum.z / normalLen)
                      : zPoint(0.0f, 1.0f, 0.0f);

  zPoint worldUp(0.0f, 0.0f, 1.0f);
  float upDotN =
      worldUp.x * normal.x + worldUp.y * normal.y + worldUp.z * normal.z;
  zPoint yAxisRaw(worldUp.x - upDotN * normal.x, worldUp.y - upDotN * normal.y,
                  worldUp.z - upDotN * normal.z);
  float yLen = std::sqrt(yAxisRaw.x * yAxisRaw.x + yAxisRaw.y * yAxisRaw.y +
                         yAxisRaw.z * yAxisRaw.z);

  if (yLen < 1e-4f) {
    zPoint worldY(0.0f, 1.0f, 0.0f);
    float yDotN =
        worldY.x * normal.x + worldY.y * normal.y + worldY.z * normal.z;
    yAxisRaw = zPoint(worldY.x - yDotN * normal.x, worldY.y - yDotN * normal.y,
                      worldY.z - yDotN * normal.z);
    yLen = std::sqrt(yAxisRaw.x * yAxisRaw.x + yAxisRaw.y * yAxisRaw.y +
                     yAxisRaw.z * yAxisRaw.z);
  }

  zPoint yAxis(yAxisRaw.x / yLen, yAxisRaw.y / yLen, yAxisRaw.z / yLen);

  // xAxis completes a right-handed (xAxis, yAxis, normal) frame
  zPoint xAxis(yAxis.y * normal.z - yAxis.z * normal.y,
               yAxis.z * normal.x - yAxis.x * normal.z,
               yAxis.x * normal.y - yAxis.y * normal.x);

  // project every point against an arbitrary reference - any
  // actual graph point works, since it's guaranteed (assuming
  // the graph really is planar) to already lie ON the fitted
  // plane - purely to find the graph's own bounding box within
  // the plane; the reference itself drops back out once minU/
  // minV are subtracted below, so its choice never shows up in
  // the final result.
  const zPoint &reference = positions[0];
  zFloatArray rawU(n), rawV(n);
  float minU = 0.0f, minV = 0.0f;

  for (int i = 0; i < n; i++) {
    zPoint d(positions[i].x - reference.x, positions[i].y - reference.y,
             positions[i].z - reference.z);
    rawU[i] = d.x * xAxis.x + d.y * xAxis.y + d.z * xAxis.z;
    rawV[i] = d.x * yAxis.x + d.y * yAxis.y + d.z * yAxis.z;
    if (i == 0 || rawU[i] < minU)
      minU = rawU[i];
    if (i == 0 || rawV[i] < minV)
      minV = rawV[i];
  }

  // the local origin sits at the bounding box's own bottom-left
  // corner, offset from `reference` (an actual point on the
  // plane) by purely in-plane vectors - so the resulting origin
  // is guaranteed to sit exactly on the plane too, with no drift
  // along the normal.
  self._Plane.origin = zPoint(reference.x + minU * xAxis.x + minV * yAxis.x,
                              reference.y + minU * xAxis.y + minV * yAxis.y,
                              reference.z + minU * xAxis.z + minV * yAxis.z);
  self._Plane.xAxis = xAxis;
  self._Plane.yAxis = yAxis;
  self._Plane.normal = normal;

  self._LocalGraphPositions.resize(n);
  for (int i = 0; i < n; i++)
    self._LocalGraphPositions[i] = zPoint(rawU[i] - minU, 0.0f, rawV[i] - minV);
}

// Checks the piers computePiers actually derived from the input
// graph - Height headroom (needs self._Height, only known after
// computePiers), pier spacing, ascending order, and that every real
// pier sits strictly inside (0, Width).
bool validatePiers(zTsPrototype &self) {
  PrototypeSettings &i_Settings = self.i_Settings;

  if (self._Width <= 0.0f) {
    self._Status = "Input graph's top chord has zero (or negative) width - "
                   "check i_Graph.";
    return false;
  }
  if (self._Height <= 0.0f) {
    self._Status = "Input graph's top chord has zero (or negative) height - "
                   "check i_Graph.";
    return false;
  }
  if (i_Settings.ColumnHeight >= self._Height) {
    self._Status = "PrototypeSettings.ColumnHeight must be smaller than the "
                   "input graph's Height.";
    return false;
  }

  int realPierCount = int(self._Piers.size()) - 2; // minus the two phantom ends
  if (realPierCount < 1) {
    self._Status =
        "Input graph must have at least one pier branch off its top chord.";
    return false;
  }

  float halfW = i_Settings.ColumnWidth * 0.5f;
  float leftEdge = self._Piers.front().center;
  float rightEdge = self._Piers.back().center;
  float prevX = leftEdge;

  for (int i = 1; i + 1 < int(self._Piers.size()); i++) {
    float x = self._Piers[i].center;

    if (x - halfW <= leftEdge || x + halfW >= rightEdge) {
      self._Status = "Pier " + std::to_string(i - 1) +
                     " (from the input graph) sits outside (0, Width).";
      return false;
    }
    if (x <= prevX) {
      self._Status =
          "Pier positions (from the input graph) must be strictly ascending.";
      return false;
    }
    if ((x - halfW) <=
        (self._Piers[i - 1].isPhantom ? leftEdge : (prevX + halfW))) {
      self._Status = "Pier " + std::to_string(i - 1) +
                     " overlaps the previous pier - increase spacing or reduce "
                     "ColumnWidth.";
      return false;
    }

    prevX = x;
  }

  return true;
}

// Classifies every i_Graph edge by direction in the plane's own local
// coordinates (self._LocalGraphPositions, populated by computePlane) -
// BeamEdge (near-constant local z: part of the top chord),
// ColumnEdge (near-constant local x, different z: a pier's
// own drop), or Other (neither - ignored). Written to
// self._GraphEdgeType, index-parallel to zFnGraph::getEdgeData's own
// flat vertex-id-pair list - computePiers uses this to find each real
// pier directly from its own Column edge.
void classifyGraphEdges(zTsPrototype &self) {
  const zPointArray &positions = self._LocalGraphPositions;

  zFnGraph fnGraph(self.i_Graph);
  zIntArray edgeConnects;
  fnGraph.getEdgeData(edgeConnects);
  int numEdges = int(edgeConnects.size() / 2);

  const float tol = 1e-4f;

  self._GraphEdgeType.resize(numEdges);
  for (int e = 0; e < numEdges; e++) {
    int a = edgeConnects[2 * e];
    int b = edgeConnects[2 * e + 1];

    float dz = fabs(positions[a].z - positions[b].z);
    float dx = fabs(positions[a].x - positions[b].x);

    if (dz < tol)
      self._GraphEdgeType[e] = BeamEdge;
    else if (dx < tol)
      self._GraphEdgeType[e] = ColumnEdge;
    else
      self._GraphEdgeType[e] = Other;
  }
}

// Parses self._LocalGraphPositions/self._GraphEdgeType
// (computePlane's local projection of i_Graph, and classifyGraphEdges'
// per-edge classification) into self._Width, self._Height and
// self._Piers - the top chord is the set of local vertices sitting
// at the local graph's own highest z (self._Height); its x extent
// is self._Width; its two ends become the phantom piers. Any OTHER
// top-chord vertex that is an endpoint of a ColumnEdge
// edge is a real pier, at that vertex's own x - the column's own
// shape past that first edge is never looked at.
void computePiers(zTsPrototype &self) {
  const zPointArray &positions = self._LocalGraphPositions;

  zFnGraph fnGraph(self.i_Graph);
  zIntArray edgeConnects;
  fnGraph.getEdgeData(edgeConnects);

  const zIntArray &edgeType = self._GraphEdgeType;

  int numVerts = int(positions.size());
  int numEdges = int(edgeType.size());

  zIntArray hasColumnEdge(numVerts, 0);
  for (int e = 0; e < numEdges; e++) {
    if (edgeType[e] != ColumnEdge)
      continue;
    hasColumnEdge[edgeConnects[2 * e]] = 1;
    hasColumnEdge[edgeConnects[2 * e + 1]] = 1;
  }

  float topZ = positions[0].z;
  for (const zPoint &p : positions)
    topZ = std::max(topZ, p.z);

  const float tol = 1e-4f;

  zIntArray chordVerts;
  for (int v = 0; v < numVerts; v++)
    if (fabs(positions[v].z - topZ) < tol)
      chordVerts.push_back(v);
  std::sort(chordVerts.begin(), chordVerts.end(),
            [&](int a, int b) { return positions[a].x < positions[b].x; });

  if (chordVerts.size() < 2)
    return; // validatePiers below will reject this (_Width/_Height stay 0)

  float minX = positions[chordVerts.front()].x;
  float maxX = positions[chordVerts.back()].x;

  self._Width = maxX - minX;
  self._Height = topZ;

  ColumnPier left;
  left.center = minX;
  left.halfWidth = 0.0f;
  left.isPhantom = true;
  self._Piers.push_back(left);

  float halfW = self.i_Settings.ColumnWidth * 0.5f;

  for (size_t i = 1; i + 1 < chordVerts.size(); i++) {
    int v = chordVerts[i];
    if (!hasColumnEdge[v])
      continue;

    ColumnPier pier;
    pier.center = positions[v].x;
    pier.halfWidth = halfW;
    pier.isPhantom = false;
    self._Piers.push_back(pier);
  }

  ColumnPier right;
  right.center = maxX;
  right.halfWidth = 0.0f;
  right.isPhantom = true;
  self._Piers.push_back(right);
}

zPoint pointOnArc(const zTsPrototype &self, const Span &span, float angleDeg) {
  float angleRad = angleDeg * float(DEG_TO_RAD);
  float x = span.center.x + span.radius * cos(angleRad);
  float height = span.center.z + span.archHeight * sin(angleRad);
  return zPoint(x, 0.0f, height);
}

// The local polyline's miter (angle-bisector) normal direction at
// `cur` - built from its two neighbouring segments' own (unit)
// TANGENTS, not their normals: sum the two unit tangents to get the
// angle bisector direction, then rotate it 90 degrees (rotating a
// bisector is the same as bisecting the already-rotated vectors, so
// this is exactly the miter join a polyline offset/buffer uses,
// just without a separate per-edge normal helper). prev/next are
// nullptr at an actual end of the polyline (only one neighbouring
// segment exists there) - a mitred normal isn't well-defined with
// just one edge, so that boundary point simply projects straight up
// instead. This only ever looks at already-sampled points, never at
// the curve's own parametric formula, so it gives the arch's normal
// identically whether the points were sampled off an ellipse, a
// circular arc, or any other (even freeform) curve.
zPoint miterNormalDir(const zPoint *prev, const zPoint &cur,
                      const zPoint *next) {
  if (!prev || !next)
    return zPoint(0.0f, 0.0f, 1.0f);

  zPoint tPrev(cur.x - prev->x, 0.0f, cur.z - prev->z);
  zPoint tNext(next->x - cur.x, 0.0f, next->z - cur.z);

  float lenPrev = std::sqrt(tPrev.x * tPrev.x + tPrev.z * tPrev.z);
  float lenNext = std::sqrt(tNext.x * tNext.x + tNext.z * tNext.z);

  zPoint uPrev = (lenPrev > 1e-8f)
                     ? zPoint(tPrev.x / lenPrev, 0.0f, tPrev.z / lenPrev)
                     : zPoint(0.0f, 0.0f, 0.0f);
  zPoint uNext = (lenNext > 1e-8f)
                     ? zPoint(tNext.x / lenNext, 0.0f, tNext.z / lenNext)
                     : zPoint(0.0f, 0.0f, 0.0f);

  // tangent angle bisector, then rotated 90 degrees to a normal
  float bx = uPrev.x + uNext.x;
  float bz = uPrev.z + uNext.z;
  return zPoint(-bz, 0.0f, bx);
}

// miterNormalDir at span.curveVertexIds[k], reading positions
// straight out of self._Positions.
zPoint miterNormalDirAt(const zTsPrototype &self, const Span &span, int k) {
  int n = span.panelCount;
  const zPoint *prev =
      (k > 0) ? &self._Positions[span.curveVertexIds[k - 1]] : nullptr;
  const zPoint *next =
      (k < n) ? &self._Positions[span.curveVertexIds[k + 1]] : nullptr;
  return miterNormalDir(prev, self._Positions[span.curveVertexIds[k]], next);
}

// Intersects the line through P in direction dirP with the line
// through Q in direction dirQ, in the XZ plane. Returns false (P/Q
// left untouched) if the two directions are parallel.
bool intersectRays(const zPoint &P, const zPoint &dirP, const zPoint &Q,
                   const zPoint &dirQ, zPoint &outPoint) {
  float denom = dirP.x * dirQ.z - dirP.z * dirQ.x;
  if (fabs(denom) < 1e-8f)
    return false;

  float t = ((Q.x - P.x) * dirQ.z - (Q.z - P.z) * dirQ.x) / denom;
  outPoint = zPoint(P.x + t * dirP.x, 0.0f, P.z + t * dirP.z);
  return true;
}

// Casts a ray from `from` in direction `dir` up to z = targetZ.
// Falls back to a purely vertical projection if dir is (near)
// horizontal (its own x wouldn't reach targetZ at all).
zPoint projectAlongRay(const zPoint &from, const zPoint &dir, float targetZ) {
  if (fabs(dir.z) < 1e-6f)
    return zPoint(from.x, 0.0f, targetZ);

  float t = (targetZ - from.z) / dir.z;
  return zPoint(from.x + t * dir.x, 0.0f, targetZ);
}

float thicknessForArchPoint(const zTsPrototype &self, const Span &span,
                            float z) {
  float springZ = self.i_Settings.ColumnHeight;
  float frac =
      (span.archHeight > 1e-6f) ? (z - springZ) / span.archHeight : 0.0f;
  frac = std::min(std::max(frac, 0.0f), 1.0f);
  return self.i_Settings.ArchBottomOffset +
         (self.i_Settings.ArchTopOffset - self.i_Settings.ArchBottomOffset) *
             frac;
}

float thicknessForColumnPoint(const zTsPrototype &self, float z) {
  float springZ = self.i_Settings.ColumnHeight;
  float frac = (springZ > 1e-6f) ? z / springZ : 0.0f;
  frac = std::min(std::max(frac, 0.0f), 1.0f);
  return self.i_Settings.ColumnBottomOffset +
         (self.i_Settings.ArchBottomOffset -
          self.i_Settings.ColumnBottomOffset) *
             frac;
}

void computeSpans(zTsPrototype &self) {
  // ColumnHeight is a single global input, so every span's ellipse
  // centre sits at this same height, regardless of its own rise.
  float springZ = self.i_Settings.ColumnHeight;

  self._Spans.resize(self._Piers.size() - 1);

  for (int i = 0; i < self._Spans.size(); i++) {
    ColumnPier &pierL = self._Piers[i];
    ColumnPier &pierR = self._Piers[i + 1];

    Span span;
    span.xLeft = pierL.center + pierL.halfWidth;
    span.xRight = pierR.center - pierR.halfWidth;

    float clearWidth = span.xRight - span.xLeft;

    if (!pierL.isPhantom && !pierR.isPhantom) {
      span.radius = clearWidth * 0.5f;
      float cx = (span.xLeft + span.xRight) * 0.5f;
      span.center = zPoint(cx, 0.0f, springZ);
      span.angleStart = 180.0f;
      span.angleEnd = 0.0f;
    } else if (pierL.isPhantom) {
      span.radius = clearWidth;
      span.center = zPoint(span.xLeft, 0.0f, springZ);
      span.angleStart = 90.0f;
      span.angleEnd = 0.0f;
    } else // pierR.isPhantom
    {
      span.radius = clearWidth;
      span.center = zPoint(span.xRight, 0.0f, springZ);
      span.angleStart = 180.0f;
      span.angleEnd = 90.0f;
    }

    // each arch's own rise, independent of the others
    span.archHeight = span.radius * self.i_Settings.ArchHeightCoefficient;

    // panel count = clear width (the gap between this span's two
    // piers, one of which may be a zero-width phantom end pier)
    // / BlockSize, rounded UP (never coarser than BlockSize
    // asks for). A FULL arch (both piers real) is additionally
    // rounded up to the next ODD count, so it always has a
    // single centred crown panel; a HALF arch (against an outer
    // end) is left odd or even, but always at least 2.
    bool isFullArch = !pierL.isPhantom && !pierR.isPhantom;
    int rawPanelCount =
        std::max(1, int(std::ceil(clearWidth / self.i_Settings.BlockSize)));
    if (isFullArch && rawPanelCount % 2 == 0)
      rawPanelCount += 1;
    if (!isFullArch)
      rawPanelCount = std::max(2, rawPanelCount);
    span.panelCount = rawPanelCount;

    self._Spans[i] = span;
  }
}

int addVertex(zTsPrototype &self, const zPoint &pos, float thickness) {
  int vId = int(self._Positions.size());
  self._Positions.push_back(pos);
  self._Thickness.push_back(thickness);
  return vId;
}

// Creates every span's panelCount + 1 boundary vertices (curve and
// beam-top at each voussoir seam, including both ends) exactly once
// - there is no beam-bottom vertex row at all. beamTopVertexIds[k]
// is curveVertexIds[k]'s own local polyline miter normal
// (miterNormalDirAt - the angle bisector of its two neighbouring
// segments) cast up to the beam top, for every k. The k = 1 /
// panelCount - 1 points next to a real pier don't actually use
// their own beamTopVertexIds entry in the final mesh -
// computePierPanels instead routes them through that pier's own
// normal-ray intersection point (pMid) - but the entry is still
// created here uniformly, for the rare fallback case where that
// pier-side treatment isn't possible (see computePierPanels).
// computeArchPanels and computePierPanels both read curveVertexIds
// from this array afterwards instead of independently recomputing +
// re-adding a boundary point, which is what actually keeps the mesh
// connected now that addVertex doesn't merge by position: two
// adjacent voussoirs, or a pier and the span next to it, share an
// edge only because both sides reference the same array entry.
void computeSpanBoundaryVertices(zTsPrototype &self) {
  for (int s = 0; s < self._Spans.size(); s++) {
    Span &span = self._Spans[s];
    int n = span.panelCount;
    float angleStep = (span.angleEnd - span.angleStart) / float(n);

    span.curveVertexIds.resize(n + 1);
    span.beamTopVertexIds.resize(n + 1);

    // pass 1: every curve point first, so pass 2's miter normal
    // (which looks at self._Positions via curveVertexIds) can
    // see any point's actual neighbours regardless of which k
    // they're at.
    for (int k = 0; k <= n; k++) {
      float angle = span.angleStart + k * angleStep;
      zPoint curve = pointOnArc(self, span, angle);
      span.curveVertexIds[k] =
          addVertex(self, curve, thicknessForArchPoint(self, span, curve.z));
    }

    // pass 2: beam-top points
    for (int k = 0; k <= n; k++) {
      zPoint dir = miterNormalDirAt(self, span, k);
      zPoint beamTop = projectAlongRay(self._Positions[span.curveVertexIds[k]],
                                       dir, self._Height);

      span.beamTopVertexIds[k] =
          addVertex(self, beamTop, self.i_Settings.BeamTopOffset);
    }
  }
}

// Collapses zero-length edges (two adjacent vertices landing on the
// same position despite being different ids - addVertex no longer
// merges by position, so this checks actual positions, not id
// equality) and drops degenerate faces - zFnMesh::create's
// half-edge build can hang on those otherwise.
// groupId < 0 starts a new panel group (its own face index); an id
// from an earlier call folds this face into that same group.
// blockType tags this face's structural family (0 = Beam, 1 =
// Column - see zTsPrototype::_BlockType); it is NOT part of the
// grouping/merging logic, purely a per-face classification.
// Returns the assigned group id, or -1 if the polygon was dropped.
int addPolygon(zTsPrototype &self, const zIntArray &vertexIds, int groupId = -1,
               int blockType = 0) {
  // A small, explicit epsilon - not PRECISION, whose exact meaning
  // (and units) this codebase's other call sites treat as an
  // opaque library-defined tolerance passed straight into
  // zUtilsCore functions, not necessarily a tiny pre-squared
  // distance safe to reuse in a raw squared-distance test here.
  const float degenerateEdgeEpsilon = 1e-5f;

  auto sameSpot = [&self](int a, int b) -> bool {
    const zPoint &pa = self._Positions[a];
    const zPoint &pb = self._Positions[b];
    float dx = pa.x - pb.x, dy = pa.y - pb.y, dz = pa.z - pb.z;
    return (dx * dx + dy * dy + dz * dz) <
           (degenerateEdgeEpsilon * degenerateEdgeEpsilon);
  };

  zIntArray cleaned;
  for (int id : vertexIds) {
    if (!cleaned.empty() &&
        (cleaned.back() == id || sameSpot(cleaned.back(), id)))
      continue;
    cleaned.push_back(id);
  }
  if (cleaned.size() >= 2 && (cleaned.front() == cleaned.back() ||
                              sameSpot(cleaned.front(), cleaned.back())))
    cleaned.pop_back();

  if (cleaned.size() < 3)
    return -1;

  int faceIndex = int(self._PCounts.size());
  for (int id : cleaned)
    self._PConnects.push_back(id);
  self._PCounts.push_back(int(cleaned.size()));

  int assignedGroupId = (groupId >= 0) ? groupId : faceIndex;
  self._PanelGroupId.push_back(assignedGroupId);
  self._BlockType.push_back(blockType);

  return assignedGroupId;
}

// Which of a span's panelCount + 1 boundary indices are claimed by
// an adjacent REAL pier's own construction (computePierPanels),
// leaving [startFree, endFree] for computeArchPanels' plain
// normal-ray-projected quads. A real pier with at least one interior
// point on this side (panelCount >= 2) claims k = 0 (or n) as its
// own spring corner AND k = 1 (or n - 1) as the point whose normal
// ray it intersects with its OTHER neighbouring span's own matching
// point (see computePierPanels) - "leftExtra"/"rightExtra" is a
// further claim of ONE more index (k = 2 or n - 2) for the pentagon
// that bridges that intersection point up to the beam top, only
// when there's room (panelCount >= 3) and it doesn't collide with
// whatever the far end of this same span claims. If a collision
// would occur (both ends want the very same middle index, on a very
// coarsely-panelled span), the left side wins arbitrarily and the
// right side's extra is dropped, since only one of them can own
// that point's route to the beam top.
struct SpanClaim {
  bool leftHasInterior = false;
  bool rightHasInterior = false;
  bool leftExtra = false;
  bool rightExtra = false;
  int startFree = 0;
  int endFree = 0;
};

SpanClaim computeSpanClaim(zTsPrototype &self, int s) {
  const Span &span = self._Spans[s];
  int n = span.panelCount;
  const ColumnPier &pierL = self._Piers[s];
  const ColumnPier &pierR = self._Piers[s + 1];

  SpanClaim claim;
  claim.leftHasInterior = !pierL.isPhantom && n >= 2;
  claim.rightHasInterior = !pierR.isPhantom && n >= 2;

  int leftBaseEnd = claim.leftHasInterior ? 1 : -1;
  int rightBaseStart = claim.rightHasInterior ? (n - 1) : (n + 1);

  int leftExtraIdx = leftBaseEnd + 1;
  int rightExtraIdx = rightBaseStart - 1;

  claim.leftExtra =
      claim.leftHasInterior && n >= 3 && leftExtraIdx < rightBaseStart;
  claim.rightExtra =
      claim.rightHasInterior && n >= 3 && rightExtraIdx > leftBaseEnd;
  if (claim.leftExtra && claim.rightExtra && leftExtraIdx == rightExtraIdx)
    claim.rightExtra = false;

  // Pentagon B/C don't RESERVE index 2 / n-2 and exclude it - they
  // terminate there, sharing that point's own beam-top edge
  // (curve2-beamTop2) with the very next ordinary quad, the same
  // way curve1 is shared between pentagon A and B. So the free
  // range for ordinary quads STARTS at (not after) that index.
  claim.startFree = claim.leftExtra
                        ? leftExtraIdx
                        : (claim.leftHasInterior ? (leftBaseEnd + 1) : 0);
  claim.endFree = claim.rightExtra
                      ? rightExtraIdx
                      : (claim.rightHasInterior ? (rightBaseStart - 1) : n);

  return claim;
}

// Plain normal-ray-projected quads {curve_j, curve_j+1, beamTop_j+1,
// beamTop_j} over whichever boundary range each span's neighbouring
// piers didn't claim for themselves (see computeSpanClaim) - every
// point adjacent to a real pier is instead handled as part of that
// pier's own pentagon construction in computePierPanels.
void computeArchPanels(zTsPrototype &self) {
  for (int s = 0; s < self._Spans.size(); s++) {
    Span &span = self._Spans[s];
    SpanClaim claim = computeSpanClaim(self, s);

    for (int j = claim.startFree; j < claim.endFree; j++) {
      int idCurveA = span.curveVertexIds[j];
      int idCurveB = span.curveVertexIds[j + 1];
      int idBeamTopA = span.beamTopVertexIds[j];
      int idBeamTopB = span.beamTopVertexIds[j + 1];

      addPolygon(self, {idCurveA, idCurveB, idBeamTopB, idBeamTopA});
    }
  }
}

void computePierPanels(zTsPrototype &self) {
  float springZ = self.i_Settings.ColumnHeight;

  self._PierGroupId.assign(self._Piers.size(), -1);

  for (int i = 0; i < self._Piers.size(); i++) {
    const ColumnPier &pier = self._Piers[i];
    if (pier.isPhantom)
      continue;

    float xL = pier.center - pier.halfWidth;
    float xR = pier.center + pier.halfWidth;

    Span &leftSpan = self._Spans[i - 1];
    Span &rightSpan = self._Spans[i];

    // pier i is rightSpan's OWN left pier, and leftSpan's OWN
    // right pier - computeSpanClaim already knows that (it looks
    // the piers up itself), so these two calls describe exactly
    // what this pier is entitled to claim on each side.
    SpanClaim claimRight = computeSpanClaim(self, i);
    SpanClaim claimLeft = computeSpanClaim(self, i - 1);
    bool hasRightInterior = claimRight.leftHasInterior;
    bool hasLeftInterior = claimLeft.rightHasInterior;

    // the spring-height corners are, by construction, the exact
    // same position AND thickness as each neighbouring span's
    // own k = 0 / k = panelCount boundary vertex. Reuse those
    // ids directly instead of adding new coincident vertices,
    // since addVertex no longer merges by position.
    int v2 = rightSpan.curveVertexIds[0];
    int v5 = leftSpan.curveVertexIds[leftSpan.panelCount];

    // pier column: base to spring, stacked panels - count is this
    // pier's own height (springZ) / BlockSize, rounded to the
    // nearest count, odd or even (no parity constraint, unlike
    // a full arch's panel count). The top of the last segment
    // is v2/v5 above, not a fresh vertex, and each segment's
    // own top corners are reused as the next segment's bottom
    // corners (prevCv2/prevCv5), same reasoning as everywhere
    // else in this function.
    int columnCount =
        std::max(1, int(std::ceil(springZ / self.i_Settings.BlockSize)));
    float columnStep = springZ / float(columnCount);
    int prevCv2 = -1, prevCv5 = -1;

    for (int k = 0; k < columnCount; k++) {
      float zA = k * columnStep;
      float zB = (k + 1) * columnStep;
      float tA = thicknessForColumnPoint(self, zA);
      float tB = thicknessForColumnPoint(self, zB);
      bool isTop = (k == columnCount - 1);

      int cv0 = (k == 0) ? addVertex(self, zPoint(xL, 0.0f, zA), tA) : prevCv5;
      int cv1 = (k == 0) ? addVertex(self, zPoint(xR, 0.0f, zA), tA) : prevCv2;
      int cv2 = isTop ? v2 : addVertex(self, zPoint(xR, 0.0f, zB), tB);
      int cv5 = isTop ? v5 : addVertex(self, zPoint(xL, 0.0f, zB), tB);

      addPolygon(self, {cv0, cv1, cv2, cv5}, -1, 1);

      prevCv2 = cv2;
      prevCv5 = cv5;
    }

    // the pier's own cap: v5 -> v2 along the spring, then up
    // each side. Where an interior point exists on both sides,
    // that "first" point (curve1R/curve1L) casts its own local
    // polyline miter normal (miterNormalDirAt), and the two
    // rays are intersected to get pMid - a single new point
    // straight above which the beam top is reached (apex),
    // shared with (but not merged into the same panel group as)
    // the two bridging pentagons below - this is what connects
    // the pier smoothly into its two innermost voussoirs
    // without ever letting the two sides' beam-top points cross
    // each other or the pier's own corners, the way the old
    // span-wide clamp could.
    zIntArray cap = {v5, v2};
    int pMidId = -1, apexId = -1;

    if (hasRightInterior && hasLeftInterior) {
      int curve1R = rightSpan.curveVertexIds[1];
      int curve1L = leftSpan.curveVertexIds[leftSpan.panelCount - 1];

      zPoint dirR = miterNormalDirAt(self, rightSpan, 1);
      zPoint dirL = miterNormalDirAt(self, leftSpan, leftSpan.panelCount - 1);

      zPoint pMid;
      if (!intersectRays(self._Positions[curve1R], dirR,
                         self._Positions[curve1L], dirL, pMid))
        pMid = zPoint(
            (self._Positions[curve1R].x + self._Positions[curve1L].x) * 0.5f,
            0.0f,
            std::max(self._Positions[curve1R].z, self._Positions[curve1L].z));

      // apex sits ON the beam top (z = Height), same as every ordinary
      // beamTopVertexIds entry - so it takes that same flat BeamTopOffset
      // thickness, not an interpolated one.
      zPoint apex(pMid.x, 0.0f, self._Height);
      apexId = addVertex(self, apex, self.i_Settings.BeamTopOffset);

      // pMid, unlike apex, is NOT on the beam top - it's the actual point
      // where the two arches' own miter-normal rays cross, sitting between
      // the beam above (apex) and the two arch curves it's wired to on
      // either side (curve1R/curve1L, via the cap and the two bridging
      // pentagons) - so its thickness is the average of those three real
      // neighbours instead.
      float pMidThickness =
          (self._Thickness[apexId] + self._Thickness[curve1R] +
           self._Thickness[curve1L]) /
          3.0f;

      std::cout << "[zTsPrototype] pMid thickness interpolation - apex: "
                << self._Thickness[apexId]
                << ", curve1R: " << self._Thickness[curve1R]
                << ", curve1L: " << self._Thickness[curve1L]
                << " -> pMid: " << pMidThickness << std::endl;

      pMidId = addVertex(self, pMid, pMidThickness);

      cap.push_back(curve1R);
      cap.push_back(pMidId);
      cap.push_back(curve1L);
    } else {
      // too few panels on one/both sides for an interior
      // point (no pMid/apex possible without both) - fall
      // back to that side's own plain normal-ray beam-top
      // point(s), same treatment as an ordinary
      // (non-pier-adjacent) voussoir corner. If that side's
      // own claimExtra is still true (it only depends on its
      // own span, not on the other side having an interior
      // point too), pull its k = 2 point in the same way, so
      // computeArchPanels' claimed range (which doesn't know
      // pMid failed to materialise) still matches what this
      // cap actually consumed.
      if (hasRightInterior) {
        cap.push_back(rightSpan.curveVertexIds[1]);
        cap.push_back(rightSpan.beamTopVertexIds[1]);
        if (claimRight.leftExtra) {
          cap.push_back(rightSpan.curveVertexIds[2]);
          cap.push_back(rightSpan.beamTopVertexIds[2]);
        }
      } else
        cap.push_back(rightSpan.beamTopVertexIds[0]);

      if (hasLeftInterior) {
        if (claimLeft.rightExtra) {
          cap.push_back(leftSpan.beamTopVertexIds[leftSpan.panelCount - 2]);
          cap.push_back(leftSpan.curveVertexIds[leftSpan.panelCount - 2]);
        }
        cap.push_back(leftSpan.beamTopVertexIds[leftSpan.panelCount - 1]);
        cap.push_back(leftSpan.curveVertexIds[leftSpan.panelCount - 1]);
      } else
        cap.push_back(leftSpan.beamTopVertexIds[leftSpan.panelCount]);
    }

    // the cap and its two bridging pentagons are three SEPARATE
    // blocks (each its own default -1 -> own-face-index group),
    // not one merged panel group - they only share vertices
    // (curve1R/curve1L, pMid, apex), never a groupId.
    int panelGroupId = addPolygon(self, cap);
    self._PierGroupId[i] = panelGroupId;

    // the pentagon bridging curve1R's normal-ray route (pMid,
    // apex) to the next PLAIN voussoir corner along the arch -
    // only when pMid/apex actually exist (both sides had an
    // interior point) AND computeSpanClaim granted this pier
    // the extra index; claimRight.leftExtra alone isn't enough,
    // since it only looks at rightSpan's own panelCount, not
    // whether leftSpan also had an interior point to intersect
    // with.
    if (pMidId >= 0 && claimRight.leftExtra) {
      int curve1R = rightSpan.curveVertexIds[1];
      int curve2R = rightSpan.curveVertexIds[2];
      int beamTop2R = rightSpan.beamTopVertexIds[2];
      addPolygon(self, {curve1R, curve2R, beamTop2R, apexId, pMidId});
    }

    if (pMidId >= 0 && claimLeft.rightExtra) {
      int curve1L = leftSpan.curveVertexIds[leftSpan.panelCount - 1];
      int curve2L = leftSpan.curveVertexIds[leftSpan.panelCount - 2];
      int beamTop2L = leftSpan.beamTopVertexIds[leftSpan.panelCount - 2];
      addPolygon(self, {curve1L, pMidId, apexId, beamTop2L, curve2L});
    }
  }
}

// self._Positions offset by +half of each vertex's own _Thickness,
// along local y (the thickness axis - every self._Positions y is 0,
// see computePlane) - still in the bridge's own LOCAL plane
// coordinates at this point (remapToGlobal hasn't run yet), which
// is exactly why this step, and computeBlock() right after it, both
// run before it: computeBlock() derives the -half-thickness side
// purely by mirroring this same local y, which only means anything
// before local and world axes have been mixed together.
void computeOffsetMeshes(zTsPrototype &self) {
  self._PositionsPositive = self._Positions;
  for (int i = 0; i < int(self._PositionsPositive.size()); i++)
    self._PositionsPositive[i].y = self._Thickness[i] * 0.5f;
}

// Maps a single LOCAL plane-coordinate point (see computePlane) to
// world/global 3D: local x along the plane's xAxis, local z along its
// yAxis, local y (the thickness axis) along its normal. Shared by
// remapToGlobal (self._Positions/_PositionsPositive) and computeBlock
// (its own locally-built closed-solid positions, remapped once each
// solid's local geometry - offset and mirrored - is finished).
zPoint toGlobalPoint(const zTsPrototype &self, const zPoint &local) {
  const Plane &plane = self._Plane;
  float gx = plane.origin.x + local.x * plane.xAxis.x +
             local.y * plane.normal.x + local.z * plane.yAxis.x;
  float gy = plane.origin.y + local.x * plane.xAxis.y +
             local.y * plane.normal.y + local.z * plane.yAxis.y;
  float gz = plane.origin.z + local.x * plane.xAxis.z +
             local.y * plane.normal.z + local.z * plane.yAxis.z;
  return zPoint(gx, gy, gz);
}

// Called once, right before the output meshes are built -
// everything upstream of this (computePiers through computeBlock)
// works entirely in local coordinates.
void remapToGlobal(zTsPrototype &self) {
  for (zPoint &p : self._Positions)
    p = toGlobalPoint(self, p);
  for (zPoint &p : self._PositionsPositive)
    p = toGlobalPoint(self, p);
}

void buildMesh(zTsPrototype &self) {
  zFnMesh fnMesh(self.o_Mesh);
  fnMesh.create(self._Positions, self._PCounts, self._PConnects);
}

void buildMeshPositive(zTsPrototype &self) {
  zFnMesh fnPositive(self.o_MeshPositive);
  fnPositive.create(self._PositionsPositive, self._PCounts, self._PConnects);
}

// Edges touching a non-quad face (a fused pier/arch panel). Kept
// sharp too, since the Catmull-Clark smooth edge-point rule averages
// in that face's centroid, which for a wide n-gon visibly bows an
// otherwise-straight neighbouring edge inward.
std::unordered_set<int64_t>
findIrregularFaceCreaseEdges(const zIntArray &pCounts,
                             const zIntArray &pConnects) {
  struct EdgeAdj {
    int face0 = -1, face1 = -1;
  };
  std::unordered_map<int64_t, EdgeAdj> edgeAdj;

  int offset = 0;
  for (int f = 0; f < int(pCounts.size()); f++) {
    int n = pCounts[f];
    for (int i = 0; i < n; i++) {
      int a = pConnects[offset + i];
      int b = pConnects[offset + (i + 1) % n];
      int64_t key = edgeKey(a, b);

      EdgeAdj &adj = edgeAdj[key];
      if (adj.face0 == -1)
        adj.face0 = f;
      else
        adj.face1 = f;
    }
    offset += n;
  }

  std::unordered_set<int64_t> crease;
  for (const auto &kv : edgeAdj) {
    const EdgeAdj &adj = kv.second;
    bool n0IsQuad = (adj.face0 >= 0) && (pCounts[adj.face0] == 4);
    bool n1IsQuad =
        (adj.face1 == -1) || ((adj.face1 >= 0) && (pCounts[adj.face1] == 4));

    if (!n0IsQuad || !n1IsQuad)
      crease.insert(kv.first);
  }

  return crease;
}

// One step of Catmull-Clark subdivision with creaseIn edges treated
// as infinitely-weighted sharp creases (never split/disconnected).
// faceIdIn/faceIdOut carry a per-face id through unchanged.
// groupBoundaryIn/Out separately track (bookkeeping only, no effect
// on the smoothing math) where a caller-chosen edge subset ends up.
void subdivideCatmullClarkCrease(
    const PrototypeRawMesh &in, const std::unordered_set<int64_t> &creaseIn,
    const zIntArray &faceIdIn,
    const std::unordered_set<int64_t> &groupBoundaryIn, PrototypeRawMesh &out,
    std::unordered_set<int64_t> &creaseOut, zIntArray &faceIdOut,
    std::unordered_set<int64_t> &groupBoundaryOut, bool keepCorners = true) {
  int numFaces = int(in.pCounts.size());
  int numVerts = int(in.positions.size());

  zIntArray faceStart(numFaces);
  {
    int running = 0;
    for (int f = 0; f < numFaces; f++) {
      faceStart[f] = running;
      running += in.pCounts[f];
    }
  }

  struct EdgeInfo {
    int v0 = -1, v1 = -1, faceCount = 0, face0 = -1, face1 = -1;
    bool crease = false;
  };
  std::unordered_map<int64_t, EdgeInfo> edgeMap;

  for (int f = 0; f < numFaces; f++) {
    int n = in.pCounts[f];
    int s = faceStart[f];
    for (int i = 0; i < n; i++) {
      int a = in.pConnects[s + i];
      int b = in.pConnects[s + (i + 1) % n];
      int64_t key = edgeKey(a, b);

      EdgeInfo &info = edgeMap[key];
      info.v0 = std::min(a, b);
      info.v1 = std::max(a, b);
      if (info.faceCount == 0)
        info.face0 = f;
      else
        info.face1 = f;
      info.faceCount++;
      if (creaseIn.count(key))
        info.crease = true;
    }
  }

  zPointArray facePoints(numFaces);
  for (int f = 0; f < numFaces; f++) {
    int n = in.pCounts[f];
    int s = faceStart[f];
    zPoint c(0.0f, 0.0f, 0.0f);
    for (int i = 0; i < n; i++)
      c = c + in.positions[in.pConnects[s + i]];
    facePoints[f] = c * (1.0f / float(n));
  }

  std::vector<zPoint> vFaceSum(numVerts, zPoint(0.0f, 0.0f, 0.0f));
  std::vector<int> vFaceCount(numVerts, 0);
  std::vector<zPoint> vEdgeMidSum(numVerts, zPoint(0.0f, 0.0f, 0.0f));
  std::vector<int> vEdgeCount(numVerts, 0);
  std::vector<int> vCreaseCount(numVerts, 0);
  std::vector<zPoint> vCreaseNeighborSum(numVerts, zPoint(0.0f, 0.0f, 0.0f));

  for (int f = 0; f < numFaces; f++) {
    int n = in.pCounts[f];
    int s = faceStart[f];
    for (int i = 0; i < n; i++) {
      int v = in.pConnects[s + i];
      vFaceSum[v] = vFaceSum[v] + facePoints[f];
      vFaceCount[v]++;
    }
  }

  zPointArray newPositions(numVerts);
  std::unordered_map<int64_t, int> edgePointId;

  for (auto &kv : edgeMap) {
    const EdgeInfo &info = kv.second;
    zPoint p0 = in.positions[info.v0];
    zPoint p1 = in.positions[info.v1];
    zPoint mid = (p0 + p1) * 0.5f;

    bool sharp = info.crease || (info.faceCount == 1);

    zPoint edgePoint;
    if (sharp)
      edgePoint = mid;
    else
      edgePoint =
          (p0 + p1 + facePoints[info.face0] + facePoints[info.face1]) * 0.25f;

    int newEdgePointId = numVerts + int(edgePointId.size());
    edgePointId[kv.first] = newEdgePointId;
    newPositions.push_back(edgePoint);

    vEdgeMidSum[info.v0] = vEdgeMidSum[info.v0] + mid;
    vEdgeCount[info.v0]++;
    vEdgeMidSum[info.v1] = vEdgeMidSum[info.v1] + mid;
    vEdgeCount[info.v1]++;

    if (sharp) {
      vCreaseCount[info.v0]++;
      vCreaseNeighborSum[info.v0] = vCreaseNeighborSum[info.v0] + p1;
      vCreaseCount[info.v1]++;
      vCreaseNeighborSum[info.v1] = vCreaseNeighborSum[info.v1] + p0;
    }
  }

  for (int v = 0; v < numVerts; v++) {
    zPoint P = in.positions[v];
    int k = vCreaseCount[v];

    if (k == 2) {
      // valence 2 (both edges sharp, no others) = a true
      // corner, kept sharp instead of the crease rounding rule
      if (keepCorners && vEdgeCount[v] == 2)
        newPositions[v] = P;
      else
        newPositions[v] = (vCreaseNeighborSum[v] + P * 6.0f) * (1.0f / 8.0f);
    } else if (k == 0) {
      int n = vEdgeCount[v];
      if (n == 0) {
        newPositions[v] = P;
        continue;
      }

      zPoint F = vFaceSum[v] * (1.0f / float(vFaceCount[v]));
      zPoint R = vEdgeMidSum[v] * (1.0f / float(vEdgeCount[v]));
      newPositions[v] = (F + R * 2.0f + P * float(n - 3)) * (1.0f / float(n));
    } else {
      // dart (k == 1) or corner (k >= 3): held fixed
      newPositions[v] = P;
    }
  }

  zIntArray facePointId(numFaces);
  for (int f = 0; f < numFaces; f++) {
    facePointId[f] = int(newPositions.size());
    newPositions.push_back(facePoints[f]);
  }

  zIntArray newPCounts;
  zIntArray newPConnects;
  zIntArray newFaceId;
  creaseOut.clear();
  groupBoundaryOut.clear();

  for (int f = 0; f < numFaces; f++) {
    int n = in.pCounts[f];
    int s = faceStart[f];

    for (int i = 0; i < n; i++) {
      int vPrev = in.pConnects[s + (i - 1 + n) % n];
      int vCur = in.pConnects[s + i];
      int vNext = in.pConnects[s + (i + 1) % n];

      int64_t keyNext = edgeKey(vCur, vNext);
      int64_t keyPrev = edgeKey(vPrev, vCur);

      int idVCur = vCur;
      int idENext = edgePointId.at(keyNext);
      int idFace = facePointId[f];
      int idEPrev = edgePointId.at(keyPrev);

      newPConnects.push_back(idVCur);
      newPConnects.push_back(idENext);
      newPConnects.push_back(idFace);
      newPConnects.push_back(idEPrev);
      newPCounts.push_back(4);
      newFaceId.push_back(faceIdIn[f]);

      if (edgeMap.at(keyNext).crease)
        creaseOut.insert(edgeKey(idVCur, idENext));
      if (edgeMap.at(keyPrev).crease)
        creaseOut.insert(edgeKey(idEPrev, idVCur));

      if (groupBoundaryIn.count(keyNext))
        groupBoundaryOut.insert(edgeKey(idVCur, idENext));
      if (groupBoundaryIn.count(keyPrev))
        groupBoundaryOut.insert(edgeKey(idEPrev, idVCur));
    }
  }

  out.positions = newPositions;
  out.pCounts = newPCounts;
  out.pConnects = newPConnects;
  faceIdOut = newFaceId;
}

// Edges whose two adjacent faces belong to different panel groups -
// the original (pre-subdivision) boundary between two output blocks.
std::unordered_set<int64_t>
findPanelGroupBoundaryEdges(const zIntArray &panelGroupId,
                            const zIntArray &pCounts,
                            const zIntArray &pConnects) {
  struct EdgeAdj {
    int face0 = -1, face1 = -1;
  };
  std::unordered_map<int64_t, EdgeAdj> edgeAdj;

  int offset = 0;
  for (int f = 0; f < int(pCounts.size()); f++) {
    int n = pCounts[f];
    for (int i = 0; i < n; i++) {
      int a = pConnects[offset + i];
      int b = pConnects[offset + (i + 1) % n];
      int64_t key = edgeKey(a, b);

      EdgeAdj &adj = edgeAdj[key];
      if (adj.face0 == -1)
        adj.face0 = f;
      else
        adj.face1 = f;
    }
    offset += n;
  }

  std::unordered_set<int64_t> boundary;
  for (const auto &kv : edgeAdj) {
    const EdgeAdj &adj = kv.second;
    if (adj.face0 >= 0 && adj.face1 >= 0 &&
        panelGroupId[adj.face0] != panelGroupId[adj.face1])
      boundary.insert(kv.first);
  }

  return boundary;
}

// Sequential (0-based) connected-component id per face, cutting at
// cutEdges. Always terminates: bounded, explicit-stack flood fill.
void floodFillGroups(const PrototypeRawMesh &mesh,
                     const std::unordered_set<int64_t> &cutEdges,
                     zIntArray &outFaceGroupId) {
  int numFaces = int(mesh.pCounts.size());
  zIntArray faceStart(numFaces);
  {
    int running = 0;
    for (int f = 0; f < numFaces; f++) {
      faceStart[f] = running;
      running += mesh.pCounts[f];
    }
  }

  struct EdgeAdj {
    int face0 = -1, face1 = -1;
  };
  std::unordered_map<int64_t, EdgeAdj> edgeAdj;

  for (int f = 0; f < numFaces; f++) {
    int n = mesh.pCounts[f];
    int s = faceStart[f];
    for (int i = 0; i < n; i++) {
      int a = mesh.pConnects[s + i];
      int b = mesh.pConnects[s + (i + 1) % n];
      int64_t key = edgeKey(a, b);

      EdgeAdj &adj = edgeAdj[key];
      if (adj.face0 == -1)
        adj.face0 = f;
      else
        adj.face1 = f;
    }
  }

  std::vector<zIntArray> faceAdj(numFaces);
  for (const auto &kv : edgeAdj) {
    if (cutEdges.count(kv.first))
      continue;

    const EdgeAdj &adj = kv.second;
    if (adj.face0 >= 0 && adj.face1 >= 0) {
      faceAdj[adj.face0].push_back(adj.face1);
      faceAdj[adj.face1].push_back(adj.face0);
    }
  }

  outFaceGroupId.assign(numFaces, -1);
  int nextId = 0;
  zIntArray stack;

  for (int f = 0; f < numFaces; f++) {
    if (outFaceGroupId[f] >= 0)
      continue;

    int id = nextId++;
    stack.clear();
    stack.push_back(f);
    outFaceGroupId[f] = id;

    while (!stack.empty()) {
      int cur = stack.back();
      stack.pop_back();

      for (int nb : faceAdj[cur]) {
        if (outFaceGroupId[nb] < 0) {
          outFaceGroupId[nb] = id;
          stack.push_back(nb);
        }
      }
    }
  }
}

// Re-derives a clean, dense-from-0 id per face of the FLAT
// (un-subdivided) mesh, from real connectivity rather than trusting
// the raw (sparse, arbitrary - each its own face index, except
// where computePierPanels folded a pentagon into an existing pier's
// group via addPolygon's groupId) values computeArchPanels /
// computePierPanels happened to assign: findPanelGroupBoundaryEdges
// turns those into boundary edges, then floodFillGroups walks the
// actual face adjacency to give every mutually-reachable patch one sequential
// id. Runs once, before computeBlock, so a caller can colour the
// flat mesh directly by _PanelGroupId, and computeBlock's own
// (post-subdivision) id re-derivation starts from the same verified
// grouping instead of redoing this from scratch.
void computeFaceId(zTsPrototype &self) {
  PrototypeRawMesh flat{self._Positions, self._PCounts, self._PConnects};
  std::unordered_set<int64_t> groupBoundaryEdges = findPanelGroupBoundaryEdges(
      self._PanelGroupId, self._PCounts, self._PConnects);

  zIntArray faceId;
  floodFillGroups(flat, groupBoundaryEdges, faceId);

  self._PanelGroupId = faceId;
}

// Splits the flat mesh into one small, uniformly-coloured mesh per
// _PanelGroupId id, for display purposes - vertices are NOT shared
// between two different ids (each block gets its own local vertex
// copy), unlike colouring a single shared mesh face-by-face, where a
// vertex sitting on the boundary between two ids only has ONE colour
// slot and a renderer that interpolates vertex colours across a face
// blends the two ids' colours together right at their shared edge
// instead of showing a crisp seam.
void computeFaceBlocks(zTsPrototype &self) {
  self.o_FaceBlocks.clear();

  std::unordered_map<int, zIntArray> facesById;
  for (int f = 0; f < int(self._PanelGroupId.size()); f++)
    facesById[self._PanelGroupId[f]].push_back(f);

  zIntArray faceStart;
  {
    int running = 0;
    for (int f = 0; f < int(self._PCounts.size()); f++) {
      faceStart.push_back(running);
      running += self._PCounts[f];
    }
  }

  for (const auto &idFacesKv : facesById) {
    int id = idFacesKv.first;
    const zIntArray &faces = idFacesKv.second;

    std::unordered_map<int, int> globalToLocal;
    zIntArray localToGlobal;

    auto getLocal = [&](int g) -> int {
      auto it = globalToLocal.find(g);
      if (it != globalToLocal.end())
        return it->second;
      int localId = int(localToGlobal.size());
      globalToLocal[g] = localId;
      localToGlobal.push_back(g);
      return localId;
    };

    zIntArray localPCounts, localPConnects;
    for (int f : faces) {
      int n = self._PCounts[f];
      int s = faceStart[f];
      localPCounts.push_back(n);
      for (int i = 0; i < n; i++)
        localPConnects.push_back(getLocal(self._PConnects[s + i]));
    }

    zPointArray localPositions(localToGlobal.size());
    for (int li = 0; li < int(localToGlobal.size()); li++)
      localPositions[li] = self._Positions[localToGlobal[li]];

    zObjMesh idMesh;
    zFnMesh fnId(idMesh);
    fnId.create(localPositions, localPCounts, localPConnects);
    fnId.setFaceColor(colorForId(id), true);

    self.o_FaceBlocks.push_back(idMesh);
  }
}

// Splits the positive offset mesh into one closed solid per panel
// group. _PanelGroupId (see computeFaceId) already holds a clean,
// dense-from-0 id per face of the flat mesh; rather than trusting it
// to propagate correctly through arbitrary subdivision on its own,
// grouping is RE-DERIVED from real connectivity again here:
// findPanelGroupBoundaryEdges turns it into boundary edges, those
// are creased (so a group boundary never bows/dents - see
// findIrregularFaceCreaseEdges) and tracked through exactly one
// subdivision pass (groupBoundaryIn/Out), then floodFillGroups
// re-derives the true ids on that once-subdivided topology, cutting
// only at the propagated boundary edges (not every crease).
// Remaining subdivision levels then just propagate those verified
// ids forward normally.
//
// Only the positive side is ever subdivided. The negative side is
// never built or subdivided as its own mesh - subdivideCatmullClarkCrease
// combines each new point from its neighbours' existing xyz only, so
// negating y commutes with it: subdividing then mirroring y equals
// mirroring y then subdividing. The negative-side positions are
// therefore only ever computed at the very last step, per group, by
// mirroring the (possibly subdivided) positive positions' y.
//
// For each group id: combines that id's faces (positive as-is, plus
// a mirrored-y copy) into one closed solid - front faces as-is, back
// faces reversed winding (so both halves face outward consistently),
// plus one bridge quad per boundary edge of that shared local
// topology - coloured uniformly. Runs entirely in the bridge's own
// LOCAL plane coordinates (self._PositionsPositive is still local
// when this is called, before remapToGlobal), since mirroring y to
// get the negative side only means anything before local and world
// axes have been mixed together - each solid's own combined
// positions are mapped to global (toGlobalPoint) only once they're
// otherwise finished, right before the mesh itself is built.
void computeBlock(zTsPrototype &self) {
  self.o_ClosedMeshes.clear();

  PrototypeRawMesh curPositive{self._PositionsPositive, self._PCounts,
                               self._PConnects};

  std::unordered_set<int64_t> irregularFaceCrease =
      findIrregularFaceCreaseEdges(self._PCounts, self._PConnects);
  std::unordered_set<int64_t> groupBoundaryEdges = findPanelGroupBoundaryEdges(
      self._PanelGroupId, self._PCounts, self._PConnects);

  std::unordered_set<int64_t> creasePositive = irregularFaceCrease;
  creasePositive.insert(groupBoundaryEdges.begin(), groupBoundaryEdges.end());

  zIntArray faceIdPositive;

  if (self.i_Settings.SubdivisionCount <= 0) {
    // no subdivision happens between computeFaceId() and here,
    // so re-flood-filling the same graph/cuts would just
    // reproduce self._PanelGroupId exactly - reuse it directly.
    faceIdPositive = self._PanelGroupId;
  } else {
    // level 0: subdivide once, track only where the ORIGINAL
    // group-boundary edges end up - faceIdIn/Out aren't used for
    // grouping here, just sized correctly.
    zIntArray seedFaceId(int(self._PCounts.size()), 0);

    PrototypeRawMesh nextPositive;
    std::unordered_set<int64_t> nextCreasePositive;
    std::unordered_set<int64_t> nextGroupBoundaryPositive;
    zIntArray discardFaceIdPositive;

    subdivideCatmullClarkCrease(curPositive, creasePositive, seedFaceId,
                                groupBoundaryEdges, nextPositive,
                                nextCreasePositive, discardFaceIdPositive,
                                nextGroupBoundaryPositive);

    curPositive = nextPositive;
    creasePositive = nextCreasePositive;

    floodFillGroups(curPositive, nextGroupBoundaryPositive, faceIdPositive);

    for (int iter = 1; iter < self.i_Settings.SubdivisionCount; iter++) {
      PrototypeRawMesh iterNextPositive;
      std::unordered_set<int64_t> iterNextCreasePositive;
      zIntArray iterNextFaceIdPositive;
      std::unordered_set<int64_t> unusedBoundaryPositive;

      subdivideCatmullClarkCrease(curPositive, creasePositive, faceIdPositive,
                                  {}, iterNextPositive, iterNextCreasePositive,
                                  iterNextFaceIdPositive,
                                  unusedBoundaryPositive);

      curPositive = iterNextPositive;
      creasePositive = iterNextCreasePositive;
      faceIdPositive = iterNextFaceIdPositive;
    }
  }

  // bucket by group id - ids come from floodFillGroups (sequential
  // from 0 over however many groups it found), not from an
  // original face's own index, hence a map, not a dense array.
  std::unordered_map<int, zIntArray> facesById;
  for (int f = 0; f < int(faceIdPositive.size()); f++)
    facesById[faceIdPositive[f]].push_back(f);

  zIntArray faceStart;
  {
    int running = 0;
    for (int f = 0; f < int(curPositive.pCounts.size()); f++) {
      faceStart.push_back(running);
      running += curPositive.pCounts[f];
    }
  }

  for (const auto &idFacesKv : facesById) {
    int id = idFacesKv.first;
    const zIntArray &faces = idFacesKv.second;

    // remap this id's vertex indices to a local 0-based range -
    // shared by both sides, since the mirrored negative side
    // below reuses the exact same face->vertex-index structure
    std::unordered_map<int, int> globalToLocal;
    zIntArray localToGlobal;

    auto getLocal = [&](int g) -> int {
      auto it = globalToLocal.find(g);
      if (it != globalToLocal.end())
        return it->second;
      int localId = int(localToGlobal.size());
      globalToLocal[g] = localId;
      localToGlobal.push_back(g);
      return localId;
    };

    zIntArray localPCounts, localPConnects;
    for (int f : faces) {
      int n = curPositive.pCounts[f];
      int s = faceStart[f];
      localPCounts.push_back(n);
      for (int i = 0; i < n; i++)
        localPConnects.push_back(getLocal(curPositive.pConnects[s + i]));
    }

    int localCount = int(localToGlobal.size());

    // the negative side is only ever computed here, at the very
    // last step, by mirroring the positive position's y - see
    // the comment above computeBlock for why that's valid.
    zPointArray combinedPositions(localCount * 2);
    for (int li = 0; li < localCount; li++) {
      int g = localToGlobal[li];
      const zPoint &pos = curPositive.positions[g];
      combinedPositions[li] = pos;
      combinedPositions[localCount + li] = zPoint(pos.x, -pos.y, pos.z);
    }

    // boundary edges of this id's own local patch (used by
    // exactly one face) - direction stored is the FRONT winding,
    // used to orient the bridge quad against it.
    struct BoundaryEdgeInfo {
      int count = 0;
      int a = -1, b = -1;
    };
    std::unordered_map<int64_t, BoundaryEdgeInfo> edgeInfo;
    {
      int idx = 0;
      for (int n : localPCounts) {
        for (int i = 0; i < n; i++) {
          int a = localPConnects[idx + i];
          int b = localPConnects[idx + (i + 1) % n];
          int64_t key = edgeKey(a, b);

          BoundaryEdgeInfo &info = edgeInfo[key];
          if (info.count == 0) {
            info.a = a;
            info.b = b;
          }
          info.count++;
        }
        idx += n;
      }
    }

    zIntArray combinedPCounts;
    zIntArray combinedPConnects;

    // front faces (positive side, unchanged winding)
    {
      int idx = 0;
      for (int n : localPCounts) {
        combinedPCounts.push_back(n);
        for (int i = 0; i < n; i++)
          combinedPConnects.push_back(localPConnects[idx + i]);
        idx += n;
      }
    }

    // back faces (negative side, reversed winding, so front/back
    // face outward consistently and the wall below stays
    // correctly oriented)
    {
      int idx = 0;
      for (int n : localPCounts) {
        combinedPCounts.push_back(n);
        for (int i = n - 1; i >= 0; i--)
          combinedPConnects.push_back(localCount + localPConnects[idx + i]);
        idx += n;
      }
    }

    // bridge quads: a consistently-oriented manifold needs every
    // edge shared by two faces walked in opposite directions, so
    // the wall starts at b -> a (opposite the front face's a -> b)
    for (const auto &kv : edgeInfo) {
      const BoundaryEdgeInfo &info = kv.second;
      if (info.count != 1)
        continue; // interior edge - nothing to bridge

      int a = info.a;
      int b = info.b;

      combinedPCounts.push_back(4);
      combinedPConnects.push_back(b);
      combinedPConnects.push_back(a);
      combinedPConnects.push_back(localCount + a);
      combinedPConnects.push_back(localCount + b);
    }

    for (zPoint &p : combinedPositions)
      p = toGlobalPoint(self, p);

    zObjMesh idMesh;
    zFnMesh fnId(idMesh);
    fnId.create(combinedPositions, combinedPCounts, combinedPConnects);
    fnId.setFaceColor(colorForId(id), true);

    self.o_ClosedMeshes.push_back(idMesh);
  }
}
} // anonymous namespace

//---------------------------------------------------------------------------
//---- COMPUTE METHOD
//---------------------------------------------------------------------------

ZSPACE_TOOLSETS_INLINE bool zTsPrototype::compute() {
  if (!validateSettings(*this))
    return false;

  _Width = 0.0f;
  _Height = 0.0f;
  _LocalGraphPositions.clear();
  _GraphEdgeType.clear();
  _Piers.clear();
  _Spans.clear();
  _Positions.clear();
  _PCounts.clear();
  _PConnects.clear();
  _PanelGroupId.clear();
  _PierGroupId.clear();
  _BlockType.clear();
  _Thickness.clear();

  computePlane(*this);
  classifyGraphEdges(*this);
  computePiers(*this);
  if (!validatePiers(*this))
    return false;

  computeSpans(*this);
  computeSpanBoundaryVertices(*this);
  computeArchPanels(*this);
  computePierPanels(*this);
  computeFaceId(*this);
  computeOffsetMeshes(*this);
  computeBlock(*this); // must run before remapToGlobal - see there
  remapToGlobal(*this);
  buildMesh(*this);
  buildMeshPositive(*this);
  computeFaceBlocks(*this);

  _Status = "Ready";
  return true;
}
} // namespace zSpace
