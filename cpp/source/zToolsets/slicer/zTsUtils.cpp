// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Weijie Tao <weijie.tao@zspace.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/slicer/zTsUtils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace zSpace {
namespace zTsUtils {

//---- Basic Vector Utilities

ZSPACE_INLINE zVector vecSub(const zPoint &a, const zPoint &b) {
  return zVector(a.x - b.x, a.y - b.y, a.z - b.z);
}

ZSPACE_INLINE double vecDistSq(const zPoint &a, const zPoint &b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

ZSPACE_INLINE double vecDot(const zVector &a, const zVector &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

ZSPACE_INLINE zVector vecCross(const zVector &a, const zVector &b) {
  return zVector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x);
}

ZSPACE_INLINE unsigned long long makeEdgeKey(int u, int v) {
  const int lo = (std::min)(u, v);
  const int hi = (std::max)(u, v);
  return (static_cast<unsigned long long>(lo) << 32) |
         static_cast<unsigned int>(hi);
}

//---- HalfEdge Mesh Topology Walk Utilities

ZSPACE_INLINE zItMeshHalfEdge nextEdgeRing(zItMeshHalfEdge he) {
  return he.getNext().getNext().getSym();
}

ZSPACE_INLINE zItMeshHalfEdge nextEdgeLoop(zItMeshHalfEdge he) {
  return he.getNext().getSym().getNext();
}

ZSPACE_INLINE zItMeshHalfEdge zUpEdge(zItMeshVertex v) {
  zItMeshHalfEdgeArray hEdges;
  v.getConnectedHalfEdges(hEdges);
  if (hEdges.empty())
    return zItMeshHalfEdge();

  return *std::max_element(hEdges.begin(), hEdges.end(),
                           [](zItMeshHalfEdge a, zItMeshHalfEdge b) {
                             return a.getVector().z < b.getVector().z;
                           });
}

ZSPACE_INLINE void getFaceVerticesFromHalfedge(zItMeshHalfEdge &heStart,
                                               bool forward,
                                               zIntArray &faceVerts) {
  faceVerts.clear();
  zItMeshHalfEdge he = heStart;
  do {
    faceVerts.push_back(forward ? he.getVertex().getId()
                                : he.getStartVertex().getId());
    he = forward ? he.getNext() : he.getPrev();
  } while (he != heStart);
}

//---- Combine Multiple Graph Objects

ZSPACE_INLINE void combineGraphObjects(const zObjGraphArray &sources,
                                       zObjGraph &outGraph) {
  zPointArray positions;
  zIntArray edgeConnects;
  for (const zObjGraph &source : sources) {
    zFnGraph fnSource(const_cast<zObjGraph &>(source));
    zPointArray srcPos;
    zIntArray srcEdges;
    fnSource.getVertexPositions(srcPos);
    fnSource.getEdgeData(srcEdges);
    const int baseIndex = static_cast<int>(positions.size());
    positions.insert(positions.end(), srcPos.begin(), srcPos.end());
    for (int id : srcEdges)
      edgeConnects.push_back(baseIndex + id);
  }
  zFnGraph fnOut(outGraph);
  fnOut.clear();
  if (!positions.empty())
    fnOut.create(positions, edgeConnects);
}

//---- Universal Polyline Sampling by Arc-Length with Optional Tangent Vector

ZSPACE_INLINE bool samplePolyline(const zPoint *positions,
                                  const zScalar *scalars, size_t count,
                                  double t, zPoint &outPoint,
                                  zVector *outTangent) {
  if (!positions || !scalars || count < 1)
    return false;

  if (count == 1 || t <= 0.0) {
    outPoint = positions[0];
    if (outTangent && count > 1) {
      *outTangent = vecSub(positions[1], positions[0]);
      if (outTangent->length() > 1e-6)
        outTangent->normalize();
    }
    return true;
  }
  if (t >= 1.0) {
    outPoint = positions[count - 1];
    if (outTangent && count > 1) {
      *outTangent = vecSub(positions[count - 1], positions[count - 2]);
      if (outTangent->length() > 1e-6)
        outTangent->normalize();
    }
    return true;
  }

  const double targetDist = t * scalars[count - 1];
  auto it = std::upper_bound(scalars, scalars + count, targetDist);
  size_t s = (it == scalars) ? 0 : std::distance(scalars, it) - 1;
  if (s + 1 >= count) {
    outPoint = positions[count - 1];
    return true;
  }

  const double segLen = scalars[s + 1] - scalars[s];
  const double localT =
      (segLen > 1e-6) ? (targetDist - scalars[s]) / segLen : 0.0;
  const zPoint &a = positions[s];
  const zPoint &b = positions[s + 1];
  outPoint.x = a.x + (b.x - a.x) * localT;
  outPoint.y = a.y + (b.y - a.y) * localT;
  outPoint.z = a.z + (b.z - a.z) * localT;

  if (outTangent) {
    *outTangent = vecSub(b, a);
    if (outTangent->length() > 1e-6)
      outTangent->normalize();
  }
  return true;
}

ZSPACE_INLINE bool samplePolyline(const zObjGraph &polylineGraph,
                                  const zScalarArray &scalars, double t,
                                  zPoint &outPoint, zVector *outTangent) {
  zFnGraph fnGraph(const_cast<zObjGraph &>(polylineGraph));
  return samplePolyline(fnGraph.getRawVertexPositions(), scalars.data(),
                        scalars.size(), t, outPoint, outTangent);
}

//---- Universal Mesh Barycentric Parameterizer (2D flat and 3D space meshes)

ZSPACE_INLINE bool parameterizeGraphOnMesh(const zObjGraph &inGraph,
                                           zObjMesh &inMesh,
                                           BracingGraph &outParams) {
  outParams.clear();
  zFnGraph fnGraph(const_cast<zObjGraph &>(inGraph));
  zPointArray graphPositions;
  fnGraph.getVertexPositions(graphPositions);
  fnGraph.getEdgeData(outParams.edgeConnects);
  if (graphPositions.empty() || outParams.edgeConnects.empty())
    return false;

  zFnMesh fnMesh(inMesh);
  zPointArray meshPositions;
  fnMesh.getVertexPositions(meshPositions);
  if (meshPositions.empty() || fnMesh.numPolygons() == 0)
    return false;

  struct TriCache {
    int v0, v1, v2;
    zPoint p0, p1, p2;
    zVector e0, e1;
    double d00, d01, d11, denom;
    bool valid;
  };

  std::vector<TriCache> allTris;
  allTris.reserve(fnMesh.numPolygons() * 2);

  for (zItMeshFace f(inMesh); !f.end(); f++) {
    zIntArray faceVIds;
    f.getVertices(faceVIds);
    if (faceVIds.size() < 3)
      continue;

    for (size_t tri = 1; tri + 1 < faceVIds.size(); tri++) {
      TriCache t;
      t.v0 = faceVIds[0];
      t.v1 = faceVIds[tri];
      t.v2 = faceVIds[tri + 1];
      t.p0 = meshPositions[t.v0];
      t.p1 = meshPositions[t.v1];
      t.p2 = meshPositions[t.v2];

      t.e0 = vecSub(t.p1, t.p0);
      t.e1 = vecSub(t.p2, t.p0);
      t.d00 = vecDot(t.e0, t.e0);
      t.d01 = vecDot(t.e0, t.e1);
      t.d11 = vecDot(t.e1, t.e1);
      t.denom = t.d00 * t.d11 - t.d01 * t.d01;
      t.valid = (std::abs(t.denom) > 1e-12);
      allTris.push_back(t);
    }
  }

  if (allTris.empty())
    return false;

  const size_t numGraphVerts = graphPositions.size();
  outParams.vertices.resize(numGraphVerts);

  for (size_t i = 0; i < numGraphVerts; i++) {
    const zPoint &P = graphPositions[i];
    double bestDistSq = (std::numeric_limits<double>::max)();
    int bestV0 = allTris[0].v0, bestV1 = allTris[0].v1, bestV2 = allTris[0].v2;
    double bestU = 1.0, bestV = 0.0, bestW = 0.0;
    bool foundInside = false;

    for (const auto &t : allTris) {
      if (!t.valid)
        continue;

      const zVector v2 = vecSub(P, t.p0);
      const double d20 = vecDot(v2, t.e0);
      const double d21 = vecDot(v2, t.e1);
      const double v = (t.d11 * d20 - t.d01 * d21) / t.denom;
      const double w = (t.d00 * d21 - t.d01 * d20) / t.denom;
      const double u = 1.0 - v - w;

      if (u >= -1e-4 && v >= -1e-4 && w >= -1e-4) {
        const double cu = (std::max)(0.0, u);
        const double cv = (std::max)(0.0, v);
        const double cw = (std::max)(0.0, w);
        const double csum = cu + cv + cw;
        const double nu = cu / csum, nv = cv / csum, nw = cw / csum;

        zPoint proj(t.p0.x * nu + t.p1.x * nv + t.p2.x * nw,
                    t.p0.y * nu + t.p1.y * nv + t.p2.y * nw,
                    t.p0.z * nu + t.p1.z * nv + t.p2.z * nw);
        const double distSq = vecDistSq(P, proj);
        if (distSq < bestDistSq) {
          bestDistSq = distSq;
          bestV0 = t.v0;
          bestV1 = t.v1;
          bestV2 = t.v2;
          bestU = nu;
          bestV = nv;
          bestW = nw;
          foundInside = true;
          if (distSq < 1e-8)
            break;
        }
      } else if (!foundInside) {
        const double cu = (std::max)(0.0, (std::min)(1.0, u));
        const double cv = (std::max)(0.0, (std::min)(1.0, v));
        const double cw = (std::max)(0.0, (std::min)(1.0, w));
        const double csum = cu + cv + cw;
        if (csum > 1e-6) {
          const double nu = cu / csum, nv = cv / csum, nw = cw / csum;
          zPoint proj(t.p0.x * nu + t.p1.x * nv + t.p2.x * nw,
                      t.p0.y * nu + t.p1.y * nv + t.p2.y * nw,
                      t.p0.z * nu + t.p1.z * nv + t.p2.z * nw);
          const double distSq = vecDistSq(P, proj);
          if (distSq < bestDistSq) {
            bestDistSq = distSq;
            bestV0 = t.v0;
            bestV1 = t.v1;
            bestV2 = t.v2;
            bestU = nu;
            bestV = nv;
            bestW = nw;
          }
        }
      }
    }

    outParams.vertices[i] = {bestV0,
                             bestV1,
                             bestV2,
                             static_cast<float>(bestU),
                             static_cast<float>(bestV),
                             static_cast<float>(bestW)};
  }

  return true;
}

//---- Universal Barycentric Mesh Graph Reconstructor

ZSPACE_INLINE bool reconstructGraphFromBarycentricParams(
    const BracingGraph &inParams, zObjMesh &targetMesh,
    zObjGraph &outGraph, zVectorArray *outNormals,
    const zFloatArray *vertexScalars, zFloatArray *outHeights) {
  if (inParams.empty())
    return false;

  zFnMesh fnTarget(targetMesh);
  zPointArray targetPositions;
  fnTarget.getVertexPositions(targetPositions);
  if (targetPositions.empty())
    return false;

  const size_t numVerts = inParams.vertices.size();
  zPointArray newPositions(numVerts);
  if (outNormals)
    outNormals->assign(numVerts, zVector(0, 0, 1));
  if (outHeights)
    outHeights->assign(numVerts, 0.0f);

  const size_t numMeshVerts = targetPositions.size();

  for (size_t i = 0; i < numVerts; i++) {
    const auto &bv = inParams.vertices[i];
    if (bv.v0 < 0 || bv.v0 >= (int)numMeshVerts || bv.v1 < 0 ||
        bv.v1 >= (int)numMeshVerts || bv.v2 < 0 || bv.v2 >= (int)numMeshVerts) {
      return false;
    }

    newPositions[i].x = targetPositions[bv.v0].x * bv.u +
                        targetPositions[bv.v1].x * bv.v +
                        targetPositions[bv.v2].x * bv.w;
    newPositions[i].y = targetPositions[bv.v0].y * bv.u +
                        targetPositions[bv.v1].y * bv.v +
                        targetPositions[bv.v2].y * bv.w;
    newPositions[i].z = targetPositions[bv.v0].z * bv.u +
                        targetPositions[bv.v1].z * bv.v +
                        targetPositions[bv.v2].z * bv.w;

    if (vertexScalars && outHeights) {
      float h0 =
          (bv.v0 < (int)vertexScalars->size()) ? (*vertexScalars)[bv.v0] : 0.0f;
      float h1 =
          (bv.v1 < (int)vertexScalars->size()) ? (*vertexScalars)[bv.v1] : 0.0f;
      float h2 =
          (bv.v2 < (int)vertexScalars->size()) ? (*vertexScalars)[bv.v2] : 0.0f;
      (*outHeights)[i] = bv.u * h0 + bv.v * h1 + bv.w * h2;
    }
  }

  if (outNormals) {
    for (size_t i = 0; i < numVerts; i++) {
      const auto &bv = inParams.vertices[i];
      zVector e0 = vecSub(targetPositions[bv.v1], targetPositions[bv.v0]);
      zVector e1 = vecSub(targetPositions[bv.v2], targetPositions[bv.v0]);
      zVector norm = vecCross(e0, e1);
      if (norm.length() > 1e-6)
        norm.normalize();
      else
        norm = zVector(0, 0, 1);
      (*outNormals)[i] = norm;
    }
  }

  zFnGraph fnOut(outGraph);
  fnOut.clear();
  fnOut.create(newPositions, const_cast<zIntArray &>(inParams.edgeConnects));
  return true;
}

//---- Barycentric Mesh-to-Mesh Graph Projection

ZSPACE_INLINE bool
barycentricProjectionTriMesh(zObjGraph &graph, zObjMesh &inMesh2D,
                             zObjMesh &projectionMesh3D,
                             zVectorArray *mappedNormals,
                             const zFloatArray *vertexHeights,
                             zFloatArray *mappedHeights) {
  BracingGraph tempParams;
  if (!parameterizeGraphOnMesh(graph, inMesh2D, tempParams))
    return false;
  return reconstructGraphFromBarycentricParams(tempParams, projectionMesh3D,
                                               graph, mappedNormals,
                                               vertexHeights, mappedHeights);
}

} // namespace zTsUtils
} // namespace zSpace
