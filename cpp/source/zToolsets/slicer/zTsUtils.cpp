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
#include <iostream>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

namespace zSpace {
namespace zTsUtils {

//---- Thread-Safe Console Logging

namespace {
// Private: shared mutex guarding console output.
std::mutex &_getConsoleLogMutex() {
  static std::mutex mtx;
  return mtx;
}
} // anonymous namespace

//--------------------------
//---- PUBLIC FUNCTIONS
//--------------------------

ZSPACE_INLINE void printSuccess(const std::string &msg) {
  std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
  std::cout << "\033[1;32m" << msg << "\033[0m\n" << std::flush;
}

ZSPACE_INLINE void printError(const std::string &msg) {
  std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
  std::cerr << "\033[1;31m" << msg << "\033[0m\n" << std::flush;
}

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

ZSPACE_INLINE bool turnAngleAt(const zPoint &prev, const zPoint &cur,
                               const zPoint &next, double &outAngleRad) {
  const zVector v0 = vecSub(prev, cur);
  const zVector v1 = vecSub(next, cur);
  const double len0 = v0.length();
  const double len1 = v1.length();
  if (len0 < 1e-6 || len1 < 1e-6)
    return false;

  double cosTheta = vecDot(v0, v1) / (len0 * len1);
  cosTheta = (std::max)(-1.0, (std::min)(1.0, cosTheta));
  outAngleRad = std::acos(cosTheta);
  return true;
}

ZSPACE_INLINE double signedArea2D(const zPointArray &positions) {
  double area = 0.0;
  const int n = static_cast<int>(positions.size());
  for (int i = 0; i < n; i++) {
    const zPoint &p0 = positions[i];
    const zPoint &p1 = positions[(i + 1) % n];
    area += (p0.x * p1.y - p1.x * p0.y);
  }
  return 0.5 * area;
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

//---- Single-Chain Graph Closed/Open Detection

ZSPACE_INLINE bool isClosedGraph(const zObjGraph &graph) {
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  const int N = fn.numVertices();
  if (N < 3)
    return false;
  zIntArray edges;
  fn.getEdgeData(edges);
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    if ((edges[e] == 0 && edges[e + 1] == N - 1) ||
        (edges[e + 1] == 0 && edges[e] == N - 1))
      return true;
  }
  return false;
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

ZSPACE_INLINE bool parameterizeGraph(const zObjGraph &inGraph, zObjMesh &inMesh,
                                     ParameterGraph &outParams) {
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

ZSPACE_INLINE bool realizeParameterGraph(const ParameterGraph &inParams,
                                         zObjMesh &targetMesh,
                                         zObjGraph &outGraph,
                                         zVectorArray *outNormals,
                                         const zFloatArray *vertexHeights,
                                         zFloatArray *outHeights) {
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
  const bool wantHeights = vertexHeights && outHeights;

  for (size_t i = 0; i < numVerts; i++) {
    const auto &bv = inParams.vertices[i];
    if (bv.v0 < 0 || bv.v0 >= (int)numMeshVerts || bv.v1 < 0 ||
        bv.v1 >= (int)numMeshVerts || bv.v2 < 0 || bv.v2 >= (int)numMeshVerts) {
      return false;
    }

    const zPoint &p0 = targetPositions[bv.v0];
    const zPoint &p1 = targetPositions[bv.v1];
    const zPoint &p2 = targetPositions[bv.v2];

    newPositions[i] = zPoint(p0.x * bv.u + p1.x * bv.v + p2.x * bv.w,
                             p0.y * bv.u + p1.y * bv.v + p2.y * bv.w,
                             p0.z * bv.u + p1.z * bv.v + p2.z * bv.w);

    if (wantHeights) {
      float h0 =
          (bv.v0 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v0] : 0.0f;
      float h1 =
          (bv.v1 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v1] : 0.0f;
      float h2 =
          (bv.v2 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v2] : 0.0f;
      (*outHeights)[i] = bv.u * h0 + bv.v * h1 + bv.w * h2;
    }

    if (outNormals) {
      zVector e0 = vecSub(p1, p0);
      zVector e1 = vecSub(p2, p0);
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

ZSPACE_INLINE bool remapGraph(zObjGraph &graph, zObjMesh &inMesh2D,
                              zObjMesh &projectionMesh3D,
                              zVectorArray *mappedNormals,
                              const zFloatArray *vertexHeights,
                              zFloatArray *mappedHeights) {
  ParameterGraph tempParams;
  if (!parameterizeGraph(graph, inMesh2D, tempParams))
    return false;
  return realizeParameterGraph(tempParams, projectionMesh3D, graph,
                               mappedNormals, vertexHeights, mappedHeights);
}

//---- Toolpath Pipe Mesh Construction

ZSPACE_INLINE void buildPathMesh(const PrintPath &path,
                                 float defaultLayerHeight,
                                 zObjMesh &outMesh) {
  zFnMesh fnOut(outMesh);
  fnOut.clear();

  zFnGraph fnContour(const_cast<zObjGraph &>(path.graph));
  zPointArray points;
  fnContour.getVertexPositions(points);
  const int N = int(points.size());
  if (N < 2)
    return;

  const bool closed = isClosedGraph(path.graph);

  const zFloatArray *layerWeights = &path.widths;
  const zFloatArray *layerHeights = &path.heights;
  const zVectorArray *layerNormals = &path.normals;

  zPointArray meshPositions;
  zIntArray polyCounts;
  zIntArray polyConnects;
  meshPositions.reserve(N * 4);

  for (int i = 0; i < N; i++) {
    zPoint p = points[i];

    // 1. Tangent direction
    zVector pathDir;
    if (closed) {
      int prevIdx = (i - 1 + N) % N;
      int nextIdx = (i + 1) % N;
      pathDir = zVector(points[nextIdx].x - points[prevIdx].x,
                        points[nextIdx].y - points[prevIdx].y,
                        points[nextIdx].z - points[prevIdx].z);
    } else {
      if (i == 0) {
        pathDir = zVector(points[1].x - points[0].x, points[1].y - points[0].y,
                          points[1].z - points[0].z);
      } else if (i == N - 1) {
        pathDir = zVector(points[N - 1].x - points[N - 2].x,
                          points[N - 1].y - points[N - 2].y,
                          points[N - 1].z - points[N - 2].z);
      } else {
        pathDir = zVector(points[i + 1].x - points[i - 1].x,
                          points[i + 1].y - points[i - 1].y,
                          points[i + 1].z - points[i - 1].z);
      }
    }
    if (pathDir.length() < 1e-6)
      pathDir = zVector(1, 0, 0);
    pathDir.normalize();

    // 2. 3D surface normal
    zVector norm = (layerNormals && i < int(layerNormals->size()))
                       ? (*layerNormals)[i]
                       : zVector(0, 0, 1);
    if (norm.length() < 1e-6)
      norm = zVector(0, 0, 1);
    norm.normalize();

    // 3. In-plane binormal perpendicular to tangent
    zVector binormal = norm ^ pathDir;
    if (binormal.length() < 1e-6) {
      zVector temp(0, 0, 1);
      if (std::abs(pathDir.z) > 0.9)
        temp = zVector(0, 1, 0);
      binormal = temp ^ pathDir;
    }
    binormal.normalize();

    // 4. Per-vertex width and height
    float w = (layerWeights && i < int(layerWeights->size()))
                  ? (*layerWeights)[i]
                  : 0.20f;
    if (w <= 1e-4f)
      w = 0.20f;
    float h = (layerHeights && i < int(layerHeights->size()))
                  ? (*layerHeights)[i]
                  : defaultLayerHeight;
    if (h <= 1e-4f)
      h = defaultLayerHeight;

    // 5. Generate 4 corner vertices for rectangular cross section
    const float halfW = w * 0.5f;
    zPoint c0(p.x + binormal.x * halfW, p.y + binormal.y * halfW,
              p.z + binormal.z * halfW); // Top-Right
    zPoint c1(p.x - binormal.x * halfW, p.y - binormal.y * halfW,
              p.z - binormal.z * halfW); // Top-Left
    zPoint c2(c1.x - norm.x * h, c1.y - norm.y * h,
              c1.z - norm.z * h); // Bottom-Left
    zPoint c3(c0.x - norm.x * h, c0.y - norm.y * h,
              c0.z - norm.z * h); // Bottom-Right

    meshPositions.push_back(c0);
    meshPositions.push_back(c1);
    meshPositions.push_back(c2);
    meshPositions.push_back(c3);
  }

  // 6. Construct 4 quads per segment along the toolpath tube
  const int segCount = closed ? N : (N - 1);
  for (int i = 0; i < segCount; i++) {
    int next = (i + 1) % N;
    // Top face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 0);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(i * 4 + 1);

    // Left face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 1);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(i * 4 + 2);

    // Bottom face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 2);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(i * 4 + 3);

    // Right face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 3);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(i * 4 + 0);
  }

  // Cap ends for open curves
  if (!closed && N >= 2) {
    polyCounts.push_back(4);
    polyConnects.push_back(0);
    polyConnects.push_back(3);
    polyConnects.push_back(2);
    polyConnects.push_back(1);

    int last = (N - 1) * 4;
    polyCounts.push_back(4);
    polyConnects.push_back(last + 0);
    polyConnects.push_back(last + 1);
    polyConnects.push_back(last + 2);
    polyConnects.push_back(last + 3);
  }

  fnOut.create(meshPositions, polyCounts, polyConnects);
}

} // namespace zTsUtils
} // namespace zSpace
