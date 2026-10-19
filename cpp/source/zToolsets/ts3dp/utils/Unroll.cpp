// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/ts3dp/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace zSpace {

namespace zTsUnroll {

//---- Face-by-face unfolding (TopologySlicer)

ZSPACE_INLINE bool bsfUnfoldOrder(zObjMesh &mesh, zIntPairArray &outHalfEdgePairs) {
  outHalfEdgePairs.clear();
  zFnMesh fnMesh(mesh);
  if (fnMesh.numPolygons() <= 1)
    return false;

  zObjGraph graph;
  zIntArray inEdgeDualEdge;
  zIntArray dualEdgeInEdge;
  fnMesh.getDualGraph(graph, inEdgeDualEdge, dualEdgeInEdge, true, false,
                      false);

  int maxValence = -1;
  zItGraphVertex maxVertex;
  for (zItGraphVertex v(graph); !v.end(); v++) {
    if (v.getValence() > maxValence) {
      maxValence = v.getValence();
      maxVertex = v;
    }
  }
  if (maxValence < 0)
    return false;

  zItGraphVertexArray bsfVertices;
  zIntPairArray bsfFacePairs;
  maxVertex.getBSF(bsfVertices, bsfFacePairs);
  if (static_cast<int>(bsfFacePairs.size()) + 1 != fnMesh.numPolygons())
    return false;

  outHalfEdgePairs.reserve(bsfFacePairs.size());
  for (const auto &facePair : bsfFacePairs) {
    zItMeshFace f1(mesh, facePair.first);
    zItMeshFace f2(mesh, facePair.second);
    zItMeshHalfEdgeArray f1HE, f2HE;
    f1.getHalfEdges(f1HE);
    f2.getHalfEdges(f2HE);

    zIntPair commonHE(-1, -1);
    for (auto &he1 : f1HE) {
      for (auto &he2 : f2HE) {
        if (he1.getEdge().getId() == he2.getEdge().getId()) {
          commonHE = zIntPair(he1.getId(), he2.getId());
          break;
        }
      }
      if (commonHE.first != -1)
        break;
    }

    if (commonHE.first == -1)
      return false;
    outHalfEdgePairs.push_back(commonHE);
  }

  return true;
}

ZSPACE_INLINE bool unfoldMesh(zObjMesh &sourceMesh, const zIntPairArray &bsfHalfEdgePairs,
                              zObjMesh &outUnrollMesh) {
  zFnMesh fnSource(sourceMesh);
  const int numVerts = fnSource.numVertices();
  const int numFaces = fnSource.numPolygons();
  if (numVerts == 0 || numFaces == 0)
    return false;

  zPoint *srcPositions = fnSource.getRawVertexPositions();
  if (!srcPositions)
    return false;

  zPointArray unrollPositions(numVerts, zPoint(0, 0, 0));
  std::vector<bool> isPlaced(numVerts, false);

  auto placeFirstEdge = [&](zItMeshHalfEdge he) -> bool {
    const float edgeLength = he.getLength();
    if (edgeLength <= 1e-3f)
      return false;
    const int pStartId = he.getStartVertex().getId();
    const int pEndId = he.getVertex().getId();
    unrollPositions[pStartId] = zPoint(0, 0, 0);
    unrollPositions[pEndId] = zPoint(0, edgeLength, 0);
    isPlaced[pStartId] = true;
    isPlaced[pEndId] = true;
    return true;
  };

  auto unrollFace = [&](zItMeshFace face, zItMeshHalfEdge baseHE) -> bool {
    const float edgeLength = baseHE.getLength();
    const int pStartId = baseHE.getStartVertex().getId();
    const int pEndId = baseHE.getVertex().getId();
    if (edgeLength <= 1e-3f || !isPlaced[pStartId] || !isPlaced[pEndId])
      return false;

    const zPoint A = srcPositions[pStartId], B = srcPositions[pEndId];
    const zPoint a = unrollPositions[pStartId], b = unrollPositions[pEndId];
    const float denom = edgeLength * edgeLength;
    zVector ba(B.x - A.x, B.y - A.y, B.z - A.z);

    zItMeshHalfEdge walker = baseHE;
    for (int v = 0; v < face.getNumVertices(); v++) {
      walker = walker.getNext();
      const int vid = walker.getVertex().getId();
      if (isPlaced[vid])
        continue;

      const zPoint C = srcPositions[vid];
      zVector ca(C.x - A.x, C.y - A.y, C.z - A.z);
      const float s = static_cast<float>((ba ^ ca).length() / denom);
      const float c = static_cast<float>((ba * ca) / denom);
      unrollPositions[vid] = zPoint(a.x + c * (b.x - a.x) - s * (b.y - a.y),
                                  a.y + c * (b.y - a.y) + s * (b.x - a.x), 0.0);
      isPlaced[vid] = true;
    }
    return true;
  };

  if (numFaces == 1) {
    zItMeshFace face(sourceMesh, 0);
    zItMeshHalfEdge he = face.getHalfEdge();
    if (!placeFirstEdge(he) || !unrollFace(face, he))
      return false;
  } else {
    if (bsfHalfEdgePairs.empty())
      return false;

    for (size_t pairId = 0; pairId < bsfHalfEdgePairs.size(); pairId++) {
      zItMeshHalfEdge parentHE(sourceMesh, bsfHalfEdgePairs[pairId].first);
      zItMeshHalfEdge childHE(sourceMesh, bsfHalfEdgePairs[pairId].second);

      if (pairId == 0) {
        if (!placeFirstEdge(parentHE) ||
            !unrollFace(parentHE.getFace(), parentHE))
          return false;
      }

      if (!unrollFace(childHE.getFace(), childHE))
        return false;
    }
  }

  outUnrollMesh = sourceMesh;
  zFnMesh fnUnrolled(outUnrollMesh);
  fnUnrolled.setVertexPositions(unrollPositions);
  return true;
}

//---- As-rigid-as-possible unrolling (HarmonicSlicer)

namespace {

// arapUnroll() on a triangle list (polygons already fanned).
bool arapSolveTriangles(const zPointArray &positions3D, const std::vector<std::array<int, 3>> &triangles,
                        const zPointArray &init2D, int iterations, zPointArray &out2D, ArapStats &stats) {
  const int n = static_cast<int>(positions3D.size());
  out2D = init2D;
  stats = ArapStats();
  if (n < 3 || triangles.empty() || static_cast<int>(init2D.size()) != n) return false;

  double initArea = 0.0;
  for (const auto &t : triangles) {
    const zPoint &a = init2D[t[0]], &b = init2D[t[1]], &c = init2D[t[2]];
    initArea += (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
  }
  const double handedness = (initArea < 0.0) ? -1.0 : 1.0;

  // Local isometric copy of every triangle; edge e runs corner e -> e + 1,
  // weighted by half the cotangent of the angle at corner e + 2 (clamped
  // positive so the system stays well posed on obtuse/degenerate triangles).
  struct LocalTriangle {
    bool valid = false;
    Eigen::Vector2d x[3];
    double w[3] = {0.0, 0.0, 0.0};
    double area = 0.0;
  };
  std::vector<LocalTriangle> local(triangles.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(triangles.size() * 12 + n);
  double diagonalSum = 0.0;
  for (size_t t = 0; t < triangles.size(); ++t) {
    const Eigen::Vector3d p[3] = {zTsUtils::toEigenVec(positions3D[triangles[t][0]]), zTsUtils::toEigenVec(positions3D[triangles[t][1]]),
                                  zTsUtils::toEigenVec(positions3D[triangles[t][2]])};
    const Eigen::Vector3d e01 = p[1] - p[0], e02 = p[2] - p[0];
    const double l01 = e01.norm();
    const double crossNorm = e01.cross(e02).norm();
    if (l01 < 1e-9 || crossNorm < 1e-8) continue;
    LocalTriangle &lt = local[t];
    lt.x[0] = Eigen::Vector2d(0.0, 0.0);
    lt.x[1] = Eigen::Vector2d(l01, 0.0);
    lt.x[2] = Eigen::Vector2d(e01.dot(e02) / l01, handedness * crossNorm / l01);
    lt.area = 0.5 * crossNorm;
    for (int e = 0; e < 3; ++e) {
      const int k = (e + 2) % 3;
      const Eigen::Vector3d u = p[e] - p[k], v = p[(e + 1) % 3] - p[k];
      const double cr = u.cross(v).norm();
      const double cot = (cr > 1e-8) ? u.dot(v) / cr : 0.0;
      lt.w[e] = 0.5 * std::clamp(cot, 1e-6, 1e3);
      const int i = triangles[t][e], j = triangles[t][(e + 1) % 3];
      triplets.emplace_back(i, i, lt.w[e]);
      triplets.emplace_back(j, j, lt.w[e]);
      triplets.emplace_back(i, j, -lt.w[e]);
      triplets.emplace_back(j, i, -lt.w[e]);
      diagonalSum += 2.0 * lt.w[e];
    }
    lt.valid = true;
  }
  const double anchor = 1e-8 * (std::max)(diagonalSum / n, 1e-12);
  for (int i = 0; i < n; ++i) triplets.emplace_back(i, i, anchor);

  Eigen::SparseMatrix<double> laplacian(n, n);
  laplacian.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(laplacian);
  if (solver.info() != Eigen::Success) return false;

  Eigen::MatrixXd uv(n, 2);
  for (int i = 0; i < n; ++i) {
    uv(i, 0) = init2D[i].x;
    uv(i, 1) = init2D[i].y;
  }
  const Eigen::MatrixXd anchorRhs = anchor * uv;

  std::vector<Eigen::Matrix2d> rotations(triangles.size(), Eigen::Matrix2d::Identity());
  for (int iteration = 0; iteration < iterations; ++iteration) {
    // Local step: closest rotation per triangle.
    for (size_t t = 0; t < triangles.size(); ++t) {
      const LocalTriangle &lt = local[t];
      if (!lt.valid) continue;
      Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();
      for (int e = 0; e < 3; ++e) {
        const int i = triangles[t][e], j = triangles[t][(e + 1) % 3];
        const Eigen::Vector2d du(uv(i, 0) - uv(j, 0), uv(i, 1) - uv(j, 1));
        covariance += lt.w[e] * du * (lt.x[e] - lt.x[(e + 1) % 3]).transpose();
      }
      Eigen::JacobiSVD<Eigen::Matrix2d> svd(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::Matrix2d u = svd.matrixU();
      Eigen::Matrix2d rotation = u * svd.matrixV().transpose();
      if (rotation.determinant() < 0.0) {
        u.col(1) *= -1.0;
        rotation = u * svd.matrixV().transpose();
      }
      rotations[t] = rotation;
    }

    // Global step: positions matching every rotated local copy.
    Eigen::MatrixXd rhs = anchorRhs;
    for (size_t t = 0; t < triangles.size(); ++t) {
      const LocalTriangle &lt = local[t];
      if (!lt.valid) continue;
      for (int e = 0; e < 3; ++e) {
        const int i = triangles[t][e], j = triangles[t][(e + 1) % 3];
        const Eigen::Vector2d d = lt.w[e] * (rotations[t] * (lt.x[e] - lt.x[(e + 1) % 3]));
        rhs(i, 0) += d.x();
        rhs(i, 1) += d.y();
        rhs(j, 0) -= d.x();
        rhs(j, 1) -= d.y();
      }
    }
    const Eigen::MatrixXd solved = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !solved.allFinite()) return false;
    uv = solved;
  }

  for (int i = 0; i < n; ++i) out2D[i] = zPoint(uv(i, 0), uv(i, 1), 0.0);

  // Stretch of every triangle: J maps its local isometric copy onto the result.
  double areaSum = 0.0, stretchSum = 0.0;
  for (size_t t = 0; t < triangles.size(); ++t) {
    const LocalTriangle &lt = local[t];
    if (!lt.valid) continue;
    const int a = triangles[t][0], b = triangles[t][1], c = triangles[t][2];
    Eigen::Matrix2d image, source;
    image << uv(b, 0) - uv(a, 0), uv(c, 0) - uv(a, 0), uv(b, 1) - uv(a, 1), uv(c, 1) - uv(a, 1);
    source.col(0) = lt.x[1] - lt.x[0];
    source.col(1) = lt.x[2] - lt.x[0];
    if (std::abs(source.determinant()) < 1e-12) continue;
    const Eigen::Matrix2d jacobian = image * source.inverse();
    if (jacobian.determinant() < 0.0) ++stats.flipped;
    const Eigen::Vector2d sigma = Eigen::JacobiSVD<Eigen::Matrix2d>(jacobian).singularValues();
    const double stretch = (std::max)(std::abs(sigma(0) - 1.0), std::abs(sigma(1) - 1.0));
    stats.stretchMax = (std::max)(stats.stretchMax, stretch);
    stretchSum += stretch * lt.area;
    areaSum += lt.area;
  }
  stats.stretchMean = (areaSum > 0.0) ? stretchSum / areaSum : 0.0;
  return true;
}

} // namespace

ZSPACE_INLINE bool arapUnroll(zObjMesh &mesh3D, const zPointArray &start2D, int iterations, zPointArray &out2D,
                               ArapStats &stats) {
  zFnMesh fn(mesh3D);
  zPointArray positions;
  fn.getVertexPositions(positions);
  zIntArray connects, counts;
  fn.getPolygonData(connects, counts);
  std::vector<std::array<int, 3>> triangles;
  size_t offset = 0;
  for (int count : counts) {
    for (int i = 1; i + 1 < count; ++i)
      triangles.push_back({connects[offset], connects[offset + i], connects[offset + i + 1]});
    offset += count;
  }
  return arapSolveTriangles(positions, triangles, start2D, iterations, out2D, stats);
}

ZSPACE_INLINE void planarProjection(zObjMesh &mesh3D, zPointArray &out2D) {
  zFnMesh fn(mesh3D);
  zPointArray positions;
  fn.getVertexPositions(positions);
  out2D.assign(positions.size(), zPoint(0, 0, 0));
  if (positions.empty()) return;

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const zPoint &p : positions) centroid += zTsUtils::toEigenVec(p);
  centroid /= static_cast<double>(positions.size());
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (const zPoint &p : positions) {
    const Eigen::Vector3d d = zTsUtils::toEigenVec(p) - centroid;
    covariance += d * d.transpose();
  }
  // Best-fit plane normal (smallest eigenvalue), turned to face +Z.
  Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
  if (solver.info() == Eigen::Success && solver.eigenvectors().col(0).allFinite() &&
      solver.eigenvectors().col(0).norm() > 1e-12)
    normal = solver.eigenvectors().col(0).normalized();
  if (normal.z() < 0.0) normal = -normal;

  // Smallest rotation taking the normal to +Z, about the centroid.
  const Eigen::Quaterniond rotation = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ());
  for (size_t i = 0; i < positions.size(); ++i) {
    const Eigen::Vector3d d =
        rotation * (zTsUtils::toEigenVec(positions[i]) - centroid);
    out2D[i] = zPoint(d.x() + centroid.x(), d.y() + centroid.y(), 0.0);
  }
}

ZSPACE_INLINE void alignRigid2D(zPointArray &points, const zPointArray &target) {
  const size_t n = (std::min)(points.size(), target.size());
  if (n == 0) return;
  double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
  for (size_t i = 0; i < n; ++i) {
    ax += points[i].x;
    ay += points[i].y;
    bx += target[i].x;
    by += target[i].y;
  }
  ax /= n;
  ay /= n;
  bx /= n;
  by /= n;
  double dotSum = 0.0, crossSum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double px = points[i].x - ax, py = points[i].y - ay;
    const double qx = target[i].x - bx, qy = target[i].y - by;
    dotSum += px * qx + py * qy;
    crossSum += px * qy - py * qx;
  }
  const double angle = std::atan2(crossSum, dotSum);
  const double c = std::cos(angle), s = std::sin(angle);
  for (zPoint &p : points) {
    const double px = p.x - ax, py = p.y - ay;
    p = zPoint(bx + c * px - s * py, by + s * px + c * py, 0.0);
  }
}

} // namespace zTsUnroll

} // namespace zSpace
