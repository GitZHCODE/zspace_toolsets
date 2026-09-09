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

#include <headers/zToolsets/slicer/zPrintSyn.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

namespace {

// Detects geometric corner feature points preserving sharp angles. Internal
// helper for resamplePrintPath - not part of zPrintSyn's public interface.
void _detectFeatureAngles(const zPointArray &positions,
                          const std::vector<int> &sequence,
                          float angleThresholdDeg,
                          std::vector<int> &outFeatureFlags) {
  const int N = static_cast<int>(sequence.size());
  outFeatureFlags.assign(N, 0);
  if (N < 3)
    return;

  const double thresholdRad =
      static_cast<double>(angleThresholdDeg) * (3.14159265358979323846 / 180.0);

  for (int i = 0; i < N; i++) {
    int prevId = sequence[(i - 1 + N) % N];
    int curId = sequence[i];
    int nextId = sequence[(i + 1) % N];

    zVector v0 = vecSub(positions[prevId], positions[curId]);
    zVector v1 = vecSub(positions[nextId], positions[curId]);
    const double len0 = v0.length();
    const double len1 = v1.length();
    if (len0 < 1e-6 || len1 < 1e-6)
      continue;

    double cosTheta = vecDot(v0, v1) / (len0 * len1);
    cosTheta = (std::max)(-1.0, (std::min)(1.0, cosTheta));
    const double angleRad = std::acos(cosTheta);

    if (angleRad < thresholdRad) {
      outFeatureFlags[i] = 1;
    }
  }
}

} // anonymous namespace

//----------------------------------------------------------------------
// Algorithm Implementations
//----------------------------------------------------------------------

ZSPACE_INLINE bool zPrintSyn::buildSequencePrintPath(PrintPath &ioPath,
                                                     double weldTolerance,
                                                     double spurLengthThreshold,
                                                     std::string *outReason) {
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  zIntArray edgeConnects;
  fnGraph.getVertexPositions(positions);
  fnGraph.getEdgeData(edgeConnects);

  const int numV = static_cast<int>(positions.size());
  const int numE = static_cast<int>(edgeConnects.size() / 2);
  if (numV < 3 || numE < 3) {
    if (outReason) {
      std::ostringstream oss;
      oss << "Insufficient elements (vertices=" << numV << ", edges=" << numE
          << " < 3)";
      *outReason = oss.str();
    }
    return false;
  }

  // Step 1: Build undirected adjacency list, discard self-loops and duplicate
  // edges
  std::vector<std::vector<int>> adj(numV);
  std::unordered_set<unsigned long long> edgeSet;
  for (size_t e = 0; e + 1 < edgeConnects.size(); e += 2) {
    int a = edgeConnects[e];
    int b = edgeConnects[e + 1];
    if (a < 0 || b < 0 || a >= numV || b >= numV || a == b)
      continue;
    if (vecDistSq(positions[b], positions[a]) < 1e-12)
      continue;

    const unsigned long long key = makeEdgeKey(a, b);
    if (!edgeSet.insert(key).second)
      continue;

    adj[a].push_back(b);
    adj[b].push_back(a);
  }

  // Step 2: Prune short dead-end branches / spurs
  std::vector<std::unordered_set<int>> activeAdj(numV);
  for (int i = 0; i < numV; i++) {
    for (int nb : adj[i]) {
      activeAdj[i].insert(nb);
    }
  }

  bool prunedAny = true;
  while (prunedAny) {
    prunedAny = false;
    for (int i = 0; i < numV; i++) {
      if (activeAdj[i].size() >= 3) {
        std::vector<int> neighbors(activeAdj[i].begin(), activeAdj[i].end());
        for (int startNb : neighbors) {
          std::vector<int> branchPath = {i, startNb};
          int prev = i;
          int curr = startNb;
          bool isDeadEnd = false;
          double branchLen = vecSub(positions[curr], positions[prev]).length();

          while (true) {
            if (activeAdj[curr].size() == 1) {
              isDeadEnd = true;
              break;
            }
            if (activeAdj[curr].size() != 2) {
              break;
            }
            int next = -1;
            for (int nb : activeAdj[curr]) {
              if (nb != prev) {
                next = nb;
                break;
              }
            }
            if (next < 0 || next == i)
              break;
            branchLen += vecSub(positions[next], positions[curr]).length();
            prev = curr;
            curr = next;
            branchPath.push_back(curr);
          }

          if (isDeadEnd && branchLen < spurLengthThreshold) {
            for (size_t k = 0; k + 1 < branchPath.size(); k++) {
              int u = branchPath[k];
              int v = branchPath[k + 1];
              activeAdj[u].erase(v);
              activeAdj[v].erase(u);
            }
            prunedAny = true;
          }
        }
      }
    }
  }

  // Step 3: Extract and order continuous chains
  std::vector<std::vector<int>> openChains;
  std::vector<std::vector<int>> closedLoops;
  std::unordered_set<unsigned long long> visitedEdges;

  // 3.1 Extract all open chains starting from endpoints (degree 1)
  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 1) {
      int next = *activeAdj[i].begin();
      if (visitedEdges.count(makeEdgeKey(i, next)))
        continue;

      std::vector<int> chain = {i};
      int prev = i;
      int curr = next;
      visitedEdges.insert(makeEdgeKey(prev, curr));

      while (true) {
        chain.push_back(curr);
        if (activeAdj[curr].size() != 2)
          break;

        int nxt = -1;
        for (int nb : activeAdj[curr]) {
          if (nb != prev) {
            nxt = nb;
            break;
          }
        }
        if (nxt < 0 || visitedEdges.count(makeEdgeKey(curr, nxt)))
          break;
        visitedEdges.insert(makeEdgeKey(curr, nxt));
        prev = curr;
        curr = nxt;
      }
      openChains.push_back(chain);
    }
  }

  // 3.2 Extract all closed loops with pure degree 2
  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 2) {
      for (int nb : activeAdj[i]) {
        if (visitedEdges.count(makeEdgeKey(i, nb)))
          continue;

        std::vector<int> loop = {i};
        int prev = i;
        int curr = nb;
        visitedEdges.insert(makeEdgeKey(prev, curr));
        bool loopClosed = false;

        while (true) {
          if (curr == i) {
            loopClosed = true;
            break;
          }
          loop.push_back(curr);
          if (activeAdj[curr].size() != 2)
            break;

          int nxt = -1;
          for (int n : activeAdj[curr]) {
            if (n != prev) {
              nxt = n;
              break;
            }
          }
          if (nxt < 0 || visitedEdges.count(makeEdgeKey(curr, nxt)))
            break;
          visitedEdges.insert(makeEdgeKey(curr, nxt));
          prev = curr;
          curr = nxt;
        }

        if (loopClosed && loop.size() >= 3) {
          closedLoops.push_back(loop);
        } else if (!loop.empty()) {
          openChains.push_back(loop);
        }
      }
    }
  }

  // Step 4: Endpoint gap welding and chain merging
  const double weldTolSq = weldTolerance * weldTolerance;

  // Check if any open chain is already self-closing
  for (auto it = openChains.begin(); it != openChains.end();) {
    if (it->size() >= 3) {
      if (vecDistSq(positions[it->front()], positions[it->back()]) <=
          weldTolSq) {
        closedLoops.push_back(*it);
        it = openChains.erase(it);
        continue;
      }
    }
    ++it;
  }

  // Iteratively weld chain endpoints
  bool mergedAny = true;
  while (mergedAny && openChains.size() > 1) {
    mergedAny = false;
    size_t bestI = 0, bestJ = 0;
    int bestMode = -1;
    double bestDistSq = weldTolSq;

    for (size_t i = 0; i < openChains.size(); i++) {
      for (size_t j = i + 1; j < openChains.size(); j++) {
        const auto &c1 = openChains[i];
        const auto &c2 = openChains[j];

        double d0 = vecDistSq(positions[c1.back()], positions[c2.front()]);
        if (d0 < bestDistSq) {
          bestDistSq = d0;
          bestI = i;
          bestJ = j;
          bestMode = 0;
        }

        double d1 = vecDistSq(positions[c1.back()], positions[c2.back()]);
        if (d1 < bestDistSq) {
          bestDistSq = d1;
          bestI = i;
          bestJ = j;
          bestMode = 1;
        }

        double d2 = vecDistSq(positions[c1.front()], positions[c2.back()]);
        if (d2 < bestDistSq) {
          bestDistSq = d2;
          bestI = i;
          bestJ = j;
          bestMode = 2;
        }

        double d3 = vecDistSq(positions[c1.front()], positions[c2.front()]);
        if (d3 < bestDistSq) {
          bestDistSq = d3;
          bestI = i;
          bestJ = j;
          bestMode = 3;
        }
      }
    }

    if (bestMode >= 0) {
      std::vector<int> &c1 = openChains[bestI];
      std::vector<int> c2 = openChains[bestJ];

      if (bestMode == 0) {
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 1) {
        std::reverse(c2.begin(), c2.end());
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 2) {
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      } else if (bestMode == 3) {
        std::reverse(c2.begin(), c2.end());
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      }

      openChains.erase(openChains.begin() + bestJ);
      mergedAny = true;

      if (c1.size() >= 3 &&
          vecDistSq(positions[c1.front()], positions[c1.back()]) <= weldTolSq) {
        closedLoops.push_back(c1);
        openChains.erase(openChains.begin() + bestI);
      }
    }
  }

  // Step 5: Select the longest continuous contour
  double maxLen = -1.0;
  std::vector<int> bestSequence;
  bool isBestClosed = false;

  auto computeLoopLength = [&](const std::vector<int> &seq,
                               bool isLoop) -> double {
    double total = 0.0;
    const int n = static_cast<int>(seq.size());
    for (int k = 0; k + 1 < n; k++) {
      total += vecSub(positions[seq[k + 1]], positions[seq[k]]).length();
    }
    if (isLoop && n > 2) {
      total += vecSub(positions[seq.front()], positions[seq.back()]).length();
    }
    return total;
  };

  for (const auto &loop : closedLoops) {
    double len = computeLoopLength(loop, true);
    if (len > maxLen) {
      maxLen = len;
      bestSequence = loop;
      isBestClosed = true;
    }
  }

  if (bestSequence.empty()) {
    for (const auto &chain : openChains) {
      double len = computeLoopLength(chain, false);
      if (len > maxLen) {
        maxLen = len;
        bestSequence = chain;
        isBestClosed = false;
      }
    }
  }

  if (bestSequence.size() < 3) {
    if (outReason) {
      std::ostringstream oss;
      oss << "No valid contour sequence extracted (openChains="
          << openChains.size() << ", closedLoops=" << closedLoops.size() << ")";
      *outReason = oss.str();
    }
    return false;
  }

  // Rewrite ioPath as the single clean, ordered sequence: reorder
  // widths/heights/normals to match, rebuild graph edges as a plain
  // chain (plus a wraparound edge if closed).
  const size_t M = bestSequence.size();
  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();

  zPointArray newPositions(M);
  zFloatArray newWidths(M, 0.04f);
  zFloatArray newHeights;
  zVectorArray newNormals;
  if (hasHeights)
    newHeights.assign(M, 0.0f);
  if (hasNormals)
    newNormals.assign(M, zVector(0, 0, 1));

  for (size_t i = 0; i < M; i++) {
    const int id = bestSequence[i];
    newPositions[i] = positions[id];
    if (id < static_cast<int>(ioPath.widths.size()))
      newWidths[i] = ioPath.widths[id];
    if (hasHeights && id < static_cast<int>(ioPath.heights.size()))
      newHeights[i] = ioPath.heights[id];
    if (hasNormals && id < static_cast<int>(ioPath.normals.size()))
      newNormals[i] = ioPath.normals[id];
  }

  zIntArray newEdges;
  newEdges.reserve(M * 2);
  for (size_t i = 0; i + 1 < M; i++) {
    newEdges.push_back(static_cast<int>(i));
    newEdges.push_back(static_cast<int>(i + 1));
  }
  if (isBestClosed && M >= 3) {
    newEdges.push_back(static_cast<int>(M - 1));
    newEdges.push_back(0);
  }

  fnGraph.clear();
  fnGraph.create(newPositions, newEdges);
  ioPath.widths = newWidths;
  if (hasHeights)
    ioPath.heights = newHeights;
  if (hasNormals)
    ioPath.normals = newNormals;

  return true;
}

ZSPACE_INLINE bool zPrintSyn::resamplePrintPath(PrintPath &ioPath) {
  zFnGraph fnGraph(ioPath.graph);
  const int N = fnGraph.numVertices();
  const float sampleLength =
      (i_Settings.sampleLength > 1e-4f) ? i_Settings.sampleLength : 0.005f;
  if (N < 2 || sampleLength <= 1e-6f)
    return false;

  const bool closed = isClosedGraph(ioPath.graph);

  zPointArray positions;
  fnGraph.getVertexPositions(positions);

  std::vector<int> sequence(N);
  for (int i = 0; i < N; i++)
    sequence[i] = i;

  std::vector<int> featureFlags;
  _detectFeatureAngles(positions, sequence, i_Settings.featureAngleThreshold,
                      featureFlags);

  const zFloatArray &weights = ioPath.widths;
  const zFloatArray &heights = ioPath.heights;
  const zVectorArray &normals = ioPath.normals;
  const bool hasHeights = !heights.empty();
  const bool hasNormals = !normals.empty();

  zPointArray outPoints;
  zFloatArray outWeights;
  zFloatArray outHeights;
  zVectorArray outNormals;
  zIntArray outEdges;

  auto appendPoint = [&](const zPoint &p, float w, float h,
                         const zVector &nrm) {
    if (!outPoints.empty()) {
      zVector d = vecSub(p, outPoints.back());
      if (d.length() < 1e-5)
        return;
    }
    outPoints.push_back(p);
    outWeights.push_back(w);
    if (hasHeights)
      outHeights.push_back(h);
    if (hasNormals)
      outNormals.push_back(nrm);
  };

  const int numSegments = closed ? N : (N - 1);
  for (int seg = 0; seg < numSegments; seg++) {
    int i0 = sequence[seg];
    int i1 = sequence[(seg + 1) % N];
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    float w0 = (i0 < static_cast<int>(weights.size())) ? weights[i0] : 0.04f;
    float w1 = (i1 < static_cast<int>(weights.size())) ? weights[i1] : 0.04f;
    float h0 = (hasHeights && i0 < static_cast<int>(heights.size()))
                   ? heights[i0]
                   : 0.0f;
    float h1 = (hasHeights && i1 < static_cast<int>(heights.size()))
                   ? heights[i1]
                   : 0.0f;
    zVector n0 = (hasNormals && i0 < static_cast<int>(normals.size()))
                     ? normals[i0]
                     : zVector(0, 0, 1);
    zVector n1 = (hasNormals && i1 < static_cast<int>(normals.size()))
                     ? normals[i1]
                     : zVector(0, 0, 1);

    zVector segVec = vecSub(p1, p0);
    const double segLen = segVec.length();

    appendPoint(p0, w0, h0, n0);

    if (segLen > sampleLength) {
      int steps = static_cast<int>(std::floor(segLen / sampleLength));
      for (int s = 1; s <= steps; s++) {
        double t = static_cast<double>(s * sampleLength) / segLen;
        if (t >= 1.0 - 1e-4)
          break;
        zPoint interpP(p0.x + t * segVec.x, p0.y + t * segVec.y,
                       p0.z + t * segVec.z);
        float interpW = static_cast<float>(w0 + t * (w1 - w0));
        float interpH = static_cast<float>(h0 + t * (h1 - h0));
        zVector interpN(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                        n0.z + t * (n1.z - n0.z));
        if (interpN.length() > 1e-6)
          interpN.normalize();
        else
          interpN = zVector(0, 0, 1);
        appendPoint(interpP, interpW, interpH, interpN);
      }
    }
  }

  if (!closed && N > 0) {
    int lastId = sequence[N - 1];
    float lastW =
        (lastId < static_cast<int>(weights.size())) ? weights[lastId] : 0.04f;
    float lastH = (hasHeights && lastId < static_cast<int>(heights.size()))
                      ? heights[lastId]
                      : 0.0f;
    zVector lastN = (hasNormals && lastId < static_cast<int>(normals.size()))
                        ? normals[lastId]
                        : zVector(0, 0, 1);
    appendPoint(positions[lastId], lastW, lastH, lastN);
  }

  const int numPts = static_cast<int>(outPoints.size());
  if (numPts >= 2) {
    for (int i = 0; i + 1 < numPts; i++) {
      outEdges.push_back(i);
      outEdges.push_back(i + 1);
    }
    if (closed && numPts >= 3) {
      zVector d = vecSub(outPoints.back(), outPoints.front());
      if (d.length() < 1e-4) {
        outPoints.pop_back();
        outWeights.pop_back();
        if (hasHeights)
          outHeights.pop_back();
        if (hasNormals)
          outNormals.pop_back();
        outEdges.pop_back();
        outEdges.pop_back();
      }
      const int finalN = static_cast<int>(outPoints.size());
      if (finalN >= 3) {
        outEdges.push_back(finalN - 1);
        outEdges.push_back(0);
      }
    }
  }

  if (outPoints.size() < 3 || outEdges.size() < 4)
    return false;

  fnGraph.clear();
  fnGraph.create(outPoints, outEdges);
  ioPath.widths = outWeights;
  if (hasHeights)
    ioPath.heights = outHeights;
  if (hasNormals)
    ioPath.normals = outNormals;
  return true;
}

ZSPACE_INLINE void zPrintSyn::filterGlitches(PrintPath &ioPath) {
  const bool closed = isClosedGraph(ioPath.graph);
  zFloatArray &widths = ioPath.widths;
  const int n = static_cast<int>(widths.size());
  if (n <= 2)
    return;
  for (int i = 0; i < n; i++) {
    int prevId = (i == 0) ? (closed ? n - 1 : 0) : i - 1;
    int nextId = (i == n - 1) ? (closed ? 0 : n - 1) : i + 1;
    if (prevId == i || nextId == i)
      continue;
    float prevWidth = widths[prevId];
    float currentWidth = widths[i];
    float nextWidth = widths[nextId];
    if (prevWidth <= 0.0f || nextWidth <= 0.0f)
      continue;
    if (std::fabs(prevWidth - nextWidth) > 0.001f)
      continue;
    if (std::fabs(currentWidth - prevWidth) <= 0.001f)
      continue;
    widths[i] = prevWidth;
  }
}

//----------------------------------------------------------------------
// High-Level Layer Processing Methods
//----------------------------------------------------------------------

ZSPACE_INLINE bool zPrintSyn::processPrintPath(std::vector<PrintPath> &ioPaths) {
  const int totalLayers = int(ioPaths.size());
  if (totalLayers == 0) {
    m_Status = "failed: No print paths available.";
    printError("[zPrintSyn] " + m_Status);
    return false;
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    if (buildSequencePrintPath(ioPaths[l]) && resamplePrintPath(ioPaths[l])) {
      filterGlitches(ioPaths[l]);
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Post-processed " << totalLayers << " print paths in " << std::fixed
      << std::setprecision(2) << wallTotalMs
      << " ms (sampleLen=" << i_Settings.sampleLength << "mm)";
  m_Status = out.str();
  printSuccess("[zPrintSyn] " + m_Status);
  return true;
}

} // namespace zSpace
