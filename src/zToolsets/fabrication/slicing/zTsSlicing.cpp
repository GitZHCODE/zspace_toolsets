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

#include <zspace/zToolsets/fabrication/slicing/zTsSlicing.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <memory>
#include <stdexcept>
#include <string>

namespace zSpace {

//--------------------------
//---- MOVE
//--------------------------

ZSPACE_INLINE zTsSlicing::zTsSlicing(zTsSlicing &&other) : _Status("Initialized") {
  *this = std::move(other);
}

ZSPACE_INLINE zTsSlicing &zTsSlicing::operator=(zTsSlicing &&other) {
  if (this == &other)
    return *this;
  i_BlockMesh = std::move(other.i_BlockMesh);
  i_BracingGraph = std::move(other.i_BracingGraph);
  i_FaceTags = std::move(other.i_FaceTags);
  i_FaceWidthIds = std::move(other.i_FaceWidthIds);
  i_SlotPoint = other.i_SlotPoint;
  _HasSlotPoint = other._HasSlotPoint;
  _Mode = other._Mode;
  _Settings = std::move(other._Settings);
  _Slicer = std::move(other._Slicer);
  _PrintSynSettings = std::move(other._PrintSynSettings);
  _Status = std::move(other._Status);
  // The operator still points at other's mesh/graph - bind it to ours, and
  // give other a fresh one so it stays usable.
  _Slicer->setInputMesh(i_BlockMesh, i_FaceTags, i_FaceWidthIds);
  _Slicer->setBracingGraph(i_BracingGraph);
  if (_HasSlotPoint)
    _Slicer->setSlotAnchor(i_SlotPoint);
  else
    _Slicer->clearSlotAnchor();
  // v3 geometry objects own pimpl storage; moved-from objects need fresh
  // storage before their operators can be reused.
  other.i_BlockMesh = zObjMesh();
  other.i_BracingGraph = zObjGraph();
  other.i_FaceTags.clear();
  other.i_FaceWidthIds.clear();
  other._Settings = SlicerSettings();
  other._PrintSynSettings = PrintSynSetting();
  other._Status = "Initialized";
  other._HasSlotPoint = false;
  other._createSlicer();
  return *this;
}

//--------------------------
//---- SET METHODS
//--------------------------

ZSPACE_INLINE void zTsSlicing::setInputMesh(zObjMesh &mesh,
                                        const zIntArray &faceTags,
                                        const zIntArray &faceWidthIds) {
  i_BlockMesh = mesh;
  i_FaceTags = faceTags;
  i_FaceWidthIds = faceWidthIds;
  _Slicer->setInputMesh(i_BlockMesh, i_FaceTags, i_FaceWidthIds);
  _Status = "Input mesh updated";
}

ZSPACE_INLINE void zTsSlicing::setBracingGraph(zObjGraph &graph) {
  i_BracingGraph = graph;
  _Slicer->setBracingGraph(i_BracingGraph);
  _Status = "Bracing graph updated";
}

ZSPACE_INLINE void zTsSlicing::setSlotAnchor(const zPoint &anchor) {
  i_SlotPoint = anchor;
  _HasSlotPoint = true;
  _Slicer->setSlotAnchor(i_SlotPoint);
  _Status = "Slot anchor updated";
}

ZSPACE_INLINE void zTsSlicing::clearSlotAnchor() {
  i_SlotPoint = zPoint(0, 0, 0);
  _HasSlotPoint = false;
  _Slicer->clearSlotAnchor();
  _Status = "Slot anchor cleared";
}

ZSPACE_INLINE bool zTsSlicing::supportsSlicerMode(zTsSlicerMode mode) {
  if (mode == zTsSlicerMode::Topology) return true;
#if defined(ZSPACE_TOOLSETS_WITH_TETGEN)
  return mode == zTsSlicerMode::Harmonic;
#else
  return false;
#endif
}

ZSPACE_INLINE void zTsSlicing::setSlicerMode(zTsSlicerMode mode) {
  if (mode != zTsSlicerMode::Topology && mode != zTsSlicerMode::Harmonic)
    throw std::invalid_argument("Unknown slicing mode.");
  if (!supportsSlicerMode(mode))
    throw std::invalid_argument("Requested slicing mode is unavailable; harmonic slicing requires ZSPACE_TOOLSETS_WITH_TETGEN=ON.");
  _Mode = mode;
  _createSlicer();
  _Status = (mode == zTsSlicerMode::Topology) ? "Slicer mode: Topology"
                                         : "Slicer mode: Harmonic";
}

//--------------------------
//---- OPERATOR LIFETIME
//--------------------------

ZSPACE_INLINE void zTsSlicing::_createSlicer() {
  _Slicer.reset(); // free the old operator before building the new one
#if defined(ZSPACE_TOOLSETS_WITH_TETGEN)
  if (_Mode == zTsSlicerMode::Harmonic)
    _Slicer = std::make_unique<HarmonicSlicer>();
  else
#endif
    _Slicer = std::make_unique<TopologySlicer>();
  _Slicer->setInputMesh(i_BlockMesh, i_FaceTags, i_FaceWidthIds);
  _Slicer->setBracingGraph(i_BracingGraph);
  if (_HasSlotPoint)
    _Slicer->setSlotAnchor(i_SlotPoint);
}

//--------------------------
//---- PIPELINE EXECUTION METHODS
//--------------------------

ZSPACE_INLINE bool zTsSlicing::computeBasicStructure(bool useParallel) {
  _Slicer->setSettings(_Settings);
  bool success = _Slicer->computeBasicStructure(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTsSlicing::computeSlicer(bool useParallel) {
  _Slicer->setSettings(_Settings);
  bool success = _Slicer->computeSlicer(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTsSlicing::computePrintPath(bool useParallel) {
  const bool success = _Slicer->computePrintPath(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTsSlicing::computePostProcess(bool useParallel) {
  bool success = _Slicer->computePostProcess(_PrintSynSettings, useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

//--------------------------
//---- GET METHODS : operator-specific
//--------------------------

ZSPACE_INLINE TopologySlicer *zTsSlicing::getTopologySlicer() {
  return (_Mode == zTsSlicerMode::Topology)
             ? static_cast<TopologySlicer *>(_Slicer.get())
             : nullptr;
}

ZSPACE_INLINE const TopologySlicer *zTsSlicing::getTopologySlicer() const {
  return const_cast<zTsSlicing *>(this)->getTopologySlicer();
}

ZSPACE_INLINE HarmonicSlicer *zTsSlicing::getHarmonicSlicer() {
#if defined(ZSPACE_TOOLSETS_WITH_TETGEN)
  return (_Mode == zTsSlicerMode::Harmonic)
             ? static_cast<HarmonicSlicer *>(_Slicer.get())
             : nullptr;
#else
  return nullptr;
#endif
}

ZSPACE_INLINE const HarmonicSlicer *zTsSlicing::getHarmonicSlicer() const {
  return const_cast<zTsSlicing *>(this)->getHarmonicSlicer();
}

ZSPACE_INLINE int zTsSlicing::getWidthGroupCount() const {
  return zTsUtils::widthGroupCount(i_FaceTags, i_FaceWidthIds);
}

// Without a Topology operator the getters below return a per-call-site empty
// object, so callers can always iterate the result.

ZSPACE_INLINE zObjMeshArray &zTsSlicing::getIsoMeshes() {
  static zObjMeshArray empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getIsoMeshes() : (empty.clear(), empty);
}

ZSPACE_INLINE zObjMeshArray &zTsSlicing::getIsoMeshes2D() {
  static zObjMeshArray empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getIsoMeshes2D() : (empty.clear(), empty);
}

ZSPACE_INLINE zObjMeshScalarField &zTsSlicing::getScalarField() {
  static zObjMeshScalarField empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getScalarField() : empty;
}

ZSPACE_INLINE std::vector<zScalarArray> &zTsSlicing::getPathFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getPathFields() : (empty.clear(), empty);
}

ZSPACE_INLINE std::vector<zScalarArray> &zTsSlicing::getWidthFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getWidthFields() : (empty.clear(), empty);
}

ZSPACE_INLINE std::vector<zScalarArray> &zTsSlicing::getTrackFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getTrackFields() : (empty.clear(), empty);
}

} // namespace zSpace
