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

#include <headers/zToolsets/ts3dp/zTs3DP.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <memory>
#include <string>

namespace zSpace {

//--------------------------
//---- MOVE
//--------------------------

ZSPACE_INLINE zTs3DP::zTs3DP(zTs3DP &&other) : _Status("Initialized") {
  *this = std::move(other);
}

ZSPACE_INLINE zTs3DP &zTs3DP::operator=(zTs3DP &&other) {
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
  other._HasSlotPoint = false;
  other._createSlicer();
  return *this;
}

//--------------------------
//---- SET METHODS
//--------------------------

ZSPACE_INLINE void zTs3DP::setInputMesh(zObjMesh &mesh,
                                        const zIntArray &faceTags,
                                        const zIntArray &faceWidthIds) {
  i_BlockMesh = mesh;
  i_FaceTags = faceTags;
  i_FaceWidthIds = faceWidthIds;
  _Slicer->setInputMesh(i_BlockMesh, i_FaceTags, i_FaceWidthIds);
  _Status = "Input mesh updated";
}

ZSPACE_INLINE void zTs3DP::setBracingGraph(zObjGraph &graph) {
  i_BracingGraph = graph;
  _Slicer->setBracingGraph(i_BracingGraph);
  _Status = "Bracing graph updated";
}

ZSPACE_INLINE void zTs3DP::setSlotAnchor(const zPoint &anchor) {
  i_SlotPoint = anchor;
  _HasSlotPoint = true;
  _Slicer->setSlotAnchor(i_SlotPoint);
  _Status = "Slot anchor updated";
}

ZSPACE_INLINE void zTs3DP::clearSlotAnchor() {
  i_SlotPoint = zPoint(0, 0, 0);
  _HasSlotPoint = false;
  _Slicer->clearSlotAnchor();
  _Status = "Slot anchor cleared";
}

ZSPACE_INLINE void zTs3DP::setSlicerMode(zTsSlicerMode mode) {
  _Mode = mode;
  _createSlicer();
  _Status = (mode == zTsSlicerMode::Topology) ? "Slicer mode: Topology"
                                         : "Slicer mode: Harmonic";
}

//--------------------------
//---- OPERATOR LIFETIME
//--------------------------

ZSPACE_INLINE void zTs3DP::_createSlicer() {
  _Slicer.reset(); // free the old operator before building the new one
  if (_Mode == zTsSlicerMode::Harmonic)
    _Slicer = std::make_unique<HarmonicSlicer>();
  else
    _Slicer = std::make_unique<TopologySlicer>();
  _Slicer->setInputMesh(i_BlockMesh, i_FaceTags, i_FaceWidthIds);
  _Slicer->setBracingGraph(i_BracingGraph);
  if (_HasSlotPoint)
    _Slicer->setSlotAnchor(i_SlotPoint);
}

//--------------------------
//---- PIPELINE EXECUTION METHODS
//--------------------------

ZSPACE_INLINE bool zTs3DP::computeBasicStructure(bool useParallel) {
  _Slicer->setSettings(_Settings);
  bool success = _Slicer->computeBasicStructure(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::computeSlicer(bool useParallel) {
  _Slicer->setSettings(_Settings);
  bool success = _Slicer->computeSlicer(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::computePrintPath(bool useParallel) {
  const bool success = _Slicer->computePrintPath(useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::computePostProcess(bool useParallel) {
  bool success = _Slicer->computePostProcess(_PrintSynSettings, useParallel);
  _Status = _Slicer->getStatus();
  return success;
}

//--------------------------
//---- GET METHODS : operator-specific
//--------------------------

ZSPACE_INLINE TopologySlicer *zTs3DP::getTopologySlicer() {
  return (_Mode == zTsSlicerMode::Topology)
             ? static_cast<TopologySlicer *>(_Slicer.get())
             : nullptr;
}

ZSPACE_INLINE const TopologySlicer *zTs3DP::getTopologySlicer() const {
  return const_cast<zTs3DP *>(this)->getTopologySlicer();
}

ZSPACE_INLINE HarmonicSlicer *zTs3DP::getHarmonicSlicer() {
  return (_Mode == zTsSlicerMode::Harmonic)
             ? static_cast<HarmonicSlicer *>(_Slicer.get())
             : nullptr;
}

ZSPACE_INLINE const HarmonicSlicer *zTs3DP::getHarmonicSlicer() const {
  return const_cast<zTs3DP *>(this)->getHarmonicSlicer();
}

ZSPACE_INLINE int zTs3DP::getWidthGroupCount() const {
  return zTsUtils::widthGroupCount(i_FaceTags, i_FaceWidthIds);
}

// Without a Topology operator the getters below return a per-call-site empty
// object, so callers can always iterate the result.

ZSPACE_INLINE zObjMeshArray &zTs3DP::getIsoMeshes() {
  static zObjMeshArray empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getIsoMeshes() : (empty.clear(), empty);
}

ZSPACE_INLINE zObjMeshArray &zTs3DP::getIsoMeshes2D() {
  static zObjMeshArray empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getIsoMeshes2D() : (empty.clear(), empty);
}

ZSPACE_INLINE zObjMeshScalarField &zTs3DP::getScalarField() {
  static zObjMeshScalarField empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getScalarField() : empty;
}

ZSPACE_INLINE std::vector<zScalarArray> &zTs3DP::getPathFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getPathFields() : (empty.clear(), empty);
}

ZSPACE_INLINE std::vector<zScalarArray> &zTs3DP::getWidthFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getWidthFields() : (empty.clear(), empty);
}

ZSPACE_INLINE std::vector<zScalarArray> &zTs3DP::getTrackFields() {
  static std::vector<zScalarArray> empty;
  TopologySlicer *topology = getTopologySlicer();
  return topology ? topology->getTrackFields() : (empty.clear(), empty);
}

} // namespace zSpace
