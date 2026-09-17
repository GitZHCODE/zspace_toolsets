#ifndef ZSPACE_TS_GEOMETRY_PRINTSYN_H
#define ZSPACE_TS_GEOMETRY_PRINTSYN_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zToolsets/slicer/zTsUtils.h>
#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Toolsets for geometry, fabrication and computation.
 *  @{
 */

/*! \struct PrintSynSetting
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 */
struct PrintSynSetting {
  float ResampleLength = 0.1f;
  float FeatureAngleTol = 90.0f;
  float SeamOffset = 0.0f;
  float SeamSize = 0.0f;
};

/*! \class zPrintSyn
 *	\brief Toolpath cleanup, corner-preserving resampling, and
 * width-filtering engine.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zPrintSyn {
protected:
  PrintSynSetting i_Settings;
  std::string m_Status;

public:
  zPrintSyn() { m_Status = "Initialized"; }
  ~zPrintSyn() {}

  void setSettings(const PrintSynSetting &inSettings) {
    i_Settings = inSettings;
  }
  const PrintSynSetting &getSettings() const { return i_Settings; }
  PrintSynSetting &getSettings() { return i_Settings; }

  const std::string &getStatus() const { return m_Status; }

  //--------------------------
  //---- PIPELINE STEPS (public) - everything else is a private helper
  //--------------------------

  /*! \brief Step 1: extracts each layer's isocontour (marching triangles on
   * scalarField, tagged by pathFields/widthFields/trackFields) into a single
   * oriented, inner-before-outer, CCW closed-loop PrintPath per layer.
   * Parallel across layers internally. */
  bool extractPrintPaths(zObjMeshScalarField &scalarField,
                        const std::vector<zScalarArray> &pathFields,
                        const std::vector<zScalarArray> &widthFields,
                        const std::vector<zScalarArray> &trackFields,
                        std::vector<std::vector<PrintPath>> &outPrintPaths2D);

  /*! \brief Step 2: arc-length resamples every layer's ring (carrying
   * widths/heights/normals/trackIds along), smooths localized width
   * glitches, then forces the ring CCW. Parallel across layers internally. */
  bool resamplePrintPath(std::vector<std::vector<PrintPath>> &ioPaths);

  /*! \brief Step 3: aligns every layer's ring start point to the same
   * barycentric point on the underlying mesh (tracks the model through
   * taper/twist), constrained to the inner track, then cuts each ring open -
   * walking backward SeamSize from the new start to find the new end point,
   * leaving that gap open between them. Seeds from the first usable layer by
   * walking SeamOffset from its extracted start along the inner segment;
   * every other layer's start is just the nearest-inner-point match to that
   * seed, no further offset. Parallel across layers internally (after a
   * short sequential seed step). */
  bool alignStartPoints(std::vector<std::vector<PrintPath>> &ioPaths,
                       const zObjMeshArray &sectionMeshes2D);

  /*! \brief Step 4: projects 2D PrintPaths back onto the 3D section meshes.
   * vScalar/layerSpacing come from zSlicer (getVScalar()/
   * getSettings().LayerSpacing) - the only zSlicer-owned inputs this step
   * needs. Parallel across layers internally. */
  bool projectPrintPath(const std::vector<std::vector<PrintPath>> &in2DPaths,
                       const zObjMeshArray &sectionMeshes2D,
                       const zObjMeshArray &sectionMeshes,
                       const std::vector<zScalarArray> &vScalar,
                       float layerSpacing,
                       std::vector<std::vector<PrintPath>> &out3DPaths);

  /*! \brief Step 5: concatenates every layer's single (open, post-
   * alignStartPoints) ring PrintPath end-to-start into one continuous
   * bottom-to-top PrintPath, connecting each layer's last point to the next
   * layer's first point. Sequential (order-dependent), no internal
   * parallelism. */
  bool linkPrintPath(const std::vector<std::vector<PrintPath>> &layeredPaths,
                    PrintPath &outCombined);
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/slicer/zPrintSyn.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_PRINTSYN_H
