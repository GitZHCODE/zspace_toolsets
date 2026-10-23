#pragma once
#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
namespace zSpace {
namespace zTsIO {

// Defaults in millimetres for settings omitted from the USD stage.
inline zSpace::SlicerSettings defaultSlicerSettings()
{
    zSpace::SlicerSettings s;
    s.BoundaryWidth = {60.0f};
    s.BracingWidth = 40.0f;
    s.OverlapWidth = 0.0f;
    s.DoubleBoundary = true;
    s.SlotOffset = 1.0f;
    s.BracingSlot = 0.5f;
    s.FieldScaleUp = 1.1f;
    s.FieldCellSize = 5.0f;
    s.LayerSpacing = 8.0f;
    s._HTraceStepFraction = 0.2;
    s._HTraceMaxSteps = 100000;
    s._HEqualiseIterations = 2;
    s._HArapIterations = 10;
    s._HUnrollStart = 0;
    return s;
}

inline zSpace::PrintSynSetting defaultPrintSynSettings()
{
    zSpace::PrintSynSetting s;
    s.ResampleLength = 5.0f;
    s.FeatureAngleTol = 120.0f;
    s.SeamOffset = 40.0f;
    s.SeamSize = 20.0f;
    return s;
}


// Also replace /TDP_printmesh (Mesh) and /TDP_printpath (BasisCurves).
// Optional layerPaths excludes inter-layer travel from visualization.
// Copy an input USD and replace /TDP_gcode with a Scope containing ordered path
// data. Custom arrays printNormals/printWidths/printHeight correspond 1:1 to
// points. Path coordinates/lengths are millimetres; export converts to the
// source stage units.
ZSPACE_TOOLSETS bool exportUSD(const PrintPath &path, const std::string &inputFilePath,
                        const std::string &filePath, int layerCount,
                        std::string *outError = nullptr, int blockid = 0,
                        const PrintPathArray *layerPaths = nullptr);

// Settings and PrintSettings are resolved native settings.
// Returned coordinates and dimensional settings are millimetres after
// object transforms. MetersPerUnit records the source file unit; ratios,
// angles, tags and IDs are unchanged.
struct SlicerUSDInput {
  zObjMesh Block;
  zIntArray FaceTags;
  zIntArray FaceWidthIds;    // Per-face widthid into settings BoundaryWidth.
  zIntArray BracingWidthIds; // Per-curve widthid into settings BracingWidth.
  zObjGraph Bracing;
  size_t BracingCurveCount = 0;
  zFloatArray BracingPrintWidths;
  zPoint StartPoint;
  SlicerSettings Settings;
  PrintSynSetting PrintSettings;
  std::vector<std::string> UnknownSettings;
  double MetersPerUnit = 0.01;
  std::string UpAxis;
};
// Uses the generic Core IO USD document API.
// Paths must identify the geometry prim itself; missing/wrong schemas are
// errors. Failure clears the output and returns a diagnostic through outError.
ZSPACE_TOOLSETS bool importUSD(const std::string &filePath, SlicerUSDInput &out,
                        std::string *outError = nullptr,
                        const std::string &blockPath = "/AG_block",
                        const std::string &bracingPath = "/Slicer_bracing",
                        const std::string &startPath = "/Slicer_slot");

} // namespace zTsIO
} // namespace zSpace

