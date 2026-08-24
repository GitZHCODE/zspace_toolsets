#ifndef ZSPACE_Z3DP_PRINT_SYNTHESIS_H
#define ZSPACE_Z3DP_PRINT_SYNTHESIS_H

#pragma once

#include <zspace/io.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/z3DPTypes.h>

namespace zSpace
{
    class ZSPACE_TOOLSETS z3DPPrintSynthesis
    {
    public:
        void setSettings(const z3DPSettings& settings);
        void setInputs(zObjectGraphArray& contours, zObjectMeshArray& unrolledMeshes, zObjectMeshArray& sliceMeshes);

        void resampleContours();
        void identifyFeaturePoints();
        void mapContoursToSlices();
        void sequenceToolpaths();
        void buildPrintMesh();

        zObjectGraphArray& printPaths();
        zObjectPointCloudArray& featurePoints();
        zObjectMesh& printMesh();

        const zObjectGraphArray& printPaths() const;
        const zObjectMesh& printMesh() const;

    private:
        z3DPSettings m_settings;
        zObjectGraphArray* m_contours = nullptr;
        zObjectMeshArray* m_unrolledMeshes = nullptr;
        zObjectMeshArray* m_sliceMeshes = nullptr;

        zObjectGraphArray m_resampledContours;
        zObjectGraphArray m_printPaths;
        zObjectPointCloudArray m_featurePoints;
        zObjectMesh m_printMesh;
    };
}

#endif
