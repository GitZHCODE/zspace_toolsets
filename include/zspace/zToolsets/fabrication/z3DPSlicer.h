#ifndef ZSPACE_Z3DP_SLICER_H
#define ZSPACE_Z3DP_SLICER_H

#pragma once

#include <vector>

#include <zspace/io.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/z3DPTypes.h>

namespace zSpace
{
    class ZSPACE_TOOLSETS z3DPSlicer
    {
    public:
        void setInputMesh(zObjectMesh& mesh);
        void setSettings(const z3DPSettings& settings);

        void identifySlicingFeatures();
        void computeSliceCount();
        void computeSliceMeshes();
        void triangulateSlices();
        void unrollSlices();
        void computeSDFFields();
        void extractSDFContours();

        zObjectMeshArray& sliceMeshes();
        zObjectMeshArray& triangulatedMeshes();
        zObjectMeshArray& unrolledMeshes();
        zObjectMeshScalarFieldArray& fields();
        zObjectGraphArray& contours();
        zObjectGraphArray& edgeLoops();
        zObjectPointCloud& cornerPoints();
        zObjectMesh& topMesh();
        zObjectMesh& bottomMesh();

        const zObjectMeshArray& sliceMeshes() const;
        const zObjectMeshArray& unrolledMeshes() const;
        const zObjectGraphArray& contours() const;
        const z3DPSlicingFeatures& features() const;

    private:
        zObjectMesh* m_inputMesh = nullptr;
        z3DPSettings m_settings;
        int m_sliceCount = 0;

        zObjectMeshArray m_sliceMeshes;
        zObjectMeshArray m_triangulatedMeshes;
        zObjectMeshArray m_unrolledMeshes;
        zObjectMeshScalarFieldArray m_fields;
        zObjectGraphArray m_contours;
        z3DPSlicingFeatures m_features;
        std::vector<zItMeshHalfEdgeArray> m_verticalLoops;
        zItMeshHalfEdgeArray m_bottomLoop;

        zPoint m_minBounds;
        zPoint m_maxBounds;
    };
}

#endif
