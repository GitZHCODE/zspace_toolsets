#ifndef ZSPACE_Z3DP_TYPES_H
#define ZSPACE_Z3DP_TYPES_H

#pragma once

#include <zspace/interface.h>

namespace zSpace
{
    struct z3DPSettings
    {
        float printLayerHeight = 0.01f;
        float printWidth = 0.004f;
        float printSpacing = 0.005f;
        float sdfThreshold = 0.0f;
        int fieldResolutionX = 200;
        int fieldResolutionY = 80;
        int slicingSeedFromVertexId = -1;
        int slicingSeedToVertexId = -1;
        int slicingSeedHorizontalVertexId = -1;
    };

    struct z3DPSliceFrame
    {
        int index = 0;
        zObjectMesh sliceMesh;
        zObjectMesh triangulatedMesh;
        zObjectMesh unrolledMesh;
        zObjectMeshScalarField sdfField;
        zObjectGraph rawContour;
        zObjectGraph printPath;
    };

    struct z3DPSlicingFeatures
    {
        zIntArray cornerVertexIds;
        zIntArray visitedEdgeIds;
        zIntArray outgoingCornerEdgeIds;
        zIntArray cornerToCornerEdgeIds;
        zIntArray blueSeedEdgeIds;
        zIntArray bottomStripFaceIds;
        zIntArray topStripFaceIds;
        int edgeLoopCount = 0;
        zObjectPointCloud cornerPoints;
        zObjectGraphArray edgeLoops;
        zObjectMesh topMesh;
        zObjectMesh bottomMesh;
    };
}

#endif
