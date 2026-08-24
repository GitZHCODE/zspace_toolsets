#ifndef ZSPACE_ZTS_3DP_H
#define ZSPACE_ZTS_3DP_H

#pragma once

#include <string>

#include <zspace/io.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/z3DPSlicer.h>
#include <zspace/zToolsets/fabrication/z3DPPrintSynthesis.h>

namespace zSpace
{
    class ZSPACE_TOOLSETS zTs3DP
    {
    public:
        zTs3DP();

        bool fromMesh(zObjectMesh& mesh);
        bool readMesh(const std::string& path);
        bool writeInputMesh(const std::string& path);
        bool writeSliceMeshes(const std::string& directory, const std::string& extension = ".obj");
        bool writePrintPaths(const std::string& directory, const std::string& extension = ".json");
        bool writePrintMesh(const std::string& path);

        void setPrintLayerHeight(float height);
        void setFieldResolution(int resX, int resY);
        void setPrintWidth(float width);
        void setPrintSpacing(float spacing);
        void setSDFThreshold(float threshold);
        void setSlicingSeedVertices(int fromVertexId, int toVertexId, int toHorizontalVertexId);

        void computeSlicingFeatures();
        void computeSlices();
        void computeUnrolledSlices();
        void computeUnrolledSDFs();
        void computePrintPaths();
        void computePrintMesh();
        void computeAll();
        void clearResults();

        zObjectMesh& inputMesh();
        zObjectMeshArray& sliceMeshes();
        zObjectMeshArray& unrolledSliceMeshes();
        zObjectGraphArray& sdfContours();
        zObjectGraphArray& edgeLoops();
        zObjectPointCloud& cornerPoints();
        zObjectMesh& topMesh();
        zObjectMesh& bottomMesh();
        zObjectGraphArray& printPaths();
        zObjectMesh& printMesh();

        const z3DPSettings& settings() const;
        const z3DPSlicingFeatures& slicingFeatures() const;
        z3DPSlicer& slicer();
        z3DPPrintSynthesis& synthesis();

    private:
        zObjectMesh m_inputMesh;
        zObjectMesh m_printMesh;
        z3DPSettings m_settings;
        z3DPSlicer m_slicer;
        z3DPPrintSynthesis m_synthesis;
        bool m_hasInput = false;
    };
}

#endif
