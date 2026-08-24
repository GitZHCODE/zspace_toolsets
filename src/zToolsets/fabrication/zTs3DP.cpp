#include <filesystem>

#include <zspace/zToolsets/fabrication/zTs3DP.h>

namespace zSpace
{
    namespace
    {
        std::string z3dpJoinIndexedPath(const std::string& directory, const std::string& stem, int index, const std::string& extension)
        {
            std::filesystem::path dir(directory);
            std::filesystem::create_directories(dir);
            return (dir / (stem + "_" + std::to_string(index) + extension)).string();
        }
    }

    zTs3DP::zTs3DP()
    {
        m_slicer.setSettings(m_settings);
        m_synthesis.setSettings(m_settings);
    }

    bool zTs3DP::fromMesh(zObjectMesh& mesh)
    {
        m_inputMesh = mesh;
        m_hasInput = true;
        m_slicer.setInputMesh(m_inputMesh);
        return true;
    }

    bool zTs3DP::readMesh(const std::string& path)
    {
        const zIOResult result = zIO::readMesh(path, m_inputMesh);
        m_hasInput = static_cast<bool>(result);
        if (m_hasInput) m_slicer.setInputMesh(m_inputMesh);
        return m_hasInput;
    }

    bool zTs3DP::writeInputMesh(const std::string& path)
    {
        return static_cast<bool>(zIO::writeMesh(path, m_inputMesh));
    }

    bool zTs3DP::writeSliceMeshes(const std::string& directory, const std::string& extension)
    {
        bool ok = true;
        auto& meshes = m_slicer.sliceMeshes();
        for (int i = 0; i < static_cast<int>(meshes.size()); ++i)
        {
            ok = static_cast<bool>(zIO::writeMesh(z3dpJoinIndexedPath(directory, "slice", i, extension), meshes[i])) && ok;
        }
        return ok;
    }

    bool zTs3DP::writePrintPaths(const std::string& directory, const std::string& extension)
    {
        bool ok = true;
        auto& paths = m_synthesis.printPaths();
        for (int i = 0; i < static_cast<int>(paths.size()); ++i)
        {
            ok = static_cast<bool>(zIO::writeGraph(z3dpJoinIndexedPath(directory, "print_path", i, extension), paths[i])) && ok;
        }
        return ok;
    }

    bool zTs3DP::writePrintMesh(const std::string& path)
    {
        return static_cast<bool>(zIO::writeMesh(path, m_printMesh));
    }

    void zTs3DP::setPrintLayerHeight(float height)
    {
        m_settings.printLayerHeight = height;
        m_slicer.setSettings(m_settings);
        m_synthesis.setSettings(m_settings);
    }

    void zTs3DP::setFieldResolution(int resX, int resY)
    {
        m_settings.fieldResolutionX = resX;
        m_settings.fieldResolutionY = resY;
        m_slicer.setSettings(m_settings);
    }

    void zTs3DP::setPrintWidth(float width)
    {
        m_settings.printWidth = width;
        m_slicer.setSettings(m_settings);
        m_synthesis.setSettings(m_settings);
    }

    void zTs3DP::setPrintSpacing(float spacing)
    {
        m_settings.printSpacing = spacing;
        m_synthesis.setSettings(m_settings);
    }

    void zTs3DP::setSDFThreshold(float threshold)
    {
        m_settings.sdfThreshold = threshold;
        m_slicer.setSettings(m_settings);
    }

    void zTs3DP::setSlicingSeedVertices(int fromVertexId, int toVertexId, int toHorizontalVertexId)
    {
        m_settings.slicingSeedFromVertexId = fromVertexId;
        m_settings.slicingSeedToVertexId = toVertexId;
        m_settings.slicingSeedHorizontalVertexId = toHorizontalVertexId;
        m_slicer.setSettings(m_settings);
    }

    void zTs3DP::computeSlicingFeatures()
    {
        if (!m_hasInput) return;
        m_slicer.identifySlicingFeatures();
    }

    void zTs3DP::computeSlices()
    {
        if (!m_hasInput) return;
        m_slicer.computeSliceMeshes();
        m_slicer.triangulateSlices();
    }

    void zTs3DP::computeUnrolledSlices()
    {
        m_slicer.unrollSlices();
    }

    void zTs3DP::computeUnrolledSDFs()
    {
        m_slicer.unrollSlices();
        m_slicer.computeSDFFields();
        m_slicer.extractSDFContours();
    }

    void zTs3DP::computePrintPaths()
    {
        m_synthesis.setInputs(m_slicer.contours(), m_slicer.unrolledMeshes(), m_slicer.sliceMeshes());
        m_synthesis.resampleContours();
        m_synthesis.identifyFeaturePoints();
        m_synthesis.mapContoursToSlices();
        m_synthesis.sequenceToolpaths();
    }

    void zTs3DP::computePrintMesh()
    {
        m_synthesis.buildPrintMesh();
        m_printMesh = m_synthesis.printMesh();
    }

    void zTs3DP::computeAll()
    {
        computeSlices();
        computeUnrolledSDFs();
        computePrintPaths();
        computePrintMesh();
    }

    void zTs3DP::clearResults()
    {
        m_slicer = z3DPSlicer();
        m_synthesis = z3DPPrintSynthesis();
        m_slicer.setSettings(m_settings);
        m_synthesis.setSettings(m_settings);
        if (m_hasInput) m_slicer.setInputMesh(m_inputMesh);
        m_printMesh = zObjectMesh();
    }

    zObjectMesh& zTs3DP::inputMesh() { return m_inputMesh; }
    zObjectMeshArray& zTs3DP::sliceMeshes() { return m_slicer.sliceMeshes(); }
    zObjectMeshArray& zTs3DP::unrolledSliceMeshes() { return m_slicer.unrolledMeshes(); }
    zObjectGraphArray& zTs3DP::sdfContours() { return m_slicer.contours(); }
    zObjectGraphArray& zTs3DP::edgeLoops() { return m_slicer.edgeLoops(); }
    zObjectPointCloud& zTs3DP::cornerPoints() { return m_slicer.cornerPoints(); }
    zObjectMesh& zTs3DP::topMesh() { return m_slicer.topMesh(); }
    zObjectMesh& zTs3DP::bottomMesh() { return m_slicer.bottomMesh(); }
    zObjectGraphArray& zTs3DP::printPaths() { return m_synthesis.printPaths(); }
    zObjectMesh& zTs3DP::printMesh() { return m_printMesh; }

    const z3DPSettings& zTs3DP::settings() const { return m_settings; }
    const z3DPSlicingFeatures& zTs3DP::slicingFeatures() const { return m_slicer.features(); }
    z3DPSlicer& zTs3DP::slicer() { return m_slicer; }
    z3DPPrintSynthesis& zTs3DP::synthesis() { return m_synthesis; }
}
