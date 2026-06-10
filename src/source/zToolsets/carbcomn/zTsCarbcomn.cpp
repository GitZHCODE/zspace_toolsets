// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms,
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
//


#include "zToolsets/carbcomn/zTsCarbcomn.h"
#include "zInterface/functionsets/zFnPointCloud.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>
//#include "zCore/base/zColor.h"
//#include "zCore/base/zEnumerators.h"
//#include "zCore/base/zTypeDef.h"
//#include "zInterface/functionsets/zFnMeshField.h"
//#include "zInterface/iterators/zItGraph.h"
//
//#include <stdio.h>

#if defined ZSPACE_USD_INTEROP
#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/vec3f.h>
#include <pxr/base/plug/registry.h>
#include <pxr/base/vt/array.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/sdf/valueTypeName.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/xform.h>

using namespace pxr;
#endif

namespace zSpace
{
#if defined ZSPACE_USD_INTEROP
	void zFnMesh::from(pxr::UsdPrim& usd, bool staticGeom) {}
	void zFnMesh::to(pxr::UsdPrim& usd) {}
	void zFnGraph::from(pxr::UsdPrim& usd, bool staticGeom) {}
	void zFnGraph::to(pxr::UsdPrim& usd) {}
	void zFnParticle::from(pxr::UsdPrim& usd, bool staticGeom) {}
	void zFnParticle::to(pxr::UsdPrim& usd) {}
	void zFnPointCloud::from(pxr::UsdPrim& usd, bool staticGeom) {}
	void zFnPointCloud::to(pxr::UsdPrim& usd) {}

	namespace
	{
		bool prepareUsdTextFile(const std::string& path, std::ofstream& out)
		{
			std::filesystem::path outPath(path);
			std::error_code ec;
			if (outPath.has_parent_path()) std::filesystem::create_directories(outPath.parent_path(), ec);
			if (std::filesystem::exists(outPath, ec)) std::filesystem::remove(outPath, ec);

			printf("\n writing USDA file: %s", outPath.string().c_str());
			out.open(outPath, std::ios::out | std::ios::trunc);
			if (!out.good())
			{
				printf("\n error opening USDA file: %s", outPath.string().c_str());
				return false;
			}

			out << std::setprecision(9);
			return true;
		}

		std::string tokenFromName(const std::string& name)
		{
			std::string clean = name;
			for (char& c : clean)
			{
				const bool valid = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
				if (!valid) c = '_';
			}
			if (clean.empty() || (clean[0] >= '0' && clean[0] <= '9')) clean = "_" + clean;
			return clean;
		}

		void writeUsdHeader(std::ofstream& out)
		{
			out << "#usda 1.0\n";
			out << "(\n";
			out << "    defaultPrim = \"World\"\n";
			out << "    upAxis = \"Z\"\n";
			out << "    metersPerUnit = 1\n";
			out << ")\n\n";
		}

		void writeMatrix(std::ofstream& out, const zTransform& transform)
		{
			out << "((";
			for (int i = 0; i < 4; i++)
			{
				if (i > 0) out << ", (";
				for (int j = 0; j < 4; j++)
				{
					if (j > 0) out << ", ";
					out << transform(i, j);
				}
				out << ")";
			}
			out << ")";
		}

		void writePointArray(std::ofstream& out, const zPointArray& positions)
		{
			out << "[";
			for (int i = 0; i < positions.size(); i++)
			{
				if (i > 0) out << ", ";
				out << "(" << positions[i].x << ", " << positions[i].y << ", " << positions[i].z << ")";
			}
			out << "]";
		}

		void writeIntArray(std::ofstream& out, const zIntArray& values)
		{
			out << "[";
			for (int i = 0; i < values.size(); i++)
			{
				if (i > 0) out << ", ";
				out << values[i];
			}
			out << "]";
		}

		void writeIndent(std::ofstream& out, int indent)
		{
			for (int i = 0; i < indent; i++) out << " ";
		}

		void writePrimXformAttrs(std::ofstream& out, const zTransform* frame, int indent)
		{
			if (frame)
			{
				writeIndent(out, indent);
				out << "matrix4d xformOp:transform = ";
				writeMatrix(out, *frame);
				out << "\n";
				writeIndent(out, indent);
				out << "uniform token[] xformOpOrder = [\"xformOp:transform\"]\n";
				writeIndent(out, indent);
				out << "custom matrix4d Frame = ";
				writeMatrix(out, *frame);
				out << "\n";
			}
		}

		void writeMeshPrim(std::ofstream& out, zObjMesh& meshObj, const std::string& primName, const zTransform* frame, int indent)
		{
			zFnMesh fnMesh(meshObj);
			zPointArray positions;
			zIntArray polyConnects;
			zIntArray polyCounts;
			fnMesh.getVertexPositions(positions);
			fnMesh.getPolygonData(polyConnects, polyCounts);

			writeIndent(out, indent);
			out << "def Xform \"" << tokenFromName(primName) << "\"\n";
			writeIndent(out, indent);
			out << "{\n";
			writePrimXformAttrs(out, frame, indent + 4);
			writeIndent(out, indent + 4);
			out << "def Mesh \"Mesh\"\n";
			writeIndent(out, indent + 4);
			out << "{\n";
			writeIndent(out, indent + 8);
			out << "point3f[] points = ";
			writePointArray(out, positions);
			out << "\n";
			writeIndent(out, indent + 8);
			out << "int[] faceVertexCounts = ";
			writeIntArray(out, polyCounts);
			out << "\n";
			writeIndent(out, indent + 8);
			out << "int[] faceVertexIndices = ";
			writeIntArray(out, polyConnects);
			out << "\n";
			writeIndent(out, indent + 8);
			out << "uniform token subdivisionScheme = \"none\"\n";
			writeIndent(out, indent + 4);
			out << "}\n";
			writeIndent(out, indent);
			out << "}\n";
		}

		void writeGraphPrim(std::ofstream& out, zObjGraph& graphObj, const std::string& primName, const zTransform* frame, const zIntArray* vertexSequence, int indent)
		{
			zFnGraph fnGraph(graphObj);
			if (fnGraph.numVertices() == 0) return;

			zPointArray positions;
			zIntArray edgeConnects;
			fnGraph.getVertexPositions(positions);
			fnGraph.getEdgeData(edgeConnects);

			zPointArray points;
			zIntArray curveVertexCounts;
			points.reserve(edgeConnects.size());
			curveVertexCounts.reserve(edgeConnects.size() / 2);

			for (size_t i = 0; i + 1 < edgeConnects.size(); i += 2)
			{
				const zPoint& start = positions[edgeConnects[i]];
				const zPoint& end = positions[edgeConnects[i + 1]];

				points.emplace_back((float)start.x, (float)start.y, (float)start.z);
				points.emplace_back((float)end.x, (float)end.y, (float)end.z);

				curveVertexCounts.push_back(2);
			}

			writeIndent(out, indent);
			out << "def Xform \"" << tokenFromName(primName) << "\"\n";
			writeIndent(out, indent);
			out << "{\n";
			writePrimXformAttrs(out, frame, indent + 4);
			writeIndent(out, indent + 4);
			out << "def BasisCurves \"Curves\"\n";
			writeIndent(out, indent + 4);
			out << "{\n";
			writeIndent(out, indent + 8);
			out << "uniform token type = \"linear\"\n";
			writeIndent(out, indent + 8);
			out << "uniform token wrap = \"nonperiodic\"\n";
			writeIndent(out, indent + 8);
			out << "point3f[] points = ";
			writePointArray(out, points);
			out << "\n";
			writeIndent(out, indent + 8);
			out << "int[] curveVertexCounts = ";
			writeIntArray(out, curveVertexCounts);
			out << "\n";
			writeIndent(out, indent + 8);
			out << "float[] widths = [0.01]\n";

			if (vertexSequence && !vertexSequence->empty())
			{
				writeIndent(out, indent + 8);
				out << "custom int[] VertexSequence = ";
				writeIntArray(out, *vertexSequence);
				out << "\n";
			}

			writeIndent(out, indent + 4);
			out << "}\n";
			writeIndent(out, indent);
			out << "}\n";
		}

		void writeWorldOpen(std::ofstream& out)
		{
			out << "def Xform \"World\"\n";
			out << "{\n";
			out << "    def Xform \"Geometry\"\n";
			out << "    {\n";
		}

		void writeWorldClose(std::ofstream& out)
		{
			out << "    }\n";
			out << "}\n";
		}

		void writeGroupOpen(std::ofstream& out, const std::string& groupName, int indent)
		{
			writeIndent(out, indent);
			out << "def Xform \"" << tokenFromName(groupName) << "\"\n";
			writeIndent(out, indent);
			out << "{\n";
		}

		void writeGroupClose(std::ofstream& out, int indent)
		{
			writeIndent(out, indent);
			out << "}\n";
		}

		bool exportMeshUsd(zUtilsCore& core, const std::string& path, zObjMesh& meshObj, const std::string& primName, const zTransform* frame = nullptr)
		{
			std::ofstream out;
			if (!prepareUsdTextFile(path, out)) return false;

			writeUsdHeader(out);
			writeWorldOpen(out);
			writeMeshPrim(out, meshObj, primName, frame, 8);
			writeWorldClose(out);
			return true;
		}

		bool exportGraphUsd(zUtilsCore& core, const std::string& path, zObjGraph& graphObj, const std::string& primName, const zTransform* frame = nullptr, const zIntArray* vertexSequence = nullptr)
		{
			std::ofstream out;
			if (!prepareUsdTextFile(path, out)) return false;

			writeUsdHeader(out);
			writeWorldOpen(out);
			writeGraphPrim(out, graphObj, primName, frame, vertexSequence, 8);
			writeWorldClose(out);
			return true;
		}

	}
#else
	namespace
	{
		bool exportMeshUsd(zUtilsCore&, const std::string&, zObjMesh&, const std::string&, const zTransform* = nullptr)
		{
			return false;
		}

		bool exportGraphUsd(zUtilsCore&, const std::string&, zObjGraph&, const std::string&, const zTransform* = nullptr, const zIntArray* = nullptr)
		{
			return false;
		}
	}
#endif
	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsCarbcomn::zTsCarbcomn()
	{


		printHeightDomain = zDomainFloat(0.006, 0.012);

	}


	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsCarbcomn::~zTsCarbcomn() {}

	//---- CREATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::createFieldMeshFromMeshBounds(float cellSize, float offset)
	{
		zFnMesh fnMesh(o_GuideMesh);
		fnMesh.setTransform(base_world, true, false);
		fnMesh.setTransform(base_local, true, true);

		zPoint bbMin, bbMax;
		fnMesh.getBounds(bbMin, bbMax);

		//using the new bounds, calculate the res in x and y, update the pounds to have complete number of cell size
		float lenX = bbMax.x - bbMin.x + offset;
		float lenY = bbMax.y - bbMin.y + offset;
		int resX = ceil(lenX / cellSize);
		int resY = ceil(lenY / cellSize);

		lenX += (resX * cellSize) - lenX;
		lenY += (resY * cellSize) - lenY;

		lenX /= 2;
		lenY /= 2;

		zDomain<zPoint> bb (zPoint(-lenX, -lenX, 0), zPoint(lenX, lenX, 0));

		zFnMeshScalarField fnField(o_field);
		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);
		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);

		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);





		// transform back
		fnMesh.setTransform(base_world, true, true);

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::createFieldMeshFromSectionBounds(float cellSize, float offset)
	{


		float minX = FLT_MAX;
		float minY = FLT_MAX;
		float maxX = FLT_MIN;
		float maxY = FLT_MIN;
		//iterate through all the section graphs bounds and pick the largest one
		for (int i = 0; i < o_sectionGraphs.size(); i++)
		{
			zFnGraph fnGraph(o_sectionGraphs[i]);
			zTransform t = sectionFrames[i];
			fnGraph.setTransform(t, true, false);
			// Transform
			zTransform tLocal;
			tLocal.setIdentity();
			fnGraph.setTransform(tLocal, true, true);

			zPoint pmin, pmax;
			fnGraph.getBounds(pmin, pmax);
			if (pmin.x < minX) minX = pmin.x;
			if (pmin.y < minY) minY = pmin.y;
			if (pmax.x > maxX) maxX = pmax.x;
			if (pmax.y > maxY) maxY = pmax.y;

			//cout << endl << "min" << pmin;
			//cout << endl << "max" << pmin << endl;


			fnGraph.setTransform(t, true, true);
		}


		zPoint bbMin(minX, minY, 0);
		zPoint bbMax(maxX, maxY, 0);

		/*cout << endl << "\n FINAL MIN MAX" << bbMin;
		cout << endl << "min" << bbMin;
		cout << endl << "max" << bbMax << endl;*/

		//fnMesh.getBounds(bbMin, bbMax);

		//using the new bounds, calculate the res in x and y, update the pounds to have complete number of cell size
		float lenX = bbMax.x - bbMin.x;
		float lenY = bbMax.y - bbMin.y;

		int resX = ceil((lenX + (offset*2)) / cellSize);
		int resY = ceil((lenY + (offset*2)) / cellSize);

	/*	cout << endl << "cell" << cellSize;
		cout << endl << "resX" << resX;
		cout << endl << "resY" << resY;*/

		lenX = (resX * cellSize) - lenX;
		lenY = (resY * cellSize) - lenY;

		lenX /= 2;
		lenY /= 2;


		bbMin.x -= offset - lenX;
		bbMin.y -= offset - lenX;
		bbMax.x += offset + lenX;
		bbMax.y += offset + lenX;


		zFnMeshScalarField fnField(o_field);
		//fnField.create()
		fnField.create(bbMin, bbMax, resX, resY, 1, true, false);
		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);

		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);


	}



	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::createFieldMesh(zDomain<zPoint>& bb, int resX, int resY)
	{
		zFnMeshScalarField fnField(o_field);

		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);


		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);


		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);


	}

		//--- SET METHODS

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::setFromJSON(string dir, int _blockID, bool runBothPlanes, bool runPlaneLeft)
	{
		string path = dir + "blockMesh_" + to_string(_blockID) + ".json";
		/*bool pathExist = coreUtils.fileExists(path);
		if (!pathExist)
		{
			throw std::invalid_argument(" error: invalid path. ");
			return;
		}*/

		json j;
		bool jsonCheck = core.json_read(path, j);

		if (jsonCheck) printf("\n %s exists", path.c_str());

		if (!jsonCheck)
		{
			printf("\n %s doesnt exists", path.c_str());
			return;
		}

		blockId = _blockID;
		planarBlock = false;
		printf("\n is planar %s ", to_string(planarBlock));

		//bool flip = j["IsCorner"];
		bool flip = false;

		readJSON(path, _blockID, runBothPlanes, runPlaneLeft, flip);
	}

			ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::setTransforms(bool toLocal)
	{
		if (toLocal)
		{
			zFnGraph fnGraphMedial(o_MedialGraph);
			fnGraphMedial.setTransform(base_world, true, false);
			fnGraphMedial.setTransform(base_local, true, true);

			zFnMesh fnMesh(o_GuideMesh);
			fnMesh.setTransform(base_world, true, false);
			fnMesh.setTransform(base_local, true, true);

			zFnMesh fnMesh_left(o_SliceMesh_Left);
			fnMesh_left.setTransform(base_world, true, false);
			fnMesh_left.setTransform(base_local, true, true);

			zFnMesh fnMesh_right(o_SliceMesh_Right);
			fnMesh_right.setTransform(base_world, true, false);
			fnMesh_right.setTransform(base_local, true, true);

			// section graphs
			for (auto& g : o_sectionGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			// contour graphs
			for (auto& g : o_contourGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			// trim graphs
			for (auto& g : o_trimGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			zFnPointCloud fnCritical_min(criticalMinLayer_pts);
			fnCritical_min.setTransform(base_world, true, false);
			fnCritical_min.setTransform(base_local, true, true);

			zFnPointCloud fnCritical_max(criticalMaxLayer_pts);
			fnCritical_max.setTransform(base_world, true, false);
			fnCritical_max.setTransform(base_local, true, true);

		}
		else
		{
			zFnGraph fnGraphMedial(o_MedialGraph);
			fnGraphMedial.setTransform(base_world, true, true);

			zFnMesh fnMesh(o_GuideMesh);
			fnMesh.setTransform(base_world, true, true);

			zFnMesh fnMesh_left(o_SliceMesh_Left);
			fnMesh_left.setTransform(base_world, true, true);

			zFnMesh fnMesh_right(o_SliceMesh_Right);
			fnMesh_right.setTransform(base_world, true, true);

			// section graphs
			for (auto& g : o_sectionGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			// contour graphs
			for (auto& g : o_contourGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			// trim graphs
			for (auto& g : o_trimGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			zFnPointCloud fnCritical_min(criticalMinLayer_pts);
			fnCritical_min.setTransform(base_world, true, true);

			zFnPointCloud fnCritical_max(criticalMaxLayer_pts);
			fnCritical_max.setTransform(base_world, true, true);

		}
	}

		ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::setCableGraph(string folderDir)
	{
		//read all files in directory

		string fileDir = folderDir;
		zStringArray files;
		core.getFilesFromDirectory(files, fileDir, zJSON);
		printf("\n readCableGraph %i", files.size());

		o_CableGraphs.assign(files.size(), zObjGraph());
		//read all files
		for (int i = 0; i < o_CableGraphs.size(); i++)
		{
			zFnGraph fnGraph(o_CableGraphs[i]);
			fnGraph.from(files[i], zJSON);
		}

		printf("\n readCableGraph %i", o_CableGraphs.size());


		//// Meshes

		fileDir = folderDir + "/Meshes";
		files.clear();
		core.getFilesFromDirectory(files, fileDir, zJSON);
		printf("\n readCableMesh %i", files.size());

		o_CableMeshes.assign(files.size(), zObjMesh());
		//read all files
		for (int i = 0; i < o_CableMeshes.size(); i++)
		{
			zFnMesh fnMesh(o_CableMeshes[i]);
			fnMesh.from(files[i], zJSON);
		}

		printf("\n readCableMesh %i", o_CableMeshes.size());


	}


	//---- GET METHODS

		ZSPACE_TOOLSETS_INLINE vector<zTransform> zTsCarbcomn::getBlockFrames()
	{
		return sectionFrames;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsCarbcomn::getBlockSectionGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_sectionGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_sectionGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

		ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsCarbcomn::getBlockCableProfileGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_trimGraphs_cableprofile.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_trimGraphs_cableprofile)
		{
			out.push_back(&graph);
		}

		return out;
	}

		ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsCarbcomn::getBlockContourGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_contourGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_contourGraphs)
		{
			out.push_back(&graph);
		}
		//printf("\n C++ num of graphs %i", o_contourGraphs.size());

		return out;
	}

			ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsCarbcomn::getBlockTrimGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_trimGraphs.size();
		//numGraphs = o_trimGraphs_bracing.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_trimGraphs)
		//for (auto& graph : o_trimGraphs_bracing)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_TOOLSETS_INLINE zObjPointCloud* zTsCarbcomn::getRawCriticalPoints(bool minHeight)
	{
		return (minHeight) ? &criticalMinLayer_pts : &criticalMaxLayer_pts;
	}

	ZSPACE_TOOLSETS_INLINE zObjMeshScalarField* zTsCarbcomn::getRawFieldMesh()
	{
		return &o_field;
	}

		ZSPACE_TOOLSETS_INLINE zObjMesh* zTsCarbcomn::getRawLeftMesh()
	{
		return &o_SliceMesh_Left;
	}

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsCarbcomn::getRawRightMesh()
	{
		return &o_SliceMesh_Right;
	}

				ZSPACE_TOOLSETS_INLINE bool zTsCarbcomn::isPlanarBlock()
	{
		return planarBlock;
	}

	//---- COMPUTE METHODS

	//---------SLICE MESH----------------
ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_SliceMesh_Regular(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		//get end he
		zItMeshHalfEdge heEnd = heStart;
		//int sideDivisionCount = 1;
		while (heEnd.getVertex().getId() != endVID)
		{
			heEnd = heEnd.getNext().getSym().getNext();
			//sideDivisionCount++;
		}

		//printf("\n sideDivisionCount %i ", sideDivisionCount);


		// walk along spine
		bool exit = false;
		zItMeshHalfEdge he = heStart;

		int featureStart = 0;

		zObjMesh* o_sliceMesh = &o_SliceMesh_Left;
		o_SliceMesh_Right.mesh.clear();

		zFnMesh fnMesh(*o_sliceMesh);


		zFnMesh fnGuideMesh(o_Mesh);
		fnGuideMesh.getVertexPositions(positions);
		fnGuideMesh.getPolygonData(pConnects, pCounts);
		fnMesh.create(positions, pCounts, pConnects);

		fnMesh.setFaceColor(zGREY);
		fnMesh.setVertexColor(zBLACK);



		//get start and end HE of the medial graph
		zFnMesh fnGuide(o_Mesh);
		zPointArray guidePts;
		fnGuide.getVertexPositions(guidePts);
		int vSIDSlice, vEIDSlice;
		vSIDSlice = core.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = core.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = util_getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
		zItMeshHalfEdge heEndSlice = heStartSlice;
		int sideDivisionCount = 1;
		while (heEndSlice.getVertex().getId() != vEIDSlice)
		{
			heEndSlice = heEndSlice.getNext().getSym().getNext();
			sideDivisionCount++;
		}

		heEndSlice = heEndSlice.getNext();

		//get all HE for feature edges, color them based on their side (exterior/interior)


		//get start of all feature edges
		zItMeshHalfEdgeArray hesFeatureInner, hesFeatureOuter;
		//zItMeshHalfEdge heTemp0 = heStartSlice;
		zItMeshHalfEdge heTemp = heStartSlice;
		zItMeshHalfEdge heSide;

		const int featureWalkNum = (FeaturedNumStrides.size());
		const int bridgingStride = FeaturedNumStrides[((FeaturedNumStrides.size() - 1) / 2)];
		printf("\n bridging stride %i %i %i ", (FeaturedNumStrides.size() - 1), ((FeaturedNumStrides.size() - 1) / 2), FeaturedNumStrides[((FeaturedNumStrides.size() - 1) / 2)]);
		bool isBackSide = true;
		bool isCornerEdge = false;
		bool fromInToOut = true;
		int innerCounter = 0;
		int outerCounter = 0;
		zColor innerColor = _col_in_feature;
		zColor outerColor = _col_out_feature;
		bool ignoreOuterEdge = false;


		////walk on the lower part of the mesh and color corner vertices
		//zItMeshHalfEdge heLowerTemp = heStartSlice.getSym().getNext();
		//int cornerCounter = 0;
		//while (heLowerTemp.getVertex().getId() != heStartSlice.getStartVertex().getId())
		//{
		//	zItMeshHalfEdge heS1, heS2;
		//	heS1 = heLowerTemp;
		//	heLowerTemp = heLowerTemp.getNext().getSym().getNext();
		//	heS2 = heLowerTemp;
		//	//get the angle between the two sides heS1 and heS2
		//	float angle = heS1.getVector().angle(heS2.getVector());
		//	bool chk = angle > 20;
		//	if (chk)
		//	{
		//		printf("\n angle0 %1.4f", angle);
		//		heLowerTemp.getStartVertex().setColor(zORANGE);
		//		cornerCounter++;
		//	}
		//}

		//get start inner corner vertex
		zItMeshVertex vCornerStart(*o_sliceMesh, StartCornerVID);
		zItMeshVertex vCornerEnd;
		vCornerStart.setColor(_col_in_corner_st);
		zItMeshHalfEdge startCornerHE;
		//get HE by finding all connected HE and choosing the one  with end valence of 3
		zItMeshHalfEdgeArray hesTemp;
		vCornerStart.getConnectedHalfEdges(hesTemp);
		for (zItMeshHalfEdge e : hesTemp)
		{
			if (e.getVertex().checkValency(3))
			{
				vCornerEnd = e.getVertex();
				break;
			}
		}


		//color vertex on the other end
		int heCornerEdgeTopId = -1;
		zItMeshHalfEdge heCornerStart, heCornerEnd;
		//walk on the strides. color feature ones. then walk again to color the corner edges
		for (int k = 0; k < featureWalkNum * 2; k++)
		{
			int featureCounter = k % FeaturedNumStrides.size();
			int blockStride = FeaturedNumStrides[featureCounter];
			if (blockStride == bridgingStride)
			{
				isCornerEdge = true;
			}
			innerColor = isCornerEdge ? _col_in_corner : _col_in_feature;
			outerColor = isCornerEdge ? _col_out_corner : _col_out_feature;

			bool chk0 = isBackSide && heTemp.getStartVertex().getId() == vCornerStart.getId();
			if (chk0)
			{
				innerColor = vCornerStart.getColor();
				startCornerHE = heTemp;
			}
			bool chk1 = !isBackSide && heTemp.getStartVertex().getId() == vCornerEnd.getId();
			if (chk1)
			{
				outerColor = _col_out_corner_st;
			}
			zColor edgeColor = isBackSide ? innerColor : outerColor;
			heTemp.getEdge().setColor(edgeColor);
			heTemp.getStartVertex().setColor(edgeColor);
			heTemp.getVertex().setColor(edgeColor);

			zItMeshHalfEdge heRest = heTemp;
			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				heRest = heRest.getNext().getSym().getNext();
				heRest.getEdge().setColor(edgeColor);
				heRest.getVertex().setColor(edgeColor);
			}
			if (chk0) heCornerEnd = heRest;
			for (int i = 0; i < blockStride; i++)
			{
				heTemp = heTemp.getSym().getNext().getNext();
			}

			if (blockStride == bridgingStride)
			{
				isBackSide = !isBackSide;
				isCornerEdge = true;
			}
			else
			{
				isCornerEdge = false;
			}
		}

		//walk on the corner edges till you reach the next corner edge. Ignore edges that are colored

		printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());
	}

	//Slice mesh: helper methods
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_MedialGraph(zObjMesh& o_Mesh, int startVID, int endVID)
	{

		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		//printf("\n hestart %i %i ", heStart.getStartVertex().getId(), heStart.getVertex().getId());

		zItMeshHalfEdge heStart_bottom = heStart;
		heStart_bottom = heStart_bottom.getPrev().getSym().getPrev();
		heStart_bottom = heStart_bottom.getPrev().getSym().getPrev();

		// compute graph
		zPointArray positions;
		zIntArray eConnects;

		zItMeshHalfEdge he = heStart;
		zItMeshHalfEdge he_bottom = heStart_bottom;

		//zPoint po = (startPoint + tmpPositions[heStart_bottom.getVertex().getId()]) * 0.5;
		zPoint po = startPoint;

		positions.push_back(po);
		bool exit = false;

		do
		{

			if (he.getVertex().getId() == endVID) exit = true;

			eConnects.push_back(positions.size() - 1);
			eConnects.push_back(positions.size());

			//zPoint p1 = (tmpPositions[he.getVertex().getId()] + tmpPositions[he_bottom.getStartVertex().getId()]) * 0.5;
			zPoint p1 = tmpPositions[he.getVertex().getId()];

			positions.push_back(p1);

			if (!exit)
			{
				he = he.getNext().getSym().getNext();
				he_bottom = he_bottom.getPrev().getSym().getPrev();
			}

		} while (!exit);

		zFnGraph fnMedial(o_MedialGraph);
		fnMedial.create(positions, eConnects);

		fnMedial.setEdgeWeight(5);
		fnMedial.setEdgeColor(zGREEN, false);
	}
		//----------PRINT BLOCKS----------
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_PrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, bool compFrames, bool compSDF)
	{


		bool frameCHECKS = false;
		bool geomCHECKS = true;
		bool sdfCHECKS = true;

		int minCriticalPtsCount = INT_MAX;
		float bestPlaneSpacing = FLT_MAX;
		printHeightDomain = _printHeightDomain;

		if (compFrames)
		{
			zVector norm(0, 0, 1);
			vector<zItMeshHalfEdgeArray> vLoops;
			zObjMesh oMesh_top, oMesh_bottom;

			computeVLoops(o_SliceMesh_Left, medialIDS, FeaturedNumStrides, norm, vLoops, oMesh_top, oMesh_bottom);

			zScalarArray scalars;
			computeGeodesicScalars(o_SliceMesh_Left, vLoops, scalars, true);

			o_sectionMeshes.clear();

			computeGeodesicContours(vLoops, scalars, 0.008, oMesh_top, oMesh_bottom, o_sectionMeshes);;
			createSectionGraphs(o_sectionMeshes, o_sectionGraphs);
			o_sectionMeshesPar.clear();
			o_sectionMeshesPar.assign(o_sectionMeshes.size(), zObjMesh());

			sectionFrames.clear();
			sectionFrames.assign(o_sectionGraphs.size(), zTransform());
			compute_PrintBlock_ComputeTrimGraphs();

			o_contourHeightLines.clear();
			o_contourHeightLines.assign(o_sectionGraphs.size(), zObjGraph());
		}

		if (compSDF)
		{
			printf("\n \n  SDF \n \n");

			compute_SDF(allSDFLayers, numSDFlayers, funcNum, numSmooth, printLayerWidth);
		}

	}

	//Print blocks: helper methods
		//Print blocks: frame methods
		//Print blocks: section methods
		//Print blocks: check methods
		ZSPACE_TOOLSETS_INLINE bool zTsCarbcomn::check_SDF_LayerHeights()
	{
		float minLayerHeight = 10;
		float maxLayerHeight = 0;

		//float minLayerHeight = 0.005;
		//float maxLayerHeight = 0.015;
		int minHeightGraphID = -1;
		//bool
		checkOranges = true;
		//bool
		checkMagentas = true;

		zFnPointCloud fnCritical_min(criticalMinLayer_pts);
		zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

		fnCritical_min.clear();
		fnCritical_max.clear();

		int r0 = 0;
		int r1 = floor(o_sectionGraphs.size() * 0.5) - 1;
		int r2 = floor(o_sectionGraphs.size() * 0.5);
		int r3 = (o_sectionGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		bool deckBlock = true;

		int end = floor(o_sectionGraphs.size() * 0.5);

		float printLength = 0;

		// PRINT LAYERS
		for (int j = 0; j < end; j++)
		{
			int kStart = 0;
			int kEnd = 2;

			for (int k = kStart; k < kEnd; k++)
			{
				int i = (k == 0) ? j : j + end;

				if (!deckBlock) i = j;

				if (i == r0) continue;

				if (deckBlock && i == r2) continue;


				zVector norm(sectionFrames[i](2, 0), sectionFrames[i](2, 1), sectionFrames[i](2, 2));
				norm *= -1;
				norm.normalize();

				zVector prevNorm(sectionFrames[i - 1](2, 0), sectionFrames[i - 1](2, 1), sectionFrames[i - 1](2, 2));
				zVector prevOrigin(sectionFrames[i - 1](3, 0), sectionFrames[i - 1](3, 1), sectionFrames[i - 1](3, 2));

				float layerWidth = 0.025;

				int orangeCounter = 0;
				int magentaCounter = 0;

				int innerSideCounter = 0;
				int outerSideCounter = 0;

				//The following check to see if there is a problem with the geometry/section graph
				//01-check if the graph is closed by checking if any vertex is boundary (val (1))


				for (zItGraphVertex v(o_contourGraphs[i]); !v.end(); v++)
				{

					zPoint p = v.getPosition();


					zPoint p1 = p + norm * 1.0;

					zPoint intPt;
					bool check = core.line_PlaneIntersection(p, p1, prevNorm, prevOrigin, intPt);

					float layerHeight = intPt.distanceTo(p);

					if (i == 1 && v.getId() == 0)
					{
						printf("\n //////  \n %i | %i | %1.4f ", i, v.getId(), layerHeight);
						cout << "\n CP " << sectionFrames[i];
						cout << "\n PP " << sectionFrames[i -1];
					}

					maxLayerHeight = (layerHeight > maxLayerHeight) ? layerHeight : maxLayerHeight;
					minHeightGraphID = (layerHeight < minLayerHeight) ? i : minHeightGraphID;
					minLayerHeight = (layerHeight < minLayerHeight) ? layerHeight : minLayerHeight;

					if (layerHeight < printHeightDomain.min)fnCritical_min.addPosition(p);
					if (layerHeight > printHeightDomain.max)fnCritical_max.addPosition(p);

				}




				for (zItGraphEdge e(o_sectionGraphs[i]); !e.end(); e++)
				{
					printLength += e.getLength();
				}


			}
		}

		fnCritical_min.setVertexColor(zORANGE);
		fnCritical_max.setVertexColor(zMAGENTA);


		bool out = true;

		if (out)
		{
			if (minLayerHeight < printHeightDomain.min) out = false;
			if (minLayerHeight > printHeightDomain.max) out = false;

			if (maxLayerHeight < printHeightDomain.min) out = false;
			if (maxLayerHeight > printHeightDomain.max) out = false;
		}

		actualPrintHeightDomain.min = minLayerHeight;
		actualPrintHeightDomain.max = maxLayerHeight;

		printf("\n block| %1.4f %1.4f| %1.1f | < min ht intersectionPts %i  | > max ht intersectionPts %i ", minLayerHeight, maxLayerHeight, printLength,  fnCritical_min.numVertices(), fnCritical_max.numVertices());


		return out;
	}
		ZSPACE_TOOLSETS_INLINE bool zTsCarbcomn::check_sectionGraphGeomCheck(zObjGraph& graph)
	{
		//checks if the graph pass all geometry checks
		//Check 1 : check if the graph is closed
		//Check 2 : check if all four corners are found
		bool chkClosed = true;
		bool corner0 = false;
		bool corner1 = false;
		bool corner2 = false;
		bool corner3 = false;
		bool featureInner = false;
		bool featureOuter = false;
		bool chkCorners = false;
		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (v.checkValency(1))
			{
				chkClosed = false;
			}
			if (v.getColor() == _col_in_corner_st) corner0 = true;
			if (v.getColor() == _col_out_corner_st) corner1 = true;
			if (v.getColor() == _col_in_corner) corner2 = true;
			if (v.getColor() == _col_out_corner) corner3 = true;
			if(v.getColor() == _col_in_feature) featureInner = true;
			if(v.getColor() == _col_out_feature) featureOuter = true;

			chkCorners = corner0 && corner1 && corner2 && corner3;
			/*if (chkCorners && chkClosed)
			{
				break
			}*/
		}
		if (chkClosed && chkCorners)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	//Print blocks: trim methods

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_PrintBlock_ComputeTrimGraphs()
	{
		o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing.clear();
		o_trimGraphs_bracing.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing_slots.clear();
		o_trimGraphs_bracing_slots.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_features_hard.clear();
		o_trimGraphs_features_hard.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_features_soft.clear();
		o_trimGraphs_features_soft.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_SlotSide.clear();
		o_trimGraphs_SlotSide.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_seamAlignment.clear();
		o_trimGraphs_seamAlignment.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing_flat.clear();
		o_trimGraphs_bracing_flat.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing_slots_flat.clear();
		o_trimGraphs_bracing_slots_flat.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_cableprofile.clear();
		o_trimGraphs_cableprofile.assign(o_sectionGraphs.size(), zObjGraph());


		for (int i = 0; i < o_sectionGraphs.size(); i++)
		{
			//The following check to see if there is a problem with the geometry/section graph
				//01-check if the graph is closed by checking if any vertex is boundary (val (1))

			if (i == 0) continue;

			for (zItGraphVertex v(o_sectionGraphs[i]); !v.end(); v++)
			{
				if (v.checkValency(1))
				{
					printf("\n section[%i] is not closed!", i);
				}
			}

			zFnGraph fng;
			zObjGraphArray allTrimGraphs;

			compute_TrimGraphs_BoundaryFeature(i, o_trimGraphs_features_hard[i], o_trimGraphs_features_soft[i]);
			compute_TrimGraphs_SlotSide(i, o_trimGraphs_SlotSide[i]);
			//compute_TrimGraphs_SeamAlignment(i, o_trimGraphs_seamAlignment[i]);

			bool chk = (i > (o_sectionGraphs.size() * 0.9));
			compute_TrimGraphs_BracingWall(o_sectionGraphs[i], o_trimGraphs_bracing[i], chk);
			fng = zFnGraph(o_trimGraphs_bracing[i]);
			fng.setEdgeColor(zBLUE);
			allTrimGraphs.push_back(o_trimGraphs_bracing[i]);

			fng = zFnGraph(o_trimGraphs_features_hard[i]);
			fng.setEdgeColor(zRED);
			fng = zFnGraph(o_trimGraphs_features_soft[i]);
			fng.setEdgeColor(zGREEN);
			fng = zFnGraph(o_trimGraphs_SlotSide[i]);
			fng.setEdgeColor(zGREEN);


			allTrimGraphs.push_back(o_trimGraphs_features_hard[i]);
			allTrimGraphs.push_back(o_trimGraphs_features_soft[i]);


			util_combineMultipleGraphs(allTrimGraphs, o_trimGraphs[i]);
		}
		printf("\n trims finished");

	}

		ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_TrimGraphs_BoundaryFeature(int graphId, zObjGraph& outGraph_hardFeature, zObjGraph& outGraph_softFeature)
	{
		//create a trim graph at each trim points. Trim points are all feature curves and corner
		//The graph will be made by avg the two outgoing vector from the feature vertex

		zPrintParamSDF _printParameters;

		zPointArray posittions;
		zIntArray eConnect;
		float length = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior + (_printParameters.printWidthInterior * 4.5);
		float angleThreshold = 30;
		float lengthTolerance = 0.05;
		zItGraphVertexArray features_hard;
		zItGraphVertexArray features_soft;
		//check if there is another point that has already been added within the tolerance
		//check the next vertex, if it is within the same tolerance, combine the two together and skip the next vertex
		//To do that, probably since we don't know if there are some points in between, better to iterate through all the vertices first, and then iterate through the ones that passes that check.

		for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
		{

			bool cornerByTopology = v.getColor() == _col_in_corner_st ||
									v.getColor() == _col_out_corner_st ||
									v.getColor() == _col_in_corner ||
									v.getColor() == _col_out_corner;

			bool otherFeature = v.getColor() == _col_in_feature ||
								v.getColor() == _col_out_feature;

			if (cornerByTopology) features_hard.push_back(v);
			else
			{
				zItGraphHalfEdgeArray hes;
				v.getConnectedHalfEdges(hes);
				zVector v0 = hes[0].getVector();
				zVector v1 = hes[1].getVector();
				v0.normalize();
				v1.normalize();

				float angle = v0.angle(v1 * -1);
				bool angleChk = angle >= angleThreshold;
				if (angleChk)
				{
					if (otherFeature)
						features_hard.push_back(v);

					else
						features_soft.push_back(v);
				}
			}
		}

		for (int index = 0; index < features_hard.size(); index++)
		{
			zVector v0, v1, result;
			zPoint pt;

			int nextIndex = (index + 1) % features_hard.size();
			//check the distance to the next index, if it is within a threshold, get middle vertex as feature and skip the next feature
			float dist = features_hard[index].getPosition().distanceTo(features_hard[nextIndex].getPosition());
			if (dist <= lengthTolerance)
			{
				//get the two vectors, average them twice. They might be sharing the same vector, but not necessarily if the have points in between
				v0 = util_averageVectorsAtGraphVertex(features_hard[index]);
				v1 = util_averageVectorsAtGraphVertex(features_hard[nextIndex]);
				result = (v0 + v1) / 2.0;
				result.normalize();
				pt = (features_hard[index].getPosition() + features_hard[nextIndex].getPosition()) / 2.0;
				//skip the next one
				index++;
			}
			else
			{
				result = util_averageVectorsAtGraphVertex(features_hard[index]);
				pt = features_hard[index].getPosition();
			}

			result *= length;

			zPoint p0 = pt + result;
			zPoint p1 = pt - result;
			posittions.push_back(p0);
			eConnect.push_back(posittions.size() - 1);
			posittions.push_back(p1);
			eConnect.push_back(posittions.size() - 1);

		}
		zFnGraph fng(outGraph_hardFeature);
		fng.create(posittions, eConnect);


		posittions.clear();
		eConnect.clear();
		for (int index = 0; index < features_soft.size(); index++)
		{
			zVector v0, v1, result;
			zPoint pt;

			int nextIndex = (index + 1) % features_soft.size();
			//check the distance to the next index, if it is within a threshold, get middle vertex as feature and skip the next feature
			float dist = features_soft[index].getPosition().distanceTo(features_soft[nextIndex].getPosition());
			if (dist <= lengthTolerance)
			{
				//get the two vectors, average them twice. They might be sharing the same vector, but not necessarily if the have points in between
				v0 = util_averageVectorsAtGraphVertex(features_soft[index]);
				v1 = util_averageVectorsAtGraphVertex(features_soft[nextIndex]);
				result = (v0 + v1) / 2.0;
				result.normalize();
				pt = (features_soft[index].getPosition() + features_soft[nextIndex].getPosition()) / 2.0;
				//skip the next one
				index++;
			}
			else
			{
				result = util_averageVectorsAtGraphVertex(features_soft[index]);
				pt = features_soft[index].getPosition();
			}

			result *= length;
			zPoint p0 = pt + result;
			zPoint p1 = pt - result;
			posittions.push_back(p0);
			eConnect.push_back(posittions.size() - 1);
			posittions.push_back(p1);
			eConnect.push_back(posittions.size() - 1);

		}
		 fng = zFnGraph(outGraph_softFeature);
		fng.create(posittions, eConnect);

		//for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
		//{
		//	bool cornerByTopology = v.getColor() == _colorCornersStart ||
		//		v.getColor() == _colorCornersEnd ||
		//		v.getColor() == _colorCornersInner ||
		//		v.getColor() == _colorCornersOuter;
		//	bool otherFeature = v.getColor() == _colorFeatureInner ||
		//		v.getColor() == _colorFeatureOuter;
		//	zItGraphHalfEdgeArray hes;
		//	v.getConnectedHalfEdges(hes);
		//	zVector v0 = hes[0].getVector();
		//	zVector v1 = hes[1].getVector();
		//	v0.normalize();
		//	v1.normalize();
		//	float angle = v0.angle(v1 * -1);
		//	bool angleChk = angle >= angleThreshold;
		//	if (cornerByTopology || angleChk)
		//	{
		//		zVector result = (v0 + v1) / 2.0;
		//		result.normalize();
		//		result *= length;
		//		zPoint p0 = v.getPosition() + result;
		//		zPoint p1 = v.getPosition() - result;
		//		posittions.push_back(p0);
		//		eConnect.push_back(posittions.size() - 1);
		//		posittions.push_back(p1);
		//		eConnect.push_back(posittions.size() - 1);
		//	}
		//
		//}
		//zFnGraph fng(o_trimGraphs[graphId]);
		//fng.create(posittions, eConnect);
		////printf("\n Trim1 nV nE - %i | %i", fng.numVertices(), fng.numEdges());


	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_TrimGraphs_SlotSide(zObjGraph& sectionGraph, zObjGraph& outGraph)
	{
		zFnGraph inFnGraph(sectionGraph);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = _col_in_corner_st;
		zColor endColor = _col_out_corner_st;
		for (zItGraphVertex v(sectionGraph); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV);
		gPts.push_back(endV);
		gEdges.push_back(0);
		gEdges.push_back(1);
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);

	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_TrimGraphs_SlotSide(int graphId, zObjGraph& outGraph)
	{
		compute_TrimGraphs_SlotSide(o_sectionGraphs[graphId], outGraph);

	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_TrimGraphs_BracingWall(zObjGraph& sectionGraph, zObjGraph& outGraph, bool remove_firstLast)
	{
		zFnGraph fnSectionG(sectionGraph);
		printf("\n remove_firstLast %i \n", remove_firstLast);
		//walk on the graph till u reach a feature vertex and store it

		zItGraphVertexArray innerVertx, outerVertx;


		util_innerOuter(sectionGraph,false , innerVertx, outerVertx);

		//make graphs between verticies
		zPointArray gPositions;
		zIntArray gEdgeCOnnects;
		zColorArray gColors;
		int startID = (remove_firstLast) ? 1 : 0;
		int endID = (remove_firstLast) ? innerVertx.size() - 1 : innerVertx.size();

		for (int i = startID; i < endID; i++)
		{

			gPositions.push_back(innerVertx[i].getPosition());
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			gColors.push_back(innerVertx[i].getColor());

			gPositions.push_back(outerVertx[i].getPosition());
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			gColors.push_back(outerVertx[i].getColor());

		}
		zFnGraph fnG(outGraph);
		fnG.create(gPositions, gEdgeCOnnects);
		fnG.setVertexColors(gColors);
		//printf("\n graph[%i] %i | %i", graphId, gPositions.size(), gEdgeCOnnects.size());

	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_TrimGraphs_BracingWall(int graphId, zObjGraph& outGraph)
	{
		compute_TrimGraphs_BracingWall(o_sectionGraphs[graphId], outGraph);
	}


	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_innerOuter(zObjGraph& sectionGraph, bool addEndStart, zItGraphVertexArray& innerVertx, zItGraphVertexArray& outerVertx)
	{
		innerVertx.clear();
		outerVertx.clear();
		zFnGraph fnSectionG(sectionGraph);
		//get inner and outer edges
		bool found = false;
		zItGraphHalfEdgeArray hesInner, hesOuter;
		found = util_getShortestHEsBetweenColors(sectionGraph, _col_in_corner_st, _col_in_corner, hesInner);
		found = util_getShortestHEsBetweenColors(sectionGraph, _col_out_corner_st, _col_out_corner, hesOuter);

		//walk on the graph till u reach a feature vertex and store it

		zItGraphHalfEdge heStart, heEnd, heTemp;

		heStart = hesInner[0];
		heEnd = hesInner[hesInner.size() - 1];
		heTemp = heStart;
		if (addEndStart) innerVertx.push_back(heTemp.getStartVertex());
		while (heTemp != heEnd)
		{
			if (heTemp.getStartVertex().getColor() == _col_in_feature)
				innerVertx.push_back(heTemp.getStartVertex());
			heTemp = heTemp.getNext();
		}
		if (addEndStart) innerVertx.push_back(heEnd.getVertex());

		heStart = hesOuter[0];
		heEnd = hesOuter[hesOuter.size() - 1];
		heTemp = heStart;
		if (addEndStart) outerVertx.push_back(heTemp.getStartVertex());
		while (heTemp != heEnd)
		{
			if (heTemp.getStartVertex().getColor() == _col_out_feature)
				outerVertx.push_back(heTemp.getStartVertex());
			heTemp = heTemp.getNext();
		}
		if (addEndStart) outerVertx.push_back(heEnd.getVertex());

		//check if the two arrays are not the same size, return
		if (innerVertx.size() != outerVertx.size())
		{
			printf("\n ERROR!  inner and outer vertices are not the same size! inner | outer  %i | %i", innerVertx.size(), outerVertx.size());
			return;
		}


	}

	//SDF MAIN method
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_SDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth)
	{

		o_contourGraphs.clear();
		o_contourGraphs.assign(o_sectionGraphs.size(), zObjGraph());
		o_contourGraphs_flatten.clear();
		o_contourGraphs_flatten.assign(o_sectionGraphs.size(), zObjGraph());

		/*o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());*/

		o_raftGraphs.clear();
		o_raftGraphs.assign(1, zObjGraph());

		printf("\n num frames : %i ", o_sectionGraphs.size());


		int r0 = 0;
		int r1 = floor(o_sectionGraphs.size() * 0.5) - 1;
		int r2 = floor(o_sectionGraphs.size() * 0.5);
		int r3 = (o_sectionGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		int end = o_sectionGraphs.size();
		numSDFlayers = (numSDFlayers > end) ? end : numSDFlayers;
		numSDFlayers = (allSDFLayers) ? end : numSDFlayers;

		int j = 1;
		int max_layers = numSDFlayers;

		printf("\n Start: %d | Max: %d\n", j, max_layers);
		for (;j < max_layers; ++j)
		{
			compute_BlockSDF_NonPlanar(funcNum, numSmooth, j, (j % 2 == 0));
		}
	}

	//SDF sub methods
	// Planar SDF slicing variants were removed; Carbcomn now uses the non-planar SDF path.

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::compute_BlockSDF_NonPlanar(int funcNum, int numSmooth, int graphId, bool alternate )
	{
		zPrintParamSDF _printParameters;

		if (graphId >= o_sectionGraphs.size())return;

		printf("\n fREP graphID %i | funcNum %i  ", graphId, funcNum);


		zFnGraph fnTmpGraph(o_sectionGraphs[graphId]);
		zPointArray positions;
		fnTmpGraph.getVertexPositions(positions);

		zPoint minBB, maxxBB;
		fnTmpGraph.getBounds(minBB, maxxBB);

		zPoint refPt = minBB;/*positions[0];*/ // update to medial point

		zObjMesh oUnrolledMesh;// = o_sectionMeshes[graphId];


		// Unroll
		int num_BSF = 0;

		zInt2DArray oriVertex_UnrollVertex_map;
		unordered_map<zIntPair, int, zPair_hash> oriFaceVertex_UnrollVertex;
		zItGraphVertexArray bsf_Vertices;
		zIntPairArray bsf_vertexPairs;

		zObjGraph oDualGraph;
		zObjMesh o_projectionMesh = o_sectionMeshes[graphId];
		//getPokeMesh(o_sectionMeshes[graphId], o_projectionMesh);
		creatUnrollMesh(o_projectionMesh, oUnrolledMesh, oDualGraph, oriVertex_UnrollVertex_map, oriFaceVertex_UnrollVertex, bsf_Vertices, bsf_vertexPairs);
		num_BSF = bsf_Vertices.size();
		zVector unitZ = zVector(0, 0, 1);
		zTransform newFrame;// = core.getPlaneFromOrigin_Normal(refPt, unitZ);
		unrollMesh(o_projectionMesh, oUnrolledMesh, oDualGraph, oriVertex_UnrollVertex_map, oriFaceVertex_UnrollVertex, bsf_vertexPairs, newFrame);
		mergeMesh(oUnrolledMesh);


		zObjGraph oFlatGraph;// = o_sectionGraphs[graphId];
		//setPtGraph(oFlatGraph, refPt, false, false, true);

		createBoundaryEdgeGraph(oUnrolledMesh, true, oFlatGraph);

		//sectionFrames[graphId] = newFrame;
		zTransform t = newFrame;

		// Transform
		zTransform tLocal;
		tLocal.setIdentity();

		//transform all graphs to local
		zFnGraph fng;
		/*fng = zFnGraph(oFlatGraph);
		fng.setTransform(t, true, false);
		fng.setTransform(tLocal, true, true);*/

		//make flatten version of the trimGraph


		//transformAllGraphs(graphId,t, true);
		zPoint refPt2(0, 0, 0);
		zObjGraph trimGraphs_bracing_flat;
		zObjGraph trimGraphs_bracing_slots_flat;

		bool chk = (graphId > o_sectionGraphs.size() * 0.9f) ? true : false;
		//bool chk = false;
		//printf("\n graphId : %i (o_sectionGraphs.size() * 0.9f) : %1.2f\n", graphId, (o_sectionGraphs.size() * 0.9f));
		compute_TrimGraphs_BracingWall(oFlatGraph, trimGraphs_bracing_flat, chk);
		o_trimGraphs_bracing_flat[graphId] = trimGraphs_bracing_flat;

		compute_TrimGraphs_BracingWall(oFlatGraph, trimGraphs_bracing_slots_flat, chk);
		o_trimGraphs_bracing_slots_flat[graphId] = trimGraphs_bracing_slots_flat;

		zObjGraph o_trimGraphs_slotSide_flat;
		compute_TrimGraphs_SlotSide(oFlatGraph, o_trimGraphs_slotSide_flat);

		zFnMeshScalarField fnField(o_field);

		zFnGraph fnGraph(oFlatGraph);
		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));

		//Polygon and offset
		zScalarArray polyField, scalar_offset_outer, scalar_offset_inner;
		getScalars_offset(oFlatGraph, numSmooth, polyField, scalar_offset_outer, scalar_offset_inner);

		//transfer color to sectiongraph
		zColorArray eColors;
		fnGraph.getEdgeColors(eColors);
		fnTmpGraph.setEdgeColors(eColors, false);

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;
		zPoint oo(0, 0, 0);
		planeXY = core.getTransformFromOrigin_Normal(oo, zAxis);

		zObjGraph slotGraph, splitGraph;
		float graphLength = 0.10;
		util_computeSlotGraph(planeXY, oFlatGraph, graphLength, graphId % 2 == 0, slotGraph);
		//splitGraph_1(planeXY, oFlatGraph, _printParameters.offset_2nd_exterior, pWidth * 1.5, splitGraph);
		util_computeSplitGraph_xy(oFlatGraph, splitGraph);

		zScalarArray scalar_slot1;
		if (funcNum >= 2)
		{
			fnField.getScalarsAsEdgeDistance(scalar_slot1, slotGraph, _printParameters.slotStartWidth, false);
		}

		zScalarArray scalar_interiorBracing;
		zScalarArray scalar_bracing;
		zScalarArray scalar_bracingSlots;


		if (funcNum >= 3)
		{

			getScalars_3dp_wall_bracing(oFlatGraph, trimGraphs_bracing_flat, trimGraphs_bracing_slots_flat ,_printParameters.slotIterating, graphId % 2 == 0, scalar_interiorBracing, scalar_bracing, scalar_bracingSlots);
		}


		zScalarArray scalar_triangles;
		zScalarArray scalar_boolean_trianglesInner;
		zScalarArray booleanField_0;

		if (funcNum >= 4)
		{

			getScalars_3dp_wall_triangles(oFlatGraph, scalar_triangles, chk);
			fnField.boolean_subtract(scalar_offset_inner, scalar_triangles, scalar_boolean_trianglesInner, false);
			if (numSmooth > 0) fnField.smoothField(scalar_boolean_trianglesInner, numSmooth); // smooth field
			fnField.boolean_subtract(scalar_boolean_trianglesInner, scalar_interiorBracing, booleanField_0, false);
		}

		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(scalar_offset_outer, booleanField_0, booleanField_1, false);


		zScalarArray scalar_booleanSlot;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_1, scalar_slot1, scalar_booleanSlot, false);

		float sdfWidth = _printParameters.printWidthInterior / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(scalar_offset_inner, zFieldSDF, sdfWidth);
			break;

		case 2:
			fnField.setFieldValues(scalar_slot1, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(scalar_interiorBracing, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(scalar_triangles, zFieldSDF, sdfWidth);
			break;
		case 5:
			fnField.setFieldValues(scalar_boolean_trianglesInner, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 7:
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, sdfWidth);
			break;


		case 8:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, _printParameters.printWidthInterior / 2.0);
			break;
		}

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		int pres = 3;
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, zVector(0, 0, 1), pres, 0.001);
		util_merge_graph(o_contourGraphs[graphId], 0.005);
		cleanContourGraph(graphId);

		zFnGraph fngraph(o_contourGraphs[graphId]);
		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);


		// project to  section Mesh
		zFnGraph fnContour(o_contourGraphs[graphId]);
		zPointArray contourPositions, projectedPositions;
		zIntArray faceIDs;

		zVectorArray pNorms, pNormsTemp;

		//project contour back to section mesh
		barycentericProjection_triMesh(o_contourGraphs[graphId], oUnrolledMesh, o_projectionMesh, pNorms);

		auto project_slot = [this](zObjGraph& graph, zObjMesh& inMesh, zObjMesh& projMesh)
		{
			zFnGraph fnGraph(graph);

			zPointArray positions;
			fnGraph.getVertexPositions(positions);

			bool done = false;

			for (int i = 0; i < 2; ++i)
			{
				zPoint& pos = positions[i];

				for (zItMeshFace face(inMesh); !face.end(); face++)
				{
					zPointArray fVerts;
					face.getVertexPositions(fVerts);
					if (core.pointInTriangle(pos, fVerts[0], fVerts[1], fVerts[2]))
					{
						zPoint pos1_bary, pos2_bary;
						getBaryCentricCoordinates_triangle(pos, fVerts[0], fVerts[1], fVerts[2], pos1_bary);
						getBaryCentricCoordinates_triangle(positions[i^1], fVerts[0], fVerts[1], fVerts[2], pos2_bary);

						zItMeshFace fProjection(projMesh, face.getId());

						zPointArray fVerts_projection;
						fProjection.getVertexPositions(fVerts_projection);

						zPoint projectionPt1, projectionPt2;
						getProjectionPoint_triangle(pos1_bary, fVerts_projection[0], fVerts_projection[1], fVerts_projection[2], projectionPt1);
						getProjectionPoint_triangle(pos2_bary, fVerts_projection[0], fVerts_projection[1], fVerts_projection[2], projectionPt2);

						positions[i] = projectionPt1;
						positions[i ^ 1] = projectionPt2;

						done = true;

						break;
					}
				}

				if (done)
					break;
			}

			fnGraph.setVertexPositions(positions);
		};

		//project splitGraph back to section mesh
		project_slot(slotGraph, oUnrolledMesh, o_projectionMesh);
		o_trimGraphs_SlotSide[graphId] = slotGraph;

		//fnContour.setVertexPositions(projectedPositions);
		fnContour.getVertexPositions(contourPositions);

		fnContour.setEdgeColor(zBLUE);
		fnContour.setEdgeWeight(3);

		o_sectionMeshesPar[graphId] = oUnrolledMesh;
		o_contourGraphs_flatten[graphId] = oFlatGraph;
	}

			//EXPORT MAIN method
	ZSPACE_TOOLSETS_INLINE bool zTsCarbcomn::exportUSD_update(string pathCurrent, string dir)
	{
		string folderName = dir + "/" + to_string(blockId);

		zObjMesh oMesh;
		zFnMesh fn(oMesh);
		json j;
		string blockPath = pathCurrent + "blockMesh_" + to_string(blockId) + ".json";
		bool fileChk = fn.json_read(blockPath, j);
		if (!fileChk) return false;

		std::filesystem::create_directories(folderName);
		for (const auto& entry : std::filesystem::directory_iterator(folderName)) std::filesystem::remove_all(entry.path());

		fn.from(j);

		string outName = folderName + "/block_" + to_string(blockId) + "_carbcomn.usda";
		std::ofstream out;
		if (!prepareUsdTextFile(outName, out)) return false;

		writeUsdHeader(out);
		writeWorldOpen(out);

		writeGroupOpen(out, "block_meshes", 8);
		writeMeshPrim(out, oMesh, "block", nullptr, 12);
		writeMeshPrim(out, o_SliceMesh_Left, "block_left", nullptr, 12);
		writeMeshPrim(out, o_SliceMesh_Right, "block_right", nullptr, 12);
		writeGroupClose(out, 8);

		auto writeContour = [&](int graphId, const string& primName, int indent)
		{
			if (graphId < 0 || graphId >= o_contourGraphs.size()) return;

			zFnGraph fnGraph(o_contourGraphs[graphId]);
			if (fnGraph.numVertices() == 0) return;

			zIntArray vSequence;
			zItGraphVertexArray vArray;

			for (zItGraphVertex v(o_contourGraphs[graphId]); !v.end(); v++)
			{
				if (!v.checkValency(2)) vArray.push_back(v);
			}

			if (vArray.size() > 0) printf("\n contour[%i] - valence != 2 verts %i", graphId, vArray.size());

			if (vArray.size() == 2)
			{
				zItGraphHalfEdge he = vArray[0].getHalfEdge();
				vSequence.push_back(vArray[0].getId());
				int safetyCounter = 0;
				do
				{
					vSequence.push_back(he.getVertex().getId());
					he = he.getNext();
				} while (he.getVertex() != vArray[1] && safetyCounter++ < 10000);

				vSequence.push_back(vArray[1].getId());
				vSequence.push_back(vArray[0].getId());
			}

			if (vArray.size() == 0)
			{
				zItGraphHalfEdge he(o_contourGraphs[graphId], 0);
				zItGraphVertex startV = he.getStartVertex();
				zItGraphHalfEdge startHe = he;
				vSequence.push_back(he.getStartVertex().getId());

				do
				{
					vSequence.push_back(he.getVertex().getId());
					he = he.getNext();
				} while (he != startHe);
			}

			writeGraphPrim(out, o_contourGraphs[graphId], primName, nullptr, &vSequence, indent);
		};

		int end = o_sectionGraphs.size();
		auto writeGraphGroup = [&](const string& groupName, zObjGraphArray& graphs, const string& childPrefix, const zTransform* frameSource)
		{
			writeGroupOpen(out, groupName, 8);
			for (int i = 0; i < end; i++)
			{
				string id = core.getPaddedIndexString(i, 3);
				const zTransform* frame = frameSource && i < sectionFrames.size() ? &sectionFrames[i] : nullptr;
				if (i >= 0 && i < graphs.size()) writeGraphPrim(out, graphs[i], childPrefix + "_" + id, frame, nullptr, 12);
			}
			writeGroupClose(out, 8);
		};

		writeGraphGroup("sections", o_sectionGraphs, "section", sectionFrames.data());

		writeGroupOpen(out, "section_meshes", 8);
		for (int i = 0; i < end; i++)
		{
			string id = core.getPaddedIndexString(i, 3);
			const zTransform* frame = i < sectionFrames.size() ? &sectionFrames[i] : nullptr;
			if (i >= 0 && i < o_sectionMeshes.size()) writeMeshPrim(out, o_sectionMeshes[i], "sectionMesh_" + id, frame, 12);
		}
		writeGroupClose(out, 8);

		writeGraphGroup("trim_bracing", o_trimGraphs_bracing, "trim_bracing", nullptr);
		writeGraphGroup("bracing_slots", o_trimGraphs_bracing_slots, "bracing_slots", nullptr);
		writeGraphGroup("trim_features_hard", o_trimGraphs_features_hard, "trim_features_hard", nullptr);
		writeGraphGroup("trim_features_soft", o_trimGraphs_features_soft, "trim_features_soft", nullptr);
		writeGraphGroup("trim_interiorSplit", o_trimGraphs_SlotSide, "trim_interiorSplit", nullptr);
		writeGraphGroup("trim_seamAlignment", o_trimGraphs_seamAlignment, "trim_seamAlignment", nullptr);

		writeGroupOpen(out, "contours", 8);
		for (int i = 0; i < end; i++)
		{
			string id = core.getPaddedIndexString(i, 3);
			writeContour(i, "contours_" + id, 12);
		}
		writeGroupClose(out, 8);

		writeWorldClose(out);
		return out.good();
	}
	//EXPORT sub method


	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::check_PrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset, bool runBothPlanes, bool runPlaneLeft)
	{
		printHeightDomain = _printHeightDomain;
		neopreneOffset = zDomainFloat(0.0f, 0.0f);

		zStringArray files;
		core.getFilesFromDirectory(files, folderDir, zJSON);

		for (auto& s : files)
		{
			printf("\n All File: %s ", s.c_str());
		}

		string outFileName = folderDir + "printLayerHeights_allBlocks.csv";

		ofstream myfile;
		myfile.open(outFileName.c_str());

		if (myfile.fail())
		{
			cout << " error in opening outPath  " << outFileName.c_str() << endl;
			return;
		}

		printf("\n numFiles %i ", files.size());
		zBoolArray Block_visitied;
		Block_visitied.assign(files.size(), false);

		myfile << "blockID" << ","
			<< "frameCHECKS" << ","
			<< "minHeight" << ","
			<< "maxHeight" << ","
			<< "sdfCHECK" << ","
			<< "geometryCHECK" << ","
			<< "criticalPtsMin" << ","
			<< "criticalPtsMax" << ","
			<< "planeSpacing" << endl;


		for (auto& s : files)
		{
			//printf("\n File: %s ", s.c_str());

			zStringArray split_0 = core.splitString(s, ".");
			zStringArray split_1 = core.splitString(split_0[split_0.size() - 2], "_");

			printf("\n File: %s ", s.c_str());
			int _blockID = atoi(split_1[split_1.size() - 1].c_str());

			setFromJSON(folderDir, _blockID, runBothPlanes, runPlaneLeft);
			Block_visitied[_blockID] = true;

			bool frameCHECKS = false;
			bool geomCHECKS = true;
			bool sdfCHECKS = true;

			bool maxHeightCheck = true;
			bool minHeightCheck = true;

			float bestPlaneSpacing = 0.01f;
			zVector norm(0, 0, 1);
			vector<zItMeshHalfEdgeArray> vLoops;
			zObjMesh oMesh_top, oMesh_bottom;

			computeVLoops(o_SliceMesh_Left, medialIDS, FeaturedNumStrides, norm, vLoops, oMesh_top, oMesh_bottom);

			zScalarArray scalars;
			computeGeodesicScalars(o_SliceMesh_Left, vLoops, scalars, true);

			o_sectionMeshes.clear();
			computeGeodesicContours(vLoops, scalars, bestPlaneSpacing, oMesh_top, oMesh_bottom, o_sectionMeshes);;
			createSectionGraphs(o_sectionMeshes, o_sectionGraphs);

			zFnPointCloud fnCritical_min(criticalMinLayer_pts);
			zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

			float minHeight = FLT_MAX;
			float maxHeight = FLT_MIN;

			for (int m = 0; m < o_sectionGraphs.size(); m++)
			{
				zPointArray secPts;
				zFnGraph fnGraph(o_sectionGraphs[m]);
				fnGraph.getVertexPositions(secPts);

				zVectorArray pNorms;
				zObjGraph outPrintHeightLines;
				for (auto& p : secPts)
				{
					for (zItMeshFace f(o_sectionMeshes[m]); !f.end(); f++)
					{
						zPointArray fVerts;
						f.getVertexPositions(fVerts);
						bool ptIntTriangle = core.pointInTriangle(p, fVerts[0], fVerts[1], fVerts[2]);
						if (ptIntTriangle)
						{
							zVector n = f.getNormal();
							pNorms.push_back(n);
							break;
						}
					}
				}
				zFloatArray pHeights;
				getPrintHeight(secPts, pNorms, o_sectionMeshes[m], pHeights, outPrintHeightLines);

				for (int h = 0; h < pHeights.size(); h++)
				{
					if (pHeights[h] < printHeightDomain.min)fnCritical_min.addPosition(secPts[h]);
					if (pHeights[h] > printHeightDomain.max)fnCritical_max.addPosition(secPts[h]);

					if (pHeights[h] < minHeight) minHeight = pHeights[h];
					if (pHeights[h] > maxHeight) maxHeight = pHeights[h];
				}
			}
			actualPrintHeightDomain.min = minHeight;
			actualPrintHeightDomain.max = maxHeight;
			frameCHECKS = fnCritical_min.numVertices() == 0 && fnCritical_max.numVertices() == 0;



			printf("\n ----------- \n BlockID %i | %s | %1.4f %1.4f \n", _blockID, (frameCHECKS) ? "True" : "False", actualPrintHeightDomain.min, actualPrintHeightDomain.max);

			int minPts, maxPts;
			zFnPointCloud fnptCloudMin, fnptCloudMax;
			fnptCloudMin = zFnPointCloud(criticalMinLayer_pts);
			fnptCloudMax = zFnPointCloud(criticalMaxLayer_pts);

			minPts = fnptCloudMin.numVertices();
			maxPts = fnptCloudMax.numVertices();


			zFloatArray ptsMin, ptsMax;

			if (fnptCloudMin.numVertices() > 0)
			{
				zPointArray pts;
				fnptCloudMin.getVertexPositions(pts);
				for (auto& p : pts)
				{
					ptsMin.push_back(p.x);
					ptsMin.push_back(p.y);
					ptsMin.push_back(p.z);
				}
			}
			if (fnptCloudMax.numVertices() > 0)
			{
				zPointArray pts;
				fnptCloudMax.getVertexPositions(pts);
				for (auto& p : pts)
				{
					ptsMax.push_back(p.x);
					ptsMax.push_back(p.y);
					ptsMax.push_back(p.z);
				}
			}

			string path = folderDir + "blockMesh_" + to_string(_blockID) + ".json";
			json j;
			core.json_read(path, j);
			//output json for critical points
			//json j;
			string outDir = folderDir + "/criticalPointsCheck_nonPlanar";
			string outPath = outDir + "/outBlock_" + to_string(_blockID) + ".json";
			if (!filesystem::is_directory(outDir) || !filesystem::exists(outDir)) filesystem::create_directory(outDir);
			j["CP_Min"] = ptsMin;
			j["CP_Max"] = ptsMax;
			j["Layer_Height_Min"] = actualPrintHeightDomain.min;
			j["Layer_Height_Max"] = actualPrintHeightDomain.max;
			j["PlaneSpacing"] = bestPlaneSpacing;

			if (ptsMin.size() > 0 || ptsMax.size() > 0)
			{

				core.json_write(outPath, j);
			}


			myfile << _blockID << ","
				<< ((frameCHECKS) ? "True" : "False") << ","
				<< actualPrintHeightDomain.min * 1000 << ","
				<< actualPrintHeightDomain.max * 1000 << ","
				<< ((sdfCHECKS) ? "True" : "False") << ","
				<< ((geomCHECKS) ? "True" : "False") << ","
				<< minPts << ","
				<< maxPts << ","
				<< bestPlaneSpacing * 1000 << endl;


			string outDir2 = folderDir + "/outSections";
			if (!filesystem::is_directory(outDir2) || !filesystem::exists(outDir2)) filesystem::create_directory(outDir2);

			outDir2 += "/block_" + to_string(_blockID);
			if (!filesystem::is_directory(outDir2) || !filesystem::exists(outDir2)) filesystem::create_directory(outDir2);
			for (const auto& entry : std::filesystem::directory_iterator(outDir2)) std::filesystem::remove_all(entry.path());

			if (true)
			{
				zFnGraph fnGraph;
				for (int ss = 0; ss < o_sectionGraphs.size(); ss++)
				{
					string outPath2 = outDir2 + "/outSection_" + to_string(ss) + ".json";
					fnGraph = zFnGraph(o_sectionGraphs[ss]);
					fnGraph.to(outPath2, zJSON);
				}
			}
			printf("\n Finished block %i ", _blockID);


		}

		myfile.close();

		cout << " \n outPath exported : " << outFileName.c_str() << endl;

		for (int i = 0; i < Block_visitied.size(); i++)
		{
			if (!Block_visitied[i]) printf("\n %i ", i);
		}
	}




			ZSPACE_TOOLSETS_INLINE zPoint zTsCarbcomn::util_getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh)
	{
		float scaleVal = core.ofMap(threshold, thresholdLow, thresholdHigh, 0.0000f, 1.0000f);
		zVector e = vertex_higher - vertex_lower;
		double edgeLen = e.length();
		e.normalize();
		return (vertex_lower + (e * edgeLen * scaleVal));
	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints)
	{
		zFnGraph fnGraph(o_graph);
		zPoint* vPositions = fnGraph.getRawVertexPositions();
		contourPoints.clear();
		for (zItGraphEdge e(o_graph); !e.end(); e++)
		{
			zIntArray eVerts;
			e.getVertices(eVerts);
			float s0 = vertexScalars[eVerts[0]];
			float s1 = vertexScalars[eVerts[1]];
			bool contour = false;
			if (s0 <= threshold && s1 >= threshold)contour = true;
			if (s0 >= threshold && s1 <= threshold)contour = true;
			if (!contour) continue;
			zPoint v0 = vPositions[eVerts[0]];
			zPoint v1 = vPositions[eVerts[1]];
			zPoint pos1 = util_getContourPosition(threshold, v1, v0, s1, s0);
			contourPoints.push_back(pos1);
		}
	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints, float threshold)
	{
		outPoints.clear();
		zScalarArray vertexScalars;
		zPoint O(inPlane(3, 0), inPlane(3, 1), inPlane(3, 2));
		zVector N(inPlane(2, 0), inPlane(2, 1), inPlane(2, 2));
		for (zItGraphVertex v(o_graph); !v.end(); v++)
		{
			zPoint P = v.getPosition();
			float minDist_Plane = core.minDist_Point_Plane(P, O, N);
			vertexScalars.push_back(minDist_Plane);
		}
		zPointArray contourPoints;
		util_isoContour(o_graph, vertexScalars, threshold, contourPoints);
		if (closestPoint)
		{
			float dist = 1000000;
			zPoint cPoint;
			for (auto& p : contourPoints)
			{
				if (p.distanceTo(O) < dist)
				{
					dist = p.distanceTo(O);
					cPoint = p;
				}
			}
			outPoints.push_back(cPoint);
		}
		else outPoints = contourPoints;
	}





	//---- PROTECTED UTILITY METHODS


	ZSPACE_TOOLSETS_INLINE zItMeshHalfEdge zTsCarbcomn::util_getStartHalfEdge(zObjMesh& o_Mesh, int startVID, int endVID)
	{

		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		zVector dir = endPoint - startPoint;
		dir.normalize();

		zItMeshVertex vStart(o_Mesh, startVID);

		//compute start half edge
		zItMeshHalfEdge heStart;

		zItMeshHalfEdgeArray cHEdges;
		vStart.getConnectedHalfEdges(cHEdges);

		float val = 10000;
		for (zItMeshHalfEdge& he : cHEdges)
		{
			zVector heVec = he.getVector();
			heVec.normalize();

			if (1 - (heVec * dir) < val)
			{
				val = 1 - (heVec * dir);
				heStart = he;
			}
		}

		return heStart;
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_createGraphFromHEArray(zItGraphHalfEdgeArray& heArray, zObjGraph& outGraph)
	{
		zPointArray positions;
		zIntArray eConnects;


		for (int i = 0; i < heArray.size(); i++)
		{
			int indexStart = -1;
			int indexEnd = -1;
			zPoint ptSt = heArray[i].getStartVertex().getPosition();
			zPoint ptEnd = heArray[i].getVertex().getPosition();
			bool existSt = core.checkRepeatVector(ptSt, positions, indexStart);
			bool existEnd = core.checkRepeatVector(ptEnd, positions, indexEnd);

			if (!existSt)
			{
				positions.push_back(ptSt);
				indexStart = positions.size() - 1;
			}
			if (!existEnd)
			{
				positions.push_back(ptEnd);
				indexEnd = positions.size() - 1;
			}
			eConnects.push_back(indexStart);
			eConnects.push_back(indexEnd);
		}
		//for (int i = 0; i < heArray.size(); i++)
		//{
		//	positions.push_back(heArray[i].getStartVertex().getPosition());

		//	if (positions.size() > 1)
		//	{
		//		eConnects.push_back(positions.size() - 2);
		//		eConnects.push_back(positions.size() - 1);
		//	}

		//	if (i == heArray.size() - 1)
		//	{
		//		int index = -1;
		//		zPoint pt = heArray[i].getVertex().getPosition();
		//		bool exist = core.checkRepeatVector(pt, positions, index);
		//		if (!exist)
		//		{
		//			positions.push_back(heArray[i].getVertex().getPosition());

		//			eConnects.push_back(positions.size() - 2);
		//			eConnects.push_back(positions.size() - 1);
		//		}
		//		else
		//		{
		//			eConnects.push_back(positions.size() - 1);
		//			eConnects.push_back(0);
		//		}
		//	}
		//}

		zFnGraph fnGraph(outGraph);
		fnGraph.create(positions, eConnects);
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_getPerpendicularVector(zPlane& plane, zVector edgeVector, zPoint midPoint, float graphLength, zObjGraph& outGraph)
	{
		zVector planeNormal(plane(2, 0), plane(2, 1), plane(2, 2));
		zVector vector = planeNormal ^ edgeVector;
		vector.normalize();

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(midPoint - (vector * (graphLength/2)));
		gPts.push_back(midPoint + (vector * (graphLength/2)));
		gEdges.push_back(0);
		gEdges.push_back(1);

		zObjGraph oGraph;
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);
	}

	ZSPACE_TOOLSETS_INLINE zVector zTsCarbcomn::util_averageVectorsAtGraphVertex(zItGraphVertex& v)
	{
		zItGraphHalfEdgeArray hes;
		v.getConnectedHalfEdges(hes);

		zVector result(0, 0, 0);
		for (zItGraphHalfEdge he : hes)
		{
			zVector vec = he.getVector();
			vec.normalize();
			result += vec;
		}
		result /= hes.size();
		result.normalize();
		return result;
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_combineMultipleGraphs(zObjGraphArray& inGraphs, zObjGraph& outGraph)
	{
		zPointArray positions;
		zIntArray eConnects;
		zColorArray eColors;

		zFnGraph fng;
		for (zObjGraph& g : inGraphs)
		{
			zPointArray pts;
			zIntArray es;
			zColorArray colors;
			fng = zFnGraph(g);
			fng.getVertexPositions(pts);
			fng.getEdgeData(es);
			fng.getEdgeColors(colors);

			// store the mapping of old indices to new indices for the current graph
			zIntArray indexMapping(pts.size());

			// Add points to the combined positions array and update the map
			for (int i = 0; i < pts.size(); ++i)
			{
				const zPoint& p = pts[i];
				int tempIndex = -1;
				bool exist = core.checkRepeatVector(pts[i], positions, tempIndex);
				if (!exist) positions.push_back(p);
				indexMapping[i] = exist? tempIndex : positions.size() - 1;
			}

			// Adjust edge connectivity based on new vertex indices
			for (int i = 0; i < es.size(); i += 2)
			{
				int newStartIdx = indexMapping[es[i]];
				int newEndIdx = indexMapping[es[i + 1]];
				eConnects.push_back(newStartIdx);
				eConnects.push_back(newEndIdx);
			}

			// Combine edge colors
			eColors.insert(eColors.end(), colors.begin(), colors.end());
		}

		// Create the combined graph
		fng = zFnGraph(outGraph);
		fng.create(positions, eConnects);
		fng.setEdgeColors(eColors, false);
	}



	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_computeSlotGraph(zPlane plane, zObjGraph& inPoly, float graphLength, bool iterate, zObjGraph& outGraph)
	{
		zPrintParamSDF _printParameters;

		//poly is sectionGraph -> we don't have to specify right/left
		zFnGraph inFnGraph(inPoly);

		zVector Y(0, 1, 0);

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = _col_out_corner;
		zColor endColor = _col_in_corner;

		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		/*float ptOffset = edgeLength / 2.0 ;
		if (iterate) ptOffset += graphLength;*/

		float ptOffset = isCableBlock ? 0.3f : _printParameters.slotStart;
		if(blockId == 52)
			ptOffset = 0.125f;
		float iterOffset = iterate ? -(_printParameters.slotIterating * 0.5f) : _printParameters.slotIterating * 0.5f;

		//startV += (edgeVector * ptOffset);
		//move slot lower to avoid too close to corner pt

			startV += (edgeVector * ((edgeLength * ptOffset) + iterOffset));
			startV += (edgeVector * (_printParameters.slotIterating * 0.5f));

		util_getPerpendicularVector(plane, edgeVector, startV, graphLength*1.5, outGraph);
	}

			ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_computeSplitGraph_xy(zObjGraph& inPoly, zObjGraph& outGraph)
	{
		zPrintParamSDF _printParameters;

		zFnGraph fnGraph(inPoly);

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zVector* inPositions = fnGraph.getRawVertexPositions();
		zColor* inColors = fnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = _col_in_corner;
		zColor endColor = _col_out_corner;

		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		//Shorten the edge
		startV += (edgeVector * _printParameters.splitTrimTimming);
		endV -= (edgeVector * _printParameters.splitTrimTimming);

		//move the edge so it is in the middle between 1st and 2nd offset
		//The vector of the offset has the following:
		////1. rotate along the normal of the plane by 90 degrees and -90 degrees (vec1 and vec2)
		////2. get a vector between the center of the edge with the center of the section graph (vecC)
		////3. find the angle between the vec1/vec2 in step1 and vecC, and choose the one with the smallest angle (vec)
		////4. move using the vector in step 3

		zVector n = zVector(planeXY(2, 0), planeXY(2, 1), planeXY(2, 2));
		zVector vec1 = edgeVector.rotateAboutAxis(n, 90.0f);
		zVector vec2 = edgeVector.rotateAboutAxis(n, -90.0f);

		zVector vecC = fnGraph.getCenter() - ((startV+endV)/2.0);
		vecC.normalize();

		float angle1 = abs( vec1.angle(vecC));
		float angle2 = abs(vec2.angle(vecC));

		zVector vec = (angle1 < angle2) ? vec1 : vec2;
		vec.normalize();
		//changed for pentagon split
		vec *= _printParameters.splitTrimOffset;
		//vec *= _printParameters.targetInteriorGap;
		startV += vec;
		endV += vec;


		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV);
		gPts.push_back(endV);
		gEdges.push_back(0);
		gEdges.push_back(1);
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);
	}


	ZSPACE_TOOLSETS_INLINE bool zTsCarbcomn::util_getShortestHEsBetweenColors(zObjGraph& graph, zColor startColor, zColor endColor, zItGraphHalfEdgeArray& outHEs)
	{
		outHEs.clear();
		vector<zItGraphHalfEdgeArray> tempHEsArray;
		zFloatArray lengths;
		zFnGraph fng(graph);


		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				//walk until you reach a vertex with color _colorFeatureOuter
				zItGraphHalfEdgeArray hes;
				v.getConnectedHalfEdges(hes);
				//traverse both way
				for (auto& he : hes)
				{
					float length = 0;

					zItGraphHalfEdge heStart = he;
					zItGraphHalfEdge he = heStart;
					zItGraphHalfEdgeArray innerHE;
					int safetyCounter = 0;
					while (safetyCounter < fng.numEdges())
					{
						innerHE.push_back(he);
						length += he.getLength();
						if (he.getVertex().getColor() == endColor) break;
						//if (he.getVertex().getColor() == startColor)
						//{
						//	innerHE.clear();
						//	length = 0;
						//}
						he = he.getNext();
						safetyCounter++;
					}

					//if (safetyCounter >= fng.numEdges())
					//{
					//	printf("\n heColor was not found");
					//}

					tempHEsArray.push_back(innerHE);
					lengths.push_back(length);
					//printf("\n slotGraph_Arch index-size %i | %i", innerHEs.size(), innerHE.size());
					//innerLengths.push_back(length);
				}


			}
		}
		if (tempHEsArray.size() == 0)
		{
			cout << startColor.r << " " << startColor.g << " " << startColor.b << endl;
			cout << endColor.r << " " << endColor.g << " " << endColor.b << endl;

			printf("\n getShotestHEsBetweenColors no tempHEsArray found! RETURN");
			return false;
		}

		int index = 0;
		int minCount = INT_MAX;
		float minLength = FLT_MAX;
		for (int i = 0; i < tempHEsArray.size(); i++)
		{
			//if (tempHEsArray[i].size() < minCount)
			if (lengths[i] < minLength)
			{
				minCount = tempHEsArray[i].size();
				minLength = lengths[i];
				index = i;
			}
		}
		outHEs = tempHEsArray[index];
		return true;
	}


	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::util_getHEsColorLen(zObjGraph& graph, zColor& startCol, zColor& endCol, float len, zItGraphHalfEdgeArray& out)
	{
		out.clear();
		zItGraphVertex start_vert(graph);

		// Find the start vertex
		for (; !start_vert.end(); start_vert++)
		{
			if (start_vert.getColor() == startCol)
				break;
		}

		if (start_vert.end())
		{
			//Could not find the start color
			printf("\n ERROR! Could not find start color in graph.");
			return;
		}

		zItGraphHalfEdgeArray hes, copy;
		start_vert.getConnectedHalfEdges(hes);
		copy = hes;


		std::array<zItGraphHalfEdgeArray, 2> half_edges;
		std::array<float, 2> lengths;

		// Gather halfedges
		for (int i = 0; i < 2; ++i)
		{
			for (int counter = 0; counter < graph.graph.n_e; ++counter)
			{
				lengths[i] += hes[i].getLength();
				half_edges[i].push_back(hes[i]);

				if (hes[i].getVertex().getColor() == endCol)
					break;

				hes[i] = hes[i].getNext();
			}
		}

		//Figure out correct direction
		int index = lengths[0] < lengths[1] ? 0 : 1;

		lengths[index] = 0.0f;
		while (lengths[index] < len)
		{
			lengths[index] += copy[index].getLength();
			out.push_back(copy[index]);

			copy[index] = copy[index].getNext();
		}
	}

	ZSPACE_TOOLSETS_INLINE int zTsCarbcomn::util_getGraphClosestPoint(zObjGraph& graph, zPoint& samplePoint, zPoint& outPoint, float& dist)
	{
		int index = -1;
		double minDist = DBL_MAX;
		zPoint cP;
		for (zItGraphEdge e(graph); !e.end(); e++)
		{
			zPointArray pts;
			e.getVertexPositions(pts);
			zPoint p;
			double d = core.minDist_Edge_Point(samplePoint, pts[0], pts[1], p);
			if (d < minDist)
			{
				minDist = d;
				index = e.getId();
				cP = p;
			}
		}
		dist = (float)minDist;
		outPoint = cP;
		return index;
	}

	ZSPACE_TOOLSETS_INLINE int zTsCarbcomn::util_getHeArrayClosestPoint(zItGraphHalfEdgeArray& hes, zPoint& samplePoint, zPoint& outPoint, float& dist)
	{
		int index = -1;
		double minDist = DBL_MAX;
		zPoint cP;
		for (int i = 0; i < hes.size(); i++)
		{

			zPoint tempCp;
			zPoint p0 = hes[i].getStartVertex().getPosition();
			zPoint p1 = hes[i].getVertex().getPosition();
			double d = core.minDist_Edge_Point(samplePoint, p0,p1, tempCp);
			if (d < minDist)
			{
				minDist = d;
				index = i;
				cP = tempCp;
			}
		}

		dist = (float)minDist;
		outPoint = cP;
		return index;
	}

		ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getScalars_3dp_wall_bracing(zObjGraph& sectionGraph, zObjGraph& bracingGraph, zObjGraph& bracing_slotsGraph, float iterateOffset, bool iterateChk, zScalarArray & outScalar_interiorBracing, zScalarArray & outScalar_bracing, zScalarArray & outScalar_bracingSlots)
	{
		/// This method create the bracing
		/// 1. get the scalar of the bracing graph
		/// 2. get the scalar of the bracing slots
		/// 3. subtract the bracing slots from the bracing graph
		/// 4. get the scalar of the interior bracing

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(sectionGraph);
		zFnMeshScalarField fnField(o_field);






		fnField.getScalarsAsEdgeDistance(outScalar_bracing, bracingGraph, _printParameters.bracingEdgeWidth, false);
		zFnGraph fnTrimG(bracingGraph);
		int bracingCount = fnTrimG.numEdges();

		zItGraphHalfEdgeArray hesTemp;
		fnTrimG.numHalfEdges();
		for (int i = 0; i < fnTrimG.numHalfEdges(); i++)
		{
			zItGraphHalfEdge he(bracingGraph, i);
			if (he.getStartVertex().getColor() == _col_in_feature)
			{
				hesTemp.push_back(he);
			}
		}

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zObjGraph o_bracingSlots;
		zObjGraphArray bracingSlotsArray;
		bracingSlotsArray.assign(hesTemp.size(), zObjGraph());
		int counter = 0;
		for (zItGraphHalfEdge& he : hesTemp)
		{
			//to get the slot offset, we have the following steps
			//1. get the length of the bracing edge
			//2. get the offset of the triangle point (using the factor value)
			//3. get the offset of the 1st offset and the 2nd offset
			//4. find the middle point between step 2 and step 3
			//5. iterate
			float triangleStart = he.getLength() * _printParameters.wall_triangleOffsetFactor;
			float exteriorStart = _printParameters.offset_1st_exterior + _printParameters.offset_2nd_exterior;
			float slotOffset = triangleStart + ((he.getLength() - exteriorStart - triangleStart) / 2);

			if (iterateChk) slotOffset -= _printParameters.slotIterating_in;
			zVector vec = he.getVector();
			vec.normalize();
			vec *= slotOffset;
			zPoint midP = he.getStartVertex().getPosition() + vec;
			float graphLength = _printParameters.bracingEdgeWidth * 4;
			util_getPerpendicularVector(planeXY, vec, midP, graphLength, bracingSlotsArray[counter]);

			counter++;

		}
		util_combineMultipleGraphs(bracingSlotsArray, o_bracingSlots);

		fnField.getScalarsAsEdgeDistance(outScalar_bracingSlots, o_bracingSlots, _printParameters.bracingEdgeSlotWidth, false);
		fnField.boolean_subtract(outScalar_bracing, outScalar_bracingSlots, outScalar_interiorBracing, false);

		bracing_slotsGraph = o_bracingSlots;
		o_debug_bracingslotsgraph = o_bracingSlots;
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getScalars_3dp_wall_triangles(zObjGraph& sectionGraph, zScalarArray& outScalar_triangles, bool remove_firstLast)
	{

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(sectionGraph);
		zFnMeshScalarField fnField(o_field);

		//inner vertex are the vertex of the sectionGraph on the inner side (_col_in_feature - magenta)
		//and outer vertex are the vertex of the sectionGraph on the outer side (_col_out_feature - orange)
		zItGraphVertexArray innerVertx, outerVertx;
		util_innerOuter(sectionGraph, true, innerVertx, outerVertx);

		if (remove_firstLast)
		{
			innerVertx.erase(innerVertx.begin() + (innerVertx.size() - 2));
			innerVertx.erase(innerVertx.begin() + 1);

			outerVertx.erase(outerVertx.begin() + (outerVertx.size() - 2));
			outerVertx.erase(outerVertx.begin() + 1);

		}

		zPointArray gPts;
		zIntArray eConnect;
		int startId = (remove_firstLast) ? /*1*/ 0 : 0;
		int endId = (remove_firstLast) ? innerVertx.size() /*- 1*/ : innerVertx.size();
		printf("\n remove_firstLast in triangle %i", remove_firstLast);

		for (int i = startId; i < endId; i++)
		{
			//this is the middle point between the inner and outer vertex - this value can be changed to any value
			zVector edgeDir = outerVertx[i].getPosition() - innerVertx[i].getPosition();
			edgeDir *= _printParameters.wall_triangleOffsetFactor;
			zPoint Pt0 = innerVertx[i].getPosition() + edgeDir;

			gPts.push_back(Pt0);
			if (i < endId - 1)
			{
				zPoint p1 = (innerVertx[i].getPosition() + innerVertx[i + 1].getPosition()) / 2;
				zPoint p2 = (outerVertx[i].getPosition() + outerVertx[i + 1].getPosition()) / 2;
				zVector v = p2 - p1;
				v.normalize();
				zPoint Pt1 = p1 + (v * 0.001); //< this is the peak of the triangle (average) , small offset so the full polygon is not through the same top edge
				gPts.push_back(Pt1);
			}

		}
		//create the triangles sides
		for (int i = 0; i < gPts.size() - 1; i++)
		{
			eConnect.push_back(i);
			eConnect.push_back(i + 1);
		}
		//create the top part of the triangle - this is the full of the inner edge
		for (int i = endId - 1; i >= startId; i--)
		{
			eConnect.push_back(gPts.size() - 1);
			//printf("\n [%i] innerVertx %i  ", i, innerVertx.size() - 1);
			gPts.push_back(innerVertx[i].getPosition());
			eConnect.push_back(gPts.size() - 1);
		}
		eConnect.push_back(gPts.size() - 1);
		eConnect.push_back(0);
		zObjGraph o_triangles;
		zFnGraph fnG(o_triangles);
		fnG.create(gPts, eConnect);


		fnField.getScalars_Polygon(outScalar_triangles, o_triangles, false);

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getScalars_offset(zObjGraph& sectionGraph, int numSmooth, zScalarArray& outScalar_polygon, zScalarArray& outScalar_offset_outer, zScalarArray & outScalar_offset_inner)
	{
		zPrintParamSDF _printParameters;
		zFnMeshScalarField fnField(o_field);
		zFnGraph fnGraph(sectionGraph);


		// Profile polygon field
		zIntArray edgeId;
		fnField.getScalars_Polygon(outScalar_polygon, sectionGraph, edgeId, false);

		//create a map of edge offset based on the edgeId
		zFloatArray outerOffsetArray, innerOffsetArray;

		outerOffsetArray.assign(fnGraph.numEdges(), _printParameters.offset_1st_exterior);//20
		//innerOffsetArray.assign(fnGraph.numEdges(), _printParameters.offset_1st_exterior + _printParameters.offset_2nd_exterior);
		innerOffsetArray.assign(fnGraph.numEdges(), _printParameters.offset_1st_exterior + _printParameters.offset_2nd_exterior);

		//color the section based on the offset color
		zItGraphHalfEdgeArray hesInterior;
		util_getShortestHEsBetweenColors(sectionGraph, _col_in_corner_st, _col_in_corner, hesInterior);

		// Get specific half-edges of the section for wall blocks.
		zItGraphHalfEdgeArray cyan_red, red_green, green_yellow, yellow_orange, cyan_orange;
		util_getShortestHEsBetweenColors(sectionGraph, zCYAN, zRED, cyan_red);
		util_getShortestHEsBetweenColors(sectionGraph, zRED, zGREEN, red_green);
		util_getShortestHEsBetweenColors(sectionGraph, zGREEN, zYELLOW, green_yellow);
		util_getShortestHEsBetweenColors(sectionGraph, zYELLOW, zORANGE, yellow_orange);
		util_getShortestHEsBetweenColors(sectionGraph, zCYAN, zORANGE, cyan_orange);

		// Default colour
		fnGraph.setEdgeColor(zBLUE, false);

		for (zItGraphHalfEdge he : cyan_red)
		{
			outerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior;
			innerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior;
			he.getEdge().setColor(zMAGENTA);
		}

		for (zItGraphHalfEdge he : red_green)
		{
			outerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior;
			innerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior;
			he.getEdge().setColor(zMAGENTA);
		}

		for (zItGraphHalfEdge he : green_yellow)
		{
			outerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior;
			innerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior;
			he.getEdge().setColor(zMAGENTA);
		}

		outScalar_offset_outer = outScalar_polygon;
		outScalar_offset_inner = outScalar_polygon;

		//update the scalar offset based on the edgeId
		for (int sf = 0; sf < outScalar_offset_outer.size(); sf++)
		{
			outScalar_offset_outer[sf] += outerOffsetArray[edgeId[sf]];
			outScalar_offset_inner[sf] += outerOffsetArray[edgeId[sf]] + innerOffsetArray[edgeId[sf]];
		}

		if (!isCableBlock)
		{
			//smooth fields
			fnField.smoothField(outScalar_offset_inner, numSmooth);
			fnField.smoothField(outScalar_offset_outer, numSmooth);
		}

	}

	ZSPACE_TOOLSETS_INLINE int zTsCarbcomn:: util_get_corrected_id(std::unordered_map<int, int>& map, int id_to_check, int id_to_set)
		{
			if (map.count(id_to_check))
				return map[id_to_check];
			else
			{
				map[id_to_check] = id_to_set;
				return -1;
			}
		};

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn:: util_merge_graph(zObjGraph& oGraph, double tol)
		{
			zIntArray new_Connects, eConnects;
			zPointArray new_Positions, vPositions;

			std::unordered_map<int, int> vertex_map;

			zFnGraph fnGraph(oGraph);
			fnGraph.getVertexPositions(vPositions);
			fnGraph.getEdgeData(eConnects);

			// Graph Edge, Number of v1 vertices
			std::vector<std::tuple<zItGraphEdge, int>> v1_edges;

			zItGraphVertexArray f_edge_verts, s_edge_verts;

			for (zItGraphEdge e(oGraph); !e.end(); e++)
			{
				f_edge_verts.clear();
				e.getVertices(f_edge_verts);

				int v1_num = 0;
				for (auto& v : f_edge_verts)
				{
					if (v.checkValency(1))
						++v1_num;
				}

				// Has v1 vertex
				if (v1_num)
				{
					v1_edges.push_back(std::make_tuple(e, v1_num));
				}
				else
				{
					// Add existing edge
					for (auto& v : f_edge_verts)
					{
						int id = util_get_corrected_id(vertex_map, v.getId(), new_Positions.size());

						if (id == -1)
						{
							new_Positions.push_back(v.getPosition());
							new_Connects.push_back(new_Positions.size() - 1);
						}
						else
						{
							new_Connects.push_back(id);
						}
					}
				}
			}

			zPoint f_pos, s_pos, n_pos;

			// Generate the new vertex positions and fill in the vertex_map;
			for (auto& [f_edge, f_valen] : v1_edges)
			{
				// Still has v1 vertices
				if (f_valen)
				{
					f_edge_verts.clear();
					f_edge.getVertices(f_edge_verts);
					// Find another v1 edge within tolerance
					for (auto& [s_edge, s_valen] : v1_edges)
					{
						if (s_valen && (f_edge != s_edge))
						{
							s_edge_verts.clear();
							s_edge.getVertices(s_edge_verts);

							for (auto& f_vert : f_edge_verts)
							{
								if (f_vert.checkValency(1))
								{
									f_pos = f_vert.getPosition();

									for (auto& s_vert : s_edge_verts)
									{
										s_pos = s_vert.getPosition();

										if (f_pos.distanceTo(s_pos) < tol)
										{
											n_pos = (f_pos + s_pos) * 0.5f;
											new_Positions.push_back(n_pos);

											// Both of them now point to the same id;
											util_get_corrected_id(vertex_map, f_vert.getId(), new_Positions.size() - 1);
											util_get_corrected_id(vertex_map, s_vert.getId(), new_Positions.size() - 1);

											--f_valen;
											--s_valen;

											break;
										}
									}
								}
							}
						}
					}
				}
			}

			// Fill in the connectivity based on the vertex_map
			for (auto& [edge, _] : v1_edges)
			{
				f_edge_verts.clear();
				edge.getVertices(f_edge_verts);

				for (auto& v : f_edge_verts)
				{
					int id = util_get_corrected_id(vertex_map, v.getId(), new_Positions.size());

					if (id == -1)
					{
						new_Positions.push_back(v.getPosition());
						new_Connects.push_back(new_Positions.size() - 1);
					}
					else
					{
						new_Connects.push_back(id);
					}
				}
			}

			int num_merged = vPositions.size() - new_Positions.size();
			if (num_merged)
			{
				printf("\033[0;32m");
				printf("\n Merged %d vert%s\n", num_merged, num_merged == 1 ? "ex" : "ices");
				printf("\033[0m");
			}

			fnGraph.create(new_Positions, new_Connects);
		};

		//-------------------------- Non planar UTILS ------------------------------

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getPokeMesh(zObjMesh& o_mesh, zObjMesh& o_TriMesh)
	{
		zFnMesh fnMesh(o_mesh);

		zColorArray vColors;
		fnMesh.getVertexColors(vColors);

		zPointArray vPositions;
		fnMesh.getVertexPositions(vPositions);

		zPointArray fCens;
		fnMesh.getCenters(zHEData::zFaceData, fCens);

		zPointArray positions;
		zIntArray pCounts, pConnects;

		positions.insert(positions.end(), vPositions.begin(), vPositions.end());
		positions.insert(positions.end(), fCens.begin(), fCens.end());

		int numOriginalVerts = vPositions.size();

		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zIntArray fVerts;
			f.getVertices(fVerts);

			int fID = f.getId();

			for (int i = 0; i < fVerts.size(); i++)
			{
				int nextID = (i + 1) % fVerts.size();

				pConnects.push_back(fVerts[i]);
				pConnects.push_back(fVerts[nextID]);
				pConnects.push_back(numOriginalVerts + fID);

				pCounts.push_back(3);
			}

			vColors.push_back(zBLACK);
		}

		zFnMesh fnPokeMesh(o_TriMesh);
		fnPokeMesh.create(positions, pCounts, pConnects);

		fnPokeMesh.setVertexColors(vColors);
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getLoop(zItMeshHalfEdge& heStart, bool forward, bool corner, int vCounter, vector<zItMeshHalfEdgeArray>& v_Loops)
	{
		zItMeshHalfEdge he_U = (forward) ? heStart.getNext() : heStart.getPrev();
		if (corner) he_U = heStart;

		bool exit_1 = false;

		zItMeshHalfEdge he_V = (forward) ? he_U.getSym().getNext() : he_U.getSym().getPrev();

		zItMeshHalfEdgeArray tempV;

		bool exit_2 = false;

		for (int i = 0; i < vCounter; i++)
		{
			if (forward) tempV.push_back(he_V.getSym());
			else tempV.push_back(he_V);

			//he_V.getEdge().setColor(zBLUE);


			if (!exit_2) he_V = (forward) ? he_V.getNext().getSym().getNext() : he_V.getPrev().getSym().getPrev();
		}


		v_Loops.push_back(tempV);

		//he_U.getEdge().setColor(zMAGENTA);



	}


	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zPointArray& fVerts, zColorArray& fVColors)
	{
		fVerts.clear();
		fVColors.clear();

		zItMeshHalfEdge he = heStart;

		do
		{
			if (forward)
			{
				fVerts.push_back(he.getVertex().getPosition());
				fVColors.push_back(he.getVertex().getColor());
				he = he.getNext();
			}
			else
			{
				fVerts.push_back(he.getStartVertex().getPosition());
				fVColors.push_back(he.getStartVertex().getColor());
				he = he.getPrev();
			}

		} while (he != heStart);

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zIntArray& fVerts)
	{
		fVerts.clear();

		zItMeshHalfEdge he = heStart;

		do
		{
			if (forward)
			{
				fVerts.push_back(he.getVertex().getId());
				he = he.getNext();
			}
			else
			{
				fVerts.push_back(he.getStartVertex().getId());
				he = he.getPrev();
			}

		} while (he != heStart);

	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::createBoundaryEdgeGraph(zObjMesh& o_mesh, bool closeGraph, zObjGraph& o_Graph)
	{
		zPointArray positions;
		zIntArray eConnects;
		zColorArray vColors;

		zItMeshHalfEdge he;

		//Find boundary edge
		for (zItMeshHalfEdge tmpHE(o_mesh); !tmpHE.end(); tmpHE++)
		{
			if (tmpHE.onBoundary())
			{
				he = tmpHE;
				break;
			}
		}

		zItMeshHalfEdge startHE = he;
		positions.push_back(he.getStartVertex().getPosition());
		vColors.push_back(he.getStartVertex().getColor());

		do
		{
			positions.push_back(he.getVertex().getPosition());
			vColors.push_back(he.getVertex().getColor());

			eConnects.push_back(positions.size() - 2);
			eConnects.push_back(positions.size() - 1);

			he = he.getNext();

		} while (he != startHE);

		if (closeGraph)
		{
			eConnects.push_back(positions.size() - 1);
			eConnects.push_back(0);
		}

		zFnGraph fnGraph(o_Graph);
		fnGraph.create(positions, eConnects);

		fnGraph.setVertexColors(vColors);
		;
	}

		ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::setPtGraph(zObjGraph& o_Graph, zPoint& refPt, bool setX, bool setY, bool setZ)
	{
		zFnGraph fnGraph(o_Graph);
		zPoint* positions = fnGraph.getRawVertexPositions();

		for (int i = 0; i < fnGraph.numVertices(); i++)
		{
			if (setX) positions[i].x = refPt.x;
			if (setY) positions[i].y = refPt.y;
			if (setZ) positions[i].z = refPt.z;
		}
	}

				ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getPrintHeight(zPointArray& pPoints, zVectorArray& pNorms, zObjMesh& o_Mesh, zFloatArray& pHeights, zObjGraph& outPrintHeightLines)
	{
		zFnGraph fnGraph(outPrintHeightLines);
		zIntArray eConnects;
		zPointArray gPositions;
		zFnMesh fnMesh(o_Mesh);
		printf("\n pPoint count %i pNorms %i numPolygons %i", pPoints.size(), pNorms.size(), fnMesh.numPolygons());
		for (int i = 0; i < pPoints.size(); i++)
		{
			float d = 10000;

			zPoint closestPt;
			bool check = false;
			for (zItMeshFace f(o_Mesh); !f.end(); f++)
			{
				zPointArray fVPositions;
				f.getVertexPositions(fVPositions);
				zPoint cP;
				//cout << endl << fVPositions[0] << "| " << fVPositions[1] << " | " << fVPositions[2] << " || " << pNorms[i] << " || " << pPoints[i];
				bool chk = core.ray_triangleIntersection(fVPositions[0], fVPositions[1], fVPositions[2], pNorms[i], pPoints[i], cP);

				if (chk)
				{
					if (cP.distanceTo(pPoints[i]) < d)
					{
						d = cP.distanceTo(pPoints[i]);
						closestPt = cP;
						check = true;
					}
				}

			}

			if (check)
			{
				pHeights.push_back(d);
				gPositions.push_back(pPoints[i]);
				gPositions.push_back(closestPt);
				eConnects.push_back(gPositions.size() - 2);
				eConnects.push_back(gPositions.size() - 1);
			}
			if (!check)
			{
				printf("\n no intersection found for point %i ", i);
			}

		}
		fnGraph.create(gPositions, eConnects);
		printf("\n printHeightGraph v %i - e %i ", fnGraph.numVertices(), fnGraph.numEdges());
	}

			ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getBaryCentricCoordinates_triangle(zPoint& pt, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& baryCoordinates)
	{
		zVector v0 = t1 - t0;
		zVector v1 = t2 - t0;
		zVector v2 = pt - t0;

		float d00 = v0 * (v0);
		float d01 = v0 * (v1);
		float d11 = v1 * (v1);
		float d20 = v2 * (v0);
		float d21 = v2 * (v1);

		float denom = d00 * d11 - d01 * d01;

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0 - v - w;

		baryCoordinates = zPoint(u, v, w);

	}
	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::getProjectionPoint_triangle(zPoint& baryCoordinates, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& projectionPt)
	{
		projectionPt = t0 * baryCoordinates.x + t1 * baryCoordinates.y + t2 * baryCoordinates.z;
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::barycentericProjection_triMesh(zObjGraph& o_graph, zObjMesh& o_inMesh, zObjMesh& o_projectionMesh, zVectorArray& outNotmals)
	{
		zFnGraph fnGraph(o_graph);

		zPointArray positions;
		fnGraph.getVertexPositions(positions);
		outNotmals.clear();
		bool check = false;
		zFnMesh fnMesh(o_inMesh);

		for (auto& p : positions)
		{
			for (zItMeshFace f(o_inMesh); !f.end(); f++)
			{

				zPointArray fVerts;
				f.getVertexPositions(fVerts);
				bool ptIntTriangle = core.pointInTriangle(p, fVerts[0], fVerts[1], fVerts[2]);
				if (ptIntTriangle)
				{
					zPoint baryCoordinates;
					getBaryCentricCoordinates_triangle(p, fVerts[0], fVerts[1], fVerts[2], baryCoordinates);


					zItMeshFace fProjection(o_projectionMesh, f.getId());

					zPointArray fVerts_projection;
					fProjection.getVertexPositions(fVerts_projection);

					zPoint projectionPt;
					getProjectionPoint_triangle(baryCoordinates, fVerts_projection[0], fVerts_projection[1], fVerts_projection[2], projectionPt);

					p = projectionPt;
					//zItMeshFace fProjection(o_projectionMesh, f.getId());
					zVector n = fProjection.getNormal();
					outNotmals.push_back(n);
					check = true;
					break;
				}
			}
			if (!check)
			{

				printf("\n barycentericProjection_triMesh check %i", check);
			}
		}

		fnGraph.setVertexPositions(positions);
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::unrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zIntPairArray& bsf_vertexPairs, zTransform& outTransformStart)
	{
		zFnMesh fnMesh(o_mesh);
		zPoint* vPositions = fnMesh.getRawVertexPositions();

		zFnMesh fnMesh_unroll(o_mesh_unroll);
		zPoint* vPositions_unroll = fnMesh_unroll.getRawVertexPositions();
		//zTransform outTransformStart;
		/*for (int i = 0; i < oriVertex_UnrollVertex_map.size(); i++)
		{
			zPoint p = vPositions[i];

			for (auto vID : oriVertex_UnrollVertex_map[i])
			{
				vPositions_unroll[vID] = p;
			}
		}*/

		// unroll
		// https://computergraphics.stackexchange.com/questions/8774/unfold-a-3d-mesh-to-a-2d-plane

		for (int i = 0; i < bsf_vertexPairs.size(); i++)
		{
			zItMeshFace f1(o_mesh, bsf_vertexPairs[i].first);
			zItMeshFace f2(o_mesh, bsf_vertexPairs[i].second);

			zItMeshFace f1_unroll(o_mesh_unroll, bsf_vertexPairs[i].first);
			zItMeshFace f2_unroll(o_mesh_unroll, bsf_vertexPairs[i].second);

			zIntPair hePair = getCommonEdge(f1, f2);

			zItMeshHalfEdge he_1(o_mesh, hePair.first);
			zItMeshHalfEdge he_2(o_mesh, hePair.second);

			zPoint A = vPositions[he_2.getStartVertex().getId()];
			zPoint B = vPositions[he_2.getVertex().getId()];

			// unroll first face
			if (i == 0)
			{
				zItMeshHalfEdge he_walker_1 = he_1;
				int f1_numV = f1.getNumVertices();

				float l_ab = he_1.getLength();

				//set the transform
				zVector y = B - A;
				zVector x = zVector(0, 0, 1) ^ y;
				zVector z = zVector(0, 0, 1);
				outTransformStart = core.getTransformFromVectors(A, x, y, z);

				zPoint a(0, 0, 0);
				zPoint b(a.x, a.y + l_ab, a.z);

				// update  postions of corresponding a & b in unroll mesh
				zIntPair hashKey_a(f1.getId(), he_1.getStartVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_a = oriFaceVertex_UnrollVertex.find(hashKey_a);
				if (got_a != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_a->second] = a;

				zIntPair hashKey_b(f1.getId(), he_1.getVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_b = oriFaceVertex_UnrollVertex.find(hashKey_b);
				if (got_b != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_b->second] = b;

				for (int j = 0; j < f1_numV; j++)
				{
					he_walker_1 = he_walker_1.getNext();
					zPoint C = vPositions[he_walker_1.getVertex().getId()];

					zVector ca = C - A;
					zVector ba = B - A;;

					float s = ((ba ^ ca).length()) / (l_ab * l_ab);
					float c = (ba * ca) / (l_ab * l_ab);

					// alternate point
					/*zPoint c1;
					c1.x = a.x + c * (b.x - a.x) - s * (b.y - a.y);
					c1.y = a.y + c * (b.y - a.y) + s * (b.x - a.x);
					c1.z = 0;*/

					zPoint c1;
					c1.x = a.x + c * (b.x - a.x) + s * (b.y - a.y);
					c1.y = a.y + c * (b.y - a.y) - s * (b.x - a.x);
					c1.z = 0;

					// update  postions of corresponding a & b in unroll mesh
					zIntPair hashKey_c(f1.getId(), he_walker_1.getVertex().getId());
					std::unordered_map<zIntPair, int>::const_iterator got_c = oriFaceVertex_UnrollVertex.find(hashKey_c);
					if (got_c != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_c->second] = c1;
				}

			}





			zItMeshHalfEdge he_walker_2 = he_2;
			int f2_numV = f2.getNumVertices();

			// get positions of the prev edge unrolled.
			zIntPair hashKey_a_prev(f1.getId(), he_2.getStartVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_a_prev = oriFaceVertex_UnrollVertex.find(hashKey_a_prev);
			zPoint a = vPositions_unroll[got_a_prev->second];

			zIntPair hashKey_b_prev(f1.getId(), he_2.getVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_b_prev = oriFaceVertex_UnrollVertex.find(hashKey_b_prev);
			zPoint b = vPositions_unroll[got_b_prev->second];

			// update  postions of corresponding a & b in unroll mesh
			zIntPair hashKey_a(f2.getId(), he_2.getStartVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_a = oriFaceVertex_UnrollVertex.find(hashKey_a);
			if (got_a != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_a->second] = a;

			zIntPair hashKey_b(f2.getId(), he_2.getVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_b = oriFaceVertex_UnrollVertex.find(hashKey_b);
			if (got_b != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_b->second] = b;


			for (int j = 0; j < f2_numV; j++)
			{
				he_walker_2 = he_walker_2.getNext();
				zPoint C = vPositions[he_walker_2.getVertex().getId()];

				float l_ab = he_2.getLength();

				zVector ca = C - A;
				zVector ba = B - A;;

				float s = ((ba ^ ca).length()) / (l_ab * l_ab);
				float c = (ba * ca) / (l_ab * l_ab);

				zPoint c1;
				c1.x = a.x + c * (b.x - a.x) - s * (b.y - a.y);
				c1.y = a.y + c * (b.y - a.y) + s * (b.x - a.x);
				c1.z = 0;

				// update  postions of corresponding a & b in unroll mesh
				zIntPair hashKey_c(f2.getId(), he_walker_2.getVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_c = oriFaceVertex_UnrollVertex.find(hashKey_c);
				if (got_c != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_c->second] = c1;

			}


		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::creatUnrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs)
	{
		zFnMesh fnMesh(o_mesh);

		computeDualGraph_BST(o_mesh, o_dualgraph, bsf_Vertices, bsf_vertexPairs);

		zPoint* vPositions = fnMesh.getRawVertexPositions();
		zColor* vColors = fnMesh.getRawVertexColors();
		zPointArray positions;
		zIntArray pConnects, pCounts;
		zColorArray colors;
		oriVertex_UnrollVertex_map.clear();
		oriVertex_UnrollVertex_map.assign(fnMesh.numVertices(), zIntArray());

		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zIntArray fVerts;
			f.getVertices(fVerts);

			for (auto fV : fVerts)
			{
				int numVerts = positions.size();

				pConnects.push_back(numVerts);
				oriVertex_UnrollVertex_map[fV].push_back(numVerts);

				zIntPair hashKey(f.getId(), fV);
				oriFaceVertex_UnrollVertex[hashKey] = numVerts;

				positions.push_back(vPositions[fV]);
				colors.push_back(vColors[fV]);

			}

			pCounts.push_back(fVerts.size());
		}


		zFnMesh fnMesh_unroll(o_mesh_unroll);
		fnMesh_unroll.create(positions, pCounts, pConnects);
		fnMesh_unroll.setVertexColors(colors);

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::computeDualGraph_BST(zObjMesh& o_mesh, zObjGraph& o_graph, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs)
	{
		zFnMesh fnMesh(o_mesh);

		zIntArray inEdge_dualEdge;
		zIntArray dualEdge_inEdge;
		fnMesh.getDualGraph(o_graph, inEdge_dualEdge, dualEdge_inEdge, true, false, false);

		zFnGraph fnGraph(o_graph);
		fnGraph.setEdgeColor(zColor(1, 1, 0, 1));

		zItGraphVertex v_MaxValence;
		int maxValence = 0;;

		for (zItGraphVertex v(o_graph); !v.end(); v++)
		{
			if (v.getValence() > maxValence)
			{
				v_MaxValence = v;
				maxValence = v.getValence();
			}
		}

		maxValence += 1;

		// breadth search first sorting

		v_MaxValence.getBSF(bsf_Vertices, bsf_vertexPairs);
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::mergeMesh(zObjMesh& o_mesh)
	{
		zObjMesh oTmpMesh = o_mesh;

		zPointArray positions;
		zIntArray pCounts, pConnects;
		zColorArray colors;

		//For every face
		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zPointArray fVPositions;
			f.getVertexPositions(fVPositions);
			zItMeshVertexArray fVertices;
			f.getVertices(fVertices);

			float avg_edge_length_square = 0;
			//Find Avg Edge Length
			for (int i = 0; i < fVPositions.size()-1; i+=2)
			{
				avg_edge_length_square += (fVPositions[i] - fVPositions[i + 1]).length2();
			}
			//Add last one
			avg_edge_length_square += (fVPositions[fVPositions.size() - 1] - fVPositions[0]).length2();

			//Get Avg
			avg_edge_length_square /= fVPositions.size();

			//Iterate over the vertices of that face
			for (auto& p : fVertices)
			{
				int id = -1;
				zPoint pPos = p.getPosition();

				//Check if the position is duplicated
				if (!core.FindDuplicateVector(pPos, positions, id, avg_edge_length_square * 0.02))
				{
					//If unique add to positions
					id = positions.size();
					positions.push_back(p.getPosition());
					colors.push_back(p.getColor());
				}

				pConnects.push_back(id);
			}

			pCounts.push_back(fVPositions.size());
		}

		zFnMesh fnMesh(o_mesh);
		fnMesh.clear();
		fnMesh.create(positions, pCounts, pConnects);
		fnMesh.setVertexColors(colors);
	}

	ZSPACE_TOOLSETS_INLINE zIntPair zTsCarbcomn::getCommonEdge(zItMeshFace& f1, zItMeshFace& f2)
	{
		zIntPair out;

		zItMeshHalfEdgeArray f1_HEdges;
		f1.getHalfEdges(f1_HEdges);

		zItMeshHalfEdgeArray f2_HEdges;
		f2.getHalfEdges(f2_HEdges);


		for (auto& f1HE : f1_HEdges)
		{
			for (auto& f2HE : f2_HEdges)
			{
				if (f1HE.getEdge().getId() == f2HE.getEdge().getId())
				{
					out = zIntPair(f1HE.getId(), f2HE.getId());
					break;
				}
			}
		}

		return out;
	}



	///NON-PLANAR BLOCKS

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::computeVLoops(zObjMesh& oMesh, zIntArray& medialIDS, zIntArray& featuredNumStrides, zVector& norm, vector<zItMeshHalfEdgeArray>& v_Loops, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom)
	{
		int stride = 0;

		for (int i = 0; i < (featuredNumStrides.size() - 1)/2; i++)
		{
			stride += featuredNumStrides[i];
		}

		int startVID = medialIDS[0];
		int endVID = medialIDS[1];

		zItMeshVertex vStart(oMesh, startVID);
		zItMeshVertex vEnd(oMesh, endVID);

		zVector dir = vEnd.getPosition() - vStart.getPosition();

		//numFrames = ceil((vEnd.getPosition().z - vStart.getPosition().z) / spacing) + 1;
		//printf("\n numFrames %i", numFrames);

		zItMeshHalfEdgeArray hEdges_Start;
		vStart.getConnectedHalfEdges(hEdges_Start);

		float ang = FLT_MAX;
		zItMeshHalfEdge heStart;

		for (auto& he : hEdges_Start)
		{
			if (he.getVector().angle(dir) < ang)
			{
				ang = he.getVector().angle(dir);
				heStart = he;
			}
		}

		zItMeshHalfEdge he = heStart;
		norm.normalize();

		zItMeshHalfEdge he_Bottom, he_Top;
		int VCounter = 0;
		int tempCounter = 0;
		do
		{
			zVector fNorm = he.getFace().getNormal();
			fNorm.normalize();

			if (norm * fNorm > 0.98)
			{
				VCounter = tempCounter;
				he_Top = he;
			}

			if (norm * fNorm < -0.98)
			{
				he_Bottom = he;
			}

			he = he.getNext().getSym().getNext();
			tempCounter++;

		} while (he != heStart);


		printf("\n VCounter %i ", VCounter);

		zPointArray positions_top, positions_bottom;
		zIntArray pCounts_top, pCounts_bottom;
		zIntArray pConnects_top, pConnects_bottom;

		zIntArray pMap_bottom, pMap_top;

		zFnMesh fnMesh_in(oMesh);


		pMap_bottom.assign(fnMesh_in.numVertices(), -1);
		pMap_top.assign(fnMesh_in.numVertices(), -1);

		bool corner = true;

		for (int i = 0; i < stride; i++)
		{
			he_Top = he_Top.getNext().getNext();
			he_Bottom = he_Bottom.getNext().getNext();

			if ((i + 1) % stride != 0)
			{
				he_Bottom = he_Bottom.getSym();
				he_Top = he_Top.getSym();
			}
		}

		zColorArray vColor_top, vColors_bottom;

		zItMeshHalfEdge heWalk_Bottom = he_Bottom;
		int walkCounter = 0;
		int loopCounter;
		do
		{
			if (corner)
			{
				loopCounter = v_Loops.size();
				getLoop(heWalk_Bottom, true, corner, VCounter, v_Loops);
				corner = false;

				pMap_bottom[v_Loops[loopCounter][0].getVertex().getId()] = positions_bottom.size();
				positions_bottom.push_back(v_Loops[loopCounter][0].getVertex().getPosition());
				vColors_bottom.push_back(v_Loops[loopCounter][0].getVertex().getColor());

				pMap_top[v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getId()] = positions_top.size();
				positions_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getPosition());
				vColor_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getColor());


			}

			loopCounter = v_Loops.size();

			//zItMeshHalfEdge he = heWalk.getNext();
			getLoop(heWalk_Bottom, true, corner, VCounter, v_Loops);
			//he.getEdge().setColor(zBLUE);

			pMap_bottom[v_Loops[loopCounter][0].getVertex().getId()] = positions_bottom.size();
			positions_bottom.push_back(v_Loops[loopCounter][0].getVertex().getPosition());
			vColors_bottom.push_back(v_Loops[loopCounter][0].getVertex().getColor());

			pMap_top[v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getId()] = positions_top.size();
			positions_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getPosition());
			vColor_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getColor());

			heWalk_Bottom = heWalk_Bottom.getNext().getNext();

			if ((walkCounter + 1) % (2 * stride) != 0) heWalk_Bottom = heWalk_Bottom.getSym();
			else corner = true;
			walkCounter++;

		} while (heWalk_Bottom != he_Bottom);


		// create meshes

		for (int i = 0; i < stride * 2; i++)
		{

			zIntArray fVerts;
			getFaceVerticesFromHalfedge(he_Bottom, true, fVerts);
			for (auto& id : fVerts) pConnects_bottom.push_back(pMap_bottom[id]);
			pCounts_bottom.push_back(fVerts.size());


			fVerts.clear();
			getFaceVerticesFromHalfedge(he_Top, true, fVerts);
			for (auto& id : fVerts) pConnects_top.push_back(pMap_top[id]);
			pCounts_top.push_back(fVerts.size());

			he_Bottom = he_Bottom.getNext().getNext().getSym();
			he_Top = he_Top.getNext().getNext().getSym();
		}



		zFnMesh fnTop(oMesh_top);
		fnTop.create(positions_top, pCounts_top, pConnects_top);
		fnTop.setVertexColors(vColor_top);

		zFnMesh fnBottom(oMesh_bottom);
		fnBottom.create(positions_bottom, pCounts_bottom, pConnects_bottom);
		fnBottom.setVertexColors(vColors_bottom);

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::computeGeodesicScalars(zObjMesh& oMesh, vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, bool normalise)
	{
		zFnMesh fnMesh(oMesh);

		scalars.clear();
		scalars.assign(fnMesh.numVertices(), -1);

		float minMaxDist = 10000;
		vector<zDomainFloat> loopDomains;
		loopDomains.assign(v_Loops.size(), zDomainFloat(10000, -10000));

		for (int l = 0; l < v_Loops.size(); l++)
		{
			float length = 0;
			for (int j = 0; j < v_Loops[l].size(); j++)
			{
				if (j == 0)
				{
					scalars[v_Loops[l][j].getVertex().getId()] = length;
					loopDomains[l].min = length;
				}

				length += v_Loops[l][j].getLength();
				scalars[v_Loops[l][j].getStartVertex().getId()] = length;

				if (j == v_Loops[l].size() - 1 && length < minMaxDist) minMaxDist = length;

				if (length > loopDomains[l].max) loopDomains[l].max = length;
			}
		}

		if (normalise)
		{
			zDomainFloat outDomain(0, minMaxDist);
			for (int l = 0; l < v_Loops.size(); l++)
			{
				for (int j = 0; j < v_Loops[l].size(); j++)
				{
					scalars[v_Loops[l][j].getStartVertex().getId()] = core.ofMap(scalars[v_Loops[l][j].getStartVertex().getId()], loopDomains[l], outDomain);
				}
			}
		}

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		zColor* mesh_vColors = fnMesh.getRawVertexColors();

		zDomainFloat distanceDomain(minScalar, maxScalar);
		printf("\n scalar domain %1.4f %1.4f | minMaxDist %1.4f ", minScalar, maxScalar, minMaxDist);

		zDomainColor colDomain(zColor(1, 0, 0, 1), zColor(0, 1, 0, 1));

		//for (int i = 0; i < fnMesh.numVertices(); i++)
		//{
		//	//mesh_vColors[i] = core.blendColor(scalars[i], distanceDomain, colDomain, zRGB);
		//}

		//fnMesh.computeFaceColorfromVertexColor();
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::computeGeodesicContours(vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, float spacing, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom, zObjMeshArray& oMeshes)
	{

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		int totalContours = ceil((maxScalar - minScalar) / spacing) ;
		float increments = (maxScalar - minScalar) / totalContours;

		printf("\n totalContours %i increments %1.2f ", totalContours, increments);
		zObjMeshArray oTmpMeshes;

		oTmpMeshes.clear();
		oTmpMeshes.assign(totalContours + 1, oMesh_bottom);



		for (int l = 0; l < oTmpMeshes.size(); l++)
		{
			float threshold = l * increments;
			printf("\n threshold %1.4f ", threshold);

			zFnMesh fnMesh(oTmpMeshes[l]);
			zPoint* points = fnMesh.getRawVertexPositions();

			for (int i = 0; i < v_Loops.size(); i++)
			{
				for (int j = 0; j < v_Loops[i].size(); j++)
				{
					float s0 = scalars[v_Loops[i][j].getStartVertex().getId()];
					float s1 = scalars[v_Loops[i][j].getVertex().getId()];

					bool contour = false;
					if (s0 <= threshold && s1 >= threshold)contour = true;
					if (s0 >= threshold && s1 <= threshold)contour = true;

					if (!contour) continue;

					zPoint v0 = v_Loops[i][j].getStartVertex().getPosition();
					zPoint v1 = v_Loops[i][j].getVertex().getPosition();


					zPoint pos1 = util_getContourPosition(threshold, v1, v0, s1, s0);
					points[i] = (pos1);
				}
			}




		}
		/*zFnMesh fnMesh_top(oMesh_top);
		oTmpMeshes.push_back(oMesh_top);*/

		// poke mesh
		oMeshes.clear();
		oMeshes.assign(oTmpMeshes.size(), zObjMesh());

		for (int l = 0; l < oTmpMeshes.size(); l++)
		{
			getPokeMesh(oTmpMeshes[l], oMeshes[l]);
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::computeGeodesicContours(zObjMesh& o_mesh, zFloatArray& scalars, float spacing, zObjGraphArray& o_contourGraphs)
	{

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		int totalContours = ceil((maxScalar - minScalar) / spacing);


		float increments = (maxScalar - minScalar) / totalContours;

		//printf("\n totalContours %i increments %1.2f ", totalContours, increments);

		o_contourGraphs.clear();
		o_contourGraphs.assign(totalContours, zObjGraph());



		for (int i = 0; i < totalContours; i++)
		{
			// Generate the isocontour using the threshold value
			zPointArray positions;
			zIntArray edgeConnects;
			zColorArray vColors;
			int pres = 3;
			zFnMesh fnMesh(o_mesh);
			fnMesh.getIsoContour(scalars, i * increments, positions, edgeConnects, vColors, pres, pow(10, -1 * pres));

			// Create graph from the isocontour
			zFnGraph tempFn(o_contourGraphs[i]);
			tempFn.create(positions, edgeConnects);
			tempFn.setEdgeColor(zColor(1, 1, 1, 1));
			tempFn.setEdgeWeight(2);
			tempFn.setVertexColors(vColors, false);
		}



	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::createSectionGraphs(zObjMeshArray& oMeshes, zObjGraphArray& o_sectionsGraphs)
	{
		o_sectionsGraphs.clear();
		o_sectionsGraphs.assign(oMeshes.size(), zObjGraph());

		int counter = 0;
		for (auto& oMesh : oMeshes)
		{
			createBoundaryEdgeGraph(oMesh, true, o_sectionsGraphs[counter]);

			zFnGraph fnGraph(o_sectionsGraphs[counter]);
			fnGraph.setEdgeColor(zGREEN);
			fnGraph.setEdgeWeight(3);

			counter++;
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::transformAllGraphs(int graphId, zTransform t , bool toLocal)
	{
		zFnGraph fng;
		if (toLocal)
		{
			// Transform
			zTransform tLocal;
			tLocal.setIdentity();

			//transform all graphs to local


			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);

			fng = zFnGraph(o_trimGraphs_cableprofile[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);

		}
		else
		{
			// transform back

			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, true);

			fng = zFnGraph(o_trimGraphs_cableprofile[graphId]);
			fng.setTransform(t, true, true);
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::cleanContourGraph(int graphId)
	{
		zFnGraph fnGraph(o_contourGraphs[graphId]);
		//check for bourndary verticies
		int boundaryCount = 0;
		zItGraphVertexArray vArray;
		for (zItGraphVertex vIt(o_contourGraphs[graphId]); !vIt.end(); vIt++)
		{
			if (vIt.getValence() == 1)
			{
				boundaryCount++;
				vArray.push_back(vIt);
			}
		}

		zPointArray positions;
		zIntArray eConnects;

		fnGraph.getVertexPositions(positions);
		fnGraph.getEdgeData(eConnects);

		if (vArray.size() == 2) //two boundary verticies - we have one loop that is open - still better to create new graph with ordered positions
		{
			eConnects.push_back(vArray[1].getId());
			eConnects.push_back(vArray[0].getId());
			fnGraph.create(positions, eConnects);
			printf("\n cleaning contour graph [%i]. valence = 2 ", graphId);
			/*zItGraphHalfEdge heStart(o_contourGraphs[graphId], vArray[0].getId());
			zItGraphHalfEdge he = heStart;
			int safetyCounter = 0;
			zItGraphHalfEdgeArray loop;
			while (safetyCounter <= fnGraph.numEdges())
			{
				loop.push_back(he);
				he = he.getNext();
				if (he.getStartVertex().checkValency(1)) break;
				safetyCounter++;
			}
			createGraphFromHEArray(loop, o_contourGraphs[graphId]);*/
		}
		else if (vArray.size() == 0) //no boundary verticies - check for how many closed loops
		{
			printf("\n cleaning contour graph [%i]. valence = 0 ", graphId);


			//check for closed loops
			vector<zItGraphHalfEdgeArray> loops;

			zBoolArray visited;
			visited.assign(fnGraph.numEdges(), false);

			int edgeCounter = 0;
			zItGraphHalfEdgeArray longestLoop;
			while (edgeCounter < fnGraph.numEdges())
			{
				zItGraphHalfEdge heStart(o_contourGraphs[graphId], edgeCounter);
				zItGraphHalfEdge he = heStart;
				int safetyCounter = 0;
				zItGraphHalfEdgeArray loop;

				while (safetyCounter <= fnGraph.numEdges())
				{
					visited[he.getEdge().getId()] = true;
					loop.push_back(he);
					he = he.getNext();
					if (he == heStart) break;
					safetyCounter++;
				}
				edgeCounter += loop.size();
				if (loop.size() > longestLoop.size())
				{
					longestLoop = loop;
				}
			}
			util_createGraphFromHEArray(longestLoop, o_contourGraphs[graphId]);
		}
		else
		{
			printf("\n cleaning contour graph [%i]. other - valence = %i ", graphId, vArray.size());

			//more than two boundary verticies - we have multiple loops
			//check for closed loops
			vector<zItGraphHalfEdgeArray> loops;
			zBoolArray visited;
			visited.assign(fnGraph.numEdges(), false);
			int edgeCounter = 0;
			zItGraphHalfEdgeArray longestLoop;
			bool loopClosed = false;
			//get the loop from all boundary verticies
			for (auto& v : vArray)
			{
				zItGraphHalfEdge heStart = v.getHalfEdge();
				if (visited[heStart.getEdge().getId()])
				{
					continue;
				}
				zItGraphHalfEdge he = heStart;
				zItGraphHalfEdgeArray loop;
				int safetyCounter = 0;

				while (safetyCounter <= fnGraph.numEdges())
				{
					visited[he.getEdge().getId()] = true;
					loop.push_back(he);
					he = he.getNext();
					if (he == heStart)
					{
						loopClosed = true;
						break;
					}
					if (he.getStartVertex().checkValency(1))
					{
						loopClosed = false;
						break;
					}
					safetyCounter++;
				}
				if (loop.size() > longestLoop.size())
				{
					longestLoop = loop;
				}
			}
			util_createGraphFromHEArray(longestLoop, o_contourGraphs[graphId]);


			if (!loopClosed)
			{
				printf("\n			cleaning contour graph [%i]. valence = %i - not closed loop", graphId, vArray.size());

				fnGraph.getVertexPositions(positions);
				fnGraph.getEdgeData(eConnects);
				zItGraphVertexArray vs;
				for (zItGraphVertex vIt(o_contourGraphs[graphId]); !vIt.end(); vIt++)
				{
					if (vIt.checkValency(1))
					{
						vs.push_back(vIt);
					}
				}
				//makesure that u have two vertex
				if (vs.size() == 2)
				{
					eConnects.push_back(vs[1].getId());
					eConnects.push_back(vs[0].getId());
					fnGraph.create(positions, eConnects);
				}
				else
				{
					printf	("\n error: invalid cleaning contour graph. valence = %i ", vs.size());
					//throw std::invalid_argument(" error: invalid cleaning contour graph.");
				}



			}


		}
	}


	ZSPACE_TOOLSETS_INLINE void zTsCarbcomn::readJSON(string path, int _blockID, bool runBothPlanes, bool runPlaneLeft, bool flip)
	{
		printf("\n readJSON 0");

		json j;
		zFnMesh fnMesh(o_GuideMesh);
		bool fileChk = fnMesh.json_read(path, j);

		if (!fileChk)
		{
			printf("\n error: invalid outPath %s", path.c_str());
			return;
		}

		fnMesh.clear();
		fnMesh.from(path, zJSON);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();

		blockId = _blockID;
		int sID = flip ? j["MedialStartEnd"][1] : j["MedialStartEnd"][0];
		int eID = flip ? j["MedialStartEnd"][0] : j["MedialStartEnd"][1];

		isCorner = false;

		runningType = runBothPlanes ? 0 : runPlaneLeft ? 1 : 2;

		base_local.setIdentity();
		base_world.setIdentity();
		//medial axis
		compute_MedialGraph(o_GuideMesh, sID, eID);


		//get feature stride
		FeaturedNumStrides.clear();
		core.json_readAttribute(j, "FeaturedNumStrides", FeaturedNumStrides);
		medialIDS.clear();
		core.json_readAttribute(j, "MedialStartEnd", medialIDS);
		StartCornerVID = j["StartCornerVID"];

		isRegular = true;
		compute_SliceMesh_Regular(o_GuideMesh, sID, eID, FeaturedNumStrides);
		printf("FeaturedNumStrides %i", FeaturedNumStrides[0]);
		printf("\n slice mesh regular");
	}


	}
