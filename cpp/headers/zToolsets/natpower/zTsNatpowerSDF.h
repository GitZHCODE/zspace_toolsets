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

#ifndef ZSPACE_TS_NATPOWER_SDF_H
#define ZSPACE_TS_NATPOWER_SDF_H



#pragma once

#include <headers/base/zSpace_Toolsets.h>

#include <headers/zCore/base/zExtern.h>

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zInterface/functionsets/zFnMeshField.h>

#include <headers/zInterOp/functionSets/zFnPlane.h>
#include <headers/zInterOp/functionSets/zFnNurbsCurve.h>


#include <igl/point_mesh_squared_distance.h>
#include <igl/boundary_loop.h>
#include <igl/lscm.h>



//#include <igl/avg_edge_length.h>
//#include <igl/cotmatrix.h>
//#include <igl/invert_diag.h>
//#include <igl/massmatrix.h>
//#include <igl/parula.h>
//#include <igl/per_corner_normals.h>
//#include <igl/per_face_normals.h>
//#include <igl/per_vertex_normals.h>
//#include <igl/principal_curvature.h>
//#include <igl/gaussian_curvature.h>
//#include <igl/read_triangle_mesh.h>
//
//#include <igl/point_mesh_squared_distance.h>



#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include<execution>

//#include<atlsafe.h>
//#include <comdef.h>
//#include<oaidl.h>

//using namespace std;


namespace zSpace
{
	enum zBlockType
	{
		Bottom,
		Top,
		Wall,
		Arch
	};

	struct zPair_hash
	{
		template <class T1, class T2>
		size_t operator()(const pair<T1, T2>& p) const
		{
			auto hash1 = hash<T1>{}(p.first);
			auto hash2 = hash<T2>{}(p.second);

			if (hash1 != hash2)
			{
				return hash1 ^ hash2;
			}

			// If hash1 == hash2, their XOR is zero.
			return hash1;
		}
	};
	struct zPrintParamSDF
	{
		const float printWidthInterior = 0.048;
		const float printWidthExterior = 0.036;
		const float printOverlap = 0.002; //2mm overlap total
		const float offset_outer_interior = printWidthInterior / 2.0; //total 48mm layer width
		const float offset_outer_exterior = printWidthExterior / 2.0; //total 36mm layer width
		const float targetInteriorGap = 0.044; //layer width - overlapping 
		const float targetExteriorGap = 0.034; //layer width - overlapping
		const float offset_inner_interior = (targetInteriorGap - offset_outer_interior);
		const float offset_inner_exterior = (targetExteriorGap - offset_outer_exterior);
		const float bracingEdgeWidth = 0.012;
		//slots width
		const float slotStartWidth = 0.022 / 2.0;
		const float slotBracingWidth = 0.018 / 2.0;
		//cable
		const float cableWidth = 0.045;



	};
	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/


	/*! \class zTsSDFSlicer
	*	\brief A toolset for SDF based slicing.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLSETS zTsNatpowerSDF
	{
		//protected:
	public:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief core utilities Object  */
		zUtilsCore core;
		//zPrintParamSDF _printParameters;

		/*!	\brief input guide mesh object  */
		zObjMesh o_GuideMesh;

		/*!	\brief left mesh object  */
		zObjMesh o_SliceMesh_Left;

		/*!	\brief right mesh object  */
		zObjMesh o_SliceMesh_Right;

		/*!	\brief medial graph object  */
		zObjGraph o_MedialGraph;

		/*!	\brief left start & end planes  */
		zTransform leftPlanes[2];

		/*!	\brief right start & end planes  */
		zTransform rightPlanes[2];

		/*!	\brief container of section frames  */
		vector<zTransform> sectionFrames;

		/*!	\brief container of section graph objects  */
		zObjGraphArray o_sectionGraphs;

		/*!	\brief container of contour graph objects  */
		zObjGraphArray o_contourGraphs;

		zObjGraphArray o_contourGraphs_flatten;

		/*!	\brief container of raft graph objects  */
		zObjGraphArray o_raftGraphs;

		/*!	\brief container of trim graph objects  */
		zObjGraphArray o_trimGraphs;
		
		/// <summary>
		/// Used to identify hard feature points in the post-processing
		/// </summary>
		zObjGraphArray o_trimGraphs_features_hard;
		/// <summary>
		/// Used to identify feature points in the post-processing
		/// </summary>
		zObjGraphArray o_trimGraphs_features_soft;

		/// <summary>
		/// Used in the SDF and in the post-processing
		/// </summary>
		zObjGraphArray o_trimGraphs_bracing;
		zObjGraphArray o_trimGraphs_bracing_flat;

		/// <summary>
		/// Used to align seam in the post-processing
		/// </summary>
		zObjGraphArray o_trimGraphs_seamAlignment;
		/// <summary>
		/// Used to separate inner path from outer path in the post-processing
		/// </summary>
		zObjGraphArray o_trimGraphs_SlotSide;

		zObjMeshArray o_sectionMeshes;
		zObjMeshArray o_sectionMeshesPar;


		zObjGraphArray o_CableGraphs;

		int cableGraphId_left = -1;
		int cableGraphId_right = -1;

		bool checkOranges = true;
		bool checkMagentas = true;

		int runningType = 0; //< 0 running both planes, 1 running left planes only, 2 running right planes only


		vector<zVectorArray> o_contourNormals;

		vector<zObjGraph> o_contourHeightLines;

		//--------------------------
		//---- PRINT ATTRIBUTES
		//--------------------------

		/*!	\brief minimum layer height  in the print */
		float minLayerHeight = 0;

		/*!	\brief maximum layer height in the print  */
		float maxLayerHeight = 0;

		/*!	\brief total length of print  */
		float totalLength = 0;

		zObjPointCloud criticalMinLayer_pts;

		zObjPointCloud criticalMaxLayer_pts;

		zDomainFloat printHeightDomain;

		zDomainFloat actualPrintHeightDomain;

		zDomainFloat neopreneOffset;

		zObjGraphArray o_allCable;

		//--------------------------
		//---- SDF ATTRIBUTES
		//--------------------------

		zObjMeshScalarField o_field;
		zObjGraph o_isoContour;

		//--------------------------
		//---- GRADIENT ATTRIBUTES
		//--------------------------

		zObjMesh o_gradientTriMesh;
		MatrixXd gradientTriMesh_V;
		MatrixXi gradientTriMesh_FTris;
		zDomainFloat offsetDomain;

		//--------------------------
		//---- COLOR ATTRIBUTES
		//--------------------------

		zColorArray blockColors;

		//bool leftPlaneExists = false;
		//bool rightPlaneExists = false;

		int blockId = 1;

		int numMagentaLoops = 1;

		zTransform base_world, base_local;

		bool planarBlock;
		bool isCorner;
		bool isRegular;
		bool isFront;

		int StartCornerVID = -1;
		bool _interpolateFramesOrigins = false;

		zBlockType blockType;
		zIntArray FeaturedNumStrides;
		zIntArray medialIDS;
	private:
		//color settings 
		zColor _col_in_corner_st = zRED;
		zColor _col_out_corner_st = zCYAN;
		zColor _col_in_corner = zGREEN;
		zColor _col_out_corner = zYELLOW;

		zColor _col_in_feature = zMAGENTA;
		zColor _col_out_feature = zORANGE;
		zColor _colorPattern = zBLUE;

		float _printOverlap = 0.001; //2mm overlap (1mm on each side)

		//bool checkWall = false;
	public:

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsNatpowerSDF();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsNatpowerSDF();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------

		void createFieldMeshFromMeshBounds(float cellSize, float offset);
		void createFieldMeshFromSectionBounds(float cellSize, float offset);

		/*! \brief This method creates the field mesh.
		*
		*	\param		[in]	bb			- input domain of bounds.
		*	\param		[in]	resX		- input resolution of field in X.
		*	\param		[in]	resY		- input resolution of field in Y.
		*	\since version 0.0.4
		*/
		void createFieldMesh(zDomain<zPoint>& bb, int resX, int resY);

		void createFieldMeshCellSize(zDomain<zPoint>& bb, float cellSizeX, float cellSizeY);

		//--------------------------
		//--- SET METHODS 
		//--------------------------

		/*! \brief This method sets mesh from json.
		*
		*	\param		[in]	dir		    - input JSON file director.
		**	\param		[in]	_blockID	- input block ID.
		*	\since version 0.0.4
		*/
		void setFromJSON(string dir, int _blockID, bool runBothPlanes = true, bool runPlaneLeft = false);

		void setCableGraph(string dir);

		/*! \brief This method sets the slice mesh object.
		*
		*	\param		[in]	_o_SliceMesh		- input mesh object.
		* 	\param		[in]	left				- input boolean indicating if the planes for the left or right side meshes.
		*	\since version 0.0.4
		*/
		void setSliceMesh(zObjMesh& _o_SliceMesh, bool left);

		/*! \brief This method sets the medial graph object.
		*
		*	\param		[in]	_o_MedialGraph			- input graph object.
		*	\since version 0.0.4
		*/
		void setMedialGraph(zObjGraph& _o_MedialGraph);

		/*! \brief This method sets the start and end plane.
		*
		*	\param		[in]	_sTransform			- input start plane.
		*	\param		[in]	_eTransform			- input end plane.
		* 	\param		[in]	left				- input boolean indicating if the planes for the left or right side meshes.
		*	\since version 0.0.4
		*/
		void setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane, bool left);

		void setGradientTriMesh(zObjMesh& _o_gradientTriMesh);

		void setOffsetDomain(zDomainFloat& _offsetDomain);

		void setTransforms(bool toLocal);

		/*! \brief This method sets the section frames from the input container of planes.
		*
		*	\param		[in]	_sectionFrames			- input container of planes.
		*	\since version 0.0.4
		*/
		void setFrames(vector<zPlane>& _sectionFrames);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method gets the block frames.
		*
		*	\param		[in]	blockId					- input block index.
		*	\return				vector<zTransform>	    - cantainer of transforms if they exist.
		*	\since version 0.0.2
		*/
		zTransform* getRawBlockStartEnd(bool left);


		/*! \brief This method gets the block frames.
		*
		*	\param		[in]	blockId					- input block index.
		*	\return				vector<zTransform>	    - cantainer of transforms if they exist.
		*	\since version 0.0.2
		*/
		vector<zTransform> getBlockFrames();

		/*! \brief This method gets the block section graphs
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	-  pointer conatiner of graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockSectionGraphs(int& numGraphs);

		/*! \brief This method gets the block section graphs
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	-  pointer conatiner of graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockRaftGraphs(int& numGraphs);


		/*! \brief This method gets the block SDF contour graphs
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer conatiner of graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockContourGraphs(int& numGraphs);

		void getBlockContourGraphsSequence(int& numGraphs, int* vSequence);

		static zIntArray getGraphSequence(zObjGraph graph);

		/*! \brief This method gets the block trim graphs
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer conatiner of graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockTrimGraphs(int& numGraphs);

		/*! \brief This method gets the critical points of the section graphs
		*
		*	\param		[in]	blockId					- input block index.
		*	\return				zPoint*					- pointer conatiner of points if they exist.
		*	\since version 0.0.2
		*/
		zObjPointCloud* getRawCriticalPoints(bool minHeight);

		/*! \brief This method gets pointer to the internal field object.
		*
		*	\return				zObjMeshScalarField*					- pointer to internal field object.
		*	\since version 0.0.4
		*/
		zObjMeshScalarField* getRawFieldMesh();

		/*! \brief This method gets pointer to the internal medial graph object.
		*
		*	\return				zObjGraph*					- pointer to internal graph object.
		*	\since version 0.0.4
		*/
		zObjGraph* getRawMedialGraph();

		/*! \brief This method gets pointer to the internal left mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawLeftMesh();

		/*! \brief This method gets pointer to the internal right mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawRightMesh();

		/*! \brief This method gets pointer to the internal guide mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawGuideMesh();

		/*! \brief This method gets pointer to the internal guide mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawGradientMesh();

		/*! \brief This method gets pointer to the internal mesh scalar field  object.
		*
		*	\return				zObjMeshScalarField*					- pointer to internal mesh scalar field object.
		*	\since version 0.0.4
		*/
		zObjMeshScalarField* getRawMeshScalarField();

		//--------------------------
		//---- COMPUTE METHODS 
		//--------------------------

		bool isPlanarBlock();



		///----------SLICE MESH METHODS

		/*! \brief This method computes the left or right slice mesh from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial spine.
		*	\param		[in]	endVID				- input end vertex id of medial spine.
		*	\param		[in]	blockStride			- input stride of edges for left and right blocks.
		* 	\param		[in]	left				- input boolean indicating if the planes for the left or right side meshes.
		*	\since version 0.0.4
		*/
		void compute_SliceMesh(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides, bool left);

		void compute_SliceMesh_Pentagon(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides, bool left);
		void compute_SliceMesh_Regular(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides);
		void compute_SliceMesh_Top(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides);
		//Slice mesh: helper methods

		/*! \brief This method computes the medial graph from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial graph.
		*	\param		[in]	endVID				- input end vertex id of medial graph.
		*	\since version 0.0.4
		*/
		void compute_MedialGraph(zObjMesh& o_Mesh, int startVID, int endVID);

		/*! \brief This method computes the medial graph from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial spine.
		*	\param		[in]	endVID				- input end vertex id of medial spine.
		*	\param		[in]	blockStride			- input stride of edges for left and right blocks.
		*	\param		[in]	braceStride			- input stride for edges for braces.
		*	\since version 0.0.4
		*/
		void compute_Medial_BraceEdges(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, int braceStride);


		///----------PRINT BLOCKS METHODS


		/*! \brief This method computes the.
		*
		* 	\param		[in]	printLayerDepth				- input print layer depth.
		*	\since version 0.0.4
		*/
		void compute_PrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum = 0, int numSmooth = 0, zDomainFloat _neopreneOffset = zDomainFloat(0, 0), bool compFrames = true, bool compSDF = true);

		//Print blocks: helper methods

		/*! \brief This method computes the.
		*
		* 	\param		[in]	printLayerDepth				- input print layer depth.
		*	\since version 0.0.4
		*/

		void compute_PrintSectionFromPlaneSpacing(float printPlaneSpacing, zDomainFloat& _printHeightDomain, zDomainFloat _neopreneOffset, bool& frameCHECKS, bool& sdfCHECKS, bool& geomCHECKS);


		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	printLayerDepth				- input print layer depth.
		*	\param		[in]	guideMesh_vertex			- input guide mesh vertex.
		*	\since version 0.0.4
		*/
		void compute_PrintBlock_Frames(float printPlaneSpacing, float neopreneOffset_start, float neopreneOffset_end, bool left);

		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void compute_PrintBlock_Sections(bool left, bool& outGeomChk);



		/*! \brief This method checks the layer heights of the block
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		bool check_PrintLayerHeights(bool& checkSDF, bool& checkGeometry);

		/*! \brief This method check the block interfaces are planar.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		bool check_InterfacePoints(bool left, float distTolerance);

		bool check_sectionGraphGeomCheck(zObjGraph& graph);



		//Print blocks: trim methods
		void compute_PrintBlock_ComputeTrimGraphs();


		void slotGraph_Arch(int graphId, float graphLength, zObjGraph& innerHE);
		void compute_TrimGraphs_Boundary_Arch(int graphId, zObjGraph& outGraph);
		void compute_TrimGraphs_BracingCable(int graphId, zObjGraph& outGraph);
		void compute_TrimGraphs_CableBracing_2(int graphId, zObjGraph& outGraph);
		void compute_TrimGraphs_BoundaryFeature(int graphId, zObjGraph & outGraph_hardFeature, zObjGraph& outGraph_softFeature);
		void compute_TrimGraphs_SlotSide(int graphId, zObjGraph& outGraph_splitGraph);
		void compute_TrimGraphs_BracingWall(int graphId, zObjGraph& outGraph);		
		
		void compute_TrimGraphs_SlotSide(zObjGraph& sectionGraph, zObjGraph& outGraph_splitGraph);
		void compute_TrimGraphs_BracingWall(zObjGraph& sectionGraph, zObjGraph& outGraph);

		void compute_topAndBottom(int graphId, zItGraphVertexArray& innerVertx, zItGraphVertexArray& outerVertx);
		void compute_topAndBottom(zObjGraph& sectionGraph, zItGraphVertexArray& innerVertx, zItGraphVertexArray& outerVertx);


		//--------------------------
		//---- UTILITY METHODS
		//--------------------------
		void getPerpendicularVector(zPlane& plane, zVector edgeVector, zPoint midPoint, float graphLength, zObjGraph& outGraph);
		zVector averageVectorsAtGraphVertex(zItGraphVertex& v);
		void combineMultipleGraphs(zObjGraphArray& inGraphs, zObjGraph& outGraph);
		zPoint getGraphPointAtParameter(zObjGraph& inGraph, float normalizedPar, int& outEdgeIndex);
		zPoint getPointAtParameterHalfEdge(zItGraphHalfEdge& he, float normalizedPar);
		
		void cleanContourGraph(zObjGraph& refGraph);

		double normalise(double value, double min, double max);
		double denormalise(double value, double min, double max);
		/*! \brief This method computes the SDF for the blocks.
		*
		*	\since version 0.0.4
		*/
		void compute_SDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth, float neopreneOffset, float raftWidth);


		/*! \brief This method compute the block frames.
		*
		*	\param		[out]	outGraph						- output trim graph.
		*	\since version 0.0.4
		*/
		void compute_PrintBlock_TrimGraphs(zObjGraph& inPolyObj, zObjGraph& o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE);
		/*! \brief This method compute the block SDF for the deck.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void compute_BlockSDF_Planar_wall(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		void compute_BlockSDF_Planar_bracing(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);
		void compute_BlockSDF_Planar_pentagon(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void compute_BlockSDF_NonPlanar(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		void compute_CableSectionPoints(int graphId, zObjGraph& o_cableGraph, zPointArray& intersectionPts, float threshold = 0.0);
		int get_CableGraphIndexPerGraph(int graphId);

		void compute_GraphEdgesForSlot_Arch();

		zPoint getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh);
		void isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints);
		void intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints, float threshold = 0.0);


		bool exportJSON(string pathCurrent, string dir, string filename, float printLyerWidth, float raftLayerWidth);
		bool exportJSON_update(string pathCurrent, string dir);
		bool exportJSON_sliceMesh(string pathCurrent, string dir);


		bool exportJSON_graphID(string dir, int graphId, bool left);
		bool exportJSON_graphID_trims(string folderName, string extName, int graphId);
		bool exportJSON_graphID_contours(string folderName, string extName, int graphId);
		bool exportJSON_graphID_section(string folderName, string extName, int graphId);




		void computeClosestPointToGradientMesh(zPointArray& inPoints, zIntArray& faceID, zPointArray& closestPoints);

		float computeWeightedGradientValue(int& faceID, zPoint& closestPt);


		/*! \brief This method checks the layer heights for all the blocks in the input directory.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void check_PrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset, bool runBothPlanes = true, bool runPlaneLeft = false);



		//--------------------------
		void getPokeMesh(zObjMesh& o_mesh, zObjMesh& o_TriMesh);
		void getLoop(zItMeshHalfEdge& heStart, bool forward, bool corner, int vCounter, vector<zItMeshHalfEdgeArray>& v_Loops);
		void getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zPointArray& fVerts, zColorArray& fVColors);
		void getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zIntArray& fVerts);
		void createBoundaryEdgeGraph(zObjMesh& o_mesh, bool closeGraph, zObjGraph& o_Graph);
		void colorMesh(zObjMesh& o_mesh, zFloatArray& scalars);
		void setPtGraph(zObjGraph& o_Graph, zPoint& refPt, bool setX, bool setY, bool setZ);
		void setPtMesh(zObjMesh& o_Mesh, zPoint& refPt, bool setX, bool setY, bool setZ);
		void getBoundaryOffset(zObjMesh& _oMesh, bool keepExistingFaces, float offset, zObjMesh& outMesh);
		void closestPointsToMesh(zPointArray& inPoints, zObjMesh oMesh, zIntArray& faceIDs, zPointArray& closestPoints, zVectorArray& printNorms);
		void getPrintHeight(zPointArray& pPoints, zVectorArray& pNorms, zObjMesh& o_Mesh, zFloatArray pHeights, zObjGraph& outPrintHeightLines);
		void projectToMesh(zPointArray& pPoints, zObjMesh& o_Mesh, zPointArray& updatePts, zVectorArray& pNorm);
		void UVParametrisation(zObjMesh& oMesh, zObjMesh& oParamMesh); 
		void getBaryCentricCoordinates_triangle(zPoint& pt, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& baryCoordinates);
		void getProjectionPoint_triangle(zPoint& baryCoordinates, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& projectionPt);
		void barycentericProjection_triMesh(zObjGraph& o_graph, zObjMesh& o_inMesh, zObjMesh& o_projectionMesh, zVectorArray& outNotmals);

		void unrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zIntPairArray& bsf_vertexPairs, zTransform& outTransformStart);
		void creatUnrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs);
		void computeDualGraph_BST(zObjMesh& o_mesh, zObjGraph& o_graph, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs);
		void mergeMesh(zObjMesh& o_mesh);
		zIntPair getCommonEdge(zItMeshFace& f1, zItMeshFace& f2);
		void computeVLoops(zObjMesh& oMesh, zIntArray& medialIDS, zIntArray& featuredNumStrides, zVector& norm, vector<zItMeshHalfEdgeArray>& v_Loops, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom);

		void computeGeodesicScalars(zObjMesh& oMesh, vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, bool normalise);
		void computeGeodesicContours(vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, float spacing, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom, zObjMeshArray& oMeshes);
		void computeGeodesicContours(zObjMesh& o_mesh, zFloatArray& scalars, float spacing, zObjGraphArray& o_contourGraphs);
		void createSectionGraphs(zObjMeshArray& oMeshes, zObjGraphArray& o_sectionsGraphs);

		void transformAllGraphs_planar(int graphId, bool toLocal);
		void transformAllGraphs(int graphId, zTransform t,  bool toLocal);


		//--------------------------
		//---- POST-PROCESSING METHODS
		//--------------------------
		void cleanContourGraph(int graphId);
		 

		void scalar_getOffset(int graphId, zScalarArray& polyField, zIntArray& edgeId, zScalarArray & out_scalar_offset_outer, zScalarArray & out_scalar_offset_inner);

		void graphIntersection(zObjGraph& graph, zObjGraphArray& trims, zObjGraph& outGraph);



		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:

		zItMeshHalfEdge getStartHalfEdge(zObjMesh& o_mesh, int startVID, int endVID);
		void createGraphFromHEArray(zItGraphHalfEdgeArray& heArray, zObjGraph& outGraph);
		void polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength);
		void slotGraph_0(zObjGraph& inPoly, zObjGraph& innerHE, float& innerLength);
		/// <summary>
		/// 
		/// </summary>
		/// <param name="inPoly"></param>
		/// <param name="innerHE"></param>
		/// <param name="innerLength"></param>
		void slotGraph_1(zPlane plane, zObjGraph& inPoly, float graphLength, bool iterate, zObjGraph& outGraph);
		void splitGraph_1(zPlane plane, zObjGraph& inPoly, float offset, float trim, zObjGraph& outGraph);
		
		


		/// <summary>
		/// 
		/// </summary>
		/// <param name="graph"></param>
		/// <param name="startColor"></param>
		/// <param name="endColor"></param>
		bool getShortestHEsBetweenColors(zObjGraph& graph, zColor startColor, zColor endColor, zItGraphHalfEdgeArray& outHEs);
		/// <summary>
		/// 
		/// </summary>
		/// <param name="graph"></param>
		/// <param name="samplePoint"></param>
		/// <param name="outPoint"></param>
		/// <param name="dist"></param>
		/// <returns>Edge index if closest point</returns>
		int getGraphClosestPoint(zObjGraph& graph, zPoint& samplePoint, zPoint& outPoint, float& dist);

		int getHeArrayClosestPoint(zItGraphHalfEdgeArray& hes, zPoint& samplePoint, zPoint& outPoint, float& dist);
		
		
		
		
		
		void getScalars_3dp_slot(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset);

		void getScalars_3dp_pattern(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, bool alternate, bool alternateCheck);
		void getScalars_3dp_pattern_2(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, float offset2, bool alternate, bool alternateCheck);
		void getScalars_3dp_pattern_3(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, float offset2, bool alternate, bool alternateCheck);
		void getScalars_3dp_patternMesh(zScalarArray& scalars, zObjGraph& o_sectionGraph, zIntArray& faceIDs, zPointArray& closestPoints, zDomainFloat& offsetDomain);


		void getScalars_3dp_brace(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset, bool alternate);
		void getScalars_3dp_brace2(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset, bool alternate);

		void getScalars_3dp_trim(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate);
		void getScalars_3dp_trimStart(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate);
		void getScalars_3dp_cable(zScalarArray& scalars, zPoint cablePoint, float radius);


		void getScalars_wall_all(zScalarArray& scalars, zObjGraph& o_sectionGraph);


		void readJSON_planarBlock(string path, int _blockID, bool runBothPlanes = true, bool runPlaneLeft = false);
		void get2DArrayFromTransform(zTransform& transform, vector<zDoubleArray>& arr);
	};

}





#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/natpower/zTsNatpowerSDF.cpp>
#endif

#endif