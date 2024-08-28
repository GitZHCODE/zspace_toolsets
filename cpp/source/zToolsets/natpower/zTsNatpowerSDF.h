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
		zUtilsCore coreUtils;

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

		/*!	\brief container of raft graph objects  */
		zObjGraphArray o_raftGraphs;

		/*!	\brief container of trim graph objects  */
		zObjGraphArray o_trimGraphs;

		zObjGraphArray o_CableGraphs;

		int cableGraphId_left = -1;
		int cableGraphId_right = -1;

		bool checkOranges = true;
		bool checkMagentas = true;

		int runningType = 0; //< 0 running both planes, 1 running left planes only, 2 running right planes only


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
	private:
		//color settings 
		zColor _colorCornersStart = zRED;
		zColor _colorCornersEnd = zCYAN;
		zColor _colorCornersOuter = zYELLOW;
		zColor _colorCornersInner = zGREEN;

		zColor _colorFeatureInner = zMAGENTA;
		zColor _colorFeatureOuter = zORANGE;
		zColor _colorPattern = zBLUE;

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


		/*! \brief This method computes the.
		*
		* 	\param		[in]	printLayerDepth				- input print layer depth.
		*	\since version 0.0.4
		*/

		void compute_FrontBackHE(zObjGraph& graph, zItGraphHalfEdgeArray& outHeBack, zItGraphHalfEdgeArray& outHEFront);


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

		void slotGraph_Arch(int graphId, float graphLength, zObjGraph& innerHE);
		void compute_TrimGraphsHEs_Boundary_Arch(int graphId, zItGraphHalfEdgeArray& outHEs);
		void compute_TrimGraphsHEs_CableBracing(int graphId, zItGraphHalfEdgeArray& outHEs);
		void compute_TrimGraphsHEs_BoundaryFeature(int graphId, zItGraphHalfEdgeArray& outHEs);



		//--------------------------
		//---- UTILITY METHODS
		//--------------------------



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

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void compute_BlockSDF_NonPlanar(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		void compute_CableSectionPoints(int graphId, zObjGraph& o_cableGraph, zPointArray& intersectionPts);
		int get_CableGraphIndexPerGraph(int graphId);

		void compute_GraphEdgesForSlot_Arch();

		zPoint getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh);
		void isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints);
		void intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints);


		bool exportJSON(string pathCurrent, string dir, string filename, float printLyerWidth, float raftLayerWidth);





		void computeClosestPointToGradientMesh(zPointArray& inPoints, zIntArray& faceID, zPointArray& closestPoints);

		float computeWeightedGradientValue(int& faceID, zPoint& closestPt);


		/*! \brief This method checks the layer heights for all the blocks in the input directory.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void check_PrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset, bool runBothPlanes = true, bool runPlaneLeft = false);


		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:

		zItMeshHalfEdge getStartHalfEdge(zObjMesh& o_mesh, int startVID, int endVID);
		void createGraphFromHEArray(zItGraphHalfEdgeArray& heArray, zObjGraph& outGraph);
		void getPerpendicularVector(zPlane& plane, zVector edgeVector, zPoint midPoint, float graphLength, zObjGraph& outGraph);
		void polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength);
		void slotGraph_0(zObjGraph& inPoly, zObjGraph& innerHE, float& innerLength);
		/// <summary>
		/// 
		/// </summary>
		/// <param name="inPoly"></param>
		/// <param name="innerHE"></param>
		/// <param name="innerLength"></param>
		void slotGraph_1(zPlane plane, zObjGraph& inPoly, float graphLength, zObjGraph& innerHE);
		/// <summary>
		/// 
		/// </summary>
		/// <param name="plane"></param>
		/// <param name="inPoly"></param>
		/// <param name="graphLength"></param>
		/// <param name="innerHE"></param>
		

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



		void readJSON_planarBlock(string path, int _blockID, bool runBothPlanes = true, bool runPlaneLeft = false);

	};

}





#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/natpower/zTsNatpowerSDF.cpp>
#endif

#endif