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
//

#ifndef ZSPACE_TS_GEOMETRY_SDFSLICER_H
#define ZSPACE_TS_GEOMETRY_SDFSLICER_H



#pragma once

#include "base/zSpace_Toolsets.h"

#include <zCore/base/zExtern.h>

#include <zInterface/functionsets/zFnMesh.h>
#include <zInterface/functionsets/zFnGraph.h>
#include <zInterface/functionsets/zFnParticle.h>

#include <zInterface/functionsets/zFnMeshField.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

//#include<atlsafe.h>
//#include <comdef.h>
//#include<oaidl.h>

using namespace std;


namespace zSpace
{
	
		
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

	class ZSPACE_TOOLSETS zTsSDFSlicer
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
				
		//--------------------------
		//---- SDF ATTRIBUTES
		//--------------------------

		zObjMeshScalarField o_field;
		zObjGraph o_isoContour;

		

		//--------------------------
		//---- COLOR ATTRIBUTES
		//--------------------------
		
		zColor red, yellow, green, cyan, blue, magenta, grey , orange;

		zColorArray blockColors;

		bool leftPlaneExists, rightPlaneExists;

		int blockId = 1;

		zTransform base_world, base_local;

		bool deckBlock;

	public:

		

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsSDFSlicer();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsSDFSlicer();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------


		/*! \brief This method creates the field mesh.
		*
		*	\param		[in]	bb			- input domain of bounds.
		*	\param		[in]	resX		- input resolution of field in X.
		*	\param		[in]	resY		- input resolution of field in Y.
		*	\since version 0.0.4
		*/
		void createFieldMesh(zDomain<zPoint> &bb,  int resX , int resY);

		//--------------------------
		//--- SET METHODS 
		//--------------------------

		/*! \brief This method creates the filed mesh.
		*
		*	\param		[in]	bb			- input domain of bounds.
		*	\param		[in]	resX		- input resolution of field in X.
		*	\param		[in]	resY		- input resolution of field in Y.
		*	\since version 0.0.4
		*/
		void setFromJSON(string path, int blockStride, int braceStride);

		/*! \brief This method sets mesh from json.
		*
		*	\param		[in]	dir		    - input JSON file director.
		**	\param		[in]	_blockID	- input block ID.
		*	\since version 0.0.4
		*/
		void setFromJSON(string dir, int _blockID);

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


		void setTransforms(bool toLocal);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		
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
		zObjGraphPointerArray getBlockSectionGraphs(int &numGraphs);

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

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		bool onDeckBlock();

		/*! \brief This method computes the.
		*
		* 	\param		[in]	printLayerDepth				- input print layer depth.
		*	\since version 0.0.4
		*/
		void computePrintBlocks(float printPlaneSpacing, float printLayerWidth , float raftLayerWidth, bool allSDFLayers , int & numSDFlayers, int funcNum = 0, zDomainFloat neopreneOffset = zDomainFloat(0,0), bool compFrames = true, bool compSDF = true);

		/*! \brief This method computes the medial graph from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial graph.
		*	\param		[in]	endVID				- input end vertex id of medial graph.
		*	\since version 0.0.4
		*/
		void computeMedialGraph(zObjMesh& o_Mesh, int startVID, int endVID);

		/*! \brief This method computes the medial graph from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial spine.
		*	\param		[in]	endVID				- input end vertex id of medial spine.
		*	\param		[in]	blockStride			- input stride of edges for left and right blocks.
		*	\param		[in]	braceStride			- input stride for edges for braces.
		*	\since version 0.0.4
		*/
		void computeMedial_BraceEdges(zObjMesh& o_Mesh,int startVID, int endVID, int blockStride, int braceStride);
	
		/*! \brief This method computes the left or right slice mesh from input mesh.
		*
		* 	\param		[in]	o_Mesh				- input guide mesh object.
		*	\param		[in]	startVID			- input start vertex id of medial spine.
		*	\param		[in]	endVID				- input end vertex id of medial spine.
		*	\param		[in]	blockStride			- input stride of edges for left and right blocks.
		* 	\param		[in]	left				- input boolean indicating if the planes for the left or right side meshes.
		*	\since version 0.0.4
		*/
		void computeSliceMesh(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, bool left);


		//--------------------------
		//---- UTILITY METHODS
		//--------------------------

		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	printLayerDepth				- input print layer depth.
		*	\param		[in]	guideMesh_vertex			- input guide mesh vertex.
		*	\since version 0.0.4
		*/
		void computePrintBlockFrames( float printPlaneSpacing, float neopreneOffset_start , float neopreneOffset_end, bool left);

		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void computePrintBlockSections(bool left);
		
		
		/*! \brief This method computes the SDF for the blocks.
		*
		*	\since version 0.0.4
		*/
		void computeSDF(bool allSDFLayers, int& numSDFlayers, int funcNum, float printWidth, float neopreneOffset, float raftWidth);

		/*! \brief This method compute the block frames.
		*
		*	\param		[out]	outGraph						- output trim graph.
		*	\since version 0.0.4
		*/
		void computePrintBlockTrimGraphs(zObjGraph& inPolyObj, zObjGraph &o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE );

		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		bool checkPrintLayerHeights();

		/*! \brief This method compute the block frames for thickned mesh.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void computePrintBlock_bounds();

		/*! \brief This method compute the legth of  medial graph of the block
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	leftBlock					- input booealn indicating if its a left block or a right
		*	\since version 0.0.4
		*/
		void computePrintBlockLength( );

		/*! \brief This method compute the block SDF for the deck.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Deck(int funcNum, int graphId, bool alternate, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Balustrade(int funcNum, int graphId, bool alternate, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);


		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Boundary( int graphId, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);
			
		string format_number(int num, size_t size = 3, char fillChar = '0');

		zIntArray shiftArray(int startID, int length);

		int closestIndex(zPointArray pos, zPoint p);

		bool exportJSON(string pathCurrent, string dir, string filename, float printLyerWidth, float raftLayerWidth);


		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
		protected:

		/*! \brief This method compute the transform from input Vectors.
		*
		*	\param		[in]	O							- input origin point.
		*	\param		[in]	X							- input X axis vector.
		* 	\param		[in]	Y							- input Y axis vector.
		*	\param		[in]	Z							- input Z axis vector.
		*	\return				zTransform					- output transform.
		*	\since version 0.0.4
		*/
		zTransform setTransformFromVectors(zPoint& O, zVector& X, zVector& Y, zVector& Z);

		/*! \brief This method compute the transform from input Vectors.
		*
		* 	\param		[in]	O							- input origin point.
		*	\param		[in]	Z							- input Z axis vector.
		*	\param		[in]	Basis						- input Basis vector.
		*	\return				zTransform					- output transform.
		*	\since version 0.0.4
		*/
		zTransform setTransformFromOrigin_Normal(zPoint& O, zVector& Z, zVector Basis = zVector(0,1,0));
		
		zItMeshHalfEdge getStartHalfEdge(zObjMesh& o_mesh, int startVID, int endVID);

		void polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength);

		void getScalars_3dp_slot(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset );

		void getScalars_3dp_brace(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset , bool alternate);

		void getScalars_3dp_trim(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate);

		void addVertexToPositionMap(unordered_map<string, int> &positionVertex, zPoint& pos, int index, int precisionfactor = 3);

		bool vertexExistsinPositionMap(unordered_map<string, int>& positionVertex, zPoint &pos, int& outVertexId, int precisionfactor = 3);

		bool readJSON(string path, json& j);

		bool fileExists(string& path);

		//Post Processing Methods

		void GCodeGenerator();

		void seamPrint();

		zIntArray getGraphStartSequence();

		void generateOpening();
		

	};

}





namespace zSpace
{
	struct ZSPACE_TOOLSETS EXT_zGraph
	{
		/*!	\brief stores number of vertices */
		int vCount;
		/*!	\brief stores number of  edges */
		int eCount;
		/*!	\brief vertex container as 1D array (num of vertices x 3) */
		float* vPositions;
		/*!	\brief vertex colors container as 1D array (num of vertices/colors x 4) */
		float* vColors;
		/*!	\brief edge vertices container as 1D array (num of edges x 2) */
		int* ePairs;
		/*!	\brief edge colors container as 1D array (num of edges/colors x 4) */
		float* eColors;
		
	};
	struct ZSPACE_TOOLSETS EXT_Test
	{
		/*!	\brief stores number of  edges */
		int eCount;
		/*!	\brief edge vertices container as 1D array (num of edges x 2) */
		int* ePairs;
	};


	ZSPACE_TOOLSETS_EXT{
		ZSPACE_TOOLSETS void EXT_SetFromJSON(zTsSDFSlicer*& slicer, char* path, int count, int blockID, bool& isDeckBlock);
		ZSPACE_TOOLSETS void EXT_SetFromJSON2(zTsSDFSlicer*& slicer, char* path, int blockStride, int braceStride);

		ZSPACE_TOOLSETS void EXT_createFieldMesh(zTsSDFSlicer*& slicer, float* domainMin, float* domainMax, int resX, int resY);
		ZSPACE_TOOLSETS int EXT_checkLayerHeight(zTsSDFSlicer*& slicer, bool& check);
		
		
		ZSPACE_TOOLSETS int EXT_computePrintBlocks(zTsSDFSlicer*& slicer, float printPlaneSpacing, float printLayerWidth, float raftLayerWidth, bool allSDFLayer, int& computeGraphCount, int funcNum,  float neopreneOffsetMin, float neopreneOffsetMax, bool compFrames, bool compSDF);
		
		//Get Methods
		ZSPACE_TOOLSETS int EXT_getRawMesh(zTsSDFSlicer* slicer, zObjMesh*& objMesh, bool rightMesh);
		ZSPACE_TOOLSETS int EXT_getRawFieldMesh(zTsSDFSlicer* slicer, zObjMesh*& objMesh);
		ZSPACE_TOOLSETS int EXT_getBlockSectionGraphs(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int EXT_getBlockContourGraphs(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int EXT_getRawMedialGraph(zTsSDFSlicer* slicer, zObjGraph*& graph);
		ZSPACE_TOOLSETS int EXT_getBlockFrames(zTsSDFSlicer* slicer, vector<zTransform>*& planes, int& Count);
		ZSPACE_TOOLSETS int EXT_getTrimGraph(zTsSDFSlicer* slicer, zObjGraphArray*& graphs, int& Count);

		//Plane Data
		//ZSPACE_TOOLSETS void EXT_getPlanesData(vector<zTransform>* graph, float* outOrigin, float* outNormal, float* outXAxis, float* outYAxis);
		ZSPACE_TOOLSETS void EXT_getPlanesData(vector<zTransform>* graph, float* matrix);


		//Graph Data
		ZSPACE_TOOLSETS void EXT_getGraphsSetFromPointersVector(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray);
		ZSPACE_TOOLSETS void EXT_getGraphsSetFromVector(zObjGraphArray* graphs, zObjGraph** outGraphArray);

		ZSPACE_TOOLSETS void EXT_getGraphCounts(zObjGraph* graph, int& outvCount, int& outeCount);
		ZSPACE_TOOLSETS void EXT_getGraphData(zObjGraph* graph, float* outVPostions, float* outvColors, int* outePair, float* outeColors);
		
		//Mesh Data
		ZSPACE_TOOLSETS void EXT_getMeshCounts(zObjMesh* objMesh, int& out_vCount, int& out_fCount);
		ZSPACE_TOOLSETS void EXT_getMeshPosition(zObjMesh* objMesh, float* outVPostions, float* outVColors);
		ZSPACE_TOOLSETS void EXT_getMeshFaceCount(zObjMesh* objMesh, int* outfCounts);
		ZSPACE_TOOLSETS void EXT_getMeshFaceConnect(zObjMesh* objMesh, int* outfConnects);

		//Export JSON
		ZSPACE_TOOLSETS void EXT_ExportJSON(zTsSDFSlicer* slicer, char* fileCurrent, char* fileNew, char* fileName, float printLayerWidth, float raftLayerWidth);
		ZSPACE_TOOLSETS void EXT_ExportJSON2(zTsSDFSlicer* slicer, char* fileCurrent, int fileCount, char* fileNew, int newCount, char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth);


		//Check
		

		//Iterate
		ZSPACE_TOOLSETS void EXT_CheckFolder(char* folderDirectoryChar, float printPlaneSpace, float printLayerWidth, float raftLayerWidth, int SDFFunc_Num, float neopreneOffsetMin, float neopreneOffsetMax);


		//Draw Mesh
		ZSPACE_TOOLSETS void EXT_DrawMesh(zObjMesh* objMesh);






		//ZSPACE_TOOLSETS void EXT_testSafeArray(SAFEARRAY* arr);

		
	}
}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSDFSlicer.cpp>
#endif

#endif