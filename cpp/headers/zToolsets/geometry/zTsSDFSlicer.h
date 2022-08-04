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

#include <headers/zCore/base/zExtern.h>

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zInterface/functionsets/zFnMeshField.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
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

	class ZSPACE_TOOLS zTsSDFSlicer
	{
	protected:
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

		zDomainFloat actualPrintHeightDomain;

		zDomainFloat neopreneOffset;
				
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

		int numMagentaLoops = 1;

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
		void computePrintBlocks(zDomainFloat &_printHeightDomain, float printLayerWidth , float raftLayerWidth, bool allSDFLayers , int & numSDFlayers, int funcNum = 0, int numSmooth = 0, zDomainFloat _neopreneOffset = zDomainFloat(0,0),  bool compFrames = true, bool compSDF = true);

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
		void computeSDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth, float neopreneOffset, float raftWidth);

		/*! \brief This method compute the block frames.
		*
		*	\param		[out]	outGraph						- output trim graph.
		*	\since version 0.0.4
		*/
		void computePrintBlockTrimGraphs(zObjGraph& inPolyObj, zObjGraph &o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE );

		/*! \brief This method checks the layer heights of the block
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		bool checkPrintLayerHeights(bool& checkGeometry);

		/*! \brief This method check the block interfaces are planar.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		bool checkInterfacePoints(bool left);

		/*! \brief This method checks the layer heights for all the blocks in the input directory.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void checkPrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset);

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
		void computeBlockSDF_Deck(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Balustrade(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);


		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Boundary( int graphId, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);
			
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

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSDFSlicer.cpp>
#endif

#endif