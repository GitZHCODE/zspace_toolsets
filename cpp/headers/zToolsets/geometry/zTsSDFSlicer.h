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

		/*!	\brief core utilities object  */
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

		/*!	\brief minimum layer height in the print */
		float minLayerHeight = 0;

		/*!	\brief maximum layer height in the print  */
		float maxLayerHeight = 0;

		/*!	\brief total length of print  */
		float totalLength = 0;

		/*!	\brief critical minimum layer points  */
		zObjPointCloud criticalMinLayer_pts;

		/*!	\brief critical maximum layer points  */
		zObjPointCloud criticalMaxLayer_pts;

		/*!	\brief print height domain  */
		zDomainFloat printHeightDomain;

		/*!	\brief actual print height domain  */
		zDomainFloat actualPrintHeightDomain;

		/*!	\brief neoprene material offset  */
		zDomainFloat neopreneOffset;
				
		//--------------------------
		//---- SDF ATTRIBUTES
		//--------------------------

		/*!	\brief scalar field  */
		zObjMeshScalarField o_field;

		//--------------------------
		//---- COLOR ATTRIBUTES
		//--------------------------
		
		zColor red, yellow, green, cyan, blue, magenta, grey , orange;

		zColorArray blockColors;

		bool leftPlaneExists, rightPlaneExists;

		/*!	\brief current block ID */
		int blockId = 1;

		/*!	\brief number of magenta loops*/
		int numMagentaLoops = 1;

		/*!	\brief transformation for the world and local position */
		zTransform base_world, base_local;

		/*!	\brief boolean indicates if current block is deck or balustrade */
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

		/*! \brief This method sets mesh from json.
		*
		*	\param		[in]	dir		            - input JSON file directory
		*	\param		[in]	blockStride			- input stride of edges for left and right blocks.
		*	\param		[in]	braceStride			- input stride for edges for braces.
		*	\since version 0.0.4
		*/
		void setFromJSON(string dir, int blockStride, int braceStride);

		/*! \brief This method sets mesh from json.
		*
		*	\param		[in]	dir		             - input JSON file directory.
		**	\param		[in]	_blockID	         - input block ID.
		*	\since version 0.0.4
		*/
		void setFromJSON(string dir, int _blockID);

		/*! \brief This method sets the left or right slice mesh object.
		*
		*	\param		[in]	_o_SliceMesh		- input mesh object.
		* 	\param		[in]	left				- input boolean indicating if it is the left or right side mesh.
		*	\since version 0.0.4
		*/
		void setSliceMesh(zObjMesh& _o_SliceMesh, bool left);

		/*! \brief This method sets the medial graph object.
		*
		*	\param		[in]	_o_MedialGraph		- input graph object.
		*	\since version 0.0.4
		*/
		void setMedialGraph(zObjGraph& _o_MedialGraph);

		/*! \brief This method sets the start and end plane.
		*
		*	\param		[in]	_sPlane		    	- input start plane.
		*	\param		[in]	_ePlane		    	- input end plane.
		* 	\param		[in]	left			    - input boolean indicating if the planes are for the left or right side meshes.
		*	\since version 0.0.4
		*/
		void setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane, bool left);

		/*! \brief This method transforms geometry.
		*
		*	\param		[in]	toLocal		    	- input boolean indicationg transformation to local or world position.
		*	\since version 0.0.4
		*/
		void setTransforms(bool toLocal);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method gets the block start and end.
		*
		*	\param		[in]	left			  - input boolean indicating if it is the left or right side mesh.
		*	\since version 0.0.2
		*/
		zTransform* getRawBlockStartEnd(bool left);

		
		/*! \brief This method gets the block frames.
		*
		*	\return		    vector<zTransform>	    - container of transforms if they exist.
		*	\since version 0.0.2
		*/
		vector<zTransform> getBlockFrames();

		/*! \brief This method gets the block section graphs.
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer container of section graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockSectionGraphs(int& numGraphs);

		/*! \brief This method gets the block raft graphs.
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer container of raft graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockRaftGraphs(int& numGraphs);

		/*! \brief This method gets the block SDF contour graphs.
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer container of contour graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockContourGraphs(int& numGraphs);

		/*! \brief This method gets the block trim graphs
		*
		*	\param		[out]	numGraphs				- output number of graphs.
		*	\return				zObjGraphPointerArray	- pointer container of graphs if they exist.
		*	\since version 0.0.2
		*/
		zObjGraphPointerArray getBlockTrimGraphs(int& numGraphs);

		/*! \brief This method gets the critical points of the section graphs.
		*
		*	\param		[in]	minHeight				- input boolean indicating if it is min or max critical points.
		*	\return				zObjPointCloud*			- pointer to container of points if they exist.
		*	\since version 0.0.2
		*/
		zObjPointCloud* getRawCriticalPoints(bool minHeight);

		/*! \brief This method gets pointer to the internal field object.
		*
		*	\return			zObjMeshScalarField*		 - pointer to internal field object.
		*	\since version 0.0.4
		*/
		zObjMeshScalarField* getRawFieldMesh();

		/*! \brief This method gets pointer to the internal medial graph object.
		*
		*	\return				zObjGraph*		- pointer to internal graph object.
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

		/*! \brief This method checks if current block is deck or balustrade.
		*
		*	\return				bool					- boolean that indicates if the block is on deck or balustrade.
		*	\since version 0.0.4
		*/
		bool onDeckBlock();

		/*! \brief This method computes the block print planes and SDF.
		*
		* 	\param		[in]	_printHeightDomain			- input print height domain.
		*  	\param		[in]	printLayerWidth				- input print layer width.
		*  	\param		[in]	raftLayerWidth				- input raft layer width.
		*  	\param		[in]	allSDFLayers				- input boolean indicating if all SDF layers will be calculated.
		*  	\param		[in]	numSDFlayers				- input number of SDF layers to compute if not all.
		*  	\param		[in]	funcNum				        - input number of functions to execute.
		* 	\param		[in]	numSmooth				    - input smooth factor.
		*  	\param		[in]	_neopreneOffset				- input neoprene material offset.
		*  	\param		[in]	compFrames				    - input boolean indicating compute frames calcualation.
		*  	\param		[in]	compSDF				        - input boolean indicating compute SDF calcualation.
		*	\since version 0.0.4
		*/
		void computePrintBlocks(zDomainFloat &_printHeightDomain, float printLayerWidth , float raftLayerWidth, bool allSDFLayers , int & numSDFlayers, int funcNum = 0, int numSmooth = 0, zDomainFloat _neopreneOffset = zDomainFloat(0,0),  bool compFrames = true, bool compSDF = true);

		/*! \brief This method computes the block print planes and SDF.
		*
		* 	\param		[in]	_printHeightDomain			- input print height domain.
		*  	\param		[in]	printLayerWidth				- input print layer width.
		*  	\param		[in]	allSDFLayers				- input boolean indicating if all SDF layers will be calculated.
		*  	\param		[in]	numSDFlayers				- input number of SDF layers to compute if not all.
		*  	\param		[in]	funcNum				        - input number of functions to execute.
		* 	\param		[in]	numSmooth				    - input smooth factor.
		*  	\param		[in]	_neopreneOffset				- input neoprene material offset.
		*  	\param		[in]	compFrames				    - input boolean indicating compute frames calcualation.
		*  	\param		[in]	compSDF				        - input boolean indicating compute SDF calcualation.
		*	\since version 0.0.4
		*/
		void computePrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum = 0, int numSmooth = 0, zDomainFloat _neopreneOffset = zDomainFloat(0, 0), bool compFrames = true, bool compSDF = true);

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
		* 	\param		[in]	left				- input boolean indicating if it computes the left or right slice mesh.
		*	\since version 0.0.4
		*/
		void computeSliceMesh(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, bool left);


		//--------------------------
		//---- UTILITY METHODS
		//--------------------------

		/*! \brief This method computes the print block frames.
		*
		*	\param		[in]	printPlaneSpacing			- input space between print planes.
		*	\param		[in]	neopreneOffset_start		- input neoprene material offset at the begining.
		*	\param		[in]	neopreneOffset_end			- input neoprene material offset at the end.
		*   \param		[in]	left            			- input boolean indicating if the planes are for the left or right side mesh.
		*	\since version 0.0.4
		*/
		void computePrintBlockFrames( float printPlaneSpacing, float neopreneOffset_start , float neopreneOffset_end, bool left);

		/*! \brief This method computes print block sections.
		*
		*	\param		[in]	left						- input boolean indicating if it computes left or right side.
		*	\since version 0.0.4
		*/
		void computePrintBlockSections(bool left);		
		
		/*! \brief This method computes the SDF for the blocks.
		*
		*	\param		[in]	allSDFLayers	                - input boolean indicating if all SDF layers will be computed.
		*	\param		[in]	numSDFlayers		            - input number of SDF layers to compute if not all.
		*	\param		[in]	funcNum			                - input number of functions to execute.
		*   \param		[in]	numSmooth            			- input smooth factor.
		*   \param		[in]	printWidth            			- input print width.
		*   \param		[in]	neopreneOffset            		- input neoprene material offset.
		*   \param		[in]	raftWidth            			- input raft width.
		*	\since version 0.0.4
		*/
		void computeSDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth, float neopreneOffset, float raftWidth);

		/*! \brief This method computes the SDF for the blocks.
		*
		*	\param		[in]	allSDFLayers	                - input boolean indicating if all SDF layers will be computed.
		*	\param		[in]	numSDFlayers		            - input number of SDF layers to compute if not all.
		*	\param		[in]	funcNum			                - input number of functions to execute.
		*   \param		[in]	numSmooth            			- input smooth factor.
		*   \param		[in]	printWidth            			- input print width.
		*   \param		[in]	neopreneOffset            		- input neoprene material offset.
		*	\since version 0.0.4
		*/
		void computeSDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth, float neopreneOffset);

		/*! \brief This method computes the print block trim graphs.
		*
		* 	\param		[in]	inPolyObj	                    - input polygon graph.
		* 	\param		[in]	topHE	                        - input top HE.
		* 	\param		[in]	bottomHE	                    - input bottom HE.
		*	\param		[out]	o_outGraph						- output trim graph.
		*	\since version 0.0.4
		*/
		void computePrintBlockTrimGraphs(zObjGraph& inPolyObj, zObjGraph& o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE );

		/*! \brief This method checks the layer heights and planarity at intersection areas of the block.
		* 
		*   \param		[in]	checkSDF				     	- input boolean indicating if SDF is correct.
		*	\param		[in]	checkGeometry					- input boolean indicating if geometry is correct (planar intersection areas).
		*   \return		        bool                            - boolean indicating if block geometry can be printed.
		*	\since version 0.0.4
		*/
		bool checkPrintBlockGeometry(bool& checkSDF, bool& checkGeometry);

		/*! \brief This method checks the layer heights and intersection planarity for all the blocks in the input directory and exports data in csv.
		*
		*   \param		[in]	folderDir						  - input folder directory for input and output files.
		*   \param		[in]	_printHeightDomain				  - input print height domain.
		*	\param		[in]	_neopreneOffset					  - input neoprene material offset domain.
		*	\since version 0.0.4
		*/
		void checkPrintBlocks_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset);

		/*! \brief This method checks if the block interfaces are planar.
		*
		*	\param		[in]	left					         - input boolean indicating if the planes for the left or right side meshes exist.
		*   \param		[in]	distTolerance					 - input distance tolerance.
		*   \return		        bool                             - boolean indicating if interface points are on the same plane.
		*	\since version 0.0.4
		*/
		bool checkInterfacePoints(bool left, float distTolerance);

		/*! \brief This method computes the block SDF for the deck.
		*
		*	\param		[in]	funcNum						- input number of functions to execute.
		*	\param		[in]	numSmooth					- input smooth factor.
		*   \param		[in]	graphId						- input current graph id.
		*	\param		[in]	alternate					- input boolean indicating the brace side.
		*   \param		[in]	printWidth					- input print width.
		*	\param		[in]	neopreneOffset				- input neoprene material offset.
		*   \param		[in]	addRaft						- input boolean indicating if raft exists.
		*	\param		[in]	raftId						- input raft id.
		*   \param		[in]	raftWidth					- input raft width.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Deck(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		/*! \brief This method computes the block SDF for the deck.
		*
		*	\param		[in]	funcNum						- input number of functions to execute.
		*	\param		[in]	numSmooth					- input smooth factor.
		*   \param		[in]	graphId						- input current graph id.
		*	\param		[in]	alternate					- input boolean indicating the brace side.
		*   \param		[in]	printWidth					- input print width.
		*	\param		[in]	neopreneOffset				- input neoprene material offset.  
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Deck(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		/*! \brief This method computes the block SDF for the deck with interpolation method.
		*
		*	\param		[in]	funcNum						- input number of functions to execute.
		*	\param		[in]	numSmooth					- input smooth number.
		*   \param		[in]	graphId						- input current graph id.
		*	\param		[in]	alternate					- input booealn indicating the brace side.
		*   \param		[in]	printWidth					- input print width.
		*	\param		[in]	neopreneOffset				- input neoprene material offset.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Deck_Interpolation(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		void computeBlockSDF_Deck_Interpolation_ss(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		/*! \brief This method computes the block SDF for the balustrade.
		*
		*	\param		[in]	funcNum						- input number of functions to execute.
		*	\param		[in]	numSmooth					- input smooth factor.
		*   \param		[in]	graphId						- input current graph id.
		*	\param		[in]	alternate					- input boolean indicating the brace side.
		*   \param		[in]	printWidth					- input print width.
		*	\param		[in]	neopreneOffset				- input neoprene material offset.
		*   \param		[in]	addRaft						- input boolean indicating if raft exists.
		*	\param		[in]	raftId						- input raft id.
		*   \param		[in]	raftWidth					- input raft width.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Balustrade(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth);

		/*! \brief This method computes the block SDF for the balustrade.
	*
	*	\param		[in]	funcNum						- input number of functions to execute.
	*	\param		[in]	numSmooth					- input smooth factor.
	*   \param		[in]	graphId						- input current graph id.
	*	\param		[in]	alternate					- input boolean indicating the brace side.
	*   \param		[in]	printWidth					- input print width.
	*	\param		[in]	neopreneOffset				- input neoprene material offset.
	*	\since version 0.0.4
	*/
		void computeBlockSDF_Balustrade(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		/*! \brief This method computes the block SDF for the balustrade with interpolation method.
	    *
    	*	\param		[in]	funcNum						- input number of functions to execute.
    	*	\param		[in]	numSmooth					- input smooth factor.
	    *   \param		[in]	graphId						- input current graph id.
    	*	\param		[in]	alternate					- input boolean indicating the brace side.
    	*   \param		[in]	printWidth					- input print width.
    	*	\param		[in]	neopreneOffset				- input neoprene material offset.
     	*	\since version 0.0.4
    	*/
		void computeBlockSDF_Balustrade_Interpolation(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		void computeBlockSDF_Balustrade_Interpolation_ss(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset);

		/*! \brief This method computes the block SDF for the boundary.
		*
		*	param		[in]	graphId						- input current graph id.
		*   \param		[in]	printWidth					- input print width.
		*	\param		[in]	neopreneOffset				- input neoprene material offset.
		*   \param		[in]	addRaft						- input boolean indicating if raft exists.
		*	\param		[in]	raftId						- input raft id.
		*   \param		[in]	raftWidth					- input raft width.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Boundary( int graphId, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);
			
		/*! \brief This method computes the start and end layer brace SDF.
		*
		*   \param		[in]	printWidth					    - input print width.
		*	\param		[in]	left				            - input boolean indicating if the planes for the left or right side meshes.
		*   \param		[in]	alternate						- input boolean indicating the brace side.
		*	\param		[in]	scalars_start					- input first layer scalar field.
		*   \param		[in]	scalars_end				    	- input last layer scalar field.
		*	\since version 0.0.4
		*/
		void computeBraceSDF_StartEnd(float printWidth, bool left, bool alternate, zFloatArray& scalars_start, zFloatArray& scalars_end);

		/*! \brief This method computes the start and end layer brace SDF.
		*
		*   \param		[in]	printWidth					    - input print width.
		*	\param		[in]	left				            - input boolean indicating if the planes for the left or right side meshes.
		*   \param		[in]	alternate						- input boolean indicating the brace side.
		*	\param		[in]	scalars_start					- input first layer scalar field.
		*   \param		[in]	scalars_end				    	- input last layer scalar field.
		*	\since version 0.0.4
		*/
		void computeBraceSDF_StartEnd_ss(float printWidth, bool left, bool alternate, zFloatArray& scalars_start, zFloatArray& scalars_end);

		/*! \brief This method computes the start and end layer trim SDF.
		*
		*   \param		[in]	printWidth					    - input print width.
		*	\param		[in]	left				            - input boolean indicating if the planes for the left or right side meshes.
		*   \param		[in]	alternate						- input boolean indicating the trim side.
		*	\param		[in]	scalars_start					- input first layer scalar field.
		*   \param		[in]	scalars_end				    	- input last layer scalar field.
		*	\since version 0.0.4
		*/
		void computeTrimSDF_StartEnd(float printWidth, bool left, bool alternate, zFloatArray& scalars_start, zFloatArray& scalars_end);

		/*! \brief This method computes the start and end layer trim SDF.
		*
		*   \param		[in]	printWidth					    - input print width.
		*	\param		[in]	left				            - input boolean indicating if the planes for the left or right side meshes.
		*   \param		[in]	alternate						- input boolean indicating the trim side.
		*	\param		[in]	scalars_start					- input first layer scalar field.
		*   \param		[in]	scalars_end				    	- input last layer scalar field.
		*	\since version 0.0.4
		*/
		void computeTrimSDF_StartEnd_ss(float printWidth, bool left, bool alternate, zFloatArray& scalars_start, zFloatArray& scalars_end);

		/*! \brief This method exports json file of a block geometry and slicing information.
		*
		*   \param		[in]	pathCurrent					        - input current path.
		*	\param		[in]	dir				                    - input folder directory.
		*   \param		[in]	filename						    - input export filename.
		*	\param		[in]	printLayerWidth					    - input print layer width.
		*   \param		[in]	raftLayerWidth				    	- input raft layer width.
		*	\since version 0.0.4
		*/
		bool exportJSON(string pathCurrent, string dir, string filename, float printLayerWidth, float raftLayerWidth);

		/*! \brief This method exports json file of a block geometry and slicing information.
         *
         *   \param		[in]	pathCurrent					        - input current path.
		 *	\param		[in]	dir				                    - input folder directory.
		 *   \param		[in]	filename						    - input export filename.
		 *	\param		[in]	printLayerWidth					    - input print layer width.
		 *	\since version 0.0.4
		 */
		bool exportJSON(string pathCurrent, string dir, string filename, float printLayerWidth);

		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
		protected:

		/*! \brief This method
		*
		*	\param		[in]	o_mesh				    - input mesh object.
		*   \param		[in]	startVID                - input start vetrex Id.
		*   \param		[in]	endVID                  - input end vetrex Id.
		*	\return				zItMeshHalfEdge     	- start half edge mesh.
		*	\since version 0.0.4
		*/
		zItMeshHalfEdge getStartHalfEdge(zObjMesh& o_mesh, int startVID, int endVID);

		/*! \brief This method computes the top and bottom edge graphs.
		*
		*	\param		[in]	inPoly				       - input graph object.
		*   \param		[out]	topHE                      - output top half edge graph.
		*   \param		[out]	bottomHE                   - output bottom half edge graph.
		*   \param		[out]	topLength                  - output top graph length.
		*   \param		[out]	bottomLength               - output bottom graph length.
		* 
		*	\since version 0.0.4
		*/
		void polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength);

		/*! \brief This method computes the slot graph.
		*
		*	\param		[in]	scalars				      - input scalars array.
		*   \param		[out]	o_trimGraph               - output trim graph object.
		*   \param		[in]	offset                    - input offset number.
		*
		*	\since version 0.0.4
		*/
		void getScalars_3dp_slot(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset);

		/*! \brief This method computes the brace graph.
		*
		*	\param		[in]	scalars				      - input scalars array.
		*   \param		[out]	o_trimGraph               - output trim graph object.
		*   \param		[in]	outer_printWidth          - input the outer print width.
		*   \param		[in]	offset                    - input offset number.
		*   \param		[in]	alternate                 - input boolean indicating the brace side.
		*
		*	\since version 0.0.4
		*/
		void getScalars_3dp_brace(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset , bool alternate);

		/*! \brief This method computes the trim graph.
		*
		*	\param		[in]	scalars				         - input scalars array.
		*   \param		[out]	o_trimGraph                  - output trim graph object.
		*   \param		[in]	offset                       - input offset number.
		* 	\param		[in]	alternate                    - input boolean indicating the trim side.
		*
		*	\since version 0.0.4
		*/
		void getScalars_3dp_trim(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate);

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSDFSlicer.cpp>
#endif

#endif