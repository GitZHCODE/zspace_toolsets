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

		/*!	\brief pointer to input mesh Object  */
		zObjMesh *o_SliceMesh;			

		/*!	\brief pointer to input mesh Object  */
		zObjGraph* o_GuideGraph;

		/*!	\brief start print plane  */
		zTransform startPlane;

		/*!	\brief end print plane  */
		zTransform endPlane;

		/*!	\brief container of section frames  */
		vector<zTransform> sectionFrames;

		/*!	\brief container of section graph objects  */
		zObjGraphArray o_sectionGraphs;

		/*!	\brief container of contour graph objects  */
		zObjGraphArray o_contourGraphs;

		/*!	\brief container of raft graph objects  */
		zObjGraphArray o_raftGraphs;

		//--------------------------
		//---- PRINT ATTRIBUTES
		//--------------------------

		/*!	\brief minimum layer height  in the print */
		float minLayerHeight = 0;

		/*!	\brief maximum layer height in the print  */
		float maxLayerHeight = 0;

		/*!	\brief total length of print  */
		float totalLength = 0;
				
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

		/*! \brief This method creates the filed mesh.
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

		/*! \brief This method sets the slice mesh object.
		*
		*	\param		[in]	_o_SliceMesh			- input slice mesh object.
		*	\since version 0.0.4
		*/
		void setSliceMesh(zObjMesh& _o_SliceMesh);

		/*! \brief This method sets the slice mesh object.
		*
		*	\param		[in]	_o_GuideGraph			- input guide graph object.
		*	\since version 0.0.4
		*/
		void setGuideGraph(zObjGraph& _o_GuideGraph);
		

		/*! \brief This method sets the start and end plane.
		*
		*	\param		[in]	_sTransform			- input start plane.
		*	\param		[in]	_eTransform			- input end plane.
		*	\since version 0.0.4
		*/
		void setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane);

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

		/*! \brief This method gets pointer to the internal field object.
		*
		*	\return				zObjMeshScalarField*					- pointer to internal field object.
		*	\since version 0.0.4
		*/
		zObjMeshScalarField* getRawFieldMesh();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the.
		*
		* 	\param		[in]	printLayerDepth				- input print layer depth.
		*	\since version 0.0.4
		*/
		void computePrintBlocks(float printPlaneSpacing, float printLayerWidth , float raftLayerWidth, zDomainFloat neopreneOffset = zDomainFloat(0,0), bool compFrames = true, bool compSDF = true);

	
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
		void computePrintBlockFrames( float printPlaneSpacing, float neopreneOffset_start , float neopreneOffset_end );

		/*! \brief This method compute the block frames.
		*
		*	\param		[in]	_block						- input block.
		*	\since version 0.0.4
		*/
		void computePrintBlockSections();	
		
		
		/*! \brief This method computes the SDF for the blocks.
		*
		*	\since version 0.0.4
		*/
		void computeSDF(float printWidth, float neopreneOffset, float raftWidth);



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

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Internal( int graphId, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);

		/*! \brief This method compute the block SDF for the balustrade.
		*
		*	\param		[in]	_block						- input block.
		*	\param		[in]	graphId						- input index of section graph.
		*	\since version 0.0.4
		*/
		void computeBlockSDF_Boundary( int graphId, float printWidth = 0.020, float neopreneOffset = 0.005, bool addRaft = false, int raftId = 0, float raftWidth = 0.030);
			
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
		

		//WIP
		void CreatePaths(JSON Path, int FieldResolution[2], float FieldBounds[6], float* outLeftPoints, float* outRightPoints, int* leftEdgeCOnnects, int* rightEdgeConnects);



	};


	extern "C" void extCreatePaths(JSON Path, int FieldResolution[2], float FieldBounds[6], float* outLeftPoints, float* outRightPoints, int* leftEdgeCOnnects, int* rightEdgeConnects)


}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSDFSlicer.cpp>
#endif

#endif