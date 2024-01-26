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

#ifndef ZSPACE_TS_GEOMETRY_SPATIALGRAPH_H
#define ZSPACE_TS_GEOMETRY_SPATIALGRAPH_H



#pragma once

#include <zCore/base/zExtern.h>

#include <zInterface/functionsets/zFnMesh.h>
#include <zInterface/functionsets/zFnGraph.h>
#include <zInterface/functionsets/zFnParticle.h>

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

	
	/*! \class zTsSpatialGraph
	*	\brief A toolset for Spatial Graphs.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLS zTsSpatialGraph
	{
	protected:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;			

		/*!	\brief input guide mesh object  */
		zObjGraph o_SpatialGraph;

		/*!	\brief left mesh object  */
		zObjMeshArray o_ConvexHulls;	
		
		

	public:
		

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsSpatialGraph();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsSpatialGraph();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------


		//--------------------------
		//--- SET METHODS 
		//--------------------------

		/*! \brief This method creates the spatial graph from input file - JSON or OBJ.
		*
		*	\param		[in]	path		- input file.
		*	\param		[in]	fileType	- input file type.
		*	\since version 0.0.4
		*/
		void setFromFile(string path, zFileTpye fileType);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method gets the raw spatial graph object.
		*
		*	\return				zObjGraph	    - pointer to spatial graph object.
		*	\since version 0.0.4
		*/
		zObjGraph* getRawSpatialGraph();			
		
		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method gets the raw spatial graph object.
		*
		*	\return				zObjGraph	    - pointer to spatial graph object.
		*	\since version 0.0.4
		*/
		void computeConvexHulls(float edgeFactor0, float edgeFactor1, int profileVerts, float profileEdgeLen);


		//--------------------------
		//---- UTILITY METHODS
		//--------------------------

		

		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
		protected:

			zPointArray computeProfilePoints(int numVerts, float edgeLen);
		

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSpatialGraph.cpp>
#endif

#endif