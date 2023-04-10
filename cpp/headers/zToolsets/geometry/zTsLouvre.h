// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
//

#ifndef ZSPACE_TS_GEOMETRY_LOUVRE_H
#define ZSPACE_TS_GEOMETRY_LOUVRE_H



#pragma once

#include <headers/base/zSpace_Toolsets.h>

#include <headers/zCore/base/zExtern.h>


#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zInterface/functionsets/zFnMeshField.h>

using namespace std;

namespace zSpace
{

	class ZSPACE_TOOLSETS zLouvre
	{
	public:

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;

		zObjMesh louvre;

		zTransform bestFitPlane;

		zLouvre();

		~zLouvre();

		double rotation;

		double height;

		double length;




	};


	class ZSPACE_TOOLSETS zTsLouvre
	{
	protected:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;

		/*!	\brief input guide mesh object  */
		zObjMesh o_GuideMesh;

		/*!	\brief topology mesh object  */
		zObjMesh o_topologyMesh;

		/*!	\brief rotation mesh object  */
		zObjMesh o_rotationGuideMesh;

		/*!	\brief height mesh object  */
		zObjMesh o_heightGuideMesh;

		/*!	\brief topology mesh functionset  */
		zFnMesh fn_topologyMesh;

		/*!	\brief container of graph for each louvre object  */
		zObjGraph o_TopologyGraph;



	public:

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsLouvre();


		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsLouvre();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------


		/*! \brief This method creates topology mesh from graph.
		*
		*	\param		[in]	graphs			- input domain of bounds.
		*	\since version 0.0.4
		*/
		void createTopologyMesh(string path, zObjGraph graphs);


		//--------------------------
		//--- SET METHODS 
		//--------------------------


		void setTopologyMeshFromOBJ(string path);
		void setRotationMeshFromOBJ(string path);

		



	};
}





#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/geometry/zTsLouvre.cpp>
#endif

#endif

