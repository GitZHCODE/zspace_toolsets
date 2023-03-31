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

#ifndef ZSPACE_TS_ROBOT_FAB_H
#define ZSPACE_TS_ROBOT_FAB_H


#pragma once

#include <headers/zCore/base/zExtern.h>
#include<headers/zToolsets/digiFab/zTsRobot.h>
#include <headers/base/zSpace_Toolsets.h>


#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

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

	/*! \class zTsRobotFab
	*	\brief A toolset for panels analysis on freeform mesh.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLSETS zTsRHWC : public zTsRobot
	{
	protected:


	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------
		
		/*!	\brief form function set  */

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------
		/*! \brief Default constructor.
		*
		*	\since version 0.0.2
		*/
		zTsRHWC();

		zTsRHWC(zObjGraph& _jointGraphObj, vector<zObjMesh>& _jointMeshObjs);


		//--------------------------
		//---- DESTRUCTOR
		//--------------------------
		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsRHWC();

		//--------------------------
        //--- SET METHODS 
        //--------------------------


		//--------------------------
		//--- GET METHODS 
		//--------------------------
		
		/*! \brief This method returns the fab mesh bounding box.
		*
		*	\since version 0.0.4
		*/
		zPointArray getFabBbox();

		/*! \brief This method returns the fab mesh.
		*
		*	\since version 0.0.4
		*/
		vector<zObjMesh> getFabMesh();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the robot targets of input mesh array.
		*
		*	\since version 0.0.4
		*/
		void computeTargets() override;


		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:
	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/digiFab/zTsRobotFab.cpp>
#endif

#endif