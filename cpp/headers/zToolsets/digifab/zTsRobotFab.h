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
				

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------
		using zTsRobot::zTsRobot;


		zTsRHWC(zTsRobot& obj);
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
		

		//--------------------------
		//---- CREATE METHODS
		//--------------------------

		//--------------------------
		//---- OVERRIDE METHODS
		//--------------------------

		/*! \brief This method creates the targets from the fabrication mesh.
		*
		*	\since version 0.0.2
		*/
		void computeTargets();

		void computeGcode();

		void FramesFromTarget(zTransform& target, zVector& _position, zVector& _rotationX, zVector& _rotationY, zVector _rotationZ);
		void getOVAngle(zTransform& target, float& angle, zObjMesh& cutMesh);

		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:

		void computeTargetsOnStrip(zObjMesh& cutMesh, vector<zTransform>& targets_strip, vector<float>& ov_angle_strip);
		void computeTargetsOnStripPerpToMesh(zObjMesh& cutMesh, vector<zTransform>& targets_strip, vector<float>& ov_angle_strip);

		void checkTargetNormal(vector<zTransform>& targets_strip);

		zTransform targetFromFrames(zVector& _position, zVector& _rotationX, zVector& _rotationY, zVector _rotationZ);

		void addSafeTargets(vector<zTransform>& targets_strip, vector<zTransform>& targets_stripPerp, vector<float>& ov_angle_strip, float multiplication, bool perp);
		void addSafeTargetsStart(vector<zTransform>& targets_strip, float multiplication, bool perp);

		vector<zTransform> computeSafeTargets(zTransform& target, float multiplication);
		vector<zTransform> computeSafeTargetsPerp(zTransform& target, float multiplication);
		
	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/digiFab/zTsRobotFab.cpp>
#endif

#endif