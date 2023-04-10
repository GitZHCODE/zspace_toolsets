//// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
//// data analysis & visualization framework.
////
//// Copyright (C) 2019 ZSPACE 
//// 
//// This Source Code Form is subject to the terms of the MIT License 
//// If a copy of the MIT License was not distributed with this file, You can 
//// obtain one at https://opensource.org/licenses/MIT.
////
//// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
////
//
//#ifndef ZSPACE_EXT_TS_ROBOT_H
//#define ZSPACE_EXT_TS_ROBOT_H
//
//
//
//#pragma once
//#include <headers/base/zSpace_Toolsets.h>
//
//#include <headers/zCore/base/zExtern.h>
//#include <headers/zInterface/functionsets/zFnMesh.h>
//#include <headers/zInterface/functionsets/zFnGraph.h>
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <iostream>
//#include <sstream>
//
//#include<execution>
//
////#include<headers/zToolsets/digifab/zTsRobot.h>
//#include<headers/zToolsets/digifab/zTsRobotFab.h>
//#include<headers/zToolsets/externalMethods/zExtTransform.h>
//#include<headers/zToolsets/externalMethods/zExtMesh.h>
//
//using namespace std;
//
//
//namespace zSpace
//{
//
//	//zTsRobot ROBOTGLOBAL;
//
//	//ZSPACE_TOOLSETS_EXT
//	//{
//
//
//	//	ZSPACE_TOOLSETS void ext_zTsRobot_createFromFile(char* path);
//	//	ZSPACE_TOOLSETS void ext_zTsRobot_setEndEffector(float* matrix);
//	// 
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_setRobotTarget(float* matrix);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_inverseKinematics(char* path, float* outJointRotation, float* outJointTransformation ); // out Joint float* 6*16*x
//	// 
//	//	ZSPACE_TOOLSETS void ext_zTsRobot_forwardKinematics(float* jointRotation, float* outJointTransformation); //<<in float[6],  out robot joint transformation 6*16
//
//
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_createRobotJointMeshesfromFile(zExtRobot& extRobot, char* directory, bool endeffector = false);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_setJointMeshTransform(zExtRobot& extRobot,  bool updatePositions = true);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_setEndEffector(zExtRobot& extRobot);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_setTarget(zExtRobot& extRobot, float* target);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_updateKinematics(zExtRobot& extRobot, bool inverseKinematics);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_forwardKinematics(zExtRobot& extRobot);
//	//	//ZSPACE_TOOLSETS void ext_zTsRobot_inverseKinematics(zExtRobot& extRobot, float* target);
//
//	//}
//
//
//	struct zExtRobot
//	{
//		zTsRHWC* robot;
//		float robotJointTransforms[96];
//		float robotJointRotation[6]; 
//		float robotJointRotationMax[6];
//		float robotJointRotationMin[6];
//		float robotJointRotationHome[6];
//		float robotJointRotationMask[6];
//		float robotJointRotationPulse[6];
//		float robotJointRotationOffset[6];
//
//		/////*!	\brief robot base matrix  */
//		////float robot_base_matrix[16];
//		/////*!	\brief robot target matrix  */
//		////float robot_target_matrix[16];
//		/////*!	\brief robot end effector matrix  */
//		////float robot_endEffector_matrix[16];
//
//		zExtRobot(zTsRHWC* r);
//		void updateAttributes();
//	};
//
//
//	ZSPACE_TOOLSETS_EXT
//	{
//		
//
//		//Read Robot JSON file <MUST-Step 0>
//		ZSPACE_TOOLSETS void ext_zTsRobot_createFromFile(zExtRobot & extRobot, char* path);
//		
//		//Read Robot JSON file <MUST-Step 1>
//		ZSPACE_TOOLSETS void ext_zTsRobot_setEndEffector(zExtRobot& extRobot, float* eeTransform);
//
//		//FK - Read joint rotation <optionA>
//		ZSPACE_TOOLSETS void ext_zTsRobot_forwardKinematics(zExtRobot& extRobot, float* jointRotation);
//
//		//IK - Set target <optionB>
//		ZSPACE_TOOLSETS void ext_zTsRobot_inverseKinematics(zExtRobot& extRobot, float* targetTransformation);
//
//
//		//Set fabrication mesh, fab Base PLane, work base Plane, and compute targets
//		ZSPACE_TOOLSETS void ext_zTsRobot_setFabricationMesh(zExtRobot& extRobot, zExtMesh* fabMesh, int fabMeshCount, float* fabBasePlane, float* workBasePlane,   int outTargetCount);
//
//		//Get computed targets
//		ZSPACE_TOOLSETS void ext_zTsRobot_getComputedTargets(zExtRobot & extRobot, float* outTargets);
//
//
//
//
//
//
//		ZSPACE_TOOLSETS void ext_zTsRobot_Test1(zExtRobot& extRobot);
//		ZSPACE_TOOLSETS void ext_zTsRobot_Test2(char* path);
//		ZSPACE_TOOLSETS void ext_zTsRobot_setDirectory(zExtRobot& extRobot, char* path);
//
//
//	}
//
//}
//
//
//
//
//#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
//// All defined OK so do nothing
//#else
//#include<source/zToolsets/externalMethods/zExtRobot.cpp>
//#endif
//
//#endif