// //This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// //data analysis & visualization framework.
//
// //Copyright (C) 2019 ZSPACE 
// //
// //This Source Code Form is subject to the terms of the MIT License 
// //If a copy of the MIT License was not distributed with this file, You can 
// //obtain one at https://opensource.org/licenses/MIT.
//
// //Author : Heba Eiz <heba.eiz@zaha-hadid.com>
//
//
//
//#include<headers/zToolsets/externalMethods/zExtRobot.h>
//
//
//
//namespace zSpace
//{
//	
//	ZSPACE_TOOLSETS_INLINE zExtRobot::zExtRobot(zTsRHWC* r)
//	{
//		robot = r;
//		updateAttributes();
//	}
//
//	ZSPACE_TOOLSETS_INLINE void zExtRobot::updateAttributes()
//	{
//		//targetCount = robot->robotTargets.size();
//		vector<zTransform> rj = robot->robotJointTransforms;
//
//		//for (int j = 0; j < 6; j++)
//		//{
//		//	Matrix4f out = rj[j];
//
//
//		//	for (int row = 0; row <4; row++)
//		//	{
//		//		for (int col = 0; col < 4; col++)
//		//		{
//		//			robotJointTransforms[j * 6 + row * 4 + col] = out(row, col);
//		//			cout << robotJointTransforms[j * 6 + row * 4 + col] << " || ";
//		//		}
//		//		cout << endl;
//		//	}
//
//		//	//cout << endl<< rj[j] << endl;
//		//	cout << endl << endl;
//
//		//}
//
//
//		for (int i = 0; i < 6; i++)
//		{
//			float* data = rj.at(i).data();
//			for (int j = 0; j < 16; j++)
//			{
//				robotJointTransforms[i * 16 + j] = data[j];
//				//printf("Update %f : ", data[j]);
//			}
//		}
//
//		for (int i = 0; i < 6; i++)
//		{
//			robotJointRotation[i] = robot->jointRotations[i].rotation;
//			robotJointRotationMax[i] = robot->jointRotations[i].maximum;
//			robotJointRotationMin[i] = robot->jointRotations[i].minimum;
//			robotJointRotationHome[i] = robot->jointRotations[i].home;
//			robotJointRotationOffset[i] = robot->jointRotations[i].offset;
//			robotJointRotationPulse[i] = robot->jointRotations[i].pulse;
//			robotJointRotationMask[i] = robot->jointRotations[i].mask;
//		}
//	}
//	
//
//	
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_Test2(char* path)
//	{
//		printf("\n Test2 - char* input");
//
//		cout << "\n Test2"; 
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_Test1(zExtRobot& extRobot)
//	{
//		printf("\n Test1 - zExtRobot&");
//		cout << "\n Test1";
//
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_setDirectory(zExtRobot& extRobot, char* dirPath)
//	{
//		printf("\n ext_zTsRobot_setDirectory - 0");
//
//		vector<string> robotFile;
//
//		zUtilsCore core;
//		zObjMeshArray r_meshObjs;
//		zObjGraph r_graphObj;
//		printf("\n ext_zTsRobot_setDirectory - 1");
//
//		core.getFilesFromDirectory(robotFile, dirPath, zJSON);
//
//		printf("\n ext_zTsRobot_setDirectory - 2");
//
//
//		int nF = core.getNumfiles_Type(dirPath, zOBJ);
//		if (nF < 8) nF = 8;
//
//		r_meshObjs.assign(nF, zObjMesh());
//
//		extRobot.robot = new zTsRHWC();
//
//		printf("\n ext_zTsRobot_setDirectory - 3");
//
//
//		cout << endl << robotFile[0];
//
//		extRobot.robot->createRobotfromFile(robotFile[0], zJSON);
//
//		printf("\n ext_zTsRobot_setDirectory - 4");
//
//		
//
//		extRobot.robot->createRobotJointMeshesfromFile(dirPath, zOBJ, true);
//
//
//		zTransform robotTarget;
//		zTransform robotEE;
//
//		// set target transform
//
//		robotTarget.setIdentity();
//		robotTarget(0, 0) = -1; robotTarget(0, 1) = 0; robotTarget(0, 2) = 0;
//		robotTarget(1, 0) = 0; robotTarget(1, 1) = 1; robotTarget(1, 2) = 0;
//		robotTarget(2, 0) = 0; robotTarget(2, 1) = 0; robotTarget(2, 2) = -1;
//		robotTarget(3, 0) = 2.0; robotTarget(3, 1) = 0; robotTarget(3, 2) = 0;
//
//		// set EE transform
//		robotEE.setIdentity();
//
//		robotEE(0, 0) = 0; robotEE(0, 1) = 0; robotEE(0, 2) = -1;
//		robotEE(1, 0) = 0; robotEE(1, 1) = 1; robotEE(1, 2) = 0;
//		robotEE(2, 0) = 1; robotEE(2, 1) = 0; robotEE(2, 2) = 0;
//		robotEE(3, 0) = -0.2; robotEE(3, 1) = 0; robotEE(3, 2) = -0.346;
//
//
//		extRobot.robot->setEndEffector(robotEE);
//
//		extRobot.updateAttributes();
//
//
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_createFromFile(zExtRobot& extRobot, char* path)
//	{
//		printf("\n ext_creatFromFile");
//
//		string pathSt(path);
//		
//
//		extRobot.robot = new zTsRHWC();
//
//		extRobot.robot->createRobotfromFile(pathSt, zJSON);
//		printf("\n createRobotfromFile success");
//
//		/*extRobot.robot->computeTargets();
//		int count = extRobot.robot->robotTargets.size();
//		printf("\n targetCount %i", count);
//		printf("\n computeTargets success");*/
//
//		extRobot.updateAttributes();
//		printf("\n ext_creatFromFile END \n ");
//
//
//	}
//	ZSPACE_TOOLSETS void ext_zTsRobot_getComputedTargets(zExtRobot& extRobot, float* outTargets)
//	{
//		printf("\n ext_zTsRobot_getComputedTargets");
//
//		vector<zTransform> tt = extRobot.robot->robotTargets;
//		for (int i = 0; i < 6; i++)
//		{
//			float* data = tt.at(i).data();
//			for (int j = 0; j < 16; j++)
//			{
//				outTargets[i * 16 + j] = data[j];
//			}
//		}
//	}
//	ZSPACE_TOOLSETS void ext_zTsRobot_setEndEffector(zExtRobot& extRobot, float* eeTransform)
//	{
//		printf("\n ext_zTsRobot_setEndEffector");
//
//		zTransform t;
//		t.setIdentity();
//
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				//t(i, j) = extRobot.EETransform[i * 4 + j];
//				t(i, j) = eeTransform[i * 4 + j];
//			}
//		}
//
//		extRobot.robot->setEndEffector(t);
//		extRobot.updateAttributes();
//
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_setFabricationMesh(zExtRobot& extRobot, zExtMesh* fabMesh, int fabMeshCount, float* fabBasePlane, float* workBasePlane, int outTargetCount)
//	{
//
//		zObjMeshArray o_FabricationMeshes;
//		o_FabricationMeshes.assign(fabMeshCount, zObjMesh());
//
//		for (int i = 0; i < fabMeshCount; i++)
//		{
//			o_FabricationMeshes[i] = *fabMesh[i].mesh;
//		}
//
//		extRobot.robot->setFabMeshes(o_FabricationMeshes);
//
//
//		zTransform fabBase;
//		fabBase.setIdentity();
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				fabBase(i, j) = fabBasePlane[i * 4 + j];
//			}
//		}
//		extRobot.robot->setFabricationWorkbase(fabBase);
//
//		zTransform workBase;
//		workBase.setIdentity();
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				workBase(i, j) = workBasePlane[i * 4 + j];
//			}
//		}
//		extRobot.robot->setFabricationWorkbase(workBase);
//
//		extRobot.robot->computeTargets();
//
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_forwardKinematics(zExtRobot& extRobot, float* jointRotation)
//	{
//		printf("\n ext_zTsRobot_forwardKinematics");
//		
//		for (int i = 0; i < 6; i++)
//		{
//			extRobot.robot->jointRotations[i].rotation = jointRotation[i];
//		}
//
//		extRobot.robot->forwardKinematics();
//
//		extRobot.updateAttributes();
//
//	}
//
//	ZSPACE_TOOLSETS void ext_zTsRobot_inverseKinematics(zExtRobot& extRobot, float* targetTransformation)
//	{
//		printf("\n ext_zTsRobot_inverseKinematics \n");
//
//		zTransform t;
//		t.setIdentity();
//
//		for (int i = 0; i < 4; i++)
//		{
//			for (int j = 0; j < 4; j++)
//			{
//				t(i, j) = targetTransformation[i * 4 + j];
//				
//			}
//
//		}
//		cout << t << endl;
//
//		zTransform tt = t;
//		extRobot.robot->setTarget(tt);
//
//		cout << endl << "inside inverseKinematics" << endl << extRobot.robot->robot_target_matrix << endl;
//
//		extRobot.robot->inverseKinematics();
//
//		extRobot.updateAttributes();
//	}
//
//
//
//	//ZSPACE_TOOLSETS void ext_zTsRobot_createRobotJointMeshesfromFile(zExtRobot& extRobot, char* directory, bool endeffector)
//	//{
//	//	std::string pathSt(directory);
//
//	//	extRobot.robot->createRobotJointMeshesfromFile(pathSt, zJSON, endeffector); 
//	//	extRobot.updateAttributes();
//
//	//	
//	//}
//
//	//ZSPACE_TOOLSETS void ext_zTsRobot_setJointMeshTransform(zExtRobot& extRobot, bool updatePositions)
//	//{
//	//	extRobot.robot->setJointMeshTransform(updatePositions);
//	//	extRobot.updateAttributes();
//
//	//	
//	//}
//
//	//ZSPACE_TOOLSETS void ext_zTsRobot_setTarget(zExtRobot& extRobot, float* target)
//	//{
//	//}
//
//
//	//ZSPACE_TOOLSETS void ext_zTsRobot_updateKinematics(zExtRobot& extRobot, bool inverseKinematics)
//	//{
//	//	if (inverseKinematics)
//	//	{
//
//	//		zTransform t;
//	//		t.setIdentity();
//
//	//		for (int i = 0; i < 4; i++)
//	//		{
//	//			for (int j = 0; j < 4; j++)
//	//			{
//	//				t(i, j) = extRobot.robot_target_matrix[i * 4 + j];
//	//			}
//	//		}
//
//	//		extRobot.robot->setTarget(t);
//	//		extRobot.robot->inverseKinematics();
//
//	//	}
//	//	else
//	//	{
//	//		zTransform t;
//	//		t.setIdentity();
//
//	//		for (int i = 0; i < 4; i++)
//	//		{
//	//			for (int j = 0; j < 4; j++)
//	//			{
//	//				t(i, j) = extRobot.robot_target_matrix[i * 4 + j];
//	//			}
//	//		}
//	//		extRobot.robot->forwardKinematics();
//	//	}
//	//}
//
//	//ZSPACE_TOOLSETS void ext_zTsRobot_inverseKinematics(zExtRobot& extRobot, float* target)
//	//{
//	//	zTransform t;
//	//	t.setIdentity();
//
//	//	for (int i = 0; i < 4; i++)
//	//	{
//	//		for (int j = 0; j < 4; j++)
//	//		{
//	//			t(i, j) = extRobot.EETransform[i * 4 + j];
//	//		}
//	//	}
//
//	//}
//
//}