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


#include<headers/zToolsets/digiFab/zTsRobotFab.h>


namespace zSpace
{
	//--zTsHWC Class
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsRHWC::zTsRHWC()
		:zTsRobot()
	{
		robot_scale = 1.0;

		for (int i = 0; i < DOF; i++)
		{
			zDHparameter DH;
			DH.alpha = DH.d = DH.theta = DH.a = 0.0;
			robot_DH.push_back(DH);

			zJointRotation jointRot;
			jointRot.home = jointRot.minimum = jointRot.maximum = jointRot.rotation = 0;
			jointRot.pulse = jointRot.mask = jointRot.offset = 0.0;
			jointRotations.push_back(jointRot);

			robotJointTransforms.push_back(zTransform(4, 4));
			robotMesh_transforms.push_back(zTransform(4, 4));
		}

		for (int i = 0; i < DOF; i++) jointMeshObjs.push_back(nullptr);

		jointGraphObj = nullptr;
	}

	ZSPACE_INLINE zTsRHWC::zTsRHWC(zObjGraph& _jointGraphObj, vector<zObjMesh>& _jointMeshObjs) 
		:zTsRobot(_jointGraphObj, _jointMeshObjs)
	{
		jointGraphObj = &_jointGraphObj;
		fnGraphJoint = zFnGraph(_jointGraphObj);

		jointMeshObjs.clear();
		fnMeshJoints.clear();

		for (int i = 0; i < _jointMeshObjs.size(); i++)
		{
			jointMeshObjs.push_back(&_jointMeshObjs[i]);
			fnMeshJoints.push_back(_jointMeshObjs[i]);
		}
	}

	//---- DESTRUCTOR
	ZSPACE_INLINE zTsRHWC::~zTsRHWC() {}

	//--- SET METHODS

	//--- GET METHODS 

	ZSPACE_INLINE zPointArray zTsRHWC::getFabBbox()
	{
		//zPointArray bbox;
		//bbox.assign(2, zPoint());

		return fabMeshBbox;
	}

	ZSPACE_INLINE vector<zObjMesh> zTsRHWC::getFabMesh()
	{
		return fabMeshObjs;
	}

	//--- CREATE METHODS 	


	//--- COMPUTE METHODS 

	ZSPACE_INLINE void zTsRHWC::computeTargets()
	{

		for (auto& cutMesh : fabMeshObjs)
		{
			//add home pos
			addTarget(robotHome);
			//myRobot.addTarget(zVector(1.5, -0.25, 1), zVector(1, 0, 0), zVector(0, -1, 0), zVector(0, 0, -1));

			zItMeshHalfEdge he(cutMesh, 0);
			for (he.begin(); !he.end(); he++)
			{
				if (he.onBoundary() && he.getVertex().checkValency(2) && he.getSym().getVertex().checkValency(2))
					break;
			}

			zVector frame_Z(0, 0, -1);

			//check orientation
			zVector check = he.getCenter() - zVector(0, 0, he.getCenter().z);
			zVector vec = he.getVector() ^ frame_Z;
			vec.normalize();
			check.normalize();

			double dot = vec * check;
			he = (dot < 0) ? he.getSym() : he;

			cout << endl << "heID" << he.getId();

			//all targets from a strip
			do
			{
				zVector frame_Y = he.getVector();
				zVector frame_X = frame_Y ^ frame_Z;
				zVector frame_O = he.getCenter();

				frame_X.normalize();
				frame_Y.normalize();
				addTarget(frame_O, frame_X, frame_Y, frame_Z);

				he = (dot < 0) ? he.getNext().getNext().getSym() : he.getSym().getNext().getNext();

			} while (!he.getVertex().checkValency(2) && !he.getSym().getVertex().checkValency(2));

			zVector frame_Y = he.getVector();
			zVector frame_X = frame_Y ^ frame_Z;
			zVector frame_O = he.getCenter();

			frame_X.normalize();
			frame_Y.normalize();
			addTarget(frame_O, frame_X, frame_Y, frame_Z);
		}
	}


	//--- EXPORT METHODS 
	
	//---- PROTECTED UTILITY METHODS
}
