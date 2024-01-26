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


#include "headers/zToolsets/digiFab/zTsRobotFab.h"


namespace zSpace
{
	//--zTsHWC Class
	

	//---- DESTRUCTOR
	ZSPACE_TOOLSETS_INLINE zTsRHWC::~zTsRHWC() {}

	zTsRHWC::zTsRHWC(zTsRobot& obj) : zTsRobot(obj) {}

	//--- SET METHODS

	//--- GET METHODS 




	//--- CREATE METHODS 	


	//--- COMPUTE METHODS 

	ZSPACE_TOOLSETS_INLINE void zTsRHWC::computeTargets()
	{

		toWorkBase();
		robotTargets.clear();
		OV_angles.clear();

		//add home pos
		addTarget(o_fabObj.robot_home);
		robotTargetTypes.push_back(1);
		OV_angles.push_back(0.0f);

		for (auto& cutMesh : o_fabObj.fabMeshes)
		{
			vector<zTransform> targets_strip;
			vector<float> ov_angle_strip;

			computeTargetsOnStrip(cutMesh, targets_strip, ov_angle_strip);

			addSafeTargets(targets_strip, ov_angle_strip,1.2);
			//checkTargetNormal(targets_strip);
			addTargets(targets_strip);
			OV_angles.insert(OV_angles.end(), ov_angle_strip.begin(), ov_angle_strip.end());


			// 0 - cut vel
			// 1 - travel vel
			robotTargetTypes.push_back(1);
			robotTargetTypes.push_back(1);
			for (int i = 0; i < targets_strip.size() - 4; i++)
			{
				robotTargetTypes.push_back(0);
			}
			robotTargetTypes.push_back(0);
			robotTargetTypes.push_back(1);

		}

		//add home pos
		addTarget(o_fabObj.robot_home);
		robotTargetTypes.push_back(1);
		OV_angles.push_back(0.0f);

		//o_fabObj.targets = robotTargets;

		cout << endl << "ov angles" << endl;
		for (auto& a : OV_angles) cout << a << ",";
		cout << endl;
	}

	ZSPACE_TOOLSETS_INLINE void zTsRHWC::computeGcode()
	{
		robotTargetReachabilities.clear();
		robotTargetSpeeds.clear();

		for (int i = 0; i < robotTargets.size(); i++)
		{
			zPoint pos = zVector(robotTargets[i](3, 0), robotTargets[i](3, 1), robotTargets[i](3, 2));
			float vel = (robotTargetTypes[i] == 0) ? o_fabObj.vel_travel : o_fabObj.vel_work;
			zRobotMoveType moveType = (robotTargetTypes[i] == 0) ? zMoveLinear : zMoveJoint;

			//cout << endl << "target_"<<i << robotTargets [i]<<endl;

			setTarget(robotTargets[i]);
			inverseKinematics();

			gCode_store(pos, vel, moveType, zEEOn);

			robotTargetSpeeds.push_back(vel);
			robotTargetReachabilities.push_back(inReach);
			
		}
	}






	//--- PRIVATE METHODS 




	
	//--- EXPORT METHODS 
	
	//---- PROTECTED UTILITY METHODS

	ZSPACE_TOOLSETS_INLINE void zTsRHWC::computeTargetsOnStrip(zObjMesh& cutMesh, vector<zTransform>& targets_strip, vector<float>& ov_angle_strip)
	{
		zItMeshHalfEdge he(cutMesh, 0);
		for (he.begin(); !he.end(); he++)
		{
			if (he.onBoundary() && he.getVertex().checkValency(2) && he.getSym().getVertex().checkValency(2))
			{
				break;
			}
		}

		zVector frame_Z(0, 0, -1);
		zVector vNormal(0, 0, 1);
		double angle = 0;
		zVector faceNormal(0, 0, 0);

		zItMeshFaceArray tempFhe;
		he.getFaces(tempFhe);
		faceNormal = zVector(0, 0, 0);
		for (zItMeshFace& f : tempFhe)
		{
			faceNormal += f.getNormal();
		}
		zVector worldX(1, 0, 0);
		zVector worldY(0, 0, 1);
		zVector f2d = faceNormal;
		f2d.normalize();
		f2d.z = 0;
		float ang = f2d.angle(worldX);
		if (f2d*worldX > 0)
		{

		}
		//all targets from a strip
		do
		{
			//target
			zVector frame_Y = he.getVector();
			zItMeshFaceArray fhe;
			he.getFaces(fhe);
			faceNormal = zVector(0, 0, 0);
			for (zItMeshFace &f : fhe)
			{
				faceNormal += f.getNormal();
			}
		
			faceNormal /= fhe.size();
			faceNormal.normalize();
			faceNormal *= -1;
			//cout << "\n positions : " << he.getVertex().getPosition() << ", " << he.getStartVertex().getPosition();

			zVector frame_X = frame_Y ^ faceNormal;
			zVector frame_O = he.getCenter();
			frame_Y.normalize();
			frame_X.normalize();

			frame_Z = frame_X ^ frame_Y;
			frame_Z.normalize();


			frame_Z = frame_X ^ frame_Y;
			targets_strip.push_back(targetFromFrames(frame_O, frame_X, frame_Y, frame_Z));

			//ov angle
			vNormal = he.getVertex().getNormal();
			angle = frame_Z.angle360(vNormal, frame_Y);
			ov_angle_strip.push_back(angle);

			he = he.getSym().getNext().getNext();

		} while (!he.getVertex().checkValency(2) && !he.getSym().getVertex().checkValency(2));

		//last target on strip


		zVector frame_Y = he.getVector();
		zItMeshFaceArray fhe;
		he.getFaces(fhe);
		faceNormal = zVector(0, 0, 0);
		for (zItMeshFace& f : fhe)
		{
			faceNormal += f.getNormal();
		}

		faceNormal /= fhe.size();
		faceNormal.normalize();
		faceNormal *= -1;
		//cout << "\n positions : " << he.getVertex().getPosition() << ", " << he.getStartVertex().getPosition();

		zVector frame_X = frame_Y ^ faceNormal;
		zVector frame_O = he.getCenter();
		frame_Y.normalize();
		frame_X.normalize();

		frame_Z = frame_X ^ frame_Y;
		frame_Z.normalize();


		frame_Z = frame_X ^ frame_Y;
		targets_strip.push_back(targetFromFrames(frame_O, frame_X, frame_Y, frame_Z));

		//ov angle
		vNormal = he.getVertex().getNormal();
		angle = frame_Z.angle360(vNormal, frame_Y);
		ov_angle_strip.push_back(angle);

		


		/*frame_X.normalize();
		frame_Y.normalize();
		frame_Z = frame_X ^ frame_Y;
		zVector frame_O = he.getCenter();

		frame_X.normalize();
		frame_Y.normalize();
		frame_Z.normalize();

		targets_strip.push_back(targetFromFrames(frame_O, frame_X, frame_Y, frame_Z));
		vNormal = he.getVertex().getNormal();
		ov_angle_strip.push_back(frame_Z.angle360(vNormal, frame_Y) * -1);*/

	}

	ZSPACE_TOOLSETS_INLINE void zTsRHWC::checkTargetNormal(vector<zTransform>& targets_strip)
	{
		zPoint robotBase = zVector(robot_base_matrix(3, 0), robot_base_matrix(3, 1), robot_base_matrix(3, 2));
		zPoint first = zVector(targets_strip[0](3, 0), targets_strip[0](3, 1), targets_strip[0](3, 2));
		zVector first_normal = zVector(targets_strip[0](0, 0), targets_strip[0](0, 1), targets_strip[0](0, 2));

		zVector check = first - robotBase;
		double dot = first_normal * check;

		if (dot < 0)
			for (auto& p : targets_strip)
			{
				p(0, 0) = -p(0, 0); p(0, 1) = -p(0, 1);
				p(1, 0) = -p(1, 0); p(1, 1) = -p(1, 1);
				p(2, 0) = -p(2, 0); p(2, 1) = -p(2, 1);
			}
	}

	ZSPACE_TOOLSETS_INLINE void zTsRHWC::addSafeTargets(vector<zTransform>& targets_strip, vector<float>& ov_angle_strip, float multiplication)
	{
		int numTargets = targets_strip.size();
		if (targets_strip[0](3, 3) > targets_strip[numTargets - 1](3, 3))
		{
			reverse(targets_strip.begin(), targets_strip.end());
			reverse(ov_angle_strip.begin(), ov_angle_strip.end());
		}

		vector<zTransform> safeTargets_first = computeSafeTargets(targets_strip[0], multiplication);
		vector<zTransform> safeTargets_last = computeSafeTargets(targets_strip[numTargets - 1], multiplication);

		reverse(safeTargets_first.begin(), safeTargets_first.end());
		targets_strip.insert(targets_strip.begin(), safeTargets_first.begin(), safeTargets_first.end());
		targets_strip.insert(targets_strip.end(), safeTargets_last.begin(), safeTargets_last.end());

		ov_angle_strip.insert(ov_angle_strip.begin(), 0);
		ov_angle_strip.insert(ov_angle_strip.begin(), 0);
		ov_angle_strip.insert(ov_angle_strip.end(), 0);
		ov_angle_strip.insert(ov_angle_strip.end(),0);
	}

	ZSPACE_TOOLSETS_INLINE vector<zTransform> zTsRHWC::computeSafeTargets(zTransform& target, float multiplication)
	{
		vector<zTransform> safeTargets;
		safeTargets.assign(2, zTransform());

		zPoint startPlanePoint = zVector(target(3, 0), target(3, 1), target(3, 2));
		zPoint startPlaneNormal = zVector(target(0, 0), target(0, 1), target(0, 2));

		zPointArray vertices;
		zFnMesh temp(o_fabObj.bbox);
		temp.getVertexPositions(vertices);
		zPoint center(o_fabObj.fabrication_base(3,0), o_fabObj.fabrication_base(3, 1), o_fabObj.fabrication_base(3, 2));
		double check = startPlaneNormal * (startPlanePoint - center);

		float dot_max = -1.0f;
		int id;
		for (int i = 0; i < vertices.size(); i++)
		{
			zVector vec = vertices[i] - startPlanePoint;
			float dot = startPlaneNormal * vec;
			if (dot > dot_max)
			{
				dot_max = dot;
				id = i;
			}
		}
		zPoint move = vertices[id];
		double dist = coreUtils.minDist_Point_Plane(move, startPlanePoint, startPlaneNormal);

		zPoint pos = startPlanePoint + startPlaneNormal * dist * multiplication;
		
		//below safe target
		safeTargets[0] = target;
		safeTargets[0](3, 0) = pos.x;
		safeTargets[0](3, 1) = pos.y;
		safeTargets[0](3, 2) = pos.z;

		//above safe target
		safeTargets[1] = safeTargets[0];
		safeTargets[1](3, 0) = pos.x;
		safeTargets[1](3, 1) = pos.y;
		safeTargets[1](3, 2) = o_fabObj.robot_home(3,2);

		return safeTargets;
	}


	ZSPACE_TOOLSETS_INLINE zTransform zTsRHWC::targetFromFrames(zVector& _position, zVector& _rotationX, zVector& _rotationY, zVector _rotationZ)
	{
		zTransform target;

		target.setIdentity();

		target(0, 0) = _rotationX.x;
		target(0, 1) = _rotationX.y;
		target(0, 2) = _rotationX.z;

		target(1, 0) = _rotationY.x;
		target(1, 1) = _rotationY.y;
		target(1, 2) = _rotationY.z;

		target(2, 0) = _rotationZ.x;
		target(2, 1) = _rotationZ.y;
		target(2, 2) = _rotationZ.z;

		target(3, 0) = _position.x;
		target(3, 1) = _position.y;
		target(3, 2) = _position.z;

		return target;
	}
	ZSPACE_TOOLSETS_INLINE void zTsRHWC::FramesFromTarget(zTransform& target, zVector& _position, zVector& _rotationX, zVector& _rotationY, zVector _rotationZ)
	{
		//zTransform target;

		target.setIdentity();

		_rotationX.x = target(0, 0);
		 _rotationX.y = target(0, 1);
		 _rotationX.z = target(0, 2);
		 _rotationY.x = target(1, 0);
		 _rotationY.y = target(1, 1);
		 _rotationY.z = target(1, 2);
		 _rotationZ.x = target(2, 0);
		 _rotationZ.y = target(2, 1);
		 _rotationZ.z = target(2, 2);
		 _position.x = target(3, 0);
		 _position.y = target(3, 1);
		 _position.z = target(3, 2);

	}
	ZSPACE_TOOLSETS_INLINE void zTsRHWC::getOVAngle(zTransform& target, float& angle, zObjMesh& cutMesh)
	{
		zItMeshHalfEdge he(cutMesh, 0);
		for (he.begin(); !he.end(); he++)
		{
			if (he.onBoundary() && he.getVertex().checkValency(2) && he.getSym().getVertex().checkValency(2))
			{
				break;
			}
		}

		zVector frame_Z(0, 0, -1);
		zVector vNormal(0, 0, 1);
		angle = 0;
		//all targets from a strip
		//do
		//{
			////target
			zVector frame_Y = he.getVector();
			zVector frame_X = frame_Y ^ frame_Z;
			frame_Z = frame_X ^ frame_Y;
			zVector frame_O = he.getCenter();


			//frame_X.normalize();
			//frame_Y.normalize();
			//frame_Z.normalize();
			//targets_strip.push_back(targetFromFrames(frame_O, frame_X, frame_Y, frame_Z));
			
			FramesFromTarget(target, frame_O, frame_X, frame_Y, frame_Z);


			//ov angle
			vNormal = he.getVertex().getNormal();
			angle = frame_Z.angle360(vNormal, frame_Y);
			
			
			//ov_angle_strip.push_back(angle);
			//cout << endl << "vNormal" << vNormal << endl;
			//cout << endl << "frame_Z" << frame_Z << endl;

			//he = he.getSym().getNext().getNext();

		//} while (!he.getVertex().checkValency(2) && !he.getSym().getVertex().checkValency(2));

	}
}
