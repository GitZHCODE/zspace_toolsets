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


#include <zToolsets/streamlines/zTsStreamsMesh.h>


//---- zTsStreamszTsStreamsMesh ------------------------------------------------------------------------------

namespace zSpace
{

	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsStreamsMesh::zTsStreamsMesh()
	{
		meshObj = nullptr;

		dSep = 1.0;
		dTest = 0.8;
		minLength = 5.0;
		maxLength = 5.0;
		dT = 0.5;
		angle = 0.0;
		flipBackward = false;

		streamType = zForwardBackward;
		integrationType = zRK4;

	}

	ZSPACE_TOOLSETS_INLINE zTsStreamsMesh::zTsStreamsMesh(zObjMesh& _meshObj, vector<zVector>& _fieldVectors)
	{
		meshObj = &_meshObj;
		fieldVectors = _fieldVectors;

		fieldIndex_streamPositions.clear();

		zFnMesh fnMesh(*meshObj);
		fieldIndex_streamPositions.assign(fnMesh.numPolygons(), zPointArray());
	
		dSep = 1.0;
		dTest = 0.8;
		minLength = 1.0;
		maxLength = 5.0;
		dT = 0.5;
		angle = 0.0;
		flipBackward = false;

		streamType = zForward;
		integrationType = zRK4;

		//---- matrices

		fnMesh.getMatrices_quadmesh(mesh_V, mesh_F);
		//fnMesh.getMatrices_trimesh(mesh_V, mesh_F);

		//---- ring neighbours
		computeRingNeighbours();
	}

	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsStreamsMesh::~zTsStreamsMesh() {}

	//----  SET METHODS

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setSeperationDistance(double _dSep)
	{
		dSep = _dSep;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setTestSeperationDistance(double _dTest)
	{
		dTest = _dTest;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setMinLength(double _minLength)
	{
		minLength = _minLength;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setMaxLength(double _maxLength)
	{
		maxLength = _maxLength;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setIntegrationType(zFieldStreamType _streamType)
	{
		streamType = _streamType;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setTimeStep(double _dT)
	{
		dT = _dT;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setIntegrationType(zIntergrationType _integrationType)
	{
		integrationType = _integrationType;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setAngle(double _angle)
	{
		angle = _angle;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::setFlipBackwards(bool _flipBackward)
	{
		flipBackward = _flipBackward;
	}

	//----  2D STREAM LINES METHODS

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::createStreams(vector<zStreamLine>& streams, vector<zVector> &start_seedPoints, bool seedStreamsOnly,int maxStreams)
	{
		streams.clear();
		streams.assign(maxStreams, zStreamLine());

		int streamCounter = 0;

		// make first stream line
		if (streamCounter == 0)
		{	

			for (int i = 0; i < start_seedPoints.size(); i++)
			{		
				printf("\n stream id %i / %i", streamCounter, maxStreams);
				bool chk = createStreamGraph(streams[streamCounter].graphObj, start_seedPoints[i]);

				if (chk)
				{
					//streams.push_back(temp);

					streams[streamCounter].isValid = true;
					streamCounter++;				

				}
			}




		}


		if (seedStreamsOnly)
		{
			printf("\n %i streamLines created. ", streamCounter);
			return;
		}

		// compute other stream lines.

		int currentStreamGraphId = 0;
		int currentStreamGraphVertexId = 0;

		bool finished = false;

		while (!finished)
		{
			vector<zVector> seedPoints;

			
			zFnGraph fnGraph(streams[currentStreamGraphId].graphObj);

			if (streams[currentStreamGraphId].isValid)
			{
				getSeedPoints(streams[currentStreamGraphId], currentStreamGraphVertexId, seedPoints);


				for (int i = 0; i < seedPoints.size(); i++)
				{

					int id = streamCounter;
					//streams.push_back(zStreamLine());
					printf("\n stream id %i / %i", id, maxStreams);

					if (streamCounter >= maxStreams) continue;

					bool chk = createStreamGraph(streams[id].graphObj, seedPoints[i]);

					if (chk)
					{
						streams[currentStreamGraphId].child.push_back(id);
						streams[id].setParent(currentStreamGraphId);

						streams[id].isValid = true;
						streamCounter++;
					}
				}
			}

			

			currentStreamGraphVertexId++;

			if (currentStreamGraphVertexId >= fnGraph.numVertices()) currentStreamGraphId++, currentStreamGraphVertexId = 0;

			if (currentStreamGraphId >= streams.size()) finished = true;

			if (streamCounter >= maxStreams) finished = true;
		}

		printf("\n %i streamLines created. ", streamCounter);

	}


	//---- COMPUTE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::computeClosestPointToMesh(zPoint& inPoint, int& faceID, zPoint& closestPoint)
	{
		MatrixXd P(1, 3);
		P(0, 0) = inPoint.x;
		P(0, 1) = inPoint.y;
		P(0, 2) = inPoint.z;

		VectorXd sqrD;
		VectorXi I;
		MatrixXd C;
		igl::point_mesh_squared_distance(P, mesh_V, mesh_F, sqrD, I, C);

		faceID = I(0);

		closestPoint.x = C(0, 0);
		closestPoint.y = C(0, 1);
		closestPoint.z = C(0, 2);

		//cout << endl << inPoint << ", " << closestPoint << ", " << faceID;

	}

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::getFieldValue(zPoint& inPoint, zVector& fieldForce, zVector& axis)
	{
		int faceId = -1;
		zPoint cP;
		computeClosestPointToMesh(inPoint, faceId, cP);

		if (faceId == -1) return false;
		
		zItMeshFace f(*meshObj, faceId);

		axis = f.getNormal();

		zPointArray vPositions;
		f.getVertexPositions(vPositions);

		zIntArray vIDs;
		f.getVertices(vIDs);

		zDoubleArray weights;
		coreUtils.getDistanceWeights(cP, vPositions, 1.0, weights);

		zVector fVal;
		double w = 0;
		for (int i = 0; i < vIDs.size(); i++)
		{
			zVector val = fieldVectors[vIDs[i]];			
			fVal += (val * weights[i]);
			w += weights[i];
			
			//printf("\n i  %1.2f %1.2f %1.2f | %1.2f ", fieldVectors[vIDs[i]].x, fieldVectors[vIDs[i]].y, fieldVectors[vIDs[i]].z, weights[i]);
			
		}

		fVal /= w;

		//printf("\n %1.2f %1.2f %1.2f | %1.2f ", fVal.x, fVal.y, fVal.z, w);

		fieldForce = fVal;

		return true;
	}

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::projectToMesh(zObjGraph& inGraph, zIntArray& faceIDs)
	{
		faceIDs.clear();

		zFnGraph fnGraph(inGraph);
		zPoint* vPositions = fnGraph.getRawVertexPositions();

		for (int i = 0; i < fnGraph.numVertices(); i++)
		{
			int faceId = -1;
			zPoint cP;
			computeClosestPointToMesh(vPositions[i], faceId, cP);

			zItMeshFace f(*meshObj, faceId);

			zPoint o = f.getCenter();
			zVector n = f.getNormal();

			float d = coreUtils.minDist_Point_Plane(vPositions[i], o, n);

			vPositions[i] += (n*d);

			faceIDs.push_back(faceId);
		}

		return false;
	}

	//---- PROTECTED METHODS

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::createStreamGraph(zObjGraph &streamGraphObj, zVector &seedPoint)
	{

		vector<zVector> positions;
		vector<int> edgeConnects;


		// move forward
		if (streamType == zForward || streamType == zForwardBackward)
		{
			bool exit = false;

			zVector startForward = seedPoint;

			zObjParticle p;
			p.particle = zParticle(startForward);

			zFnParticle seedForward(p);

			//zFnParticle seedForward;
			//seedForward.create(startForward);

			double currentLength = 0.0;

			//printf("\n working!");

			while (!exit)
			{
				bool firstVertex = (startForward == seedPoint) ? true : false;

				zVector curPos = seedForward.getPosition();

				if (firstVertex)
				{

					positions.push_back(curPos);
				}

				// get field focrce
				zVector fieldForce;
				zVector axis;
				bool checkBounds = getFieldValue(curPos, fieldForce, axis);



				if (!checkBounds)
				{
					//printf("\n bounds working!");

					exit = true;
					continue;
				}

				// local minima or maxima point
				if (fieldForce.length() == 0)
				{
					//printf("\n force working!");

					exit = true;
					continue;
				}

				// update particle force
				fieldForce.normalize();
				//cout << "\n field force " << fieldForce << "| " << dSep;
				fieldForce *= (dSep * 1.0);		

				double rotateAngle = angle;
				if(rotateAngle > 0) fieldForce = fieldForce.rotateAboutAxis(axis, rotateAngle);

			
				seedForward.addForce(fieldForce);

				// update seed particle
				seedForward.integrateForces(dT, integrationType);
				seedForward.updateParticle(true,true,true);
				zVector newPos = seedForward.getPosition();
				zVector newPos_mesh;
				int faceId = -1;
				checkBounds = checkFieldBounds(newPos, newPos_mesh, faceId);

				

				if (!checkBounds)
				{
					//printf("\n bounds 2 working!");

					//printf("\n %1.2f %1.2f %1.2f ", newPos.x, newPos.y, newPos.z);

					exit = true;
					continue;



				}				
				//newPos = newPos_mesh;
				//seedForward.setPosition(&newPos);
			
				int index = -1;
				bool checkRepeat = coreUtils.checkRepeatElement(newPos, positions, index);
				

				if (checkRepeat)
				{
					//printf("\n repeat working!");

					/*for (int k = 0; k < positions.size(); k++)
					{
						cout << "\n " << newPos << ", " << positions[k];
					}*/

					exit = true;
					continue;



				}


				bool validStreamPoint = checkValidStreamPosition(newPos, dTest);

				if (!validStreamPoint)
				{
					exit = true;

					//printf("\n validity working!");
				}

				// check length
				if (currentLength + curPos.distanceTo(newPos) > maxLength)
				{
					exit = true;

					//printf("\n length working!");
				}

				// add new stream point
				if (!exit)
				{


					if (positions.size() > 0)
					{
						edgeConnects.push_back(positions.size());
						edgeConnects.push_back(positions.size() - 1);
					}

					positions.push_back(newPos);

					currentLength += curPos.distanceTo(newPos);

				}




			}
		}


		if (streamType == zBackward || streamType == zForwardBackward)
		{
			// move backwards
			bool exit = false;

			zVector startBackward = seedPoint;

			zObjParticle p;
			p.particle = zParticle(startBackward);

			zFnParticle seedBackward(p);

			//zFnParticle seedBackward;
			//seedBackward.create(startBackward);


			double currentLength = 0.0;

			while (!exit)
			{
				bool firstVertex = (startBackward == seedPoint) ? true : false;


				zVector curPos = seedBackward.getPosition();


				// insert first point if the stream is inly for backward direction.
				if (firstVertex && streamType == zBackward)
				{
					positions.push_back(curPos);
				}

				// get field focrce
				zVector fieldForce;
				zVector axis;
				bool checkBounds = getFieldValue(curPos, fieldForce, axis);

				if (!checkBounds)
				{
					exit = true;
					continue;
				}
				// local minima or maxima point
				if (fieldForce.length() == 0)
				{
					exit = true;
					continue;
				}

				// update particle force
				fieldForce.normalize();
				fieldForce *= (dSep * 1.0);

				
				fieldForce *= -1;

				double rotateAngle = angle;
				if (!flipBackward) rotateAngle = 180.0 - angle;

				fieldForce = fieldForce.rotateAboutAxis(axis, rotateAngle);

				seedBackward.addForce(fieldForce);

				// update seed particle
				seedBackward.integrateForces(dT, integrationType);
				seedBackward.updateParticle(true);
				zVector newPos = seedBackward.getPosition();
				zVector newPos_mesh;
				int faceId;
				checkBounds = checkFieldBounds(newPos, newPos_mesh, faceId);

				if (!checkBounds)
				{
					exit = true;
					continue;
				}
			
				//newPos = newPos_mesh;
				//seedBackward.setPosition(&newPos);

				int index = -1;
				bool checkRepeat = coreUtils.checkRepeatElement(newPos, positions, index);
				if (checkRepeat)
				{
					exit = true;
					continue;
				}

				bool validStreamPoint = checkValidStreamPosition(newPos, dTest);


				if (!validStreamPoint) exit = true;

				// check length
				if (currentLength + curPos.distanceTo(newPos) > maxLength) exit = true;

				// add new stream point
				if (!exit)
				{


					if (positions.size() > 0)
					{
						(firstVertex) ? edgeConnects.push_back(0) : edgeConnects.push_back(positions.size() - 1);
						edgeConnects.push_back(positions.size());

					}

					positions.push_back(newPos);

					currentLength += curPos.distanceTo(newPos);
				}

			}
		}


		/*printf("\n v: %i e:%i ", positions.size(), edgeConnects.size());*/

		// create stream graph
		bool out = false;

		if (edgeConnects.size() > 0)
		{

			zFnGraph tempFn(streamGraphObj);
			tempFn.create(positions, edgeConnects);

			vector<double> lengths;
			double length = tempFn.getEdgeLengths(lengths);

			zIntArray faceIDs;
			projectToMesh(streamGraphObj, faceIDs);

			positions.clear();
			tempFn.getVertexPositions(positions);

			tempFn.setEdgeColor(zGREEN);

			if (length > minLength)
			{
				for (int i = 0; i < positions.size(); i++)
				{
					addToFieldStreamPositions(positions[i], faceIDs[i]);

				}

				out = true;
			}



		}



		return out;

	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::getSeedPoints(zStreamLine& currentStream, int vertexId, vector<zVector> &seedPoints)
	{
		zFnGraph tempFn(currentStream.graphObj);
		if (tempFn.numEdges() == 0) return;
		
		zItGraphVertex v(currentStream.graphObj, vertexId);

		if (v.checkValency(1)) return;		
		
		zVector up(0, 0, 1);
		zVector norm;

		zVector vPos = v.getPosition();

		int faceId = -1;
		zPoint cP;
		computeClosestPointToMesh(vPos, faceId, cP);
		zItMeshFace f(*meshObj, faceId);

		up = f.getNormal();
		zItGraphHalfEdge curEdge = v.getHalfEdge();


		if (curEdge.getVertex().isActive())
		{
			zVector v1 = curEdge.getVertex().getPosition();
			zVector e1 = v1 - vPos;
			e1.normalize();

			norm += e1 ^ up;
		}

		if (curEdge.getPrev().isActive())
		{
			zVector v2 = curEdge.getPrev().getStartVertex().getPosition();
			zVector e2 = vPos - v2;
			e2.normalize();

			norm += e2 ^ up;
		}


		if (norm.length() < EPS) return;

		norm *= 0.5;
		norm.normalize();

		zVector tempSeedPoint = vPos + (norm* dSep);
		zPoint tempSeedPoint_mesh;
		computeClosestPointToMesh(tempSeedPoint, faceId, tempSeedPoint_mesh);

		/*zItMeshFace f0(*meshObj, faceId);
		zPoint o = f0.getCenter();
		zVector n = f0.getNormal();
		float d = coreUtils.minDist_Point_Plane(tempSeedPoint, o, n);
		tempSeedPoint += (n * d);
		*/

		tempSeedPoint = tempSeedPoint_mesh;

		bool out = checkValidSeedPosition(tempSeedPoint, dSep);
		//printf("\n valid 1 %s", (out) ? "T" : "F");
		if (out)  seedPoints.push_back(tempSeedPoint);

		tempSeedPoint = vPos + (norm * dSep*-1);
		computeClosestPointToMesh(tempSeedPoint, faceId, tempSeedPoint_mesh);

	/*	zItMeshFace f1(*meshObj, faceId);
		o = f1.getCenter();
		n = f1.getNormal();
		d = coreUtils.minDist_Point_Plane(tempSeedPoint, o, n);
		tempSeedPoint += (n * d);*/

		tempSeedPoint = tempSeedPoint_mesh;

		out = checkValidSeedPosition(tempSeedPoint, dSep);
		//printf("\n valid 2 %s", (out) ? "T" : "F");
		if (out)  seedPoints.push_back(tempSeedPoint);

	}


	//----  2D FIELD UTILITIES

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::checkFieldBounds(zPoint &inPoint, zPoint& cP, int& faceID)
	{

		faceID = -1;
		zPoint closeP;
		computeClosestPointToMesh(inPoint, faceID, cP);

		if (faceID == -1) return false;
	
		zItMeshFace f(*meshObj, faceID);

		/*zPoint o = f.getCenter();
		zVector n = f.getNormal();

		float d = coreUtils.minDist_Point_Plane(inPoint, o, n);
		cP = inPoint + (n * d);*/

		return !f.onBoundary();	
	}

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::checkValidStreamPosition(zVector &inPoint, double &dTest)
	{

		int newFieldIndex;
		zPoint inPoint_Mesh;
		bool checkBounds = checkFieldBounds(inPoint, inPoint_Mesh, newFieldIndex);

		bool validStreamPoint = true;

		for (int j = 0; j < fieldIndex_streamPositions[newFieldIndex].size(); j++)
		{
			zVector streamPoint = fieldIndex_streamPositions[newFieldIndex][j];


			double dist = streamPoint.distanceTo(inPoint);

			if (dist < dTest)
			{
				validStreamPoint = false;
				break;
			}
		}

		// check in neighbour if validStreamPoint is true
		if (validStreamPoint)
		{
			

			for (int i = 0; i < ringNeighbours[newFieldIndex].size(); i++)
			{
				if (ringNeighbours[newFieldIndex][i] == newFieldIndex) continue;

				for (int j = 0; j < fieldIndex_streamPositions[ringNeighbours[newFieldIndex][i]].size(); j++)
				{
					zVector streamPoint = fieldIndex_streamPositions[ringNeighbours[newFieldIndex][i]][j];

					double dist = streamPoint.distanceTo(inPoint);

					if (dist < dTest)
					{
						validStreamPoint = false;


						break;
					}
				}

			}
		}

		return validStreamPoint;
	}

	ZSPACE_TOOLSETS_INLINE bool zTsStreamsMesh::checkValidSeedPosition(zVector &inPoint, double &dSep)
	{


		int newFieldIndex;
		zPoint inPoint_Mesh;
		bool checkBounds = checkFieldBounds(inPoint, inPoint_Mesh, newFieldIndex);

		if (!checkBounds) return false;

		bool validSeedPoint = true;

		for (int j = 0; j < fieldIndex_streamPositions[newFieldIndex].size(); j++)
		{
			zVector streamPoint = fieldIndex_streamPositions[newFieldIndex][j];

			double dist = streamPoint.distanceTo(inPoint);

			if (dist < dSep)
			{
				validSeedPoint = false;
				break;
			}
		}

		// check in neighbour if validStreamPoint is true
		if (validSeedPoint)
		{
			

			for (int i = 0; i < ringNeighbours[newFieldIndex].size(); i++)
			{
				if (ringNeighbours[newFieldIndex][i] == newFieldIndex) continue;

				for (int j = 0; j < fieldIndex_streamPositions[ringNeighbours[newFieldIndex][i]].size(); j++)
				{
					zVector streamPoint = fieldIndex_streamPositions[ringNeighbours[newFieldIndex][i]][j];

					double dist = streamPoint.distanceTo(inPoint);

					if (dist < dSep)
					{
						validSeedPoint = false;
						break;
					}
				}

			}
		}

		return validSeedPoint;
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::addToFieldStreamPositions(zVector &inPoint, int faceId)
	{

		fieldIndex_streamPositions[faceId].push_back(inPoint);
	}

	ZSPACE_TOOLSETS_INLINE void zTsStreamsMesh::computeRingNeighbours()
	{
		ringNeighbours.clear();

		for (zItMeshFace f(*meshObj); !f.end(); f++)
		{
			zIntArray fNeighbours;

			zItMeshVertexArray fVerts;
			f.getVertices(fVerts);

			for (auto& v : fVerts)
			{
				zIntArray cFaces;
				v.getConnectedFaces(cFaces);

				for (auto &cF: cFaces)
				{
					int id = -1;
					bool chkRepeat = coreUtils.checkRepeatElement(cF,fNeighbours, id);

					if (!chkRepeat) fNeighbours.push_back(cF);
				}
				
			}

			ringNeighbours.push_back(fNeighbours);
		}

	}
	
}