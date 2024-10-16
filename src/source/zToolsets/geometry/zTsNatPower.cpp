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


#include "zToolsets/geometry/zTsNatPower.h"

namespace zSpace
{



	//ZSPACE_TOOLSETS_INLINE zHalfMedial::zHalfMedial()
	//{
	//}
	//ZSPACE_TOOLSETS_INLINE zHalfMedial::~zHalfMedial()
	//{
	//}

	//ZSPACE_TOOLSETS_INLINE zTsNatPower::zTsNatPower()
	//{
	//	/*medials_inner.assign(0, zHalfMedial());
	//	medials_outer.assign(0, zHalfMedial());*/
	//}

	//ZSPACE_TOOLSETS_INLINE zTsNatPower::~zTsNatPower()
	//{
	//}

	/*
	* PUBLIC METHODS
	*/
	ZSPACE_TOOLSETS_INLINE void zTsNatPower::setMeshFromPath(string& innerMeshPath, string& outerMeshPath)
	{
		o_inMesh_inner = new zObjMesh();
		o_inMesh_outer = new zObjMesh();

		int innerStride, outerStride;
		readMeshJson(innerMeshPath, o_inMesh_inner, startVertexId_inner, innerStride);

		zFnMesh fnMesh(*o_inMesh_inner);

		readMeshJson(outerMeshPath, o_inMesh_outer, startVertexId_outer, outerStride);

		if (innerStride != outerStride) printf("\n innerStride != outerStride  %i  %i ", innerStride, outerStride);
		meshStride = (innerStride != outerStride)? -1 : innerStride;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeMedials()
	{

		double total;
		std::chrono::time_point<std::chrono::system_clock> s_Time = std::chrono::system_clock::now();

		

		computeAllMedialGraphs(*o_inMesh_inner, startVertexId_inner, meshStride, 1, feature1Map_inner, medials_inner);
		computeAllMedialGraphs(*o_inMesh_outer, startVertexId_outer, meshStride , 1, feature1Map_outer, medials_outer);

		std::chrono::time_point<std::chrono::system_clock> e_Time = std::chrono::system_clock::now();
		total = std::chrono::duration_cast<std::chrono::milliseconds>(e_Time - s_Time).count();

		printf("\n medial computed! Compute time %1.4f on second", total/1000 );
		//computeMedialHE(o_inMesh_inner, startVertexId_inner, medial_inner, medialMap_inner, feature1_inner, feature1Map_inner);
		//computeMedialHE(o_inMesh_outer, startVertexId_outer, medial_outer, medialMap_outer, feature1_outer, feature1Map_outer);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeInterfacePlanes(float interfaceAngle, float dividLength)
	{

		//o_interfacePlanes.assign(medials_inner.size(), zObjPlaneArray());
		//feature_inner.assign(medials_inner.size(), zObjNurbsCurveArray());
		//feature_outer.assign(medials_inner.size(), zObjNurbsCurveArray());



		double total;
		std::chrono::time_point<std::chrono::system_clock> s_Time = std::chrono::system_clock::now();

		for (int i = 0; i < medials_inner.size(); i++)
		{
			computeInterfacePlanesOnMedial(medials_inner[i], medials_outer[i], interfaceAngle, dividLength);

			//for(auto& p : o_interfacePlanes[i]) printf("\n fnp %1.4f - %1.4f", p.origin, p.normal);


			//zObjGraphArray outGraph;
			//computePlaneIntersectionOnMedial(medials_inner[i], medials_outer[i], o_interfacePlanes[i],  outGraph);

		}
		std::chrono::time_point<std::chrono::system_clock> e_Time = std::chrono::system_clock::now();
		total = std::chrono::duration_cast<std::chrono::milliseconds>(e_Time - s_Time).count();


		printf("\n interface planes computed! Compute time %1.4f on second", total / 1000);

		s_Time = std::chrono::system_clock::now();
		for (int i = 0; i < medials_inner.size(); i++)
		{

			//for(auto& p : o_interfacePlanes[i]) printf("\n fnp %1.4f - %1.4f", p.origin, p.normal);

			//printf(" \n Medial ID %i \n ", i);
			zObjGraphArray outGraph;
			computePlaneIntersectionOnMedial(medials_inner[i], medials_outer[i], o_interfacePlanes[i],  outGraph);

		}
		e_Time = std::chrono::system_clock::now();
		total = std::chrono::duration_cast<std::chrono::milliseconds>(e_Time - s_Time).count();

		printf("\n plane intersection computed! Compute time %1.4f on second", total / 1000);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getFeatureNurbs(vector<zObjNurbsCurveArray> & outInnerFeature, vector<zObjNurbsCurveArray> & outOuterFeature)
	{
		/*
		printf("\n feature size : %i - %i ", feature1_inner.size(), feature1_outer.size());
		outInnerFeature.clear();
		outInnerFeature.assign(feature1_inner.size(), zObjNurbsCurveArray());
		for (int i = 0; i < feature1_inner.size(); i++)
		{
			outInnerFeature[i].assign(feature1_inner[i].size(), zObjNurbsCurve());

			for (int j = 0; j < feature1_inner[i].size(); j++)
			{
				createNurbsFromHEArray(feature1_inner[i][j], outInnerFeature[i][j]);

			}
			//zFnNurbsCurve fn(outInnerFeature[i]);
			///printf("\n outInnerFeature[i]  : %i ", fn.numControlVertices());


		}

		*/

		/*outOuterFeature.clear();
		outOuterFeature.assign(feature1_outer.size(), zObjNurbsCurve());
		for (int i = 0; i < outOuterFeature.size(); i++)
		{
			createNurbsFromHEArray(feature1_outer[i], outOuterFeature[i]);

			zFnNurbsCurve fn(outOuterFeature[i]);
			printf("\n outOuterFeature[i]  : %i ", fn.numControlVertices());
		}*/
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getMesh(zObjMesh& innerMesh, zObjMesh& outerMesh)
	{
		innerMesh = *o_inMesh_inner;
		outerMesh = *o_inMesh_outer;

	}

	/*
	* PRIVATE METHODS
	*/
	ZSPACE_TOOLSETS_INLINE void zTsNatPower::readMeshJson(string& meshPath, zObjMesh* outMesh, int& outStartVertexInd, int& outStride)
	{
		json oJson;
		core.json_read(meshPath, oJson);
		
		zFnMesh fnMesh(*outMesh);

		fnMesh.from(oJson);

		printf("\n mesh numV %i numE %i ", fnMesh.numVertices(), fnMesh.numEdges());


		int index, strid;
		string keyId = "startVID";
		string keyStride = "blockStride";


		bool chk = core.json_readAttribute(oJson, keyId, index);
		outStartVertexInd = chk ? index : -1;
		if (!chk) printf("\n json unable to read '%s'  ", keyId);
		//else printf("Able to read '%s' %i ", keyId, index);


		bool chk2 = core.json_readAttribute(oJson, keyStride, strid);
		outStride = chk2 ? strid : -1;
		if (!chk2) printf("\n json unable to read '%s'  ", keyStride);
		//else printf("Able to read '%s' %i ", keyStride, strid);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeAllMedialGraphs(zObjMesh& o_mesh, int startVID, int numFeatureLoops, int numStride, vector<int>& _He_Medial_Map, vector<zHalfMedial>& _medials)
	{
		_medials.clear();
		_He_Medial_Map.clear();

		zFnMesh fnMesh(o_mesh);
		_He_Medial_Map.assign(fnMesh.numHalfEdges(), -1);

		// get corner vertex
		zItMeshVertex v(o_mesh, startVID);

		zItMeshHalfEdgeArray cHEdges;
		v.getConnectedHalfEdges(cHEdges);

		zItMeshHalfEdge startHe;
		for (auto& he : cHEdges)
		{
			if (!he.getEdge().onBoundary()) startHe = he;
		}

		zItMeshHalfEdge he = startHe;

		zItMeshHalfEdge he_top_Start;
		bool compHeTopStart = true;

		int numMedials = 0;

		//////////  COMPUTE NUMBER OF MEDIALS;
		//bottom loops
		do
		{
			he.getEdge().setColor(zRED);
			zItMeshHalfEdge heSec = he;

			// store he for top loop walk
			if (compHeTopStart)
			{
				if (he.getStartVertex().getValence() > 4)
				{
					he_top_Start = he.getSym();
					he_top_Start = he_top_Start.getNext().getSym().getNext();

					compHeTopStart = false;
				}
			}

			// add
			if (he.getStartVertex().getValence() > 4) numMedials++;

			for (int i = 0; i < numFeatureLoops; i++)
			{
				for (int j = 0; j < numStride; j++) heSec = heSec.getPrev().getPrev().getSym();

				heSec.getEdge().setColor(zBLUE);
			}

			if (he.getVertex().onBoundary())
			{
				he = he.getSym();
				numMedials++;
			}
			else
			{
				he = he.getNext().getSym().getNext();
			}

		} while (he != startHe);


		// top loop
		he = he_top_Start;

		do
		{
			he.getEdge().setColor(zYELLOW);
			zItMeshHalfEdge heSec = he;

			// add
			if (he.getStartVertex().getValence() > 4) numMedials++;

			for (int i = 0; i < numFeatureLoops; i++)
			{
				for (int j = 0; j < numStride; j++) heSec = heSec.getPrev().getPrev().getSym();

				heSec.getEdge().setColor(zGREEN);

			}

			he = he.getNext().getSym().getNext();

		} while (he != he_top_Start);

		printf("\n num Medials %i ", numMedials);

		_medials.assign(numMedials, zHalfMedial());

		////////// CREATE MEDIAL GRAPHS
		he = startHe;
		int medialCounter = 0;

		zItMeshHalfEdgeArray tmp_medials;
		vector<zItMeshHalfEdgeArray> tmp_features;
		tmp_features.assign(numFeatureLoops, zItMeshHalfEdgeArray());

		do
		{
			zItMeshHalfEdge heSec = he;

			// store
			tmp_medials.push_back(he);
			_He_Medial_Map[he.getId()] = medialCounter;

			// add medial graph
			if (he.getVertex().getValence() > 4 || he.getVertex().onBoundary())
			{
				//printf("\n medial graph ");
				createGraphFromHalfEdgeLoop(_medials[medialCounter].o_medial, tmp_medials, zRED);
				_medials[medialCounter].meshHEStartEnd.assign(2, zItMeshHalfEdge());
				_medials[medialCounter].meshHEStartEnd[0] = tmp_medials[0];
				_medials[medialCounter].meshHEStartEnd[1] = tmp_medials[tmp_medials.size() - 1];

		
				_medials[medialCounter].nextMedial = &_medials[medialCounter + 1];
				_medials[medialCounter + 1].prevMedial = &_medials[medialCounter];

				_medials[medialCounter].type = zMedialType::Bottom;


				tmp_medials.clear();
				medialCounter++;
			}

			for (int i = 0; i < numFeatureLoops; i++)
			{
				for (int j = 0; j < numStride; j++) heSec = heSec.getPrev().getPrev().getSym();
				tmp_features[i].push_back(heSec);

			}

			if (he.getVertex().onBoundary())
			{
				// add feature graph
				_medials[medialCounter - 3].o_features.assign(numFeatureLoops, zObjGraph());
				_medials[medialCounter - 2].o_features.assign(numFeatureLoops, zObjGraph());
				_medials[medialCounter - 1].o_features.assign(numFeatureLoops, zObjGraph());

				for (int i = 0; i < numFeatureLoops; i++)
				{
					//printf("\n feature graph ");
					createGraphFromHalfEdgeLoop(_medials[medialCounter - 3].o_features[i], tmp_features[i], zBLUE);
					createGraphFromHalfEdgeLoop(_medials[medialCounter - 2].o_features[i], tmp_features[i], zBLUE);
					createGraphFromHalfEdgeLoop(_medials[medialCounter - 1].o_features[i], tmp_features[i], zBLUE);
					tmp_features[i].clear();
				}

				// walk
				he = he.getSym();

			}
			else
			{
				he = he.getNext().getSym().getNext();
			}

		} while (he != startHe);

		_medials[0].prevMedial = &_medials[medialCounter - 1];
		_medials[medialCounter - 1].nextMedial = &_medials[0];
		// top loop
		he = he_top_Start;

		int startTopInd = medialCounter;
		do
		{
			he.getEdge().setColor(zYELLOW);
			zItMeshHalfEdge heSec = he;

			// store
			tmp_medials.push_back(he);
			_He_Medial_Map[he.getId()] = medialCounter;

			for (int i = 0; i < numFeatureLoops; i++)
			{
				for (int j = 0; j < numStride; j++) heSec = heSec.getPrev().getPrev().getSym();
				tmp_features[i].push_back(heSec);
			}

			// add medial and feature graph
			if (he.getVertex().getValence() > 4)
			{
				//printf("\n medial graph ");
				createGraphFromHalfEdgeLoop(_medials[medialCounter].o_medial, tmp_medials, zYELLOW);

				_medials[medialCounter].o_features.assign(numFeatureLoops, zObjGraph());
				for (int i = 0; i < numFeatureLoops; i++)
				{
					createGraphFromHalfEdgeLoop(_medials[medialCounter].o_features[i], tmp_features[i], zGREEN);
					tmp_features[i].clear();
				}

				_medials[medialCounter].meshHEStartEnd.assign(2, zItMeshHalfEdge());
				_medials[medialCounter].meshHEStartEnd[0] = tmp_medials[0];
				_medials[medialCounter].meshHEStartEnd[1] = tmp_medials[tmp_medials.size()-1];


				_medials[medialCounter].nextMedial = &_medials[medialCounter + 1];
				_medials[medialCounter + 1].prevMedial = &_medials[medialCounter];

				_medials[medialCounter].type = zMedialType::Top;


				tmp_medials.clear();
				medialCounter++;
			}

			he = he.getNext().getSym().getNext();

		} while (he != he_top_Start);
		_medials[startTopInd].prevMedial = &_medials[medialCounter - 1];
		_medials[medialCounter - 1].nextMedial = &_medials[startTopInd];
		//  set sym medial indicies
		for (int i = 0; i < _He_Medial_Map.size(); i += 2)
		{
			if (_He_Medial_Map[i] != -1)
			{
				//printf("\n %i %i | %i %i ", i, _He_Medial_Map[i], i + 1, _He_Medial_Map[i + 1]);
				_medials[_He_Medial_Map[i]].symMedial = &_medials[_He_Medial_Map[i + 1]];
				_medials[_He_Medial_Map[i + 1]].symMedial = &_medials[_He_Medial_Map[i]];
			}
		}


		for (int i = 0; i < _medials.size(); i++)
		{
			zFnNurbsCurve fnN(_medials[i].o_medial_nurbs);
			fnN.create(_medials[i].o_medial, 0, 3, false, true, 100);

			_medials[i].o_features_nurbs.clear();
			_medials[i].o_features_nurbs.assign(_medials[i].o_features.size(), zObjNurbsCurve());
			for (int f = 0; f < _medials[i].o_features.size(); f++)
			{
				zFnNurbsCurve fn(_medials[i].o_features_nurbs[f]);
				fn.create(_medials[i].o_features[f],0, 3, false, true, 100);
			}
		}

		/*for (auto& m : _medials)
		{
			printf("\n medial HE %i - %i \n", m.meshHEStartEnd[0].getId(), m.meshHEStartEnd[1].getId());
			cout << m.meshHEStartEnd[0].getStartVertex().getPosition() << endl;
			cout << m.meshHEStartEnd[1].getVertex().getPosition() << endl;
		}*/

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeInterfacePlanesOnMedial(zHalfMedial& medialInner, zHalfMedial& medialOuter, float interfaceAngle, float dividLength)
	{

		zFnGraph fnGraph(medialInner.o_medial);
		zFnNurbsCurve fnMedialNurbs(medialInner.o_medial_nurbs);


		int divisionCount = ceil( fnMedialNurbs.getLength() / dividLength);
		
		divisionCount = divisionCount % 2 != 0? divisionCount : divisionCount + 1;


		medialInner.o_splitPlanes.clear();
		medialOuter.o_splitPlanes.clear();

		//outPlanes.clear();
		//outPlanes.assign(divisionCount + 1, zObjPlane());

		zFnMesh fnMesh(*o_inMesh_inner);
		zVectorArray meshNormals;
		fnMesh.getVertexNormals(meshNormals);




		zItMeshVertexArray vItStartEnd_inner = { medialInner.meshHEStartEnd[0].getStartVertex() , medialInner.meshHEStartEnd[1].getVertex() };
		zItMeshVertexArray vItStartEnd_outer = { medialOuter.meshHEStartEnd[0].getStartVertex() , medialOuter.meshHEStartEnd[1].getVertex() };
		zObjPlaneArray o_planeStartEnd;
		o_planeStartEnd.assign(2, zObjPlane());
		zVector zAxisUp(0, 0, 1);
		zVector prevVec = medialInner.prevMedial->meshHEStartEnd[1].getVector();
		zVector nextVec = medialInner.nextMedial->meshHEStartEnd[0].getVector();
		prevVec.normalize();
		nextVec.normalize();
		zFnPlane fnPlane;

		for (int i = 0; i < 2; i++)
		{
			zItMeshVertex vIt_inner = vItStartEnd_inner[i];
			zItMeshVertex vIt_outer = vItStartEnd_outer[i];
			//get bisector
			zVector heVec = medialInner.meshHEStartEnd[i].getVector();
			zVector bisector = heVec + prevVec;

			bool isHighValence = vIt_outer.getValence() > 4;
			zVector zAxis = (zAxisUp.angle(heVec) > 90) ? zVector(0, 0, -1) : zVector(0, 0, 1);
			zPoint o = vIt_outer.getPosition();
			float angle_bisectorMeshN = zAxis.angle360(bisector, meshNormals[vIt_inner.getId()]);
			bool checkAngle = (abs(angle_bisectorMeshN - interfaceAngle)) < (abs(angle_bisectorMeshN - (interfaceAngle * 2)));
			float rotateAngle = checkAngle ? -interfaceAngle : interfaceAngle;
			zVector vec1 = isHighValence ? zAxis.rotateAboutAxis(meshNormals[vIt_inner.getId()], rotateAngle) : zAxis;
			zVector vec2 = vIt_outer.getPosition() - vIt_inner.getPosition();  //vector between outer mesh and inner mesh
			zVector normal = isHighValence ? (checkAngle ? (vec1 ^ vec2) : (vec2 ^ vec1)) : (zAxis);

			fnPlane = zFnPlane(o_planeStartEnd[i]);
			bool chk = fnPlane.createFromNormal(o, normal, vec2);

			if (!chk)
			{
				printf("\n ERROR IN CREATING PLANE!");
				throw;
			}
		}

		//zFnGraph fnGraph(medialInner.o_medial);

		////get point on graph
		//zFloatArray params = { 0, 0.25, 0.5, 0.75, 1 };
		//zPointArray origins;
		//zIntArray eIndices;

		bool isArch = medialInner.meshHEStartEnd[0].getStartVertex().getValence() > 4 && medialInner.meshHEStartEnd[1].getVertex().getValence() > 4;


		//zObjNurbsCurve nurbs_inner, nurbs_outer;
		zFnNurbsCurve fnNurbs_inner(medialInner.o_medial_nurbs);
		zFnNurbsCurve fnNurbs_outer(medialOuter.o_medial_nurbs);
		//fnNurbs_inner.create(medialInner.o_medial, 0, 3, false, true, 50);
		//fnNurbs_outer.create(medialOuter.o_medial, 0, 3, false, true, 50);
		zPointArray origins_inner, origins_outer;
		zDoubleArray originsParms_inner, originsParms_outer;
		fnNurbs_inner.divideByCount(divisionCount, origins_inner, originsParms_inner);
		//fnNurbs_inner.divideByLength(fnNurbs_inner.getLength()/ divisionCount, origins_inner, originsParms_inner);
		fnNurbs_outer.divideByCount(divisionCount, origins_outer, originsParms_outer);
		//printf("\n division %i - %i - %i - %i", origins_inner.size(), origins_outer.size(), originsParms_inner.size(), originsParms_outer.size());


		medialInner.o_splitPlanes.assign(origins_outer.size(), zObjPlane());
		medialOuter.o_splitPlanes.assign(origins_outer.size(), zObjPlane());
		
		zObjPlaneArray tempPlanes;
		tempPlanes.assign(origins_outer.size(), zObjPlane());


		for (int i = 0; i < origins_outer.size(); i++)
		{
			zPoint pt0 = fnNurbs_outer.getPointAt(originsParms_outer[i]);
			zPoint pt1;

			double closestParam;
			fnNurbs_inner.closestPoint(pt0, pt1, closestParam);

			zVector tan0 = fnNurbs_outer.getTangentAt(originsParms_outer[i]);
			zVector tan1 = fnNurbs_inner.getTangentAt(closestParam);


			tan0.normalize();
			tan1.normalize();

			zVector vec1 = pt1 - pt0;
			zVector normal = isArch ? (tan0 + tan1) / 2 : zAxisUp;



			//zPoint o_inner = fnNurbs_inner.getPointAt(originsParms_inner[i]);
			//zPoint o_outer = fnNurbs_outer.getPointAt(originsParms_outer[i]);

			//
			////fnNurbs_inner.closestPoint(o_inner, o_outer, closestParam);
			////fnNurbs_outer.closestPoint(o_outer, o_inner, closestParam);


			//zVector tan_inner = fnNurbs_inner.getTangentAt(originsParms_inner[i]);
			//zVector tan_outer = fnNurbs_outer.getTangentAt(originsParms_outer[i]);
			//
			//
			//tan_outer = fnNurbs_outer.getTangentAt(closestParam);
			//
			////tan_inner = fnNurbs_inner.getTangentAt(closestParam);
			//
			//tan_inner.normalize();
			//tan_outer.normalize();

			//zVector vec1 = o_outer - o_inner;  //vector between outer mesh and inner mesh
			//zVector normal = isArch ? (tan_inner + tan_outer) / 2 : zAxisUp;

			zFnPlane fnPlane(tempPlanes[i]);
			bool chkPlane =  fnPlane.createFromNormal(pt0, normal, vec1);

			if (!chkPlane)
			{
				printf("\n ERROR IN CREATING PLANE! _2_ \n");
				cout << pt0 << endl << normal << endl << vec1 << endl;
				cout << pt0 << endl << pt1;
				throw;

			}

		}

		if (medialInner.type == zMedialType::Top)
		{
			zFnPlane fnPlane(tempPlanes[0]);
			fnPlane.create(fnPlane.getOrigin(), fnPlane.getYAxis(), zAxisUp*-1);
			//fnPlane.createFromNormal(fnPlane.getOrigin(), zAxisUp, fnPlane.getYAxis());

			fnPlane = zFnPlane(tempPlanes[tempPlanes.size() - 1]);
			fnPlane.create(fnPlane.getOrigin(), fnPlane.getYAxis(), zAxisUp*-1);
			

			//fnPlane.createFromNormal(fnPlane.getOrigin(), zAxisUp, fnPlane.getYAxis());
			
		}
		else
		{
			//printf("\n tempPlanes %i - %i ", origins_outer.size(), tempPlanes.size());
			tempPlanes[0] = o_planeStartEnd[0];
			tempPlanes[tempPlanes.size() - 1] = o_planeStartEnd[1];
		}

		for (int i = 0; i < tempPlanes.size(); i++)
		{
			medialInner.o_splitPlanes[i] = tempPlanes[i];
			medialOuter.o_splitPlanes[i] = tempPlanes[i];
		}



		




		//getPointsAtGraphParams(medialInner.o_medial, params, origins, eIndices);
		//interpolatePlanes(o_planeStartEnd[0], o_planeStartEnd[1], params, outPlanes, origins);

		//outPlanes.assign(2, zObjPlane());





	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computePlaneIntersectionOnMedial(zHalfMedial& medialInner, zHalfMedial& medialOuter, zObjPlaneArray& inPlanes, zObjGraphArray& outGraph)
	{
		double total;
		std::chrono::time_point<std::chrono::system_clock> s_Time, e_Time;

		s_Time = std::chrono::system_clock::now();

		vector<zPointArray> polylines;
		vector<zFloatArray> shatterParams_inner, shatterParams_outer;
		shatterParams_inner.assign(medialInner.o_features.size(), zFloatArray());
		shatterParams_outer.assign(medialInner.o_features.size(), zFloatArray());

		//zObjNurbsCurveArray features_inner, features_outer;
		//features_inner.assign(medialInner.o_features.size(), zObjNurbsCurve());
		//features_outer.assign(medialInner.o_features.size(), zObjNurbsCurve());

		zFnPlane fnPlane;
		zFnNurbsCurve fncrv;
		zPointArray tempPts;
		zFloatArray tempPars;

		bool isArch = medialInner.meshHEStartEnd[0].getStartVertex().getValence() > 4 && medialInner.meshHEStartEnd[1].getVertex().getValence() > 4;
		bool isTop = medialInner.type == zMedialType::Top;

		//get the parameters to shatter each feature curve
		for (int featureId = 0; featureId < medialInner.o_features.size(); featureId++)
		{

			//printf("\n isArch %d - isTop %i", isArch, medialInner.type == zMedialType::Top);
			//printf("\n o_splitPlanes %i inner-outer - %i - %i", featureId, medialInner.o_splitPlanes.size(), medialOuter.o_splitPlanes.size());


			for (int planeID = 0; planeID < medialInner.o_splitPlanes.size() ; planeID++)
			{
				if (!isTop)
				{
					if (!isArch)
					{
						if (medialInner.meshHEStartEnd[0].getStartVertex().getValence() > 4 && planeID == 0)
						{
							//printf("\n continue 1 %i \n", planeID);
							continue;
						}
						else if (planeID = medialInner.o_features.size() - 1)
						{
							//printf("\n continue 2 %i \n", planeID);

							continue;
						}
					}
				}
				else
				{
					if (planeID == 0 || planeID == medialInner.o_splitPlanes.size() - 1)
					{
						//printf("\n >>>  continue 3 %i \n", planeID);
						continue;
					}
				
				}

				


				int index = -1;

				curvePlaneIntersection(medialInner.o_features_nurbs[featureId], medialInner.o_splitPlanes[planeID], tempPts, tempPars);
				index = core.getClosest_PointCloud(medialInner.o_splitPlanes[planeID].origin, tempPts);
				//printf("\n curvePlaneIntersection_inner %i %i", index, tempPts.size());

				if (index != -1)
				{
					shatterParams_inner[featureId].push_back(tempPars[index]);
				}
				else
				{
					printf(" \n INTERSECTION FAILD! _1_ %i", featureId);

					cout << endl << medialInner.o_splitPlanes[planeID].origin;
					cout << endl << medialInner.o_splitPlanes[planeID].normal;

					zFnGraph fngraph(medialInner.o_features[featureId]);
					string path = "data/NatPower/ERROR/graph_error_inner_" + to_string(featureId) + ".json";

					

					//printf("\n INTERSECTION FAILD! _1_");
					throw;
				}

				tempPts.clear();
				tempPars.clear();

				curvePlaneIntersection(medialOuter.o_features_nurbs[featureId], medialInner.o_splitPlanes[planeID], tempPts, tempPars);
				index = core.getClosest_PointCloud(medialOuter.o_splitPlanes[planeID].origin, tempPts);
				//printf("\n curvePlaneIntersection_outer %i %i", index, tempPts.size());

				if (index != -1)
				{
					shatterParams_outer[featureId].push_back(tempPars[index]);
				}
				else
				{
					printf(" \n INTERSECTION FAILD! _2_ %i %i", featureId, planeID);
					cout << endl << medialOuter.o_splitPlanes[planeID].origin;
					cout << endl << medialOuter.o_splitPlanes[planeID].normal;

					zFnGraph fngraph(medialOuter.o_features[featureId]);
					string path = "data/NatPower/ERROR/graph_error_" + to_string(featureId);

					fngraph.to(path, zJSON);


					throw;

				}

			}

			//printf("\n shatterParams %i inner-outer - %i - %i \n",featureId, shatterParams_inner[featureId].size(), shatterParams_outer[featureId].size());


		}

		e_Time = std::chrono::system_clock::now();
		total = std::chrono::duration_cast<std::chrono::milliseconds>(e_Time - s_Time).count();

		//printf("\n medial __ intersection computed! Compute time %1.4f on millisecond",  total);

		e_Time = std::chrono::system_clock::now();
		s_Time = std::chrono::system_clock::now();

		//for (int i = 0; i < shatterParams_inner[0].size(); i++)
		//{
		//	//printf("\n params in_out %1.4f - %1.4f", shatterParams_inner[0][i], shatterParams_outer[0][i]);
		//}

		int sectionCount = shatterParams_inner[0].size() + 1;
		int featureCount = medialInner.o_features_nurbs.size();

		int subDivisionCount = 8;

		vector<zObjNurbsCurveArray> featureTrimed_inner, featureTrimed_outer;
		featureTrimed_inner.assign(sectionCount, zObjNurbsCurveArray());
		featureTrimed_outer.assign(sectionCount, zObjNurbsCurveArray());

		vector<vector<zPointArray>>blockPolylines;
		blockPolylines.assign(subDivisionCount, vector<zPointArray>());


		//shatterParams_inner organized as: for each featureCurve, list of params for shattering
		//featureTrimed_inner organized as: for each shattered section/param, we have multiple featureCurve 
		//blockPolylines organized as: for each shattered section (block), we have a list of polylines based on the division of the shatter


		for (int sectionId = 0; sectionId < sectionCount; sectionId++)
		{
			featureTrimed_inner[sectionId].assign(featureCount, zObjNurbsCurve());
			featureTrimed_outer[sectionId].assign(featureCount, zObjNurbsCurve());
			blockPolylines[sectionId].assign(subDivisionCount, zPointArray());

			for (int featureId = 0; featureId < featureCount; featureId++)
			{
				featureTrimed_inner[sectionId][featureId] = zObjNurbsCurve(medialInner.o_features_nurbs[featureId]);
				featureTrimed_outer[sectionId][featureId] = zObjNurbsCurve(medialOuter.o_features_nurbs[featureId]);

				double t0, t1;
				ON_Interval interval;
				bool trimChk;

				fncrv = zFnNurbsCurve(featureTrimed_inner[sectionId][featureId]);
				t0 = sectionId != 0 ?					shatterParams_inner[featureId][sectionId - 1]	: fncrv.getDomain().min ;
				t1 = sectionId != sectionCount - 1 ?	shatterParams_inner[featureId][sectionId]		: fncrv.getDomain().max;

				if (abs(t0-t1) < EPS)
				{
					continue;
				}

				zPointArray tempPts;
				zDoubleArray tempPtsParams;
				fncrv.divideByCount(subDivisionCount, tempPts, tempPtsParams);


				for (int sd = 0; sd < tempPts.size(); sd++)
				{
					blockPolylines[sectionId][sd];

				}





				interval.Set(t0, t1);
				trimChk = fncrv.getRawON_Curve()->Trim(interval);
				//fncrv.computeSubCurve

				if (!trimChk)
				{
					printf("\n \n TRIM FAILED! _1_");
					printf("\n featureTrimed_inner[%i][%i] = %1.4f  ", sectionId, featureId, shatterParams_inner[sectionId][featureId]);
					printf("\n domain() = %1.4f,%1.4f  ", t0, t1);
					printf("\n fncrv.getDomain() = %1.4f,%1.4f  ", fncrv.getDomain().min, fncrv.getDomain().max);

					throw;
			
				}

				fncrv = zFnNurbsCurve(featureTrimed_outer[sectionId][featureId]);
				t0 = sectionId != 0 ?					shatterParams_outer[featureId][sectionId - 1]	: fncrv.getDomain().min;
				t1 = sectionId != sectionCount - 1 ?	shatterParams_outer[featureId][sectionId]		: fncrv.getDomain().max;
				if (abs(t0 - t1) < EPS)
				{
					continue;
				}
				interval.Set(t0, t1);
				trimChk = fncrv.getRawON_Curve()->Trim(interval);

				if (!trimChk)
				{
					printf("\n \n TRIM FAILED! _2_");
					printf("\n featureTrimed_inner[%i][%i] = %1.4f  ", sectionId, featureId, shatterParams_inner[sectionId][featureId]);
					printf("\n domain() = %1.4f,%1.4f  ", t0, t1);
					printf("\n fncrv.getDomain() = %1.4f,%1.4f  ", fncrv.getDomain().min, fncrv.getDomain().max);

					throw;
				}

			}
		}

		e_Time = std::chrono::system_clock::now();
		total = std::chrono::duration_cast<std::chrono::milliseconds>(e_Time - s_Time).count();

		//printf("\n medial __ shatter computed! Compute time %1.4f on millisecond", total);
		//printf("\n featureTrimed inner-outer - %i - %i", shatterParams_inner.size(), shatterParams_outer.size());

		//divide each shattered curve




	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::createNurbsFromMedial(zHalfMedial& medial, zObjNurbsCurveArray& outCurve)
	{
		outCurve.clear();
		outCurve.assign(medial.o_features.size(), zObjNurbsCurve());
		zFnNurbsCurve fnNurbs;
		for (int i = 0; i < outCurve.size(); i++)
		{
			fnNurbs = zFnNurbsCurve(outCurve[i]);
			fnNurbs.create(medial.o_features[i], 0, 3, false, true, 100);
		}
	}

	/*
	* GENERIC METHODS
	*/

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints)
	{
		outPoints.clear();

		zScalarArray vertexScalars;

		zPoint O(inPlane(3, 0), inPlane(3, 1), inPlane(3, 2));
		zVector N(inPlane(2, 0), inPlane(2, 1), inPlane(2, 2));

		for (zItGraphVertex v(o_graph); !v.end(); v++)
		{
			zPoint P = v.getPosition();
			float minDist_Plane = core.minDist_Point_Plane(P, O, N);
			vertexScalars.push_back(minDist_Plane);
		}

		zPointArray contourPoints;
		isoContour(o_graph, vertexScalars, 0.0, contourPoints);

		if (closestPoint)
		{
			float dist = 1000000;
			zPoint cPoint;
			for (auto& p : contourPoints)
			{
				if (p.distanceTo(O) < dist)
				{
					dist = p.distanceTo(O);
					cPoint = p;
				}

			}

			outPoints.push_back(cPoint);
		}
		else outPoints = contourPoints;

	}

	ZSPACE_TOOLSETS_INLINE zPoint zTsNatPower::getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh)
	{
		float scaleVal = core.ofMap(threshold, thresholdLow, thresholdHigh, 0.0000f, 1.0000f);

		zVector e = vertex_higher - vertex_lower;
		double edgeLen = e.length();
		e.normalize();

		return (vertex_lower + (e * edgeLen * scaleVal));
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints)
	{
		zFnGraph fnGraph(o_graph);
		zPoint* vPositions = fnGraph.getRawVertexPositions();

		contourPoints.clear();

		for (zItGraphEdge e(o_graph); !e.end(); e++)
		{
			zIntArray eVerts;
			e.getVertices(eVerts);

			float s0 = vertexScalars[eVerts[0]];
			float s1 = vertexScalars[eVerts[1]];

			bool contour = false;
			if (s0 <= threshold && s1 >= threshold)contour = true;
			if (s0 >= threshold && s1 <= threshold)contour = true;

			if (!contour) continue;

			zPoint v0 = vPositions[eVerts[0]];
			zPoint v1 = vPositions[eVerts[1]];

			zPoint pos1 = getContourPosition(threshold, v1, v0, s1, s0);
			contourPoints.push_back(pos1);
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::interpolatePlanes(zObjPlane startPlane, zObjPlane endPlane, float t, zObjPlane& outPlane)
	{

		zFnPlane fnS(startPlane);
		zFnPlane fnE(endPlane);
		zFnPlane fnO(outPlane);

		zVector xs = fnS.getXAxis();
		zVector ys = fnS.getYAxis();
		zVector zs = fnS.getNormal();

		zVector xe = fnE.getXAxis();
		zVector ye = fnE.getYAxis();
		zVector ze = fnE.getNormal();

		zVector os = fnS.getOrigin();
		zVector oe = fnE.getOrigin();

		// Perform linear interpolation for each axis and the origin
		zVector xn = xs * (1 - t) + xe * t;
		zVector yn = ys * (1 - t) + ye * t;
		zVector zn = zs * (1 - t) + ze * t;
		zVector on = os * (1 - t) + oe * t; // Interpolated origin

		xn.normalize();
		yn.normalize();
		zn.normalize();
		//fnO.create(on, xn, yn);
		fnO.create(on, zn, yn);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::interpolatePlanes(zObjPlane startPlane, zObjPlane endPlane, zFloatArray& params, zObjPlaneArray& outPlanes, zPointArray originAtParams)
	{

		zFnPlane fnS(startPlane);
		zFnPlane fnE(endPlane);

		zVector xs = fnS.getXAxis();
		zVector ys = fnS.getYAxis();
		zVector zs = fnS.getNormal();

		zVector xe = fnE.getXAxis();
		zVector ye = fnE.getYAxis();
		zVector ze = fnE.getNormal();

		zVector os = fnS.getOrigin();
		zVector oe = fnE.getOrigin();

		outPlanes.clear();
		outPlanes.assign(params.size(), zObjPlane());

		bool interpolateOrigin = originAtParams.size() != params.size();
		if (interpolateOrigin)
		{
			originAtParams.clear();
			originAtParams.assign(params.size(), zPoint());
		}
		for (int i = 0; i < params.size(); i++)
		{
			float t = params[i];
			zFnPlane fnO(outPlanes[i]);

			// Perform linear interpolation for each axis and the origin
			zVector xn = xs * (1 - t) + xe * t;
			zVector yn = ys * (1 - t) + ye * t;
			zVector zn = zs * (1 - t) + ze * t;
			if (interpolateOrigin) originAtParams[i] = os * (1 - t) + oe * t;
			zVector on = originAtParams[i]; // Interpolated origin
			//zVector on = interpolateOrigin? os * (1 - t) + oe * t : originAtParams[i]; // Interpolated origin
			xn.normalize();
			yn.normalize();
			zn.normalize();
			//fnO.create(on, xn, yn);
			fnO.createFromNormal(on, zn, yn);
		}
		
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::createNurbsFromHEArray(zItMeshHalfEdgeArray& heArray, zObjNurbsCurve& outCurve)
	{
		zFnNurbsCurve fnCurve(outCurve);
		zObjGraph graph;
		createGraphFromHalfEdgeLoop(graph, heArray);
		fnCurve.create(graph, 0, 3, false, true, 100);
		printf("\n nurbsCreated %i ", fnCurve.numControlVertices());
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::createGraphFromHalfEdgeLoop(zObjGraph& o_graph, zItMeshHalfEdgeArray& heLoop, zColor eColor)
	{
		zPointArray positions;
		zIntArray eConnects;

		for (int i = 0; i < heLoop.size(); i++)
		{
			positions.push_back(heLoop[i].getStartVertex().getPosition());

			if (positions.size() > 1)
			{
				eConnects.push_back(positions.size() - 2);
				eConnects.push_back(positions.size() - 1);
			}

			if (i == heLoop.size() - 1)
			{
				positions.push_back(heLoop[i].getVertex().getPosition());

				eConnects.push_back(positions.size() - 2);
				eConnects.push_back(positions.size() - 1);
			}
		}

		zFnGraph fnGraph(o_graph);
		fnGraph.create(positions, eConnects);

		//fnGraph.setEdgeColor(eColor);
		//printf("\n g: %i %i ", fnGraph.numVertices(), fnGraph.numEdges());
	}

	ZSPACE_TOOLSETS_INLINE float zTsNatPower::getGraphCumulativeLengths(zObjGraph& graph, zFloatArray& cumulativeLengths, zFloatArray& cumulativeLengthsReparms)
	{
		zFnGraph fn(graph);

		zDoubleArray edgeLengths ;
		double totalLength = fn.getEdgeLengths(edgeLengths);

		cumulativeLengths.clear();
		cumulativeLengthsReparms.clear();

		cumulativeLengths.assign(edgeLengths.size(), -1);
		cumulativeLengthsReparms.assign(edgeLengths.size(), -1);

		float cumulativeSum = 0.0f;

		for (int i = 0; i < edgeLengths.size(); i++)
		{
			cumulativeSum += edgeLengths[i];
			cumulativeLengths[i] = cumulativeSum;
			cumulativeLengthsReparms[i] = cumulativeSum / totalLength;
		}

		return cumulativeSum;


	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getPointsAtGraphParams(zObjGraph& graph, zFloatArray& params, zPointArray& outPoints, zIntArray& outEdgeIndex)
	{
		//The following only works if the graph edges are sorted, and all in the same direction. i.e. if the end of the first edge is the start of the second
		//The following only works for graphs with maximum 2 vertex valence

		zFnGraph fn(graph);

		zFloatArray accLengeths, accLengthsRep;

		float totalLength = getGraphCumulativeLengths(graph, accLengeths, accLengthsRep);

		outPoints.clear();
		outPoints.assign(params.size(), zPoint());

		outEdgeIndex.clear();
		outEdgeIndex.assign(params.size(), -1);

		for (int i = 0; i < params.size(); i++)
		{

			if (params[i] <= 1)
			{
				auto it = std::lower_bound(accLengthsRep.begin(), accLengthsRep.end(), params[i]);
				int index = std::distance(accLengthsRep.begin(), it);

				if (index >= fn.numEdges())
				{
					printf("\n index out of range!! skipping!!");
					continue;
				}
				

				float moveLength = (accLengthsRep[index] - params[i]) * totalLength;

				
				zItGraphEdge edge(graph, index);
				float eLength = edge.getLength();
				zVector eDir = edge.getVector();
				eDir.normalize();
				eDir *= moveLength;

				zPointArray pts;
				edge.getVertexPositions(pts);

				outPoints[i] = pts[1] - eDir;
				outEdgeIndex[i] = index;
			}
			else
			{
				printf("\n failed to get point! param[%i] = %1.4f > 1", i, params[i]);
			}

		}




	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getNurbsMidPointFromGraphVertex(zObjNurbsCurve& crv, zObjGraph& graph, double& outPar, zPoint& outPoint)
	{
		zFnGraph fnGraph(graph);
		zFnNurbsCurve fnNurbs(crv);
		int verIndex = ceil (fnGraph.numVertices()/2);
		outPoint = zItGraphVertex(graph, verIndex).getPosition();
		fnNurbs.closestPoint(outPoint, outPoint, outPar);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::curvePlaneIntersection(zObjNurbsCurve& o_curve, zObjPlane& o_plane, zPointArray& intersectionPoint, zFloatArray& crvParameters)
	{
		zFnPlane fnp(o_plane);
		ON_Plane* plane = fnp.getRawON_Plane();
		ON_SimpleArray<ON_X_EVENT> intersectionEvents;
		zFnNurbsCurve fnCurve(o_curve);
		fnCurve.setDomain();
		ON_NurbsCurve* onCurve = fnCurve.getRawON_Curve();
		//cout << "\n fnp" << o_plane.origin << endl;
		//cout << "\n plane" << plane->origin.x <<"," << plane->origin.y <<","<< plane->origin.z << endl;

		int intersection_result = onCurve->IntersectPlane(plane->plane_equation, intersectionEvents, 0.001);

		//ON_SimpleArray<ON_X_EVENT>* rc = NULL;
		//if (onCurve != NULL && plane != NULL)
		//{
		//	rc = new ON_SimpleArray<ON_X_EVENT>();
		//	if (onCurve->IntersectPlane(plane.plane_equation, *rc) < 1)
		//	{
		//		// no intersections found. No need to create a list of intersections
		//		delete rc;
		//		rc = NULL;
		//	}
		//}


		intersectionPoint.clear();
		crvParameters.clear();

		//printf("\n intersection %i %i", intersection_result,  intersectionEvents.Count());
		for (int i = 0; i < intersectionEvents.Count(); i++)
		{
			ON_3dPoint onpt = intersectionEvents[i].m_A[0];
			intersectionPoint.push_back(zPoint(onpt.x, onpt.y, onpt.z));
			crvParameters.push_back(intersectionEvents[i].m_a[0]);

		}

	}


}