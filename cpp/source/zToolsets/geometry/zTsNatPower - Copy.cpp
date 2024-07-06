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


#include<headers/zToolsets/geometry/zTsNatPower.h>

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


	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getMedialHE()
	{

		computeMedialGraph(*o_inMesh_inner, startVertexId_inner, meshStride, 1, feature1Map_inner, medials_inner);
		computeMedialGraph(*o_inMesh_outer, startVertexId_outer, meshStride , 1, feature1Map_outer, medials_outer);

		printf("\n medial computed!");

		o_interfacePlanes.assign(medials_inner.size(), zObjPlaneArray());
		for (int i = 0; i < medials_inner.size(); i++)
		{
			calculateInterfacePlanes(medials_inner[i], medials_outer[i], o_interfacePlanes[i]);
		}
		printf("\n interface computed!");

		//computeMedialHE(o_inMesh_inner, startVertexId_inner, medial_inner, medialMap_inner, feature1_inner, feature1Map_inner);
		//computeMedialHE(o_inMesh_outer, startVertexId_outer, medial_outer, medialMap_outer, feature1_outer, feature1Map_outer);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getFeatureNurbs(vector<zObjNurbsCurveArray> & outInnerFeature, vector<zObjNurbsCurveArray> & outOuterFeature)
	{
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
			/*zFnNurbsCurve fn(outInnerFeature[i]);
			printf("\n outInnerFeature[i]  : %i ", fn.numControlVertices());*/


		}

		/*outOuterFeature.clear();
		outOuterFeature.assign(feature1_outer.size(), zObjNurbsCurve());
		for (int i = 0; i < outOuterFeature.size(); i++)
		{
			createNurbsFromHEArray(feature1_outer[i], outOuterFeature[i]);

			zFnNurbsCurve fn(outOuterFeature[i]);
			printf("\n outOuterFeature[i]  : %i ", fn.numControlVertices());
		}*/
	}


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
		else printf("Able to read '%s' %i ", keyId, index);


		bool chk2 = core.json_readAttribute(oJson, keyStride, strid);
		outStride = chk2 ? strid : -1;
		if (!chk2) printf("\n json unable to read '%s'  ", keyStride);
		else printf("Able to read '%s' %i ", keyStride, strid);

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
	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeMedialGraph(zObjMesh& o_mesh, int startVID, int numFeatureLoops, int numStride, vector<int>& _He_Medial_Map, vector<zHalfMedial>& _medials)
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

		for (auto& m : _medials)
		{
			printf("\n medial HE %i - %i \n", m.meshHEStartEnd[0].getId(), m.meshHEStartEnd[1].getId());
			cout << m.meshHEStartEnd[0].getStartVertex().getPosition() << endl;
			cout << m.meshHEStartEnd[1].getVertex().getPosition() << endl;
		}

	}



	ZSPACE_TOOLSETS_INLINE void zTsNatPower::getMesh(zObjMesh& innerMesh, zObjMesh& outerMesh)
	{
		innerMesh = *o_inMesh_inner;
		outerMesh = *o_inMesh_outer;

	}

	

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::computeMedialHE(
		zObjMesh* oMesh, int startVID, 
		zItMeshHalfEdgeArray& outHEArray, zIntArray& outHEMap,
		vector<vector<zItMeshHalfEdgeArray>>& outFeatureHE, zIntArray& outFeatureMap
		)
	{
		zFnMesh fnMesh(*oMesh);
		outHEMap.clear();
		outHEMap.assign(fnMesh.numHalfEdges(), -1);

		outFeatureMap.clear();
		outFeatureMap.assign(fnMesh.numHalfEdges(), -1);

		outHEArray.clear();
		outFeatureHE.clear();
		outFeatureHE.assign(baysCount*2, vector<zItMeshHalfEdgeArray>());

		zItMeshHalfEdge heStart; //Start HE


		zItMeshVertex startVertex(*oMesh, startVID);
		zItMeshHalfEdgeArray startVCHE; // connected HE of start vertex

		startVertex.getConnectedHalfEdges(startVCHE);

		bool startHEFound = false;
		for (auto& e : startVCHE)
		{
			printf("\n e Id %i", e.getId());

			if (!e.getEdge().onBoundary())
			{
				heStart = e;
				startHEFound = true;
				break;
			}
		}

		if (!startHEFound)
		{
			printf("\n No internal HE was found!!");
			throw;
		}

		int safetyMaxWalk = 10000000; //Exist while loop if this number was reached. It means the loop was not able to exit
		int counter = 0;
		zItMeshHalfEdge heWalk = heStart; // HE 
		int map = 0;
		int featureCounter = 0;
		bool topSideHE = false;
		outFeatureMap[heWalk.getId()] = featureCounter;

		outFeatureHE[featureCounter].assign(meshStride + 1, zItMeshHalfEdgeArray());
		outFeatureHE[featureCounter + baysCount].assign((meshStride / 2) + 1, zItMeshHalfEdgeArray());

		
		

		while (counter < safetyMaxWalk)
		{
			counter++;

			if(map == 0) heWalk.getEdge().setColor(zMAGENTA);
			else heWalk.getEdge().setColor(zGREEN);
			
			if (heWalk.getStartVertex().checkValency(6))
			{
				topSideHE = !topSideHE;
			}

			outHEArray.push_back(heWalk);
			outHEMap[heWalk.getId()] = map;
			outHEMap[heWalk.getSym().getId()] = map == 0 ? 1 : 0;

			//zItMeshHalfEdgeArray heFeatures1;
			//heFeatures1.assign(meshStride + 1, heWalk);



			zItMeshHalfEdge heWalkFeature = heWalk; // HE 

			for (int i = 0; i < meshStride + 1; i++)
			{
				//outFeatureHE[featureCounter][i].push_back(heWalkFeature);

				//if (i != 0)
				//{
				//	if (map == 0) heWalkFeature.getEdge().setColor(zRED);
				//	else heWalkFeature.getEdge().setColor(zBLUE);
				//}

				//if (i != outFeatureHE[featureCounter].size()-1)
				//{
				//	if (map == 0) heWalkFeature.getFace().setColor(zBLUE);
				//	else heWalkFeature.getFace().setColor(zRED);
				//}

				////heWalkFeature = heWalkFeature.getSym().getNext().getNext();
				//heWalkFeature = heWalkFeature.getPrev().getPrev().getSym();
				////heFeatures1[i] = heWalkFeature;
				
			}


			if (topSideHE)
			{
				heWalkFeature = heWalk; // HE 
				outFeatureMap[heWalkFeature.getSym().getId()] = featureCounter + baysCount;



				//get the maximum number of strides 



				//outFeatureHE[featureCounter + baysCount].assign(meshStride + 1, heWalk);
				//outFeatureHE[featureCounter + baysCount].assign((meshStride/2)+1, heWalk);

				//printf("\n INSIDE0 outFeatureHE %i", outFeatureHE.size());


				//for (int i = 0; i < meshStride + 1; i++)
				for (int i = 0; i < (meshStride/2)+1; i++)
				{
					/*outFeatureHE[featureCounter + baysCount][i].push_back(heWalkFeature);
					heWalkFeature = heWalkFeature.getSym().getNext().getNext();
					if (i != 0)
					{
						if (map == 0) heWalkFeature.getEdge().setColor(zCYAN);
						else heWalkFeature.getEdge().setColor(zYELLOW);
					}
					if (i != outFeatureHE[featureCounter + baysCount].size() -1)
					{
						if (map == 0) heWalkFeature.getFace().setColor(zCYAN);
						else heWalkFeature.getFace().setColor(zYELLOW);
					}*/


					
					//heWalkFeature = heWalkFeature.getPrev().getPrev().getSym();
					//heFeatures1[i] = heWalkFeature;
					
				}
			}




			bool reachBoundary = heWalk.getVertex().onBoundary();
			if (reachBoundary)
			{
				heWalk = heWalk.getSym();
				map = map == 0 ? 1 : 0;

				featureCounter++;
			}
			else
			{
				heWalk = heWalk.getNext().getSym().getNext();
			}
			//heWalk = !reachBoundary ? heWalk.getNext().getSym().getNext() : heWalk.getSym(); 

			if (heWalk.getId() == heStart.getId()) break;
			else
			{
				outFeatureMap[heWalk.getId()] = featureCounter;
			}

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

	ZSPACE_TOOLSETS_INLINE void zTsNatPower::createNurbsFromHEMedial(zHalfMedial& medial, zObjNurbsCurveArray & outCurve)
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


	ZSPACE_TOOLSETS_INLINE zPoint zTsNatPower:: getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh)
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


	ZSPACE_TOOLSETS_INLINE void zTsNatPower::calculateInterfacePlanes(zHalfMedial& medialInner, zHalfMedial& medialOuter, vector<zObjPlane>& outPlanes)
	{
		int divisionCount = 5;
		outPlanes.clear();
		outPlanes.assign(divisionCount + 1, zObjPlane());

		zFnMesh fnMesh(*o_inMesh_inner);
		zVectorArray meshNormals;
		fnMesh.getVertexNormals(meshNormals);
		

		//zObjNurbsCurve nurbs_inner, nurbs_outer;

		//zFnNurbsCurve fnNurbs_inner(nurbs_inner);
		//zFnNurbsCurve fnNurbs_outer(nurbs_outer);

		//fnNurbs_inner.create(medialInner.o_medial, 0, 3, false, true, 50);
		//fnNurbs_outer.create(medialOuter.o_medial, 0, 3, false, true, 50);

		//zPointArray origins_inner, origins_outer ;
		//zDoubleArray originsParms_inner, originsParms_outer;
		//fnNurbs_inner.divideByCount(divisionCount, origins_inner, originsParms_inner);
		//fnNurbs_outer.divideByCount(divisionCount, origins_outer, originsParms_outer);

		//printf("\n division %i - %i - %i - %i", origins_inner.size(), origins_outer.size(), originsParms_inner.size(), originsParms_outer.size());

		//for (int i = 0; i < origins_outer.size(); i++)
		//{
		//	zVector tan_inner = fnNurbs_inner.getTangentAt(originsParms_inner[i]);
		//	zVector tan_outer = fnNurbs_outer.getTangentAt(originsParms_outer[i]);
		//	tan_inner.normalize();
		//	tan_outer.normalize();
		//	
		//	zPoint o_inner = fnNurbs_inner.getPointAt(originsParms_inner[i]);
		//	zPoint o_outer = fnNurbs_outer.getPointAt(originsParms_outer[i]);


		//	zVector vec1 = o_outer - o_inner;  //vector between outer mesh and inner mesh
		//	zVector normal = (tan_inner + tan_outer) / 2;
		//	
		//		zFnPlane fnPlane(outPlanes[i]);
		//	fnPlane.createFromNormal(o_outer, normal, vec1);


		//}


		zItMeshVertex vItStartInner = medialInner.meshHEStartEnd[0].getStartVertex();
		zItMeshVertex vItEndInner = medialInner.meshHEStartEnd[1].getVertex();

		zItMeshVertex vItStartOuter = medialOuter.meshHEStartEnd[0].getStartVertex();
		zItMeshVertex vItEndOuter = medialOuter.meshHEStartEnd[1].getVertex();
		zObjPlane o_startPlane, o_endPlane;
		zFnPlane fnPlane;

		//get 120 angle 
		zPoint zAxis(0, 0, 1);
		bool isHighValence;

		//plane at start
		
		
		//get bisector
		zItMeshVertex v0 = vItStartInner;
		zItMeshVertex v1 = vItStartOuter;
		zVector prevVec = medialInner.prevMedial->meshHEStartEnd[1].getVector();
		zVector nextVec = medialInner.nextMedial->meshHEStartEnd[0].getVector();
		prevVec.normalize();
		nextVec.normalize();


		////////////////
		zVector heVec = medialInner.meshHEStartEnd[0].getVector();
		zVector bisector = heVec + prevVec;

		isHighValence = v1.getValence() > 4;
		zAxis = (abs(zVector(0,0,1).angle(heVec) > 90))? zVector(0,0,-1) : zVector(0, 0, 1);
		zPoint o = v1.getPosition();
		float angle = zAxis.angle360(bisector, meshNormals[v0.getId()]);
		bool checkAngle = (abs(angle - 120)) < (abs(angle - 240));
		float rotateAngle = checkAngle ? -120.0 : 120.0;
		zVector vec1 = isHighValence? zAxis.rotateAboutAxis(meshNormals[v0.getId()], rotateAngle) : zAxis;
		zVector vec2 = v1.getPosition() - v0.getPosition();  //vector between outer mesh and inner mesh
		zVector normal = isHighValence ? checkAngle? vec1 ^ vec2 : vec2 ^ vec1 : zAxis;
		
		
		zPlane startPlane = core.getTransformFromOrigin_Normal(o, normal);
		fnPlane = zFnPlane(o_startPlane);
		bool chk = fnPlane.createFromNormal(o, normal, vec2);

		cout << endl << "start o" << o  << chk;
		cout << endl << "start o" << fnPlane.getOrigin() << endl;

		//plane at end
		 v0 = vItEndInner;
		 v1 = vItEndOuter;
		heVec = medialInner.meshHEStartEnd[1].getVector();
		 bisector = heVec + nextVec;

		isHighValence = v1.getValence() > 4;
		zAxis = (abs(zVector(0, 0, 1).angle(heVec) > 90)) ? zVector(0, 0, -1) : zVector(0, 0, 1);
		 o = v1.getPosition();
		 angle = zAxis.angle360(bisector, meshNormals[v0.getId()]);
		 checkAngle = (abs(angle - 120)) < (abs(angle - 240));
		 rotateAngle = checkAngle ? -120.0 : 120.0;
		 vec1 = isHighValence ? zAxis.rotateAboutAxis(meshNormals[v0.getId()], rotateAngle) : zAxis;
		 vec2 = v1.getPosition() - v0.getPosition();  //vector between outer mesh and inner mesh
		 normal = isHighValence ? checkAngle ? vec1 ^ vec2 : vec2 ^ vec1 : zAxis;

		
		zPlane endPlane = core.getTransformFromOrigin_Normal(o, normal);
		fnPlane = zFnPlane(o_endPlane);
		chk = fnPlane.createFromNormal(o, normal, vec2);
		cout << endl << "end o" << o << chk;
		cout << endl << "end o" << fnPlane.getOrigin() << endl;


		
		////o = vItEndOuter.getPosition();
		////vec1 = zAxis.rotateAboutAxis(meshNormals[vItEndInner.getId()], -120.0);
		////vec2 = vItEndOuter.getPosition() - vItEndInner.getPosition();  //vector between outer mesh and inner mesh
		////normal = vec2 ^ vec1;
		////zPlane endPlane = core.getTransformFromOrigin_Normal(o, normal);
		////fnPlane = zFnPlane(o_endPlane);
		////fnPlane.createFromNormal(o, normal);// , vec2);


		//outPlanes.assign(5, zObjPlane());
		for (int i = 0; i < outPlanes.size(); i++)
		{
			float t = i / ((outPlanes.size() - 1)*1.0) ;
			interpolatePlanes(o_startPlane, o_endPlane, t, outPlanes[i]);
		}

		//zRED for bay bottom , zYELLOW for bay top
		if (medialInner.meshHEStartEnd[0].getColor() == zRED)
		{
			//get the arch from the bay when the start and end vertex are high valence
			if (!(vItStartInner.getValence() > 4 && vItEndInner.getValence() > 4))
			{
				for (int i = 1; i < outPlanes.size()-1; i++)
				{
					zFnPlane fn(outPlanes[i]);
					//fn.create(outPlanes[i].origin, outPlanes[i].xAxis, zAxis ^ outPlanes[i].xAxis );
					//fn.create(outPlanes[i].origin, outPlanes[i].xAxis, zAxis ^ outPlanes[i].xAxis );
				}
			}
		}
		else
		{
			zFnPlane fn(outPlanes[0]);
			//fn.create(outPlanes[0].origin, outPlanes[0].xAxis, zAxis ^ outPlanes[0].xAxis );

			fn = zFnPlane(outPlanes[1]);
			//fn.create(outPlanes[1].origin, outPlanes[1].xAxis, zAxis ^ outPlanes[1].xAxis );
		}

		outPlanes.assign(2, zObjPlane());
		outPlanes[0] = o_startPlane;
		outPlanes[1] = o_endPlane;

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

	






}