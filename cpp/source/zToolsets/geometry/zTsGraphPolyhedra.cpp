// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Leo Bieling <leo.bieling@zaha-hadid.com>
//


#include<headers/zToolsets/geometry/zTsGraphPolyhedra.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsGraphPolyhedra::zTsGraphPolyhedra() {}

	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsGraphPolyhedra::~zTsGraphPolyhedra() {}

	//---- CREATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setFormGraphFromFile(string& _path, zFileTpye _type, bool _staticGeom)
	{
		// create graph
		zFnGraph fnGraph(o_formGraph);
		fnGraph.from(_path, _type, _staticGeom);

		// allocate memory
		o_convexHullMeshes.assign(fnGraph.numVertices(), zObjMesh());
		o_forceMeshes.assign(fnGraph.numVertices(), zObjMesh());
		//c_graphHalfEdge_dualCellFace.assign(fnGraph.numHalfEdges(), zIntPairArray());

		// color vertices
		for (zItGraphVertex v(o_formGraph); !v.end(); v++)
			if (v.getValence() > 1) v.setColor(zColor(0.0, 1.0, 0.0, 1.0));

		// color edges
		colorGraphEdges();

	}


	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setFormGraphFromMesh(zObjMesh& _inMeshObj, zVector& _verticalForce)
	{
		zFnMesh fnMesh(_inMeshObj);
		zFnGraph fnGraph(o_formGraph);

		// allocate memory
		o_convexHullMeshes.assign(fnMesh.numVertices(), zObjMesh());
		o_forceMeshes.assign(fnMesh.numVertices(), zObjMesh());
		//c_graphHalfEdge_dualCellFace.assign(fnMesh.numHalfEdges(), zIntPairArray());

		vector<int>edgeConnects;
		vector<zVector> vertexPositions;

		fnMesh.getEdgeData(edgeConnects);
		fnMesh.getVertexPositions(vertexPositions);

		for (zItMeshVertex v(_inMeshObj); !v.end(); v++)
		{
			zPoint newPos = v.getPosition() + (_verticalForce * -1);

			edgeConnects.push_back(v.getId());
			edgeConnects.push_back(vertexPositions.size());

			vertexPositions.push_back(newPos);
		}

		fnGraph.create(vertexPositions, edgeConnects, true);

		// color vertices
		for (zItGraphVertex v(o_formGraph); !v.end(); v++)
			if (v.getValence() > 1) v.setColor(zColor(0.0, 1.0, 0.0, 1.0));

		// color edges
		colorGraphEdges();

	}

	//---- GET METHODS

	ZSPACE_TOOLSETS_INLINE zObjGraph* zTsGraphPolyhedra::getRawFormGraph()
	{
		return &o_formGraph;
	}

	ZSPACE_TOOLSETS_INLINE zObjMeshPointerArray zTsGraphPolyhedra::getRawForceMeshes(int& numMeshes)
	{
		zObjMeshPointerArray out;
		numMeshes = 0;

		numMeshes = o_forceMeshes.size();

		if (numMeshes == 0)return out;

		for (auto& m : o_forceMeshes)
		{
			out.push_back(&m);
		}

		return out;
	}

	//---- CREATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::createForceMeshes()
	{
		int activeNodes = 0;
		double tol = 0.002;

		zFnGraph fnGraph(o_formGraph);

		forceCellFace_formHalfedge.clear();
		formHalfedge_forceCellFace.clear();

		for (zItGraphVertex g_v(o_formGraph); !g_v.end(); g_v++)
		{
			if (!g_v.checkValency(1))
			{
				nodeId.push_back(g_v.getId());

				//cout << "\n\t\t\tNODE_ID: " << g_v.getId();

				// make clean convex hull
				cleanConvexHull(g_v);

				// dual mesh from convex hull
				createForceMesh(g_v, tol);

				activeNodes++;
			}
		}
		cout << "\n////////////////////////////////////////////////////////////////////\n";
		cout << "convexhull mesh size: " << o_convexHullMeshes.size() << endl;
		cout << "dual mesh size: " << o_forceMeshes.size() << endl;
		cout << "active nodes: " << activeNodes << endl;

		//colors.insert(colors.end(),
		//	{
		//		zColor(1, 0, 0, 1), //red
		//		zColor(1, 1, 0, 1), //yellow
		//		zColor(0, 1, 0, 1), //green
		//		zColor(0, 1, 1, 1), //cyan
		//		zColor(0, 0, 1, 1), //blue
		//		zColor(1, 0, 1, 1), //magenta
		//		zColor(0.5, 0.5, 0.5, 1), //gray
		//		zColor(1, 0.5, 0, 1) //orange
		//	});


		// set edge and face centers
		/*for (auto& m : o_convexHullMeshes)
		{
			zPointArray edgeCenters, faceCenters;

			for (zItMeshEdge e(m); !e.end(); e++)
				edgeCenters.push_back(e.getCenter());

			for (zItMeshFace f(m); !f.end(); f++)
				faceCenters.push_back(f.getCenter());

			m.setEdgeCenters(edgeCenters);
			m.setFaceCenters(faceCenters);
		}

		for (auto& m : o_forceMeshes)
		{
			zPointArray edgeCenters, faceCenters;

			for (zItMeshEdge e(m); !e.end(); e++)
				edgeCenters.push_back(e.getCenter());

			for (zItMeshFace f(m); !f.end(); f++)
				faceCenters.push_back(f.getCenter());

			m.setEdgeCenters(edgeCenters);
			m.setFaceCenters(faceCenters);
		}

		zPointArray graphEdgeCenters;
		for (zItGraphEdge e(o_formGraph); !e.end(); e++) graphEdgeCenters.push_back(e.getCenter());
		o_formGraph.setEdgeCenters(graphEdgeCenters);*/

		// color connected faces
		colorDualFaceConnectivity();

		// snap and merge dual cells
		//zItGraphVertexArray(graphCenters);
		//fnGraph.getGraphEccentricityCenter(graphCenters);
		//snapDualCells(sortedGraphVertices, graphCenters);


		/*cout << "------MAP 1-------" << endl;

		for (auto itr = formHalfedge_forceCellFace.begin();
			itr != formHalfedge_forceCellFace.end(); itr++)
		{
			cout << "HE:" << itr->first << "-->" <<
				"F(node, face):" << (itr->second).first << "," << (itr->second).second << endl;
		}

		cout << "------MAP 2-------" << endl;

		for (auto itr = forceCellFace_formHalfedge.begin();
			itr != forceCellFace_formHalfedge.end(); itr++)
		{
			cout << "F(node, face):" << itr->first << "-->" <<
				"HE:" << itr->second << endl;
		}*/

	}

	//---- UPDATE METHODS

	ZSPACE_TOOLSETS_INLINE bool zTsGraphPolyhedra::equilibrium(bool& compTargets, float formWeight, float areaScale,  double dT, zIntergrationType type, float minMax_formEdge, float minMax_forceEdge, int numIterations, double angleTolerance, double areaTolerance, bool printInfo)
	{
		if (compTargets)
		{
			printf("\n compute targets ");
			computeTargets(formWeight, areaScale);
			compTargets = !compTargets;
		}

		bool out = checkDeviations(deviations[0], angleTolerance, deviations[1], areaTolerance, printInfo);

		// update diagrams

		if (formWeight != 1.0)
		{
			updateFormDiagram(minMax_formEdge, dT, type, numIterations);
		}

		if (formWeight != 0.0)
		{
			updateForceDiagram(minMax_forceEdge, dT, type, numIterations);
		}

		// check deviations	
		out = checkDeviations(deviations[0], angleTolerance, deviations[1], areaTolerance, printInfo);

		return out;

	}

	//---- PRIVATE CREATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::createForceMesh(zItGraphVertex& _graphVertex, double _tol)
	{
		
		zFnMesh fnMesh(o_convexHullMeshes[_graphVertex.getId()]);
		zIntArray inEdge_dualEdge;
		zIntArray dualEdge_inEdge;
		fnMesh.getDualMesh(o_forceMeshes[_graphVertex.getId()], inEdge_dualEdge, dualEdge_inEdge, true);

		zFnMesh fnDual(o_forceMeshes[_graphVertex.getId()]);
		fnDual.computeMeshNormals();

		// get connected edges
		zItGraphHalfEdgeArray connectedHalfEdges;
		_graphVertex.getConnectedHalfEdges(connectedHalfEdges);

		zFnMesh fnConHull(o_convexHullMeshes[_graphVertex.getId()]);

		// map positionId -> graphVertexId
		zPoint* hullPos = fnConHull.getRawVertexPositions();

		//cout << "num cHalEdges" << connectedHalfEdges.size() << endl;
		for (auto &he : connectedHalfEdges)
		{
			//zPoint hullP = (he.getVertex().checkValency(1)) ? he.getVertex().getPosition() : he.getCenter();

			zVector heVec = he.getVector();
			heVec.normalize();

			zPoint hullP = ((_graphVertex.getPosition() + heVec));

			for (int i = 0; i < fnConHull.numVertices(); i++)
			{
				if (hullP.distanceTo(hullPos[i]) < _tol)
				{
					//old map; seems adding the symmetry he to the map?
					//c_graphHalfEdge_dualCellFace[he.getId()].push_back(make_pair(_graphVertex.getId(), i));

					//set the values to the two maps
					//string hash1 = to_string(_graphVertex.getId()) + "," + to_string(i);
					zIntPair hash1(_graphVertex.getId(), i);
					forceCellFace_formHalfedge[hash1] = he.getId();

					//string hash2 = to_string(he.getId());
					int hash2 = he.getId();
					formHalfedge_forceCellFace[hash2] = zIntPair(_graphVertex.getId(), i);

				}
			}
		}



		////scale for visulisation
		zPoint center = fnMesh.getCenter();
		zPointArray vertices;
		fnDual.getVertexPositions(vertices);
		for (auto& v : vertices)
		{
			zVector vec = v - center;
			vec.normalize();
			vec *= 1;
			v += vec;
		}
		fnDual.setVertexPositions(vertices);

		//move it to the graph node
		zVector move = _graphVertex.getPosition() - center;
		fnDual.setTranslation(move);
	}

	//---- PRIVATE COMPUTE / UPDATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::computeTargets(float formWeight, float areaScale)
	{
		zFnGraph fnGraph(o_formGraph);
		
		form_targetHalfEdges.clear();
		form_targetHalfEdges.assign(fnGraph.numHalfEdges(), zVector());


		// Normals
		force_targetNormals.clear();
		force_targetNormals.assign(o_forceMeshes.size(), zVectorArray());

		// Areas
		force_targetAreas.clear();
		force_targetAreas.assign(o_forceMeshes.size(), zDoubleArray());

		for (int i = 0; i < o_forceMeshes.size(); i++)
		{
			zFnMesh tempFn(o_forceMeshes[i]);
			force_targetNormals[i].assign(tempFn.numPolygons(), zVector());

			force_targetAreas[i].assign(tempFn.numPolygons(), 1.0);
		}

		// Compute
		for (zItGraphHalfEdge he_form(o_formGraph); !he_form.end(); he_form++)
		{
			int he_form_id = he_form.getId();
			int he_form_sym_id = he_form.getSym().getId();

			zVector he_form_Vec = he_form.getVector();
			he_form_Vec.normalize();
			
			zVector f1_force_Vec;
			zVector f2_force_Vec;

			float f1_force_Area = 0;
			float f2_force_Area = 0;

			std::unordered_map<int, zIntPair>::const_iterator gotForce_1 = formHalfedge_forceCellFace.find(he_form_id);
			if (gotForce_1 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				f1_force_Vec = f1_force.getNormal();
				f1_force_Vec.normalize();

				f1_force_Area = f1_force.getPlanarFaceArea();
			}

			std::unordered_map<int, zIntPair>::const_iterator gotForce_2 = formHalfedge_forceCellFace.find(he_form_sym_id);
			if (gotForce_2 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f2_force(o_forceMeshes[(gotForce_2->second).first], (gotForce_2->second).second);
				f2_force_Vec = f2_force.getNormal();
				f2_force_Vec.normalize();

				f2_force_Vec *= -1;

				f2_force_Area = f2_force.getPlanarFaceArea();
			}
			
			zVector median_force_Vec = f1_force_Vec * 0.5 + f2_force_Vec * 0.5;
			if (f1_force_Vec.length() == 0) median_force_Vec = f2_force_Vec;
			if (f2_force_Vec.length() == 0) median_force_Vec = f1_force_Vec;
			median_force_Vec.normalize();

			float median_force_Area = (f1_force_Area + f2_force_Area) * 0.5 * areaScale;
			if (f1_force_Vec.length() == 0) median_force_Area = f2_force_Area * areaScale;
			if (f2_force_Vec.length() == 0) median_force_Area = f1_force_Area * areaScale;

			// target edge 
			zVector target_Vec = (he_form_Vec * formWeight) + (median_force_Vec * (1 - formWeight));
			target_Vec.normalize();


			form_targetHalfEdges[he_form_id] = target_Vec;

			if (gotForce_1 != formHalfedge_forceCellFace.end())
			{
				force_targetNormals[(gotForce_2->second).first][(gotForce_2->second).second] = target_Vec;
				force_targetAreas[(gotForce_2->second).first][(gotForce_2->second).second] = median_force_Area;
			}
		}

		/*for (int i = 0; i < o_forceMeshes.size(); i++)
		{
			zFnMesh tempFn(o_forceMeshes[i]);
			dualCellFace_NormTargets[i].assign(tempFn.numPolygons(), zVector());

			dualCellFace_AreaTargets[i].assign(tempFn.numPolygons(), 1.0);
		}

		for (int i = 0; i < c_graphHalfEdge_dualCellFace.size(); i++)
		{

			for (int j = 0; j < c_graphHalfEdge_dualCellFace[i].size(); j++)
			{
				int cellId = c_graphHalfEdge_dualCellFace[i][j].first;
				int faceId = c_graphHalfEdge_dualCellFace[i][j].second;

				zItGraphHalfEdge he(o_formGraph, i);

				dualCellFace_NormTargets[cellId][faceId] = he.getVector();
				dualCellFace_NormTargets[cellId][faceId].normalize();

				dualCellFace_AreaTargets[cellId][faceId] = he.getLength();
			}
		}*/

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::updateFormDiagram(float minmax_Edge, float dT, zIntergrationType type, int numIterations)
	{
		zFnGraph fnForm(o_formGraph);

		if (fnFormParticles.size() != fnForm.numVertices())
		{
			fnFormParticles.clear();
			o_formParticles.clear();


			for (int i = 0; i < fnForm.numVertices(); i++)
			{
				bool fixed = false;

				zObjParticle p;
				p.particle = zParticle(o_formGraph.graph.vertexPositions[i], fixed);
				o_formParticles.push_back(p);

			}

			for (int i = 0; i < o_formParticles.size(); i++)
			{
				fnFormParticles.push_back(zFnParticle(o_formParticles[i]));
			}
		}

		vector<zVector> v_residual;
		v_residual.assign(fnForm.numVertices(), zVector());

		vector<double> edgelengths;
		fnForm.getEdgeLengths(edgelengths);

		double minEdgeLength, maxEdgeLength;
		minEdgeLength = coreUtils.zMin(edgelengths);
		maxEdgeLength = coreUtils.zMax(edgelengths);

		minEdgeLength = maxEdgeLength * minmax_Edge;

		zVector* positions = fnForm.getRawVertexPositions();

		for (int k = 0; k < numIterations; k++)
		{
			for (zItGraphVertex v(o_formGraph); !v.end(); v++)
			{
				int i = v.getId();
				vector<zItGraphHalfEdge> cEdges;
				v.getConnectedHalfEdges(cEdges);

				zVector v_i = positions[v.getId()];

				// compute barycenter per vertex
				zVector b_i(0, 0, 0);
				for (auto& he : cEdges)
				{

					zVector v_j = positions[he.getVertex().getId()];

					zVector e_ij = v_i - v_j;
					double len_e_ij = e_ij.length();

					if (len_e_ij < minEdgeLength) len_e_ij = minEdgeLength;
					if (len_e_ij > maxEdgeLength) len_e_ij = maxEdgeLength;

					zVector t_ij = form_targetHalfEdges[he.getSym().getId()];
					t_ij.normalize();

					b_i += (v_j + (t_ij * len_e_ij));

				}

				b_i /= cEdges.size();

				// compute residue force				
				v_residual[i] = (b_i - v_i);

				zVector forceV = v_residual[i] * 1.0;
				fnFormParticles[i].addForce(forceV);
			}

			// update positions
			for (int i = 0; i < fnFormParticles.size(); i++)
			{
				fnFormParticles[i].integrateForces(dT, type);
				fnFormParticles[i].updateParticle(true);
			}
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::updateForceDiagram(float minmax_Edge, float dT, zIntergrationType type, int numIterations)
	{
		// create particles
		if (fnForceParticles.size() != o_forceMeshes.size())
		{
			fnForceParticles.clear();
			o_forceParticles.clear();

			o_forceParticles.assign(o_forceMeshes.size(), zObjParticleArray());
			fnForceParticles.assign(o_forceMeshes.size(), vector<zFnParticle>());

			int cellId = 0;
			for (auto& cell : o_forceMeshes)
			{
				zFnMesh tmpFn(cell);
				if (tmpFn.numVertices() == 0)
				{
					cellId++;
					continue;
				}

				zPoint* pos = tmpFn.getRawVertexPositions();



				for (zItMeshVertex v(cell); !v.end(); v++)
				{
					bool fixed = false;

					int i = v.getId();

					zObjParticle p;
					p.particle = zParticle(pos[i], fixed);
					o_forceParticles[cellId].push_back(p);

				}

				for (int i = 0; i < o_forceParticles[cellId].size(); i++)
				{
					fnForceParticles[cellId].push_back(zFnParticle(o_forceParticles[cellId][i]));
				}

				cellId++;
			}

		}

		for (int k = 0; k < numIterations; k++)
		{
			int cellId = 0;
			for (auto& cell : o_forceMeshes)
			{
				zFnMesh tmpFn(cell);
				if (tmpFn.numVertices() == 0)
				{
					cellId++;
					continue;
				}

				zPoint* positions = tmpFn.getRawVertexPositions();

				
				zPointArray fCenters;
				zFnMesh fnForceMesh(cell);

				fnForceMesh.getCenters(zFaceData,fCenters);

				for (zItMeshFace f(cell); !f.end(); f++)
				{
					int i = f.getId();

					vector<int> fVerts;
					f.getVertices(fVerts);

					float currentArea = f.getPlanarFaceArea();
					float scale = currentArea / force_targetAreas[cellId][i];

					for (int j = 0; j < fVerts.size(); j++)
					{
						zPoint p = positions[fVerts[j]];

						// Normal alignment force
						double dist = coreUtils.minDist_Point_Plane(p, fCenters[i], force_targetNormals[cellId][i]);
						zVector pForce = force_targetNormals[cellId][i] * dist * -1.0;

						fnForceParticles[cellId][fVerts[j]].addForce(pForce);

						// Area  force
						zVector dir = p - fCenters[i];
						float len = dir.length();
						dir.normalize();

						zPoint proj_p = fCenters[i] + dir * len * scale;
						zVector aForce = proj_p - p;

						fnForceParticles[cellId][fVerts[j]].addForce(aForce);
					}	
					

				}

				// update Particles
				for (int i = 0; i < fnForceParticles[cellId].size(); i++)
				{
					fnForceParticles[cellId][i].integrateForces(dT, type);
					fnForceParticles[cellId][i].updateParticle(true);
				}

				tmpFn.computeMeshNormals();

				cellId++;
			}


		}


	}

	ZSPACE_TOOLSETS_INLINE bool zTsGraphPolyhedra::checkDeviations(zDomainDouble& angleDeviation, double& angleTolerance, zDomainDouble& areaDeviation, double& areaTolerance, bool& printInfo)
	{
		bool out = true;
		vector<double> deviations;
		angleDeviation = zDomainDouble(10000, -10000);
		areaDeviation = zDomainDouble(10000, -10000);

		for (zItGraphHalfEdge he_form(o_formGraph); !he_form.end(); he_form++)
		{
			int he_form_id = he_form.getId();

			zVector he_form_Vec = he_form.getVector();
			he_form_Vec.normalize();

			zVector f1_force_Vec;
			float f1_force_Area = 0;

			std::unordered_map<int, zIntPair>::const_iterator gotForce_1 = formHalfedge_forceCellFace.find(he_form_id);
			if (gotForce_1 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				f1_force_Vec = f1_force.getNormal();
				f1_force_Vec.normalize();

				//angle
				double angle_i = f1_force_Vec.angle(he_form_Vec);

				if (angle_i > angleTolerance)
				{
					out = false;
				}

				if (angle_i < angleDeviation.min) angleDeviation.min = angle_i;
				if (angle_i > angleDeviation.max) angleDeviation.max = angle_i;

				//area
				float area_i = f1_force.getPlanarFaceArea();
				float target_area_i = force_targetAreas[(gotForce_1->second).first][(gotForce_1->second).second];

				if (abs(area_i - target_area_i) > areaTolerance)
				{
					out = false;
				}

				if (area_i < areaDeviation.min) areaDeviation.min = area_i;
				if (area_i > areaDeviation.max) areaDeviation.max = area_i;
			}
		}

		if (printInfo)
		{
			printf("\n Angle tolerance : %1.4f minDeviation : %1.4f , maxDeviation: %1.4f ", angleTolerance, angleDeviation.min, angleDeviation.max);
			printf("\n Angle tolerance : %1.4f minDeviation : %1.4f , maxDeviation: %1.4f \n", areaTolerance, areaDeviation.min, areaDeviation.max);

		}

		return out;
	}


	//---- PRIVATE UTILITY METHODS


	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::cleanConvexHull(zItGraphVertex& _vIt)
	{
		zPointArray _hullPts;

		// get connected edges
		zItGraphHalfEdgeArray cHalfEdges;
		_vIt.getConnectedHalfEdges(cHalfEdges);

		for (auto he : cHalfEdges)
		{
			//if (he.getVertex().checkValency(1))	_hullPts.push_back(he.getVertex().getPosition());
			//else _hullPts.push_back((_vIt.getPosition() + he.getVertex().getPosition()) / 2);

			zVector heVec = he.getVector();
			heVec.normalize();

			_hullPts.push_back((_vIt.getPosition() + heVec));
		}

		zFnMesh fnMesh(o_convexHullMeshes[_vIt.getId()]);
		fnMesh.makeConvexHull(_hullPts);

		while (fnMesh.numPolygons() > cHalfEdges.size())
		{
			zDoubleArray hedralAngles;
			fnMesh.getEdgeDihedralAngles(hedralAngles);

			vector<std::pair<double, int>> tmpPair;

			for (int i = 0; i < hedralAngles.size(); i++)
				tmpPair.push_back(make_pair(hedralAngles[i], i));

			std::sort(tmpPair.begin(), tmpPair.end());

			zItMeshEdge e(o_convexHullMeshes[_vIt.getId()], tmpPair[0].second);

			fnMesh.deleteEdge(e, true);
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::colorDualFaceConnectivity()
	{

		for (auto itr = formHalfedge_forceCellFace.begin();	itr != formHalfedge_forceCellFace.end(); itr++)
		{
			//zItGraphHalfEdge he(o_formGraph, atoi(itr->first.c_str()));
			zItGraphHalfEdge he(o_formGraph, itr->first);
			zItMeshFace f(o_forceMeshes[(itr->second).first], (itr->second).second);

			f.setColor(he.getEdge().getColor());			 
		}
		
	}


	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::colorGraphEdges()
	{
		zFnGraph fnGraph(o_formGraph);

		zItGraphVertex v_MaxValence;
		int maxValence = 0;;

		for (zItGraphVertex v(o_formGraph); !v.end(); v++)
		{
			if (v.getValence() > maxValence)
			{
				v_MaxValence = v;
				maxValence = v.getValence();
			}
		}				
		
		maxValence += 1;

		// breadth search first sorting
		zItGraphVertexArray bsfVertices;
		v_MaxValence.getBSF(bsfVertices);

		zIntArray e_colorIndex;
		e_colorIndex.assign(fnGraph.numEdges(), -1);		

		float h_step = 360 / maxValence;

		
		// Matchings method
		//based on https://www.youtube.com/watch?v=UN36EXS4w2I

		for (int i = 0; i < maxValence; i++)
		{
			zBoolArray e_colorMatching_visited;
			e_colorMatching_visited.assign(fnGraph.numEdges(), false);

			for (auto &bsfV : bsfVertices)
			{
				zItGraphHalfEdgeArray cHEdges;
				bsfV.getConnectedHalfEdges(cHEdges);

				int currentLocalID = -1;
				for (int j = 0; j < cHEdges.size(); j++)
				{
					if (e_colorIndex[cHEdges[j].getEdge().getId()] == -1 && !e_colorMatching_visited[cHEdges[j].getEdge().getId()])
					{
						currentLocalID = j;
						e_colorIndex[cHEdges[j].getEdge().getId()] = i;

						zColor col(i * h_step, 1, 1);
						cHEdges[j].getEdge().setColor(col);

						break;
					}
				}

				if (currentLocalID != -1)
				{
					for (int j = 0; j < cHEdges.size(); j++)
					{
						e_colorMatching_visited[cHEdges[j].getEdge().getId()] = true;
					}

					// get connect edges of end vertex and mark as visited
					zItGraphVertex endV = cHEdges[currentLocalID].getVertex();

					zItGraphHalfEdgeArray cHEdges_1;
					endV.getConnectedHalfEdges(cHEdges_1);

					for (int j = 0; j < cHEdges_1.size(); j++)
					{
						e_colorMatching_visited[cHEdges_1[j].getEdge().getId()] = true;
					}
				}
			}
			
		}
		
		// multi color of vertex and applying to edge. 
		
		/*zIntArray v_colorIndex;
		v_colorIndex.assign(fnGraph.numVertices(), -1);*/

		//for (int k =0; k < bsfVertices.size(); k++)
		//{
		//	/*zItGraphVertexArray cVertices;
		//	bsfVertices[k].getConnectedVertices(cVertices);

		//	int v_colorID = -1;

		//	for (int tmpID = maxValence - 1; tmpID >= 0; tmpID--)
		//	{
		//		bool chkRepeat = false;

		//		for (int j = 0; j < cVertices.size(); j++)
		//		{

		//			if (tmpID == v_colorIndex[cVertices[j].getId()])
		//			{
		//				chkRepeat = true;
		//				break;
		//			}
		//		}

		//		if (!chkRepeat)
		//		{
		//			v_colorID = tmpID;
		//			break;
		//		}
		//	}

		//	if (v_colorID != -1)
		//	{
		//		v_colorIndex[bsfVertices[k].getId()] = v_colorID;

		//		zColor col(v_colorID * h_step, 1, 1);
		//		bsfVertices[k].setColor(col);
		//	}*/			
		//	
		//}

		
		//fnGraph.computeEdgeColorfromVertexColor();

	}


}