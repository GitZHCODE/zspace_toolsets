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

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setFormGraphFromFile(string _path, zFileTpye _type, bool _staticGeom)
	{
		zColor red(1, 0, 0, 1);
		zColor green(0, 1, 0, 1);
		zColor blue(0, 0, 1, 1);
		zColor black(0, 0, 0, 1);

		// create graph
		zFnGraph fnGraph(o_formGraph);
		fnGraph.clear();
		fnGraph.from(_path, _type, _staticGeom);

		// allocate memory
		o_convexHullMeshes.clear();
		o_forceMeshes.clear();

		o_convexHullMeshes.assign(fnGraph.numVertices(), zObjMesh());
		o_forceMeshes.assign(fnGraph.numVertices(), zObjMesh());

		// color vertices
		
		zFloatArray formVertexWeights;
		zColorArray formVertexColor;

		/*for (zItGraphVertex v(o_formGraph); !v.end(); v++)
			if (v.getValence() > 1) v.setColor(zColor(0.0, 1.0, 0.0, 1.0));*/

		for (zItGraphVertex v(o_formGraph); !v.end(); v++)
		{
			if (v.getValence() == 1) v.setColor(zColor(1.0, 0.0, 0.0, 1.0));

			if (v.getColor() == blue)
			{
				formVertexWeights.push_back(1.0);
				formVertexColor.push_back(blue);
			}
			else
			{
				formVertexWeights.push_back(0.0);
				formVertexColor.push_back(red);
			}
		}

		fnGraph.setVertexColors(formVertexColor);
		setVertexWeights(zFormDiagram, formVertexWeights);

		// color edges
		//setDiagrams_uniqueColor();
		setDiagrams_forceColor(zDomainFloat(1,1));
	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setFormGraphFromOffsetMeshes(string _path_top, string _path_bottom, zFileTpye _type)
	{
		zColor red(1, 0, 0, 1);
		zColor green(0, 1, 0, 1);
		zColor blue(0, 0, 1, 1);
		zColor black(0, 0, 0, 1);

		zObjMesh o_topMesh;
		zFnMesh fnMesh_top(o_topMesh);
		fnMesh_top.from(_path_top, _type);

		zObjMesh o_bottomMesh;
		zFnMesh fnMesh_bottom(o_bottomMesh);
		fnMesh_bottom.from(_path_bottom, _type);

		zPointArray positions;
		zIntArray eConnects;

		int numVertices = 0;
		int numTopVertices = fnMesh_top.numVertices();
		zIntArray meshToGraphVertexIndexMap;
		meshToGraphVertexIndexMap.assign(fnMesh_top.numVertices() * 2, -1);

		zFloatArray formVertexWeights;
		zColorArray formVertexColor;

		for (zItMeshVertex v(o_topMesh); !v.end(); v++)
		{
			if (!v.checkValency(2))
			{
				positions.push_back(v.getPosition());
				meshToGraphVertexIndexMap[v.getId()] = numVertices;
				numVertices++;
				
				if (v.getColor() == blue)
				{
					formVertexWeights.push_back(0.0);
					formVertexColor.push_back(blue);
				}
				else
				{
					formVertexWeights.push_back(1.0);
					formVertexColor.push_back(red);
				}
			}
		}

		
		for (zItMeshVertex v(o_bottomMesh); !v.end(); v++)
		{		
			if (!v.checkValency(2))
			{
				positions.push_back(v.getPosition());
				meshToGraphVertexIndexMap[numTopVertices + v.getId()] = numVertices;
				numVertices++;

				if (v.getColor() == blue)
				{
					formVertexWeights.push_back(0.0);
					formVertexColor.push_back(blue);
				}
				else
				{
					formVertexWeights.push_back(1.0);
					formVertexColor.push_back(red);
				}

				if (!v.onBoundary())
				{
					int vID_0 = meshToGraphVertexIndexMap[v.getId()];
					int vID_1 = meshToGraphVertexIndexMap[v.getId() + numTopVertices];

					eConnects.push_back(vID_0);
					eConnects.push_back(vID_1);

				}
				
			}
			
		}

		for (zItMeshHalfEdge he(o_topMesh); !he.end(); he++)
		{
			if (!he.getEdge().onBoundary())
			{
				int vID_0 = meshToGraphVertexIndexMap[he.getStartVertex().getId()];
				int vID_1 = meshToGraphVertexIndexMap[he.getVertex().getId()];

				eConnects.push_back(vID_0);
				eConnects.push_back(vID_1);
			}

			he++;
		}

		for (zItMeshHalfEdge he(o_bottomMesh); !he.end(); he++)
		{
			if (!he.getEdge().onBoundary())
			{
				int vID_0 = meshToGraphVertexIndexMap[he.getStartVertex().getId() + numTopVertices];
				int vID_1 = meshToGraphVertexIndexMap[he.getVertex().getId() + numTopVertices];

				eConnects.push_back(vID_0);
				eConnects.push_back(vID_1);				
			}

			he++;
		}

		zVector up(0, 0, 1);
		for (zItMeshVertex v(o_topMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				positions.push_back(v.getPosition() + up);				

				int vID_0 = meshToGraphVertexIndexMap[v.getId()];

				eConnects.push_back(vID_0);
				eConnects.push_back(numVertices);

				numVertices++;

				if (v.getColor() == blue)
				{
					formVertexWeights.push_back(0.0);
					formVertexColor.push_back(blue);
				}
				else
				{
					formVertexWeights.push_back(1.0);
					formVertexColor.push_back(red);
				}
			}
		}

		

		// create graph
		zFnGraph fnGraph(o_formGraph);
		fnGraph.clear();
		fnGraph.create(positions, eConnects);

		fnGraph.setVertexColors(formVertexColor);
		setVertexWeights(zFormDiagram, formVertexWeights);

		// allocate memory
		o_convexHullMeshes.clear();
		o_forceMeshes.clear();

		o_convexHullMeshes.assign(fnGraph.numVertices(), zObjMesh());
		o_forceMeshes.assign(fnGraph.numVertices(), zObjMesh());

		// color vertices
		for (zItGraphVertex v(o_formGraph); !v.end(); v++)
			if (v.getValence() > 1) v.setColor(zColor(0.0, 1.0, 0.0, 1.0));

		// color edges
		//setDiagrams_uniqueColor();
		setDiagrams_forceColor(zDomainFloat(1, 1));

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
		setDiagrams_uniqueColor();

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setDiagram_Colors(int type, zDomainFloat wtDomain)
	{	
		switch (type)
		{
		case 0:
			setDiagrams_forceColor(wtDomain);
			break;

		case 1:
			setDiagrams_uniqueColor();
			break;

		case 2:
			setDiagrams_deviationColor();
			break;

		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setVertexWeights(zDiagramType type, const vector<float>& vWeights)
	{
		if (type == zFormDiagram)
		{
			zFnGraph fnForm(o_formGraph);

			if (vWeights.size() == 0)
			{
				formVWeights.clear();
				formVWeights.assign(fnForm.numVertices(), 1.0);
			}
			else
			{
				if (vWeights.size() != fnForm.numVertices()) throw std::invalid_argument("size of loads contatiner is not equal to number of form vertices.");

				formVWeights = vWeights;
			}

		}

		else if (type == zForceDiagram)
		{
			forceVWeights.clear();
			forceVWeights.assign(o_forceMeshes.size(), zFloatArray());
			

			/*if (vWeights.size() == 0)
			{*/
				for (int i = 0; i < o_forceMeshes.size(); i++)
				{
					zFnMesh fnForce(o_forceMeshes[i]);
					forceVWeights[i].assign(fnForce.numVertices(), 1.0);
				}

			/*}
			else
			{
				for (int i = 0; i < o_forceMeshes.size(); i++)
				{
					zFnMesh fnForce(o_forceMeshes[i]);
					forceVWeights[i].assign(fnForce.numVertices(), 1.0);
				}
			}*/

		}
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

	ZSPACE_TOOLSETS_INLINE zDomainFloat* zTsGraphPolyhedra::getRawAngleDeviation()
	{
		return &deviations[0];
	}

	ZSPACE_TOOLSETS_INLINE zDomainFloat* zTsGraphPolyhedra::getRawAreaDeviation()
	{
		return &deviations[1];
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
			if (g_v.checkValency(1)) continue;
			if (g_v.checkValency(2)) continue;
			if (g_v.checkValency(3)) continue;

			
				nodeId.push_back(g_v.getId());

				cout << "\n NODE_ID: " << g_v.getId() << " | " << g_v.getValence();

				// make clean convex hull
				cleanConvexHull(g_v);

				// dual mesh from convex hull
				createForceMesh(g_v, tol);

				activeNodes++;
			
		}
		cout << "\n////////////////////////////////////////////////////////////////////\n";
		cout << "convexhull mesh size: " << o_convexHullMeshes.size() << endl;
		cout << "dual mesh size: " << o_forceMeshes.size() << endl;
		cout << "active nodes: " << activeNodes << endl;

		// color connected faces
		setDiagrams_forceColor(zDomainFloat(1, 1));



	}

	//---- UPDATE METHODS

	ZSPACE_TOOLSETS_INLINE bool zTsGraphPolyhedra::equilibrium(bool& compTargets, float formWeight, float areaScale,  double dT, zIntergrationType type, float minMax_formEdge, bool areaForce, int numIterations, double angleTolerance, double areaTolerance, bool printInfo)
	{
		if (compTargets)
		{
			//printf("\n compute targets ");
			computeTargets(formWeight, areaScale);

			if (formVWeights.size() == 0) setVertexWeights(zDiagramType::zFormDiagram);
			if (forceVWeights.size() == 0) setVertexWeights(zDiagramType::zForceDiagram);

			compTargets = !compTargets;
		}

		angleTol = angleTolerance;

		bool out = checkDeviations(deviations[0], angleTolerance, deviations[1], areaTolerance, printInfo);

		// update diagrams

		if (formWeight != 1.0)
		{
			updateFormDiagram(minMax_formEdge, dT, type, numIterations);
		}

		if (formWeight != 0.0)
		{
			updateForceDiagram(areaForce, dT, type, numIterations);
		}

		// check deviations	
		out = checkDeviations(deviations[0], angleTolerance, deviations[1], areaTolerance, printInfo);

		return out;

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::exportFiles(string path)
	{

		string folderName = path  + "/out";
		_mkdir(folderName.c_str());
		for (const auto& entry : std::filesystem::directory_iterator(folderName)) std::filesystem::remove_all(entry.path());

		zFnGraph fnGraph(o_formGraph);
		fnGraph.to(folderName + "/formGraph.json", zJSON);		

		int i = 0;
		for (auto& m : o_forceMeshes)
		{
			zFnMesh fnMesh(m);
			fnMesh.to(folderName + "/forceMesh_" + to_string(i) + ".json", zJSON);

			i++;
		}

		exportHEMapToTXT(folderName);

		printf("\n files exported! ");
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
			
			// weight has to be 1 for force edges ie valence 1 edges
			float wt = (he_form.getVertex().checkValency(1) || he_form.getSym().getVertex().checkValency(1)) ? 1.0 : formWeight;
			zVector target_Vec = (he_form_Vec * wt) + (median_force_Vec * (1.0 - wt));
			target_Vec.normalize();


			form_targetHalfEdges[he_form_id] = target_Vec;

			if (gotForce_1 != formHalfedge_forceCellFace.end())
			{
				force_targetNormals[(gotForce_1->second).first][(gotForce_1->second).second] = target_Vec;
				force_targetAreas[(gotForce_1->second).first][(gotForce_1->second).second] = median_force_Area;
			}

			
			//printf("\n %i | %1.2f %1.2f %1.2f ", he_form_id, target_Vec.x, target_Vec.y, target_Vec.z);

		}

		

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

				if(b_i.length2() > 0) b_i /= cEdges.size();

				// compute residue force				
				v_residual[i] = (b_i - v_i);

				zVector forceV = v_residual[i] * formVWeights[i];
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

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::updateForceDiagram(bool areaForce, float dT, zIntergrationType type, int numIterations)
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


		zFnGraph fnFormGraph(o_formGraph);
		zPoint* graphPositions = fnFormGraph.getRawVertexPositions();

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
					float scale = force_targetAreas[cellId][i] / currentArea;

					

					for (int j = 0; j < fVerts.size(); j++)
					{
						zPoint p = positions[fVerts[j]];

						// Normal alignment force
						double dist = coreUtils.minDist_Point_Plane(p, fCenters[i], force_targetNormals[cellId][i]);
						zVector pForce = force_targetNormals[cellId][i] * dist * -1.0;
						
						pForce *= forceVWeights[cellId][fVerts[j]];
						fnForceParticles[cellId][fVerts[j]].addForce(pForce);

						// Area  force
						if (areaForce)
						{
							zVector dir = p - fCenters[i];
							float len = dir.length();
							dir.normalize();

							zPoint proj_p = fCenters[i] + dir * len * scale;
							zVector aForce = proj_p - p;

							aForce *= forceVWeights[cellId][fVerts[j]];
							fnForceParticles[cellId][fVerts[j]].addForce(aForce);
						}

						

					}	
					

				}

				
				// Centroid force
				zPoint cell_center;
				for (zItMeshVertex v(cell); !v.end(); v++)
				{
					cell_center += positions[v.getId()];
				}
				cell_center /= fnForceMesh.numVertices();

				// Scale Force

				float minLen = 0.25;
				float maxLen = 2.0;
				for (zItMeshVertex v(cell); !v.end(); v++)
				{
					zVector dir = positions[v.getId()] - cell_center;
					float len = dir.length();
					dir.normalize();

					if (len < minLen )
					{
						zPoint proj_p = cell_center + dir * minLen;
						zVector sForce = proj_p - positions[v.getId()];

						sForce *= forceVWeights[cellId][v.getId()];
						fnForceParticles[cellId][v.getId()].addForce(sForce);
					}

					else if (len > maxLen)
					{
						zPoint proj_p = cell_center + dir * maxLen;
						zVector sForce = proj_p - positions[v.getId()];

						sForce *= forceVWeights[cellId][v.getId()];
						fnForceParticles[cellId][v.getId()].addForce(sForce);
					}
					

				}
				
				

				zVector centroidForce = graphPositions[cellId] - cell_center;

				for (zItMeshVertex v(cell); !v.end(); v++)
				{

					zVector cForce = centroidForce * forceVWeights[cellId][v.getId()];
					fnForceParticles[cellId][v.getId()].addForce(cForce);
				}

				// min edge length
				float minEdgeLen = 0.2;
				for (zItMeshEdge e(cell); !e.end(); e++)
				{
					int i = e.getId();

					if (e.getLength() < minEdgeLen)
					{

						zItMeshHalfEdge he0 = e.getHalfEdge(0);
						zItMeshHalfEdge he1 = e.getHalfEdge(1);

						zVector  he0_vec = he0.getVector();
						zVector  he1_vec = he1.getVector();
						he0_vec.normalize();
						he1_vec.normalize();

						zVector pForce1 = (he1.getStartVertex().getPosition() + he1_vec * minEdgeLen) - he1.getVertex().getPosition();
						
						pForce1 *= forceVWeights[cellId][he1.getVertex().getId()];
						fnForceParticles[cellId][he1.getVertex().getId()].addForce(pForce1);


						zVector pForce0 = (he0.getStartVertex().getPosition() + he0_vec * minEdgeLen) - he0.getVertex().getPosition();
						
						pForce0 *= forceVWeights[cellId][he0.getVertex().getId()];
						fnForceParticles[cellId][he0.getVertex().getId()].addForce(pForce0);


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

	ZSPACE_TOOLSETS_INLINE bool zTsGraphPolyhedra::checkDeviations(zDomainFloat& angleDeviation, double& angleTolerance, zDomainFloat& areaDeviation, double& areaTolerance, bool& printInfo)
	{
		bool out = true;
		vector<double> deviations;
		angleDeviation = zDomainFloat(10000, -10000);
		areaDeviation = zDomainFloat(10000, -10000);

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

				//if (angle_i > 90)f1_force.setColor(zColor(1, 0, 0, 1));
				//else f1_force.setColor(zColor(0.5, 0.5, 0.5, 1));

				////area
				//float area_i = f1_force.getPlanarFaceArea();
				//float target_area_i = force_targetAreas[(gotForce_1->second).first][(gotForce_1->second).second];

				//if (abs(area_i - target_area_i) > areaTolerance)
				//{
				//	out = false;
				//}

				//if (area_i < areaDeviation.min) areaDeviation.min = area_i;
				//if (area_i > areaDeviation.max) areaDeviation.max = area_i;
			}

			zVector f2_force_Vec;
			float f2_force_Area = 0;

			float diffArea = 0;
			std::unordered_map<int, zIntPair>::const_iterator gotForce_2 = formHalfedge_forceCellFace.find(he_form.getSym().getId());
			if (gotForce_1 != formHalfedge_forceCellFace.end() && gotForce_2 != formHalfedge_forceCellFace.end())
			{
				//area
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				float area_1 = f1_force.getPlanarFaceArea();

				zItMeshFace f2_force(o_forceMeshes[(gotForce_2->second).first], (gotForce_2->second).second);
				float area_2 = f2_force.getPlanarFaceArea();

				float target_area = force_targetAreas[(gotForce_1->second).first][(gotForce_1->second).second];

				diffArea = abs(area_1 - area_2);	

				//printf("\n %i %i %i %i | %1.4f %1.4f | %1.4f ", (gotForce_1->second).first, (gotForce_1->second).second, (gotForce_2->second).first, (gotForce_2->second).second, area_1, area_2, target_area);
			}

			if (diffArea < areaDeviation.min) areaDeviation.min = diffArea;
			if (diffArea > areaDeviation.max) areaDeviation.max = diffArea;
		}

		if (printInfo)
		{
			printf("\n Angle tolerance : %1.4f minDeviation : %1.4f , maxDeviation: %1.4f ", angleTolerance, angleDeviation.min, angleDeviation.max);
			printf("\n Area tolerance : %1.4f minDeviation : %1.4f , maxDeviation: %1.4f \n", areaTolerance, areaDeviation.min, areaDeviation.max);

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

			cout <<endl << _hullPts[_hullPts.size() - 1];
		}

		zFnMesh fnMesh(o_convexHullMeshes[_vIt.getId()]);
		fnMesh.makeConvexHull(_hullPts);

		/*while (fnMesh.numPolygons() > cHalfEdges.size())
				{
					zDoubleArray hedralAngles;
					fnMesh.getEdgeDihedralAngles(hedralAngles);

					vector<std::pair<double, int>> tmpPair;

					for (int i = 0; i < hedralAngles.size(); i++)
						tmpPair.push_back(make_pair(hedralAngles[i], i));

					std::sort(tmpPair.begin(), tmpPair.end());

					zItMeshEdge e(o_convexHullMeshes[_vIt.getId()], tmpPair[0].second);

					fnMesh.deleteEdge(e, true);
				}*/

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setDiagrams_uniqueColor()
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

		
		// Matchings method based on https://www.youtube.com/watch?v=UN36EXS4w2I

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
						cHEdges[j].getEdge().setWeight(1);

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
		
		for (auto itr = formHalfedge_forceCellFace.begin(); itr != formHalfedge_forceCellFace.end(); itr++)
		{
			zItGraphHalfEdge he(o_formGraph, itr->first);
			zItMeshFace f(o_forceMeshes[(itr->second).first], (itr->second).second);

			f.setColor(he.getEdge().getColor());
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setDiagrams_forceColor(zDomainFloat weightDomain)
	{

		//compute edgeWeights

		float minArea = 100000;
		float maxArea = -100000;

		int cellId = 0;
		for (auto& cell : o_forceMeshes)
		{
			zFnMesh tmpFn(cell);
			if (tmpFn.numVertices() == 0)
			{
				cellId++;
				continue;
			}

			zDoubleArray fAreas;
			tmpFn.getPlanarFaceAreas(fAreas);

			for (auto &fArea: fAreas)
			{		
				if (fArea < minArea) minArea = fArea;
				if (fArea > maxArea) maxArea = fArea;
			}

			cellId++;
		}
			

		zDomainFloat areaDomain(minArea, maxArea);
		zDomainColor colDomain(zColor(0.027, 0, 0.157, 1), zColor(0.784, 0, 0.157, 1));

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
			std::unordered_map<int, zIntPair>::const_iterator gotForce_2 = formHalfedge_forceCellFace.find(he_form_sym_id);

			if (gotForce_1 != formHalfedge_forceCellFace.end() && gotForce_2 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				f1_force_Area = f1_force.getPlanarFaceArea();

				zItMeshFace f2_force(o_forceMeshes[(gotForce_2->second).first], (gotForce_2->second).second);
				f2_force_Area = f2_force.getPlanarFaceArea();

				float fArea = (f1_force_Area + f2_force_Area) * 0.5;

				float wt = coreUtils.ofMap((float)fArea, areaDomain, weightDomain);
				zColor col = coreUtils.blendColor(fArea, areaDomain, colDomain, zRGB);

				he_form.getEdge().setWeight(wt);
				he_form.getEdge().setColor(col);

				f1_force.setColor(col);
				f2_force.setColor(col);

			}
			else
			{
				if (he_form.getVertex().checkValency(1) || he_form.getSym().getVertex().checkValency(1))
				{
					he_form.getEdge().setColor(zColor(0, 1, 0, 1));
					he_form.getEdge().setWeight(1);
				}
				else
				{
					he_form.getEdge().setColor(zColor(0, 0, 0, 1));
					he_form.getEdge().setWeight(1);
				}
								
				if (gotForce_1 != formHalfedge_forceCellFace.end())
				{
					zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
					f1_force.setColor(zColor(0, 1, 0, 1));
				}
			}
		
		}
			
	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::setDiagrams_deviationColor()
	{		
		zDomainColor colDomain(zColor(180, 1, 1), zColor(360, 1, 1));

		for (zItGraphHalfEdge he_form(o_formGraph); !he_form.end(); he_form++)
		{
			int he_form_id = he_form.getId();

			zVector he_form_Vec = he_form.getVector();
			he_form_Vec.normalize();

			zVector he_form_sym_Vec = he_form.getSym().getVector();
			he_form_sym_Vec.normalize();
					

			zVector f1_force_Vec;
			zVector f2_force_Vec;

			std::unordered_map<int, zIntPair>::const_iterator gotForce_1 = formHalfedge_forceCellFace.find(he_form_id);
			std::unordered_map<int, zIntPair>::const_iterator gotForce_2 = formHalfedge_forceCellFace.find(he_form.getSym().getId());

			
			if (gotForce_1 != formHalfedge_forceCellFace.end() && gotForce_2 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				f1_force_Vec = f1_force.getNormal();
				f1_force_Vec.normalize();

				zItMeshFace f2_force(o_forceMeshes[(gotForce_2->second).first], (gotForce_2->second).second);
				f2_force_Vec = f2_force.getNormal();
				f2_force_Vec.normalize();

				//angle
				float angle_1 = f1_force_Vec.angle(he_form_Vec);
				float angle_2 = f2_force_Vec.angle(he_form_sym_Vec);

				float angle_i = (angle_1 > angle_2) ? angle_1 : angle_2;

				zColor col = coreUtils.blendColor(angle_i, deviations[0], colDomain, zHSV);
								
				
				if(angle_i < angleTol)  col = (zColor());
				else if (angle_i > 90)  col = (zColor(1,0.65,0,1));
				else he_form.getEdge().setColor(col);
				
				he_form.getEdge().setColor(col);
				f1_force.setColor(col);
				f2_force.setColor(col);

				he_form.getEdge().setWeight(1);
			}
			else if (gotForce_1 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f1_force(o_forceMeshes[(gotForce_1->second).first], (gotForce_1->second).second);
				f1_force_Vec = f1_force.getNormal();
				f1_force_Vec.normalize();				

				//angle
				float angle_i = f1_force_Vec.angle(he_form_Vec);						

				zColor col = coreUtils.blendColor(angle_i, deviations[0], colDomain, zHSV);

				if (angle_i < angleTol)  col = (zColor());
				else if (angle_i > 90)  col = (zColor(1, 0.65, 0, 1));
				else he_form.getEdge().setColor(col);

				he_form.getEdge().setColor(col);
				f1_force.setColor(col);
				
				he_form.getEdge().setWeight(1);
			}

			else if (gotForce_2 != formHalfedge_forceCellFace.end())
			{
				zItMeshFace f2_force(o_forceMeshes[(gotForce_2->second).first], (gotForce_2->second).second);
				f2_force_Vec = f2_force.getNormal();
				f2_force_Vec.normalize();

				//angle
				float angle_i = f2_force_Vec.angle(he_form_sym_Vec);

				zColor col = coreUtils.blendColor(angle_i, deviations[0], colDomain, zHSV);

				if (angle_i < angleTol)  col = (zColor());
				else if (angle_i > 90)  col = (zColor(1, 0.65, 0, 1));
				else he_form.getEdge().setColor(col);

				he_form.getEdge().setColor(col);
				f2_force.setColor(col);

				he_form.getEdge().setWeight(1);
			}


			he_form++;
		}
	

		

	}

	ZSPACE_TOOLSETS_INLINE void zTsGraphPolyhedra::exportHEMapToTXT(string path)
	{
		string outfilename = path;
		outfilename += "/heMap.txt";
		
		ofstream myfile;
		myfile.open(outfilename.c_str());

		if (myfile.fail())
		{
			cout << " error in opening file  " << outfilename.c_str() << endl;
			return;
		}

		int counter = 0;

		for (auto itr = formHalfedge_forceCellFace.begin(); itr != formHalfedge_forceCellFace.end(); itr++)
		{
			myfile << "\n" << itr->first << "," << (itr->second).first << "," << (itr->second).second;
		}			

		myfile.close();

	}


}