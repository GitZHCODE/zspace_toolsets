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
// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
//


#include<headers/zToolsets/natpower/zTsNatpowerSDF.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsNatpowerSDF::zTsNatpowerSDF()
	{


		printHeightDomain = zDomainFloat(0.006, 0.012);

	}


	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsNatpowerSDF::~zTsNatpowerSDF() {}

	//---- CREATE METHODS

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createFieldMeshFromMeshBounds(float cellSize, float offset)
	{
		// Transform
		//zTransform t = sectionFrames[0];
		zTransform t = leftPlanes[0]; //<given that start left and right planes are the same

		zFnMesh fnMesh(o_GuideMesh);
		fnMesh.setTransform(base_world, true, false);
		fnMesh.setTransform(base_local, true, true);

		zPoint bbMin, bbMax;
		fnMesh.getBounds(bbMin, bbMax);

		//using the new bounds, calculate the res in x and y, update the pounds to have complete number of cell size
		float lenX = bbMax.x - bbMin.x + offset;
		float lenY = bbMax.y - bbMin.y + offset;
		int resX = ceil(lenX / cellSize);
		int resY = ceil(lenY / cellSize);

		lenX += (resX * cellSize) - lenX;
		lenY += (resY * cellSize) - lenY;

		lenX /= 2;
		lenY /= 2;

		zDomain<zPoint> bb (zPoint(-lenX, -lenX, 0), zPoint(lenX, lenX, 0));

		zFnMeshScalarField fnField(o_field);
		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);
		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);

		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);





		// transform back 
		fnMesh.setTransform(base_world, true, true);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createFieldMeshFromSectionBounds(float cellSize, float offset)
	{
		

		float minX = FLT_MAX;
		float minY = FLT_MAX;
		float maxX = FLT_MIN;
		float maxY = FLT_MIN;
		//iterate through all the section graphs bounds and pick the largest one
		for (int i = 0; i < o_sectionGraphs.size(); i++)
		{
			zFnGraph fnGraph(o_sectionGraphs[i]);
			zTransform t = sectionFrames[i];
			fnGraph.setTransform(t, true, false);
			// Transform
			zTransform tLocal;
			tLocal.setIdentity();
			fnGraph.setTransform(tLocal, true, true);

			zPoint pmin, pmax;
			fnGraph.getBounds(pmin, pmax);
			if (pmin.x < minX) minX = pmin.x;
			if (pmin.y < minY) minY = pmin.y;
			if (pmax.x > maxX) maxX = pmax.x;
			if (pmax.y > maxY) maxY = pmax.y;
			
			//cout << endl << "min" << pmin;
			//cout << endl << "max" << pmin << endl;


			fnGraph.setTransform(t, true, true);
		}


		zPoint bbMin(minX, minY, 0);
		zPoint bbMax(maxX, maxY, 0);
		
		/*cout << endl << "\n FINAL MIN MAX" << bbMin;
		cout << endl << "min" << bbMin;
		cout << endl << "max" << bbMax << endl;*/

		//fnMesh.getBounds(bbMin, bbMax);

		//using the new bounds, calculate the res in x and y, update the pounds to have complete number of cell size
		float lenX = bbMax.x - bbMin.x;
		float lenY = bbMax.y - bbMin.y;
		
		int resX = ceil((lenX + (offset*2)) / cellSize);
		int resY = ceil((lenY + (offset*2)) / cellSize);
		
	/*	cout << endl << "cell" << cellSize;
		cout << endl << "resX" << resX;
		cout << endl << "resY" << resY;*/

		lenX = (resX * cellSize) - lenX;
		lenY = (resY * cellSize) - lenY;

		lenX /= 2;
		lenY /= 2;


		bbMin.x -= offset - lenX;
		bbMin.y -= offset - lenX;
		bbMax.x += offset + lenX;
		bbMax.y += offset + lenX;


		zFnMeshScalarField fnField(o_field);
		//fnField.create()
		fnField.create(bbMin, bbMax, resX, resY, 1, true, false);
		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);

		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);


	}



	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createFieldMesh(zDomain<zPoint>& bb, int resX, int resY)
	{
		zFnMeshScalarField fnField(o_field);

		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);


		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);


		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);


	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createFieldMeshCellSize(zDomain<zPoint>& bb, float cellSizeX, float cellSizeY)
	{
		zFnMeshScalarField fnField(o_field);
		int resX = ceil((bb.max.x - bb.min.x) / cellSizeX);
		int resY = ceil((bb.max.y - bb.min.y) / cellSizeY);
		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);
		zDomainColor dCol(zBLUE, zRED);
		fnField.setFieldColorDomain(dCol);

		printf("\n dCol Min %1.2f , %1.2f , %1.2f", dCol.min.r, dCol.min.g, dCol.min.b);
		printf("dCol max %1.2f , %1.2f , %1.2f", dCol.max.r, dCol.max.g, dCol.max.b);
	}

	//--- SET METHODS 

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setFromJSON(string dir, int _blockID, bool runBothPlanes, bool runPlaneLeft)
	{
		string path = dir + "blockMesh_" + to_string(_blockID) + ".json";
		/*bool pathExist = coreUtils.fileExists(path);
		if (!pathExist)
		{
			throw std::invalid_argument(" error: invalid path. ");
			return;
		}*/

		json j;
		bool jsonCheck = core.json_read(path, j);

		if (jsonCheck) printf("\n %s exists", path.c_str());

		if (!jsonCheck)
		{
			printf("\n %s doesnt exists", path.c_str());
			throw std::invalid_argument(" error: invalid outPath. ");
			return;
		}

		blockId = _blockID;
		planarBlock = j["IsPlanar"];
		printf("\n is planar %s ", to_string(planarBlock));


		readJSON(path, _blockID, runBothPlanes, runPlaneLeft);


		

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setSliceMesh(zObjMesh& _o_SliceMesh, bool left)
	{
		(left) ? o_SliceMesh_Left = _o_SliceMesh : o_SliceMesh_Right = _o_SliceMesh;

		//(left) ? leftPlaneExists = true : rightPlaneExists = true;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setMedialGraph(zObjGraph& _o_MedialGraph)
	{
		o_MedialGraph = _o_MedialGraph;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane, bool left)
	{
		(left) ? leftPlanes[0] = _sPlane : rightPlanes[0] = _sPlane;
		(left) ? leftPlanes[1] = _ePlane : rightPlanes[1] = _ePlane;

		//(left) ? leftPlaneExists = true : rightPlaneExists = true;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setGradientTriMesh(zObjMesh& _o_gradientTriMesh)
	{

		zFnMesh fnGradMesh(o_gradientTriMesh);
		o_gradientTriMesh = _o_gradientTriMesh;

		//fnGradMesh.setFaceColor(zColor(), true);

		zPoint* vPositions = fnGradMesh.getRawVertexPositions();

		MatrixXd V(fnGradMesh.numVertices(), 3);
		MatrixXi FTris(fnGradMesh.numPolygons(), 4);

		fnGradMesh.getMatrices_quadmesh(V, FTris);

		// fill vertex matrix
		/*for (int i = 0; i < fnGradMesh.numVertices(); i++)
		{
			V(i, 0) = vPositions[i].x;
			V(i, 1) = vPositions[i].y;
			V(i, 2) = vPositions[i].z;
		}*/

		gradientTriMesh_V = V;

		// fill triangle matrix		

		/*int nTris = 0;
		for (zItMeshFace f(o_gradientTriMesh);!f.end(); f++)
		{
			int i = f.getId();

			zIntArray fVerts;
			f.getVertices(fVerts);


			FTris(i, 0) = fVerts[0];
			FTris(i, 1) = fVerts[1];
			FTris(i, 2) = fVerts[2];
		}*/

		/*zIntArray pCounts, pConnects;
		fnGradMesh.getPolygonData(pConnects, pCounts);

		int pConnectsIndex = 0;
		for (int i = 0; i < pCounts.size(); i++)
		{
			if (pCounts[i] == 3)
			{
				FTris(i, 0) = pConnects[pConnectsIndex + 0];
				FTris(i, 1) = pConnects[pConnectsIndex + 1];
				FTris(i, 2) = pConnects[pConnectsIndex + 2];
			}
			else
			{
				printf("\n non Tri mesh %i ", i);
			}

			pConnectsIndex += pCounts[i];
		}*/

		gradientTriMesh_FTris = FTris;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setOffsetDomain(zDomainFloat& _offsetDomain)
	{
		offsetDomain = _offsetDomain;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setTransforms(bool toLocal)
	{
		if (toLocal)
		{
			zFnGraph fnGraphMedial(o_MedialGraph);
			fnGraphMedial.setTransform(base_world, true, false);
			fnGraphMedial.setTransform(base_local, true, true);

			zFnMesh fnMesh(o_GuideMesh);
			fnMesh.setTransform(base_world, true, false);
			fnMesh.setTransform(base_local, true, true);

			zFnMesh fnMesh_left(o_SliceMesh_Left);
			fnMesh_left.setTransform(base_world, true, false);
			fnMesh_left.setTransform(base_local, true, true);

			zFnMesh fnMesh_right(o_SliceMesh_Right);
			fnMesh_right.setTransform(base_world, true, false);
			fnMesh_right.setTransform(base_local, true, true);

			// section graphs
			for (auto& g : o_sectionGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			// contour graphs
			for (auto& g : o_contourGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			// trim graphs
			for (auto& g : o_trimGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, false);
				fnGraph.setTransform(base_local, true, true);
			}

			zFnPointCloud fnCritical_min(criticalMinLayer_pts);
			fnCritical_min.setTransform(base_world, true, false);
			fnCritical_min.setTransform(base_local, true, true);

			zFnPointCloud fnCritical_max(criticalMaxLayer_pts);
			fnCritical_max.setTransform(base_world, true, false);
			fnCritical_max.setTransform(base_local, true, true);

			/*if (leftPlaneExists)
			{
				for (int i = 0; i < 2; i++)
				{
					zTransform transform = coreUtils.PlanetoPlane(leftPlanes[i], base_local);
					leftPlanes[i] *= transform;
				}
			}

			if (rightPlaneExists)
			{

					zTransformationMatrix tMat;
					tMat.setTransform(rightPlanes[0]);

					zTransformationMatrix tLocal;
					tLocal.setTransform(base_local);
					zTransform t = tMat.getToMatrix(tLocal).transpose();

					rightPlanes[0] *= t;
					rightPlanes[1] *= t;

			}*/

		}
		else
		{
			zFnGraph fnGraphMedial(o_MedialGraph);
			fnGraphMedial.setTransform(base_world, true, true);

			zFnMesh fnMesh(o_GuideMesh);
			fnMesh.setTransform(base_world, true, true);

			zFnMesh fnMesh_left(o_SliceMesh_Left);
			fnMesh_left.setTransform(base_world, true, true);

			zFnMesh fnMesh_right(o_SliceMesh_Right);
			fnMesh_right.setTransform(base_world, true, true);

			// section graphs
			for (auto& g : o_sectionGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			// contour graphs
			for (auto& g : o_contourGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			// trim graphs
			for (auto& g : o_trimGraphs)
			{
				zFnGraph fnGraph(g);
				fnGraph.setTransform(base_world, true, true);
			}

			zFnPointCloud fnCritical_min(criticalMinLayer_pts);
			fnCritical_min.setTransform(base_world, true, true);

			zFnPointCloud fnCritical_max(criticalMaxLayer_pts);
			fnCritical_max.setTransform(base_world, true, true);

		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setFrames(vector<zPlane>& _sectionFrames)
	{
		sectionFrames.clear();
		sectionFrames = _sectionFrames;

		o_sectionGraphs.assign(sectionFrames.size(), zObjGraph());
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setCableGraph(string folderDir)
	{
		//read all files in directory
		zStringArray files;
		core.getFilesFromDirectory(files, folderDir, zJSON);
		printf("\n readCableGraph %i", files.size());

		o_CableGraphs.assign(files.size(), zObjGraph());
		//read all files
		for (int i = 0; i < o_CableGraphs.size(); i++)
		{
			zFnGraph fnGraph(o_CableGraphs[i]);
			fnGraph.from(files[i], zJSON);
		}

		printf("\n readCableGraph %i", o_CableGraphs.size());

	}


	//---- GET METHODS

	ZSPACE_TOOLSETS_INLINE zTransform* zTsNatpowerSDF::getRawBlockStartEnd(bool left)
	{
		/*if (left && leftPlaneExists) return &leftPlanes[0];
		else if (!left && rightPlaneExists) return &rightPlanes[0];
		else return nullptr;*/
		if (left) return &leftPlanes[0];
		else if (!left) return &rightPlanes[0];
		else return nullptr;
	}

	ZSPACE_TOOLSETS_INLINE vector<zTransform> zTsNatpowerSDF::getBlockFrames()
	{
		return sectionFrames;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsNatpowerSDF::getBlockSectionGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_sectionGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_sectionGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsNatpowerSDF::getBlockRaftGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_raftGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_raftGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsNatpowerSDF::getBlockContourGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_contourGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_contourGraphs)
		{
			out.push_back(&graph);
		}
		//printf("\n C++ num of graphs %i", o_contourGraphs.size());

		return out;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getBlockContourGraphsSequence(int& numGraphs, int* vSequence)
	{

		float minLayerHeight = 10;
		float maxLayerHeight = 0;

		int r0 = 0;
		int r1 = floor(o_contourGraphs.size() * 0.5) - 1;
		int r2 = floor(o_contourGraphs.size() * 0.5);
		int r3 = (o_contourGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		//bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;
		bool deckBlock = true;

		//int end = (leftPlaneExists && rightPlaneExists) ? floor(o_contourGraphs.size() * 0.5) : o_contourGraphs.size();
		int end = floor(o_contourGraphs.size() * 0.5);

		float printLength = 0;

		zInt2DArray vSeq; //sequence of points in each layer
		zIntArray vSeqAll;
		int totalVerCount = 0;
		// PRINT LAYERS
		for (int j = 0; j < end; j++)
		{
			//int kStart = (!deckBlock && !rightPlaneExists) ? 1 : 0;
			//int kEnd = (!deckBlock && !leftPlaneExists) ? 1 : 2;

			int kStart = (!deckBlock) ? 1 : 0;
			int kEnd = (!deckBlock) ? 1 : 2;

			for (int k = kStart; k < kEnd; k++)
			{
				int i = (k == 0) ? j : j + end;

				if (!deckBlock) i = j;

				if (i == r0) continue;

				if (deckBlock && i == r2) continue;

				zFnGraph fnContourGraph(o_contourGraphs[i]);
				if (fnContourGraph.numVertices() == 0) continue;

				// vertex sequence
				zIntArray sequence;
				zItGraphVertexArray vArray;

				for (zItGraphVertex v(o_contourGraphs[i]); !v.end(); v++)
				{
					if (!v.checkValency(2))
					{
						vArray.push_back(v);
					}
				}


				if (vArray.size() == 2)
				{
					zItGraphHalfEdge he = vArray[0].getHalfEdge();
					sequence.push_back(vArray[0].getId());

					do
					{

						sequence.push_back(he.getVertex().getId());

						he = he.getNext();
					} while (he.getVertex() != vArray[1]);

					sequence.push_back(vArray[1].getId());
					sequence.push_back(vArray[0].getId());

				}

				if (vArray.size() == 0)
				{
					zItGraphHalfEdge he(o_contourGraphs[i], 0);

					zItGraphVertex startV = he.getStartVertex();
					zItGraphHalfEdge startHe = he;
					sequence.push_back(he.getStartVertex().getId());

					do
					{
						if (he.getVertex() == startV)
						{
							sequence.push_back(he.getVertex().getId());
						}
						else
						{
							sequence.push_back(he.getVertex().getId());
						}

						he = he.getNext();
					} while (he != startHe);

				}

				totalVerCount += sequence.size();
				vSeq.push_back(sequence);
			}
		}
		vSequence = new int[totalVerCount];
		for (int i = 0; i < vSeq.size(); i++)
		{
			for (int j = 0; j < vSeq[i].size(); j++)
			{
				vSequence[i] = vSeq[i][j];
			}

		}

	}

	ZSPACE_TOOLSETS_INLINE zIntArray zTsNatpowerSDF::getGraphSequence(zObjGraph graph)
	{
		zIntArray sequence;
		zItGraphVertexArray vArray;

		zFnGraph fnGraph(graph);
		if (fnGraph.numVertices() == 0) return sequence;

		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (!v.checkValency(2))
			{
				vArray.push_back(v);
			}
		}

		printf("\n vArray.size() - %i", vArray.size());

		if (vArray.size() == 2)
		{
			zItGraphHalfEdge he = vArray[0].getHalfEdge();
			sequence.push_back(vArray[0].getId());
			do
			{
				sequence.push_back(he.getVertex().getId());
				he = he.getNext();
			} while (he.getVertex() != vArray[1]);

			sequence.push_back(vArray[1].getId());
			sequence.push_back(vArray[0].getId());
		}
		if (vArray.size() == 0)
		{

			zItGraphHalfEdge he(graph, 0);

			zItGraphVertex startV = he.getStartVertex();

			zItGraphHalfEdge startHe = he;

			sequence.push_back(he.getStartVertex().getId());

			do
			{

				if (he.getVertex() == startV)
				{

					sequence.push_back(he.getVertex().getId());
				}
				else
				{

					sequence.push_back(he.getVertex().getId());
				}

				he = he.getNext();

			} while (he != startHe);
		}

		printf("\n num of vertices : num of sequence ---- %i : %i \n ", fnGraph.numVertices(), sequence.size());

		return sequence;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraphPointerArray zTsNatpowerSDF::getBlockTrimGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_trimGraphs.size();
		//numGraphs = o_trimGraphs_bracing.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_trimGraphs)
		//for (auto& graph : o_trimGraphs_bracing)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_TOOLSETS_INLINE zObjPointCloud* zTsNatpowerSDF::getRawCriticalPoints(bool minHeight)
	{
		return (minHeight) ? &criticalMinLayer_pts : &criticalMaxLayer_pts;
	}

	ZSPACE_TOOLSETS_INLINE zObjMeshScalarField* zTsNatpowerSDF::getRawFieldMesh()
	{
		return &o_field;
	}

	ZSPACE_TOOLSETS_INLINE zObjGraph* zTsNatpowerSDF::getRawMedialGraph()
	{
		return &o_MedialGraph;
	}

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsNatpowerSDF::getRawLeftMesh()
	{
		return &o_SliceMesh_Left;
	}

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsNatpowerSDF::getRawRightMesh()
	{
		return &o_SliceMesh_Right;
	}

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsNatpowerSDF::getRawGuideMesh()
	{
		return &o_GuideMesh;
	}

	ZSPACE_INLINE zObjMesh* zTsNatpowerSDF::getRawGradientMesh()
	{
		return &o_gradientTriMesh;
	}

	ZSPACE_INLINE zObjMeshScalarField* zTsNatpowerSDF::getRawMeshScalarField()
	{
		return &o_field;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::isPlanarBlock()
	{
		return planarBlock;
	}

	//---- COMPUTE METHODS 

	//---------SLICE MESH----------------
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_SliceMesh_Pentagon(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides, bool left)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		//get end he
		zItMeshHalfEdge heEnd = heStart;
		//int sideDivisionCount = 1;
		while (heEnd.getVertex().getId() != endVID)
		{
			heEnd = heEnd.getNext().getSym().getNext();
			//sideDivisionCount++;
		}

		//printf("\n sideDivisionCount %i ", sideDivisionCount);


		// walk along spine
		bool exit = false;
		zItMeshHalfEdge he = heStart;

		int featureWalkNum1 = (FeaturedNumStrides.size() - 1) / 2;
		int featureStart = 0;


		do
		{
			//right
			if (!left)
			{
				//for (int ff = FeaturedNumStrides.size()-1; ff >= 0; ff--)
				//{
					//int blockStride = FeaturedNumStrides[ff];
				zItMeshHalfEdge he_Left = he.getPrev();

				//for (int& blockStride : FeaturedNumStrides)
				for (int k = featureStart; k < featureWalkNum1; k++)
				{
					int blockStride = FeaturedNumStrides[k];
					for (int i = 0; i < blockStride; i++)
					{

						zItMeshHalfEdge heTmp = he_Left;
						int pCount = 0;
						do
						{
							zPoint p = heTmp.getVertex().getPosition();
							int vID;
							bool chkRepeat = core.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								core.addToPositionMap(positionVertex, p, vID, PRECISION);
								positions.push_back(p);
							}

							pConnects.push_back(vID);
							pCount++;

							heTmp = heTmp.getNext();

							//heTmp = heTmp.getPrev();

						} while (heTmp != he_Left);

						pCounts.push_back(pCount);

						he_Left = he_Left.getPrev().getSym().getPrev();


					}

				}
			}

			//sideCounter++;
			//if (sideCounter < sideDivisionCount)
			//{
			//	heTmp.getEdge().setColor(zCYAN);
			//}

			//left
			else
			{
				zItMeshHalfEdge he_Right = he.getSym().getNext();
				//for (int& blockStride : FeaturedNumStrides)
				for (int k = featureStart; k < featureWalkNum1; k++)
				{
					int blockStride = FeaturedNumStrides[k];

					for (int i = 0; i < blockStride; i++)
					{
						zItMeshHalfEdge heTmp = he_Right;
						int pCount = 0;
						do
						{
							zPoint p = heTmp.getStartVertex().getPosition();
							int vID;
							bool chkRepeat = core.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								core.addToPositionMap(positionVertex, p, vID, PRECISION);
								positions.push_back(p);
							}
							pConnects.push_back(vID);
							pCount++;

							heTmp = heTmp.getNext();

						} while (heTmp != he_Right);

						pCounts.push_back(pCount);

						he_Right = he_Right.getNext().getSym().getNext();

					}

				}

			}


			//spine walk
			he = he.getNext().getSym().getNext();

			if (he == heStart) exit = true;
			//if(he.getStartVertex().getId() == )



		} while (!exit);

		//printf("\n sliceMesh input %i %i %i ", positions.size(), pCounts.size(), pConnects.size());


		// cap faces
		//right
		if (!left)
		{
			zItMeshHalfEdge heTop = heStart;
			zItMeshHalfEdge heBottom = heStart;
			heBottom = heBottom.getPrev().getSym().getPrev();
			heBottom = heBottom.getPrev().getSym().getPrev();

			zItMeshHalfEdge heTop_corner = heTop;
			zItMeshHalfEdge heBottom_corner = heBottom;

			for (int k = featureStart; k < featureWalkNum1; k++)
			{
				int blockStride = FeaturedNumStrides[k];
				//for (int& blockStride : FeaturedNumStrides)
				//{

				for (int i = 0; i < blockStride; i++)
				{
					heTop_corner = heTop_corner.getPrev().getPrev().getSym();
					heBottom_corner = heBottom_corner.getPrev().getPrev().getSym();
				}
			}


			exit = false;
			do
			{
				if (heTop.getVertex().getId() == endVID) exit = true;

				zPoint p0 = heTop.getVertex().getPosition();
				int v0;
				core.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				core.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				core.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				core.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				core.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				core.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				core.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				core.vertexExists(positionVertex, p7, PRECISION, v7);

				pConnects.push_back(v4);
				pConnects.push_back(v5);
				pConnects.push_back(v6);
				pConnects.push_back(v7);

				pCounts.push_back(4);

				//
				heTop = heTop.getNext().getSym().getNext();
				heBottom = heBottom.getPrev().getSym().getPrev();

				heTop_corner = heTop_corner.getNext().getSym().getNext();
				heBottom_corner = heBottom_corner.getPrev().getSym().getPrev();

			} while (!exit);
		}
		else
		{
			zItMeshHalfEdge heTop = heStart.getSym();
			zItMeshHalfEdge heBottom = heStart.getSym();
			heBottom = heBottom.getNext().getSym().getNext();
			heBottom = heBottom.getNext().getSym().getNext();

			zItMeshHalfEdge heTop_corner = heTop;
			zItMeshHalfEdge heBottom_corner = heBottom;

			for (int k = 0; k < featureWalkNum1; k++)
			{
				int blockStride = FeaturedNumStrides[k];
				//for (int& blockStride : FeaturedNumStrides)
				//{

				for (int i = 0; i < blockStride; i++)
				{
					heTop_corner = heTop_corner.getPrev().getPrev().getSym();
					heBottom_corner = heBottom_corner.getPrev().getPrev().getSym();
				}
			}

			exit = false;
			do
			{
				if (heTop.getStartVertex().getId() == endVID) exit = true;

				zPoint p0 = heTop.getVertex().getPosition();
				int v0;
				core.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				core.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				core.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				core.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				core.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				core.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				core.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				core.vertexExists(positionVertex, p7, PRECISION, v7);

				pConnects.push_back(v4);
				pConnects.push_back(v5);
				pConnects.push_back(v6);
				pConnects.push_back(v7);

				pCounts.push_back(4);

				//
				heTop = heTop.getPrev().getSym().getPrev();
				heBottom = heBottom.getNext().getSym().getNext();

				heTop_corner = heTop_corner.getPrev().getSym().getPrev();
				heBottom_corner = heBottom_corner.getNext().getSym().getNext();

			} while (!exit);

		}



		zObjMesh* o_sliceMesh = (left) ? &o_SliceMesh_Left : &o_SliceMesh_Right;
		zFnMesh fnMesh(*o_sliceMesh);

		fnMesh.create(positions, pCounts, pConnects);

		(left) ? fnMesh.setFaceColor(zGREEN) : fnMesh.setFaceColor(zMAGENTA);

		fnMesh.setFaceColor(zGREY);
		fnMesh.setVertexColor(zBLACK);

		//get start and end HE of the medial graph
		zFnMesh fnGuide(o_Mesh);
		zPointArray guidePts;
		fnGuide.getVertexPositions(guidePts);
		int vSIDSlice, vEIDSlice;
		vSIDSlice = core.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = core.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = util_getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
		zItMeshHalfEdge heEndSlice = heStartSlice;
		int sideDivisionCount = 1;
		while (heEndSlice.getVertex().getId() != vEIDSlice)
		{
			heEndSlice = heEndSlice.getNext().getSym().getNext();
			sideDivisionCount++;
		}
		heEndSlice = heEndSlice.getNext();

		//get connected he for the updated cornerVertex
		//find the he which, is not heStartSlice and end valency is not 3 (given that the inner face is only one stride)
		zItMeshVertex vCornerStart(*o_sliceMesh, vSIDSlice);
		zItMeshHalfEdge heCornerStartSection, heWalkSection;
		zItMeshHalfEdgeArray hesCTemp;
		vCornerStart.getConnectedHalfEdges(hesCTemp);
		for (zItMeshHalfEdge e : hesCTemp)
		{
			if (!e.getVertex().checkValency(3) && e.getId() != heStartSlice.getId())
			{
				heCornerStartSection = e;
				break;
			}
		}
		//walk on the slice mesh to color the edges
		//color the feature edges using numFeaturesStrides
		//color the corner edges using the corner vertex
		//color the bridging edges using the bridging stride
		//the walk changes based if it is the left or right slice mesh. To go to the next loop edge
		//For the left slice mesh, the walk goes sym.next.next
		//For the right slice mesh, the walk goes prev.prev.sym
		//Lets start by adding all the feature edges, when reaching the middle stride, change the color to the outer colors
		zItMeshHalfEdge heTemp1 = heStartSlice;
		zItMeshHalfEdgeArray hesFeatures;
		const int featureWalkNum = (FeaturedNumStrides.size()) + 1;
		for (int i = 0; i < featureWalkNum; i++)
		{
			int blockStride = i < FeaturedNumStrides.size()? FeaturedNumStrides[i] : 1;
			zColor edgeColor;
			if (i == 0) edgeColor = _col_in_corner;
			else if (i == featureWalkNum - 1) edgeColor = _col_out_corner;
			else if (i == (featureWalkNum / 2) - 1) edgeColor = _col_in_corner_st;
			else if (i == (featureWalkNum / 2)) edgeColor = _col_out_corner_st;
			else if (i < featureWalkNum / 2) edgeColor = _col_in_feature;
			else edgeColor = _col_out_feature;
			heTemp1.getEdge().setColor(edgeColor);
			heTemp1.getStartVertex().setColor(edgeColor);
			heTemp1.getVertex().setColor(edgeColor);

			zItMeshHalfEdge heRest = heTemp1;
			for (int k = 0; k < sideDivisionCount - 1; k++)
			{
				heRest = heRest.getNext().getSym().getNext();
				heRest.getEdge().setColor(edgeColor);
				heRest.getStartVertex().setColor(edgeColor);
				heRest.getVertex().setColor(edgeColor);

			}

			hesFeatures.push_back(heTemp1);
			for (int j = 0; j < blockStride; j++)
			{
				if (left)
				{
					heTemp1 = heTemp1.getSym().getNext().getNext();
				}
				else
				{
					heTemp1 = heTemp1.getPrev().getPrev().getSym();
				}	
			}
		}


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_SliceMesh_Regular(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		//get end he
		zItMeshHalfEdge heEnd = heStart;
		//int sideDivisionCount = 1;
		while (heEnd.getVertex().getId() != endVID)
		{
			heEnd = heEnd.getNext().getSym().getNext();
			//sideDivisionCount++;
		}

		//printf("\n sideDivisionCount %i ", sideDivisionCount);


		// walk along spine
		bool exit = false;
		zItMeshHalfEdge he = heStart;

		int featureStart = 0;

		zObjMesh* o_sliceMesh = &o_SliceMesh_Left;
		o_SliceMesh_Right.mesh.clear();

		zFnMesh fnMesh(*o_sliceMesh);


		zFnMesh fnGuideMesh(o_Mesh);
		fnGuideMesh.getVertexPositions(positions);
		fnGuideMesh.getPolygonData(pConnects, pCounts);
		fnMesh.create(positions, pCounts, pConnects);

		fnMesh.setFaceColor(zGREY);
		fnMesh.setVertexColor(zBLACK);



		//get start and end HE of the medial graph
		zFnMesh fnGuide(o_Mesh);
		zPointArray guidePts;
		fnGuide.getVertexPositions(guidePts);
		int vSIDSlice, vEIDSlice;
		vSIDSlice = core.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = core.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = util_getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
		zItMeshHalfEdge heEndSlice = heStartSlice;
		int sideDivisionCount = 1;
		while (heEndSlice.getVertex().getId() != vEIDSlice)
		{
			heEndSlice = heEndSlice.getNext().getSym().getNext();
			sideDivisionCount++;
		}

		heEndSlice = heEndSlice.getNext();

		//get all HE for feature edges, color them based on their side (exterior/interior)


		//get start of all feature edges
		zItMeshHalfEdgeArray hesFeatureInner, hesFeatureOuter;
		//zItMeshHalfEdge heTemp0 = heStartSlice;
		zItMeshHalfEdge heTemp = heStartSlice;
		zItMeshHalfEdge heSide;

		const int featureWalkNum = (FeaturedNumStrides.size());
		const int bridgingStride = FeaturedNumStrides[((FeaturedNumStrides.size() - 1) / 2)];
		printf("\n bridging stride %i %i %i ", (FeaturedNumStrides.size() - 1), ((FeaturedNumStrides.size() - 1) / 2), FeaturedNumStrides[((FeaturedNumStrides.size() - 1) / 2)]);
		bool isBackSide = true;
		bool isCornerEdge = false;
		bool fromInToOut = true;
		int innerCounter = 0;
		int outerCounter = 0;
		zColor innerColor = _col_in_feature;
		zColor outerColor = _col_out_feature;
		bool ignoreOuterEdge = false;


		////walk on the lower part of the mesh and color corner vertices
		//zItMeshHalfEdge heLowerTemp = heStartSlice.getSym().getNext();
		//int cornerCounter = 0;
		//while (heLowerTemp.getVertex().getId() != heStartSlice.getStartVertex().getId())
		//{
		//	zItMeshHalfEdge heS1, heS2;
		//	heS1 = heLowerTemp;
		//	heLowerTemp = heLowerTemp.getNext().getSym().getNext();
		//	heS2 = heLowerTemp;
		//	//get the angle between the two sides heS1 and heS2
		//	float angle = heS1.getVector().angle(heS2.getVector());
		//	bool chk = angle > 20;
		//	if (chk)
		//	{
		//		printf("\n angle0 %1.4f", angle);
		//		heLowerTemp.getStartVertex().setColor(zORANGE);
		//		cornerCounter++;
		//	}
		//}

		//get start inner corner vertex
		zItMeshVertex vCornerStart(*o_sliceMesh, StartCornerVID);
		zItMeshVertex vCornerEnd;
		vCornerStart.setColor(_col_in_corner_st);
		zItMeshHalfEdge startCornerHE;
		//get HE by finding all connected HE and choosing the one  with end valence of 3
		zItMeshHalfEdgeArray hesTemp;
		vCornerStart.getConnectedHalfEdges(hesTemp);
		for (zItMeshHalfEdge e : hesTemp)
		{
			if (e.getVertex().checkValency(3))
			{
				vCornerEnd = e.getVertex();
				break;
			}
		}


		//color vertex on the other end 
		int heCornerEdgeTopId = -1;
		zItMeshHalfEdge heCornerStart, heCornerEnd;
		//walk on the strides. color feature ones. then walk again to color the corner edges
		for (int k = 0; k < featureWalkNum * 2; k++)
		{
			int featureCounter = k % FeaturedNumStrides.size();
			int blockStride = FeaturedNumStrides[featureCounter];
			if (blockStride == bridgingStride)
			{
				isCornerEdge = true;
			}
			innerColor = isCornerEdge ? _col_in_corner : _col_in_feature;
			outerColor = isCornerEdge ? _col_out_corner : _col_out_feature;

			bool chk0 = isBackSide && heTemp.getStartVertex().getId() == vCornerStart.getId();
			if (chk0)
			{
				innerColor = vCornerStart.getColor();
				startCornerHE = heTemp;
			}
			bool chk1 = !isBackSide && heTemp.getStartVertex().getId() == vCornerEnd.getId();
			if (chk1)
			{
				outerColor = _col_out_corner_st;
			}
			zColor edgeColor = isBackSide ? innerColor : outerColor;
			heTemp.getEdge().setColor(edgeColor);
			heTemp.getStartVertex().setColor(edgeColor);
			heTemp.getVertex().setColor(edgeColor);

			zItMeshHalfEdge heRest = heTemp;
			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				heRest = heRest.getNext().getSym().getNext();
				heRest.getEdge().setColor(edgeColor);
				heRest.getVertex().setColor(edgeColor);
			}
			if (chk0) heCornerEnd = heRest;
			for (int i = 0; i < blockStride; i++)
			{
				heTemp = heTemp.getSym().getNext().getNext();
			}

			if (blockStride == bridgingStride)
			{
				isBackSide = !isBackSide;
				isCornerEdge = true;
			}
			else
			{
				isCornerEdge = false;
			}
		}

		//walk on the corner edges till you reach the next corner edge. Ignore edges that are colored

		if (blockType == zBlockType::Top)
		{



			zItMeshHalfEdge heWalk = heCornerEnd.getNext();
			bool flip = heWalk.getVertex().checkValency(3);
			if (flip) heWalk = heWalk.getSym().getNext();
			int safetyCounter = 0;
			while (safetyCounter < 10000)
			{
				heWalk.getEdge().setColor(_col_in_corner);
				heWalk.getVertex().setColor(_col_in_corner);

				zItMeshHalfEdge heOuter = heWalk.getSym().getNext().getNext();
				if (flip) heOuter = heWalk.getNext().getNext();

				heOuter.getEdge().setColor(_col_out_corner);
				heOuter.getVertex().setColor(_col_out_corner);

				heWalk = heWalk.getNext().getSym().getNext();
				if (heWalk.getStartVertex().checkValency(3))
				{
					break;
				}
				counter++;
			}
		}



		printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());
	}

	//Slice mesh: helper methods
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_MedialGraph(zObjMesh& o_Mesh, int startVID, int endVID)
	{

		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		//printf("\n hestart %i %i ", heStart.getStartVertex().getId(), heStart.getVertex().getId());

		zItMeshHalfEdge heStart_bottom = heStart;
		heStart_bottom = heStart_bottom.getPrev().getSym().getPrev();
		heStart_bottom = heStart_bottom.getPrev().getSym().getPrev();

		// compute graph
		zPointArray positions;
		zIntArray eConnects;

		zItMeshHalfEdge he = heStart;
		zItMeshHalfEdge he_bottom = heStart_bottom;

		//zPoint po = (startPoint + tmpPositions[heStart_bottom.getVertex().getId()]) * 0.5;
		zPoint po = startPoint;

		positions.push_back(po);
		bool exit = false;

		do
		{

			if (he.getVertex().getId() == endVID) exit = true;

			eConnects.push_back(positions.size() - 1);
			eConnects.push_back(positions.size());

			//zPoint p1 = (tmpPositions[he.getVertex().getId()] + tmpPositions[he_bottom.getStartVertex().getId()]) * 0.5;
			zPoint p1 = tmpPositions[he.getVertex().getId()];

			positions.push_back(p1);

			if (!exit)
			{
				he = he.getNext().getSym().getNext();
				he_bottom = he_bottom.getPrev().getSym().getPrev();
			}

		} while (!exit);

		zFnGraph fnMedial(o_MedialGraph);
		fnMedial.create(positions, eConnects);

		fnMedial.setEdgeWeight(5);
		fnMedial.setEdgeColor(zGREEN, false);
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_Medial_BraceEdges(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, int braceStride)
	{
		//compute start half edge
		zItMeshHalfEdge heStart = util_getStartHalfEdge(o_Mesh, startVID, endVID);

		// check if associated face normal is pointing down
		zVector down(0, 0, -1);
		zVector norm = heStart.getFace().getNormal();
		norm.normalize();

		bool flipHE = false;
		if (norm * down < 0.8)
		{
			flipHE = true;
			heStart = heStart.getSym();
		}

		// compute spine & brace edges
		bool exit = false;
		zItMeshHalfEdge he = heStart;



		do
		{
			//U direction
			zItMeshHalfEdge he_U = (flipHE) ? he.getNext() : he.getPrev();
			for (int i = 0; i < blockStride; i++)
			{
				if (i % braceStride == 0)
				{
					(flipHE) ? he_U.getPrev().getEdge().setColor(zMAGENTA) : he_U.getNext().getEdge().setColor(zMAGENTA);
				}
				if (i == blockStride - 1)
				{
					(flipHE) ? he_U.getNext().getEdge().setColor(zORANGE) : he_U.getPrev().getEdge().setColor(zORANGE);
				}

				he_U = (flipHE) ? he_U.getNext().getSym().getNext() : he_U.getPrev().getSym().getPrev();
			}


			//spine
			he.getEdge().setColor(zBLUE);
			he = (flipHE) ? he.getPrev().getSym().getPrev() : he.getNext().getSym().getNext();

			if (he == heStart) exit = true;

		} while (!exit);


	}


	//----------PRINT BLOCKS----------
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_PrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, bool compFrames, bool compSDF)
	{
	

		bool frameCHECKS = false;
		bool geomCHECKS = true;
		bool sdfCHECKS = true;

		int minCriticalPtsCount = INT_MAX;
		float bestPlaneSpacing = FLT_MAX;


		if (compFrames)
		{
			if (planarBlock)
			{
				for (float printPlaneSpacing = printHeightDomain.max; printPlaneSpacing >= printHeightDomain.min; printPlaneSpacing -= 0.00025)
				{
					compute_PrintSectionFromPlaneSpacing(printPlaneSpacing, printHeightDomain, frameCHECKS, sdfCHECKS, geomCHECKS);

					if (frameCHECKS) break;

					zFnPointCloud fnCritical_min(criticalMinLayer_pts);
					zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

					if (minCriticalPtsCount > fnCritical_min.numVertices() + fnCritical_max.numVertices())
					{
						minCriticalPtsCount = fnCritical_min.numVertices() + fnCritical_max.numVertices();
						bestPlaneSpacing = printPlaneSpacing;
					}

					if (actualPrintHeightDomain.min < printHeightDomain.min && actualPrintHeightDomain.max < printHeightDomain.max)
					{
						printf("\n minimum print height was reached and no solution was found.");
						break;
					}

				}

				if (!frameCHECKS)
				{
					printf("\n Layer height check FALSE!  Choose best planeXY spacing %1.4f", bestPlaneSpacing);

					compute_PrintSectionFromPlaneSpacing(bestPlaneSpacing, printHeightDomain, frameCHECKS, sdfCHECKS, geomCHECKS);
				}

				printf("\n frameCHECKS %s ", (frameCHECKS) ? "T" : "F");
				printf("\n layerCheck %i | geomChk %i | sdCheck %i ", (frameCHECKS), geomCHECKS, sdfCHECKS);
				printf("\n sectionFrames %i | o_sectionGraphs %i", sectionFrames.size(), o_sectionGraphs.size());
				compute_PrintBlock_ComputeTrimGraphs();

			}
			else
			{
				int numBlendFrames;
				zVector norm(0, 0, 1);
				vector<zItMeshHalfEdgeArray> vLoops;
				zObjMesh oMesh_top, oMesh_bottom;
				//vector<zObjMesh> oSectionMeshes;

				computeVLoops(o_SliceMesh_Left, medialIDS, FeaturedNumStrides, norm, vLoops, oMesh_top, oMesh_bottom);

				zScalarArray scalars;
				computeGeodesicScalars(o_SliceMesh_Left, vLoops, scalars, true);

				o_sectionMeshes.clear();
				computeGeodesicContours(vLoops, scalars, 0.01, oMesh_top, oMesh_bottom, o_sectionMeshes);;
				createSectionGraphs(o_sectionMeshes, o_sectionGraphs);
				o_sectionMeshesPar.clear();
				o_sectionMeshesPar.assign(o_sectionMeshes.size(), zObjMesh());
				
				sectionFrames.clear();
				sectionFrames.assign(o_sectionGraphs.size(), zTransform());
				compute_PrintBlock_ComputeTrimGraphs();

				o_contourHeightLines.clear();
				o_contourHeightLines.assign(o_sectionGraphs.size(), zObjGraph());

			}
		}
		

		


		if (compSDF)
		{
			printf("\n \n  SDF \n \n");

			compute_SDF(allSDFLayers, numSDFlayers, funcNum, numSmooth, printLayerWidth);
		}

	}

	//Print blocks: helper methods
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_PrintSectionFromPlaneSpacing(float printPlaneSpacing, zDomainFloat& _printHeightDomain, bool& frameCHECKS, bool& sdfCHECKS, bool& geomCHECKS)
	{
		frameCHECKS = false;
		geomCHECKS = true;
		sdfCHECKS = true;

		printf("\n printPlaneSpace %1.4f ", printPlaneSpacing);

		sectionFrames.clear();
		compute_PrintBlock_Frames(printPlaneSpacing, true, neopreneOffset.min, neopreneOffset.max);
		if (!isRegular) compute_PrintBlock_Frames(printPlaneSpacing, false, neopreneOffset.min, neopreneOffset.max);

		o_sectionGraphs.clear();
		o_sectionGraphs.assign(sectionFrames.size(), zObjGraph());
		bool geomChk0 = true;
		bool geomChk1 = true;
		//geomChk0 = 
		compute_PrintBlock_Sections(true, geomChk0);
		if (!isRegular) compute_PrintBlock_Sections(false, geomChk1);
		//if (geomChk0 && geomChk1)
		//{
		frameCHECKS = check_PrintLayerHeights(sdfCHECKS, geomCHECKS);

		//}
		if (!geomChk0 && !geomChk1)
		{
			printf("\n planeSpacing did not pass geometry check. Skip!");
		}

		printf("\n");
	}


	//Print blocks: frame methods
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_PrintBlock_Frames(float printPlaneSpacing, bool leftBlock, float neopreneOffset_start, float neopreneOffset_end)
	{
		// getLength of guide graph
		zDoubleArray eLens;
		zFnGraph fnGraph(o_MedialGraph);
		float totalLength = fnGraph.getEdgeLengths(eLens);

		//zFloatArray weights = { 0.0, 0.35, 1.0 };
		//zFloatArray multVals = { 0.0, 0.5,  1.0 };

		zFloatArray weights = { 0.0, 0.35, 0.70, 1.0 };
		zFloatArray multVals = { 0.0, 0.45, 0.80, 1.0 };



		float len = totalLength - (neopreneOffset_start + neopreneOffset_end);

		int numLayers = floor(len / printPlaneSpacing);

		float equalisedPlaneSpacing = len / numLayers;

		printf("\n %i %1.2f %1.2f ", numLayers, len, equalisedPlaneSpacing);

		zVector startNorm = (leftBlock) ? zVector(leftPlanes[0](2, 0), leftPlanes[0](2, 1), leftPlanes[0](2, 2)) : zVector(rightPlanes[0](2, 0), rightPlanes[0](2, 1), rightPlanes[0](2, 2));
		zVector endNorm = (leftBlock) ? zVector(leftPlanes[1](2, 0), leftPlanes[1](2, 1), leftPlanes[1](2, 2)) : zVector(rightPlanes[1](2, 0), rightPlanes[1](2, 1), rightPlanes[1](2, 2));

		zPoint startOrig = (leftBlock) ? zVector(leftPlanes[0](3, 0), leftPlanes[0](3, 1), leftPlanes[0](3, 2)) : zVector(rightPlanes[0](3, 0), rightPlanes[0](3, 1), rightPlanes[0](3, 2));
		zPoint endOrig = (leftBlock) ? zVector(leftPlanes[1](3, 0), leftPlanes[1](3, 1), leftPlanes[1](3, 2)) : zVector(rightPlanes[1](3, 0), rightPlanes[1](3, 1), rightPlanes[1](3, 2));



		zItGraphVertex v(o_MedialGraph, 0);
		zItGraphHalfEdge startHe = v.getHalfEdge();

		zItGraphVertex vEnd(o_MedialGraph, fnGraph.numVertices() - 1);
		zItGraphHalfEdge endHe = vEnd.getHalfEdge();

		// START
		zPoint start, end;
		float dIncrement = 0.0001;

		zItGraphHalfEdge walkHe = startHe;
		walkHe = walkHe.getSym();

		zPoint pOnCurve = v.getPosition();;

		bool exit = false;

		bool left = true;
		bool right = !isRegular;

		float dStart = 0;
		float dEnd = 0;

		while (!exit)
		{
			zPoint eEndPoint = walkHe.getVertex().getPosition();
			dStart += dIncrement;
			float distance_increment = dIncrement;
			while (pOnCurve.distanceTo(eEndPoint) < distance_increment)
			{
				distance_increment = distance_increment - pOnCurve.distanceTo(eEndPoint);
				pOnCurve = eEndPoint;

				walkHe = walkHe.getNext();
				eEndPoint = walkHe.getVertex().getPosition();
			}

			zVector he_vec = walkHe.getVector();
			he_vec.normalize();

			start = pOnCurve + he_vec * distance_increment;
			pOnCurve = start;
			// check 



			if (!right)
			{
				zPoint startPlanePoint = zVector(rightPlanes[0](3, 0), rightPlanes[0](3, 1), rightPlanes[0](3, 2));
				zPoint startPlaneNormal = zVector(rightPlanes[0](2, 0), rightPlanes[0](2, 1), rightPlanes[0](2, 2));

				float dStart = core.minDist_Point_Plane(pOnCurve, startPlanePoint, startPlaneNormal);

				if (abs(dStart) >= neopreneOffset_start)  right = true;


			}

			if (!left)
			{
				zPoint startPlanePoint = zVector(leftPlanes[0](3, 0), leftPlanes[0](3, 1), leftPlanes[0](3, 2));
				zPoint startPlaneNormal = zVector(leftPlanes[0](2, 0), leftPlanes[0](2, 1), leftPlanes[0](2, 2));

				float dStart = core.minDist_Point_Plane(pOnCurve, startPlanePoint, startPlaneNormal);
				if (abs(dStart) >= neopreneOffset_start)  left = true;
			}

			if (right && left) exit = true;

		}


		// END

		walkHe = endHe;
		walkHe = walkHe.getSym();
		pOnCurve = vEnd.getPosition();;

		exit = false;
		left = true;
		right = !isRegular;

		while (!exit)
		{
			zPoint eEndPoint = walkHe.getVertex().getPosition();
			dEnd += dIncrement;
			float distance_increment = dIncrement;
			while (pOnCurve.distanceTo(eEndPoint) < distance_increment)
			{
				distance_increment = distance_increment - pOnCurve.distanceTo(eEndPoint);
				pOnCurve = eEndPoint;

				walkHe = walkHe.getNext();
				eEndPoint = walkHe.getVertex().getPosition();
			}

			zVector he_vec = walkHe.getVector();
			he_vec.normalize();

			end = pOnCurve + he_vec * distance_increment;
			pOnCurve = end;
			// check 



			if (!right)
			{
				zPoint endPlanePoint = zVector(rightPlanes[1](3, 0), rightPlanes[1](3, 1), rightPlanes[1](3, 2));
				zPoint endPlaneNormal = zVector(rightPlanes[1](2, 0), rightPlanes[1](2, 1), rightPlanes[1](2, 2));

				float dEnd = core.minDist_Point_Plane(pOnCurve, endPlanePoint, endPlaneNormal);

				if (abs(dEnd) >= neopreneOffset_end)  right = true;

			}

			if (!left)
			{
				zPoint endPlanePoint = zVector(leftPlanes[1](3, 0), leftPlanes[1](3, 1), leftPlanes[1](3, 2));
				zPoint endPlaneNormal = zVector(leftPlanes[1](2, 0), leftPlanes[1](2, 1), leftPlanes[1](2, 2));

				float dEnd = core.minDist_Point_Plane(pOnCurve, endPlanePoint, endPlaneNormal);

				if (abs(dEnd) >= neopreneOffset_end)  left = true;

			}

			if (right && left) exit = true;

		}

		bool interpolateOrigin = _interpolateFramesOrigins;

		// Start point

		zPoint O = start;
		zVector Z = startNorm;

		zVector tempZ = Z;
		tempZ.normalize();

		zPoint tempO;

		zVector X;
		zVector Y(0, 1, 0);

		float weight = 0;

		float mult;

		for (int l = 0; l < weights.size() - 1; l++)
		{
			if (weight >= weights[l] && weight <= weights[l + 1])   mult = core.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
		}

		tempZ.x = (startNorm.x * (1 - mult)) + (endNorm.x * mult);
		tempZ.y = (startNorm.y * (1 - mult)) + (endNorm.y * mult);
		tempZ.z = (startNorm.z * (1 - mult)) + (endNorm.z * mult);
		tempZ.normalize();

		X = Y ^ tempZ;
		X.normalize();

		Y = tempZ ^ X;
		Y.normalize();

		//zTransform pFrame = setTransformFromVectors(O, X, Y, tempZ);
		zTransform pFrame = core.getTransformFromOrigin_Normal(O, tempZ, zVector(1, 1, 0));
		sectionFrames.push_back(pFrame);

		pOnCurve = O;
		dStart = core.minDist_Point_Plane(pOnCurve, startOrig, startNorm);


		// in between points
		walkHe = startHe;
		for (int j = 0; j < numLayers; j++)
		{
			zPoint prevPoint = pOnCurve;

			zPoint eEndPoint = walkHe.getVertex().getPosition();

			float distance_increment = equalisedPlaneSpacing;

			while (pOnCurve.distanceTo(eEndPoint) < distance_increment)
			{
				distance_increment = distance_increment - pOnCurve.distanceTo(eEndPoint);
				pOnCurve = eEndPoint;

				walkHe = (walkHe.onBoundary()) ? walkHe.getNext() : walkHe.getNext().getSym().getNext();
				eEndPoint = walkHe.getVertex().getPosition();
			}

			zVector he_vec = walkHe.getVector();
			he_vec.normalize();


			//O
			O = pOnCurve + he_vec * distance_increment;

			if (j == numLayers - 1)
			{
				O = end;
				Z = endNorm;
			}

			if (interpolateOrigin)
			{

				tempO = (startOrig * (1 - mult)) + (endOrig * mult);
				O = tempO;
			}

			tempZ = Z;
			tempZ.normalize();

			Y = zVector(1, 1, 0);

			weight = (float)(j + 1) / numLayers;

			for (int l = 0; l < weights.size() - 1; l++)
			{
				if (weight >= weights[l] && weight <= weights[l + 1])   mult = core.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
			}

			tempZ.x = (startNorm.x * (1 - mult)) + (endNorm.x * mult);
			tempZ.y = (startNorm.y * (1 - mult)) + (endNorm.y * mult);
			tempZ.z = (startNorm.z * (1 - mult)) + (endNorm.z * mult);
			tempZ.normalize();

			X = Y ^ tempZ;
			X.normalize();

			Y = tempZ ^ X;
			Y.normalize();

			// add frame
			//pFrame = setTransformFromVectors(O, X, Y, tempZ);
			pFrame = core.getTransformFromOrigin_Normal(O, tempZ, zVector(1, 1, 0));
			sectionFrames.push_back(pFrame);

			pOnCurve = O;

			if (j == numLayers - 1)
			{
				float dEnd = core.minDist_Point_Plane(pOnCurve, endOrig, endNorm);

				if (leftBlock) printf(" \n left sD %1.4f eD %1.4f ", dStart, dEnd);
				else printf(" \n right sD %1.4f eD %1.4f ", dStart, dEnd);
			}
		}




		//supsitute frames in a new calculation way
		//zPlane startPlane = coreUtils.getPlaneFromOrigin_Normal(startOrig, startNorm);
		zPointArray origins;
		numLayers += 1;// isRegular ? sectionFrames.size() : sectionFrames.size() / 2;
		for (int i = 0; i < numLayers; i++)
		{
			origins.push_back(zVector(sectionFrames[i](3, 0), sectionFrames[i](3, 1), sectionFrames[i](3, 2)));
		}
		vector<zPlane> leftFrames, rightFrames;
		//printf("\n numLayers %i | sectionFrames %i | origins %i", numLayers, sectionFrames.size(), origins.size());
		core.interpolatePlanes_slerp(leftPlanes[0], leftPlanes[1], numLayers, origins, leftFrames);
		if (!isRegular) core.interpolatePlanes_slerp(rightPlanes[0], rightPlanes[1], numLayers, origins, rightFrames);
		sectionFrames.clear();
		sectionFrames.insert(sectionFrames.end(), leftFrames.begin(), leftFrames.end());
		if (!isRegular) sectionFrames.insert(sectionFrames.end(), rightFrames.begin(), rightFrames.end());

		printf("\n numLayers %i | sectionFrames %i | origins %i \n", numLayers, sectionFrames.size(), origins.size());




	}

	//Print blocks: section methods
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_PrintBlock_Sections(bool left, bool& outGeomChk)
	{


		zScalarArray scalars;

		int start = 0;
		int end = sectionFrames.size();

		if (!isRegular)
		{
			start = (left) ? 0 : floor(sectionFrames.size() * 0.5);
			end = (left) ? floor(sectionFrames.size() * 0.5) : sectionFrames.size();
		}

		zObjMesh* o_SliceMesh = (left) ? &o_SliceMesh_Left : &o_SliceMesh_Right;
		zFnMesh  fn_sliceMesh(*o_SliceMesh);
		outGeomChk = true;
		for (int i = start; i < end; i++)
		{
			scalars.clear();
			zPoint O(sectionFrames[i](3, 0), sectionFrames[i](3, 1), sectionFrames[i](3, 2));
			zVector N(sectionFrames[i](2, 0), sectionFrames[i](2, 1), sectionFrames[i](2, 2));

			zVector X(sectionFrames[i](0, 0), sectionFrames[i](0, 1), sectionFrames[i](0, 2));
			X.normalize();

			for (zItMeshVertex v(*o_SliceMesh); !v.end(); v++)
			{
				zPoint P = v.getPosition();
				float minDist_Plane = core.minDist_Point_Plane(P, O, N);
				scalars.push_back(minDist_Plane);
			}
			int startPres = 3;
			int endPres = 6;
			bool geomChk = false;
			for (int pres = startPres; pres <= endPres; pres++)
			{
				zPointArray positions;
				zIntArray edgeConnects;
				zColorArray vColors;
				//int pres = 3;
				fn_sliceMesh.getIsoContour(scalars, 0.0, positions, edgeConnects, vColors, pres, pow(10, -1 * pres));

				// create graphs
				zFnGraph tempFn(o_sectionGraphs[i]);
				tempFn.create(positions, edgeConnects);
				//tempFn.setEdgeColor(zGREY);
				tempFn.setVertexColors(vColors, false);
				geomChk = check_sectionGraphGeomCheck(o_sectionGraphs[i]);
				if (geomChk || i == 0)
				{
					break;
				}
				else
				{
					//printf("\n -- geom check failed for pres %i section [%i] ", pres,i);
				}

			}
			if (!geomChk)
			{
				outGeomChk = false;
				printf("\n outGeometryCheck FAILED for block[%i] section [%i]!", blockId, i);

			}
			//return true;
		}
		//outGeomChk = geom
	}

	//Print blocks: check methods
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::check_PrintLayerHeights(bool& checkSDF, bool& checkGeometry)
	{
		float minLayerHeight = 10;
		float maxLayerHeight = 0;

		//float minLayerHeight = 0.005;
		//float maxLayerHeight = 0.015;
		int minHeightGraphID = -1;
		//bool 
		checkOranges = true;
		//bool 
		checkMagentas = true;

		zFnPointCloud fnCritical_min(criticalMinLayer_pts);
		zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

		fnCritical_min.clear();
		fnCritical_max.clear();

		int r0 = 0;
		int r1 = floor(o_sectionGraphs.size() * 0.5) - 1;
		int r2 = floor(o_sectionGraphs.size() * 0.5);
		int r3 = (o_sectionGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		//bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;
		bool deckBlock = true;

		//int end = (leftPlaneExists && rightPlaneExists) ? floor(o_sectionGraphs.size() * 0.5) : o_sectionGraphs.size();
		int end = floor(o_sectionGraphs.size() * 0.5);

		float printLength = 0;

		// PRINT LAYERS
		for (int j = 0; j < end; j++)
		{
			//int kStart = (!deckBlock && !rightPlaneExists) ? 1 : 0;
			//int kEnd = (!deckBlock && !leftPlaneExists) ? 1 : 2;

			int kStart = 0;
			int kEnd = 2;

			for (int k = kStart; k < kEnd; k++)
			{
				int i = (k == 0) ? j : j + end;

				if (!deckBlock) i = j;

				if (i == r0) continue;

				if (deckBlock && i == r2) continue;


				zVector norm(sectionFrames[i](2, 0), sectionFrames[i](2, 1), sectionFrames[i](2, 2));
				norm *= -1;
				norm.normalize();

				zVector prevNorm(sectionFrames[i - 1](2, 0), sectionFrames[i - 1](2, 1), sectionFrames[i - 1](2, 2));
				zVector prevOrigin(sectionFrames[i - 1](3, 0), sectionFrames[i - 1](3, 1), sectionFrames[i - 1](3, 2));

				float layerWidth = 0.025;

				int orangeCounter = 0;
				int magentaCounter = 0;

				int innerSideCounter = 0;
				int outerSideCounter = 0;

				//The following check to see if there is a problem with the geometry/section graph
				//01-check if the graph is closed by checking if any vertex is boundary (val (1))


				for (zItGraphVertex v(o_sectionGraphs[i]); !v.end(); v++)
				{
					/*if (v.checkValency(1))
					{
						printf("\n ")
					}*/



					zPoint p = v.getPosition();

					if (v.getColor() == zORANGE) orangeCounter++;
					if (v.getColor() == zMAGENTA) magentaCounter++;

					if (v.getColor() == _col_in_feature) innerSideCounter++;
					if (v.getColor() == _col_out_feature) outerSideCounter++;


					zPoint p1 = p + norm * 1.0;

					zPoint intPt;
					bool check = core.line_PlaneIntersection(p, p1, prevNorm, prevOrigin, intPt);

					float layerHeight = intPt.distanceTo(p);
					maxLayerHeight = (layerHeight > maxLayerHeight) ? layerHeight : maxLayerHeight;
					minHeightGraphID = (layerHeight < minLayerHeight) ? i : minHeightGraphID;
					minLayerHeight = (layerHeight < minLayerHeight) ? layerHeight : minLayerHeight;

					if (layerHeight < printHeightDomain.min)fnCritical_min.addPosition(p);
					if (layerHeight > printHeightDomain.max)fnCritical_max.addPosition(p);

				}

				if (orangeCounter != 2)
				{
					//printf("\n checkOrange error Layer %i | %i ", i, orangeCounter);

					zFnGraph fnGraph(o_sectionGraphs[i]);
					//fnGraph.setEdgeColor(zORANGE);
					//checkOranges = false;


				}

				if (magentaCounter != numMagentaLoops * 2)
				{
					//printf("\n checkMagenta error Layer %i | %i ", i, magentaCounter);
					zFnGraph fnGraph(o_sectionGraphs[i]);
					//fnGraph.setEdgeColor(zMAGENTA);
					//checkMagentas = false;
				}

				if (innerSideCounter == outerSideCounter)
				{
					//printf("\n inner side and outer side has the same count %i", innerSideCounter);
				}
				else
				{
					//printf("\n inner side and outer side are NOT equal !!!  inner | outer  %i | %i", innerSideCounter, outerSideCounter);
				}

				for (zItGraphEdge e(o_sectionGraphs[i]); !e.end(); e++)
				{
					printLength += e.getLength();
				}


			}
		}

		fnCritical_min.setVertexColor(zORANGE);
		fnCritical_max.setVertexColor(zMAGENTA);


		bool out = true;

		if (out)
		{
			if (minLayerHeight < printHeightDomain.min) out = false;
			if (minLayerHeight > printHeightDomain.max) out = false;

			if (maxLayerHeight < printHeightDomain.min) out = false;
			if (maxLayerHeight > printHeightDomain.max) out = false;
		}

		actualPrintHeightDomain.min = minLayerHeight;
		actualPrintHeightDomain.max = maxLayerHeight;
		checkSDF = (!checkOranges || !checkMagentas) ? false : true;
		if (!checkSDF) out = false;

		int precision = 4;
		bool checkLeftPlanar = check_InterfacePoints(true, pow(10, -1 * precision));
		bool checkRightPlanar = isRegular ? true : check_InterfacePoints(false, pow(10, -1 * precision));

		checkGeometry = (!checkLeftPlanar || !checkRightPlanar) ? false : true;
		//if (!checkGeometry) out = false;







		printf("\n block| %1.4f %1.4f| %1.1f | chkOrange %s | chkMagenta %s | pathExist %s | < min ht intersectionPts %i  | > max ht intersectionPts %i ", minLayerHeight, maxLayerHeight, printLength, (checkOranges) ? "T" : "F", (checkMagentas) ? "T" : "F", (checkGeometry) ? "T" : "F", fnCritical_min.numVertices(), fnCritical_max.numVertices());

		/*zPointArray maxPts;
		fnCritical_max.getVertexPositions(maxPts);
		for (auto& p : maxPts)
		{
			cout << endl << p;
		}*/

		return out;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::check_InterfacePoints(bool left, float distTolerance)
	{
		zObjMesh* oMesh = (left) ? &o_SliceMesh_Left : &o_SliceMesh_Right;
		zTransform* starEnd = (left) ? &leftPlanes[0] : &rightPlanes[0];

		//if (left && !leftPlaneExists) return true;
		//if (!left && !rightPlaneExists) return true;



		zFnMesh fnMesh(*oMesh);
		float tol = 0.02;
		fnMesh.setFaceColor(zGREY);

		// get face points
		zIntArray interfaceVerts_start;
		zIntArray interfaceVerts_end;

		zPoint sPoint(starEnd[0](3, 0), starEnd[0](3, 1), starEnd[0](3, 2));
		zVector sNorm(starEnd[0](2, 0), starEnd[0](2, 1), starEnd[0](2, 2));
		sNorm.normalize();
		sNorm *= -1;

		zPoint ePoint(starEnd[1](3, 0), starEnd[1](3, 1), starEnd[1](3, 2));
		zVector eNorm(starEnd[1](2, 0), starEnd[1](2, 1), starEnd[1](2, 2));
		eNorm.normalize();

		for (zItMeshFace f(*oMesh); !f.end(); f++)
		{
			zVector fNorm = f.getNormal();
			fNorm.normalize();
			if (fNorm * sNorm >= 1 - tol)
			{
				zIntArray fVerts;
				f.getVertices(fVerts);
				for (auto& fV : fVerts) interfaceVerts_start.push_back(fV);
				f.setColor(zGREEN);
			}
			else if (fNorm * eNorm >= 1 - tol)
			{
				zIntArray fVerts;
				f.getVertices(fVerts);
				for (auto& fV : fVerts) interfaceVerts_end.push_back(fV);
				f.setColor(zBLUE);
			}
		}

		// check for in Plane
		zPoint* vPositions = fnMesh.getRawVertexPositions();
		bool chkStart = true;
		bool chkEnd = true;


		for (auto& vId : interfaceVerts_start)
		{
			float d = core.minDist_Point_Plane(vPositions[vId], sPoint, sNorm);
			//printf("\n start %1.4f ", d);
			if (d > distTolerance) chkStart = false;
			if (!chkStart) break;
		}
		for (auto& vId : interfaceVerts_end)
		{
			float d = core.minDist_Point_Plane(vPositions[vId], ePoint, eNorm);
			//printf("\n end %1.4f ", d);
			if (d > distTolerance) chkEnd = false;
			if (!chkEnd) break;
		}
		bool out = (!chkStart || !chkEnd) ? false : true;
		if (!out) printf("\n checkGeom Fail! %s | start-End %s - %s ", left ? "Left" : "Right", chkStart ? "True" : "false", chkEnd ? "True" : "false");

		return out;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::check_sectionGraphGeomCheck(zObjGraph& graph)
	{
		//checks if the graph pass all geometry checks
		//Check 1 : check if the graph is closed
		//Check 2 : check if all four corners are found
		bool chkClosed = true;
		bool corner0 = false;
		bool corner1 = false;
		bool corner2 = false;
		bool corner3 = false;
		bool featureInner = false;
		bool featureOuter = false;
		bool chkCorners = false;
		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (v.checkValency(1))
			{
				chkClosed = false;
			}
			if (v.getColor() == _col_in_corner_st) corner0 = true;
			if (v.getColor() == _col_out_corner_st) corner1 = true;
			if (v.getColor() == _col_in_corner) corner2 = true;
			if (v.getColor() == _col_out_corner) corner3 = true;
			if(v.getColor() == _col_in_feature) featureInner = true;
			if(v.getColor() == _col_out_feature) featureOuter = true;

			chkCorners = corner0 && corner1 && corner2 && corner3;
			if(blockType == zBlockType::Bottom) chkCorners = chkCorners && featureInner && featureOuter;
			/*if (chkCorners && chkClosed)
			{
				break
			}*/
		}
		if (chkClosed && chkCorners)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	

	//Print blocks: trim methods

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_PrintBlock_ComputeTrimGraphs()
	{
		o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing.clear();
		o_trimGraphs_bracing.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_features_hard.clear();
		o_trimGraphs_features_hard.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_features_soft.clear();
		o_trimGraphs_features_soft.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_SlotSide.clear();
		o_trimGraphs_SlotSide.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_seamAlignment.clear();
		o_trimGraphs_seamAlignment.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs_bracing_flat.clear();
		o_trimGraphs_bracing_flat.assign(o_sectionGraphs.size(), zObjGraph());
		

		for (int i = 0; i < o_sectionGraphs.size(); i++)
		{

			//The following check to see if there is a problem with the geometry/section graph
				//01-check if the graph is closed by checking if any vertex is boundary (val (1))

			if (i == 0) continue;
			if(!isRegular && i == floor(o_sectionGraphs.size() * 0.5)) continue;

			for (zItGraphVertex v(o_sectionGraphs[i]); !v.end(); v++)
			{
				if (v.checkValency(1))
				{
					printf("\n section[%i] is not closed!", i);
				}
			}
			zFnGraph fng;
			zObjGraphArray allTrimGraphs;

			compute_TrimGraphs_BoundaryFeature(i, o_trimGraphs_features_hard[i], o_trimGraphs_features_soft[i]);
			compute_TrimGraphs_SlotSide(i, o_trimGraphs_SlotSide[i]);
			//compute_TrimGraphs_SeamAlignment(i, o_trimGraphs_seamAlignment[i]);
			if (blockType != zBlockType::Wall)
			{
				compute_TrimGraphs_BracingCable(i, o_trimGraphs_bracing[i]);
				fng = zFnGraph(o_trimGraphs_bracing[i]);
				fng.setEdgeColor(zBLUE);
				allTrimGraphs.push_back(o_trimGraphs_bracing[i]);
			}
			else
			{
				compute_TrimGraphs_BracingWall(i, o_trimGraphs_bracing[i]);
				fng = zFnGraph(o_trimGraphs_bracing[i]);
				fng.setEdgeColor(zBLUE);
				allTrimGraphs.push_back(o_trimGraphs_bracing[i]);
			}

			fng = zFnGraph(o_trimGraphs_features_hard[i]);
			fng.setEdgeColor(zRED);
			fng = zFnGraph(o_trimGraphs_features_soft[i]);
			fng.setEdgeColor(zGREEN);
			fng = zFnGraph(o_trimGraphs_SlotSide[i]);
			fng.setEdgeColor(zGREEN);


			allTrimGraphs.push_back(o_trimGraphs_features_hard[i]);
			allTrimGraphs.push_back(o_trimGraphs_features_soft[i]);
			

			util_combineMultipleGraphs(allTrimGraphs, o_trimGraphs[i]);




		}
		printf("\n trims finished");

	}
	
	
	
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_BracingCable(int graphId, zObjGraph& outGraph)
	{
		zFnGraph fnSectionG(o_sectionGraphs[graphId]);
		//get intersection point with cable
		int cableGraphId = compute_cable_CableGraphIndexPerGraph(graphId);
		if (cableGraphId == -1) return;
		zPointArray pts;
		compute_cable_CableSectionPoints(graphId, o_CableGraphs[cableGraphId], pts);
		//cout << endl << "cableGraphId " << graphId << " | " << pts[0];
		//we should have one point (will update method later)
		zPoint cablePoint = pts[0];
		//find closest pt from p to graph section points
		//first we need to choose which graph segments we want to find closest point to
		//then we create trim graph (these graphs will be the bracing as well)

		bool found = false;
		zItGraphHalfEdgeArray hes0, hes1, hes2;
		zPointArray gPositions;
		zIntArray gEdgeCOnnects;
		zFloatArray distances;
		gPositions.push_back(cablePoint);
		//printf("\n section: ");

		//00
		found = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_in_corner_st, _col_out_corner_st, hes0);
		if (found)
		{
			zPoint p;
			float d;
			util_getHeArrayClosestPoint(hes0, cablePoint, p, d);
			if (blockType == zBlockType::Bottom)
			{
				zObjGraph tg;
				zFnGraph tgf(tg);
				zPointArray tpts;
				util_createGraphFromHEArray(hes0, tg);
				tgf.getVertexPositions(tpts);
				zVector v = tpts[tpts.size() - 1] - tpts[0];
				v *= 0.3;
				p = tpts[0] + v;
			}
			gPositions.push_back(p);
			gEdgeCOnnects.push_back(0);
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			distances.push_back(d);
			if (d < 0.005) printf("\n [%i] not found 0 %i | %1.4f", graphId, hes0.size(), d);

		}
		else printf("\n [%i] not found 0 %i", graphId, hes0.size());

		//01 
		found = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_in_corner_st, _col_in_corner, hes1);
		if(blockType == zBlockType::Bottom) found = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_in_corner, _col_out_corner, hes1);
		if (found)
		{
			zPoint p;
			float d;
			util_getHeArrayClosestPoint(hes1, cablePoint, p, d);

			/*if (blockType == zBlockType::Bottom)
			{
				zObjGraph tg;
				zFnGraph tgf(tg);
				zPointArray tpts;
				createGraphFromHEArray(hes1, tg);
				tgf.getVertexPositions(tpts);
				zVector v = tpts[tpts.size() - 1] - tpts[0];
				v *= 0.3;
				p = tpts[0] + v;

			}*/


			gPositions.push_back(p);
			gEdgeCOnnects.push_back(0);
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			distances.push_back(d);

			for (zItGraphHalfEdge e : hes1)
			{
				e.getEdge().setColor(zCYAN);
			}

			if (d < 0.005) printf("\n [%i] not found 1 %i | %1.4f", graphId, hes1.size(), d);



		}
		else printf("\n [%i] not found 1 %i", graphId, hes1.size());

		//02
		if(isCorner || blockType == zBlockType::Bottom) found = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_out_corner_st, _col_out_feature, hes2);
		else found = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_out_corner, _col_out_feature, hes2);
		if (found)
		{
			zPoint p;
			float d;
			util_getHeArrayClosestPoint(hes2, cablePoint, p, d);

			if (isCorner || blockType == zBlockType::Bottom)
			{
				zObjGraph tg;
				zFnGraph tgf(tg);
				zPointArray tpts;
				util_createGraphFromHEArray(hes2, tg);
				tgf.getVertexPositions(tpts);
				zVector v = tpts[tpts.size() - 1] - tpts[0];
				v *= 0.3;
				p = tpts[0] + v;

				//extend point
				v = p-cablePoint;
				float tempExtend = 0.0;
				zPoint outPt;
				util_getGraphClosestPoint(o_sectionGraphs[graphId], p, outPt, tempExtend);
				v.normalize();
				v *= tempExtend;
				p += v;
			
			}
			
			
			
			gPositions.push_back(p);
			gEdgeCOnnects.push_back(0);
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			distances.push_back(d);
			if (d < 0.005) printf("\n [%i] not found 2 %i | %1.4f", graphId, hes2.size(), d);


		}
		else printf("\n [%i] not found 2 %i", graphId, hes2.size());

		if (isCorner || blockType == zBlockType::Bottom)
		{
			//add extra bracing between corner point exterior and mid of interior
			
			//if (found)
			{
				zItGraphHalfEdgeArray hes3, hes4;
				//found = getShortestHEsBetweenColors(o_sectionGraphs[graphId], _colorCornersEnd, _colorFeatureOuter, hes3);
				bool found1 = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_out_corner_st, _col_out_feature, hes3);
				bool found2 = util_getShortestHEsBetweenColors(o_sectionGraphs[graphId], _col_in_corner, _col_out_corner, hes4);
				
				zPoint bracingStart, bracingEnd;

				
				zObjGraph tempGraphStart;
				zFnGraph tempFnStart(tempGraphStart);
				util_createGraphFromHEArray(hes3, tempGraphStart);

				zPointArray tempPts;
				tempFnStart.getVertexPositions(tempPts);
				zVector tempVec = tempPts[tempPts.size() - 1] - tempPts[0];
				tempVec *= 0.75;
				bracingEnd = tempPts[0] + tempVec;
				if (isCorner)
				{
					//bracingStart = tempFnStart.getCenter();
					for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
					{
						//iterate till finding corener exterior, then iterate through edges till finding the one between corner start and feature (snap to the closest point)
						if (v.getColor() == _col_in_corner)
						{
							bracingStart = v.getPosition();
							found == true;
							break;
						}

					}
				}
				else
				{
					//bracingStart = tempFnStart.getCenter();
					zObjGraph tempGraphEnd;
					zFnGraph tempFnEnd(tempGraphEnd);
					util_createGraphFromHEArray(hes4, tempGraphEnd);
					//bracingEnd = tempFnEnd.getCenter();
					zPointArray tempPts2;
					tempFnEnd.getVertexPositions(tempPts2);
					zVector tempVec2 = tempPts2[tempPts2.size() - 1] - tempPts2[0];
					tempVec2 *= 0.75;
					bracingStart = tempPts2[0] + tempVec2;

					
					//bracingEnd = (tempPts2[0] + tempPts2[tempPts2.size() - 1]) * 0.5;
				}

				//extend point
				zVector v = bracingEnd - bracingStart;
				float tempExtend = 0.0;
				zPoint outPt;
				util_getGraphClosestPoint(o_sectionGraphs[graphId], bracingEnd, outPt, tempExtend);
				v.normalize();
				v *= tempExtend;
				bracingEnd += v;

				
				printf("\n found %i", found);
				//if (found)
				//if (found)
				{
					////get the center then snap it to the closest point on that graph
					//zObjGraph tempGraph;
					//createGraphFromHEArray(hes3, tempGraph);
					//zFnGraph tempFn(tempGraph);
					//zPointArray tempPts;
					//tempFn.getVertexPositions(tempPts);
					//zPoint tempCenter = tempFn.getCenter();
					//int index = coreUtils.getClosest_PointCloud(tempCenter, tempPts);

					gPositions.push_back(bracingStart);
					gEdgeCOnnects.push_back(gPositions.size() - 1);
					//gPositions.push_back(tempPts[index]);
					gPositions.push_back(bracingEnd);
					gEdgeCOnnects.push_back(gPositions.size() - 1);


				}
			}
		}
		




		//make graph from the positions and choose the longest 2 to create perpendicular slot graph
		zFnGraph fnTrimG(outGraph);
		fnTrimG.create(gPositions, gEdgeCOnnects);


		////printf("\n graph[%i] %i | %i", graphId, gPositions.size(), gEdgeCOnnects.size());
		//for (int i = 0; i < fnTrimG.numEdges(); i++)
		//{
		//	int prevId = (i + 1) % fnTrimG.numEdges();
		//	zItGraphEdge e(o_trimGraphs[graphId], i);
		//	zItGraphEdge ePrev(o_trimGraphs[graphId], prevId);
		//	if (abs(ePrev.getLength() - zItGraphEdge(o_trimGraphs[graphId], i).getLength()) < EPS)
		//	{
		//		printf("\n graph[%i] distance %1.4f", graphId, e.getLength());
		//		zPointArray pts;
		//		e.getVertexPositions(pts);
		//		cout << endl << pts[0];
		//		cout << endl << pts[1] << endl;
		//		ePrev.getVertexPositions(pts);
		//		cout << endl << pts[0];
		//		cout << endl << pts[1] << endl;
		//	}
		//}
		//coreUtils.checkRepeatVector()

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_BoundaryFeature(int graphId, zObjGraph& outGraph_hardFeature, zObjGraph& outGraph_softFeature)
	{
		//create a trim graph at each trim points. Trim points are all feature curves and corner
		//The graph will be made by avg the two outgoing vector from the feature vertex

		zPrintParamSDF _printParameters;

		zPointArray posittions;
		zIntArray eConnect;
		float length = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior + (_printParameters.printWidthInterior * 3);
		float angleThreshold = 30;
		float lengthTolerance = 0.05;
		zItGraphVertexArray features_hard;
		zItGraphVertexArray features_soft;
		//check if there is another point that has already been added within the tolerance
		//check the next vertex, if it is within the same tolerance, combine the two together and skip the next vertex
		//To do that, probably since we don't know if there are some points in between, better to iterate through all the vertices first, and then oterate through the ones that passes that check.

		for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
		{

			bool cornerByTopology = v.getColor() == _col_in_corner_st ||
									v.getColor() == _col_out_corner_st ||
									v.getColor() == _col_in_corner ||
									v.getColor() == _col_out_corner;

			bool otherFeature = v.getColor() == _col_in_feature ||
								v.getColor() == _col_out_feature;

			/*if (blockType == zBlockType::Bottom)
			{
				cornerByTopology =	v.getColor() == _col_in_corner_st ||
									v.getColor() == _col_in_corner;
			}*/
			if (cornerByTopology) features_hard.push_back(v);
			else
			{
				zItGraphHalfEdgeArray hes;
				v.getConnectedHalfEdges(hes);
				zVector v0 = hes[0].getVector();
				zVector v1 = hes[1].getVector();
				v0.normalize();
				v1.normalize();

				float angle = v0.angle(v1 * -1);
				bool angleChk = angle >= angleThreshold;
				if (angleChk)
				{
					if (otherFeature) features_hard.push_back(v);
					else features_soft.push_back(v);
				}
			}
		}

		for (int index = 0; index < features_hard.size(); index++)
		{
			zVector v0, v1, result;
			zPoint pt;

			int nextIndex = (index + 1) % features_hard.size();
			//check the distance to the next index, if it is within a threshold, get middle vertex as feature and skip the next feature
			float dist = features_hard[index].getPosition().distanceTo(features_hard[nextIndex].getPosition());
			if (dist <= lengthTolerance)
			{
				//get the two vectors, average them twice. They might be sharing the same vector, but not necessarily if the have points in between
				v0 = util_averageVectorsAtGraphVertex(features_hard[index]);
				v1 = util_averageVectorsAtGraphVertex(features_hard[nextIndex]);
				result = (v0 + v1) / 2.0;
				result.normalize();
				pt = (features_hard[index].getPosition() + features_hard[nextIndex].getPosition()) / 2.0;
				//skip the next one
				index++;
			}
			else
			{
				result = util_averageVectorsAtGraphVertex(features_hard[index]);
				pt = features_hard[index].getPosition();
			}

			result *= length;
			zPoint p0 = pt + result;
			zPoint p1 = pt - result;
			posittions.push_back(p0);
			eConnect.push_back(posittions.size() - 1);
			posittions.push_back(p1);
			eConnect.push_back(posittions.size() - 1);

		}
		zFnGraph fng(outGraph_hardFeature);
		fng.create(posittions, eConnect);


		posittions.clear();
		eConnect.clear();
		for (int index = 0; index < features_soft.size(); index++)
		{
			zVector v0, v1, result;
			zPoint pt;

			int nextIndex = (index + 1) % features_soft.size();
			//check the distance to the next index, if it is within a threshold, get middle vertex as feature and skip the next feature
			float dist = features_soft[index].getPosition().distanceTo(features_soft[nextIndex].getPosition());
			if (dist <= lengthTolerance)
			{
				//get the two vectors, average them twice. They might be sharing the same vector, but not necessarily if the have points in between
				v0 = util_averageVectorsAtGraphVertex(features_soft[index]);
				v1 = util_averageVectorsAtGraphVertex(features_soft[nextIndex]);
				result = (v0 + v1) / 2.0;
				result.normalize();
				pt = (features_soft[index].getPosition() + features_soft[nextIndex].getPosition()) / 2.0;
				//skip the next one
				index++;
			}
			else
			{
				result = util_averageVectorsAtGraphVertex(features_soft[index]);
				pt = features_soft[index].getPosition();
			}

			result *= length;
			zPoint p0 = pt + result;
			zPoint p1 = pt - result;
			posittions.push_back(p0);
			eConnect.push_back(posittions.size() - 1);
			posittions.push_back(p1);
			eConnect.push_back(posittions.size() - 1);

		}
		 fng = zFnGraph(outGraph_softFeature);
		fng.create(posittions, eConnect);


		//for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
		//{
		//	bool cornerByTopology = v.getColor() == _colorCornersStart ||
		//		v.getColor() == _colorCornersEnd ||
		//		v.getColor() == _colorCornersInner ||
		//		v.getColor() == _colorCornersOuter;
		//	bool otherFeature = v.getColor() == _colorFeatureInner ||
		//		v.getColor() == _colorFeatureOuter;
		//	zItGraphHalfEdgeArray hes;
		//	v.getConnectedHalfEdges(hes);
		//	zVector v0 = hes[0].getVector();
		//	zVector v1 = hes[1].getVector();
		//	v0.normalize();
		//	v1.normalize();
		//	float angle = v0.angle(v1 * -1);
		//	bool angleChk = angle >= angleThreshold;
		//	if (cornerByTopology || angleChk)
		//	{
		//		zVector result = (v0 + v1) / 2.0;
		//		result.normalize();
		//		result *= length;
		//		zPoint p0 = v.getPosition() + result;
		//		zPoint p1 = v.getPosition() - result;
		//		posittions.push_back(p0);
		//		eConnect.push_back(posittions.size() - 1);
		//		posittions.push_back(p1);
		//		eConnect.push_back(posittions.size() - 1);
		//	}
		//	
		//}
		//zFnGraph fng(o_trimGraphs[graphId]);
		//fng.create(posittions, eConnect);
		////printf("\n Trim1 nV nE - %i | %i", fng.numVertices(), fng.numEdges());


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_SlotSide(zObjGraph& sectionGraph, zObjGraph& outGraph)
	{
		zFnGraph inFnGraph(sectionGraph);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = (blockType != zBlockType::Arch) ? _col_in_corner_st : _col_in_corner;
		zColor endColor = (blockType != zBlockType::Arch) ? _col_out_corner_st : _col_out_corner;
		for (zItGraphVertex v(sectionGraph); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV);
		gPts.push_back(endV);
		gEdges.push_back(0);
		gEdges.push_back(1);
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_SlotSide(int graphId, zObjGraph& outGraph)
	{
		compute_TrimGraphs_SlotSide(o_sectionGraphs[graphId], outGraph);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_BracingWall(zObjGraph& sectionGraph, zObjGraph& outGraph)
	{
		zFnGraph fnSectionG(sectionGraph);
		
		//walk on the graph till u reach a feature vertex and store it

		zItGraphVertexArray innerVertx, outerVertx;
		

		util_innerOuter(sectionGraph,false , innerVertx, outerVertx);

		//make graphs between verticies
		zPointArray gPositions;
		zIntArray gEdgeCOnnects;
		zColorArray gColors;
		for (int i = 0; i < innerVertx.size(); i++)
		{
			gPositions.push_back(innerVertx[i].getPosition());
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			gColors.push_back(innerVertx[i].getColor());

			gPositions.push_back(outerVertx[i].getPosition());
			gEdgeCOnnects.push_back(gPositions.size() - 1);
			gColors.push_back(outerVertx[i].getColor());

		}
		zFnGraph fnG(outGraph);
		fnG.create(gPositions, gEdgeCOnnects);
		fnG.setVertexColors(gColors);
		//printf("\n graph[%i] %i | %i", graphId, gPositions.size(), gEdgeCOnnects.size());

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_TrimGraphs_BracingWall(int graphId, zObjGraph& outGraph)
	{
		compute_TrimGraphs_BracingWall(o_sectionGraphs[graphId], outGraph);
	}

	
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_innerOuter(zObjGraph& sectionGraph, bool addEndStart, zItGraphVertexArray& innerVertx, zItGraphVertexArray& outerVertx)
	{
		innerVertx.clear();
		outerVertx.clear();
		zFnGraph fnSectionG(sectionGraph);
		//get inner and outer edges 
		bool found = false;
		zItGraphHalfEdgeArray hesInner, hesOuter;
		found = util_getShortestHEsBetweenColors(sectionGraph, _col_in_corner_st, _col_in_corner, hesInner);
		found = util_getShortestHEsBetweenColors(sectionGraph, _col_out_corner_st, _col_out_corner, hesOuter);

		//walk on the graph till u reach a feature vertex and store it

		zItGraphHalfEdge heStart, heEnd, heTemp;

		heStart = hesInner[0];
		heEnd = hesInner[hesInner.size() - 1];
		heTemp = heStart;
		if (addEndStart) innerVertx.push_back(heTemp.getStartVertex());
		while (heTemp != heEnd)
		{
			if (heTemp.getStartVertex().getColor() == _col_in_feature)
				innerVertx.push_back(heTemp.getStartVertex());
			heTemp = heTemp.getNext();
		}
		if (addEndStart) innerVertx.push_back(heEnd.getVertex());

		heStart = hesOuter[0];
		heEnd = hesOuter[hesOuter.size() - 1];
		heTemp = heStart;
		if (addEndStart) outerVertx.push_back(heTemp.getStartVertex());
		while (heTemp != heEnd)
		{
			if (heTemp.getStartVertex().getColor() == _col_out_feature)
				outerVertx.push_back(heTemp.getStartVertex());
			heTemp = heTemp.getNext();
		}
		if (addEndStart) outerVertx.push_back(heEnd.getVertex());

		//check if the two arrays are not the same size, return
		if (innerVertx.size() != outerVertx.size())
		{
			printf("\n ERROR!  inner and outer vertices are not the same size! inner | outer  %i | %i", innerVertx.size(), outerVertx.size());
			return;
		}


	}





	//SDF MAIN method
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_SDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth)
	{

		o_contourGraphs.clear();
		o_contourGraphs.assign(o_sectionGraphs.size(), zObjGraph());
		o_contourGraphs_flatten.clear();
		o_contourGraphs_flatten.assign(o_sectionGraphs.size(), zObjGraph());

		/*o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());*/

		o_raftGraphs.clear();
		o_raftGraphs.assign(1, zObjGraph());

		printf("\n num frames : %i ", o_sectionGraphs.size());


		int r0 = 0;
		int r1 = floor(o_sectionGraphs.size() * 0.5) - 1;
		int r2 = floor(o_sectionGraphs.size() * 0.5);
		int r3 = (o_sectionGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		//bool blockPlanarCheck = (leftPlaneExists && rightPlaneExists) ? true : false;
		bool blockPlanarCheck = planarBlock;// && !isRegular;

		int end = (isRegular) ? o_sectionGraphs.size() : floor(o_sectionGraphs.size() * 0.5);
		//int end =  floor(o_sectionGraphs.size() * 0.5);
		numSDFlayers = (numSDFlayers > end) ? end : numSDFlayers;
		numSDFlayers = (allSDFLayers) ? end : numSDFlayers;

		for (int j = 1; j < numSDFlayers; j++)
		{
			if (planarBlock)
			{
				if (isRegular)
				{
					if (blockType == zBlockType::Wall)
					{

						compute_BlockSDF_Planar_wall(funcNum, numSmooth, j, (j % 2 == 0), printWidth);

					}

					else
					{
						compute_BlockSDF_Planar_regular(funcNum, numSmooth, j, (j % 2 == 0), printWidth);
					}
				}
				else
				{
					compute_BlockSDF_Planar_pentagon(funcNum, numSmooth, j, (j % 2 == 0), printWidth);
					compute_BlockSDF_Planar_pentagon(funcNum, numSmooth, j + end, (j % 2 == 0), printWidth);
				}
			}
			else
			{
				compute_BlockSDF_NonPlanar(funcNum, numSmooth, j, (j % 2 == 0));

			}
			

		}



	}

	//SDF sub methods
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_BlockSDF_Planar_regular(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth)
	{
		zPrintParamSDF _printParameters;


		//printf("\n 0 fREP graphID %i  o_sectionGraphs.size() %i", graphId, o_sectionGraphs.size());

		if (graphId >= o_sectionGraphs.size())return;


		printf("\n fREP graphID %i | funcNum %i  ", graphId, funcNum);

		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		float pWidth =  printWidth;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];

		transformAllGraphs_planar(graphId, true);

		
		// field
		zFnMeshScalarField fnField(o_field);


		
		//Polygon and offset
		zScalarArray polyField, scalar_offset_outer, scalar_offset_inner;
		getScalars_offset(o_sectionGraphs[graphId], numSmooth, polyField, scalar_offset_outer, scalar_offset_inner);


		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zObjGraph slotGraph, splitGraph;
		float graphLength =   (_printParameters.offset_1st_interior + _printParameters.offset_2nd_interior + _printParameters.printWidthInterior) ;
		util_computeSlotGraph(planeXY, o_sectionGraphs[graphId], graphLength,graphId %2 == 0 , slotGraph);
		//splitGraph_1(planeXY, o_sectionGraphs[graphId], _printParameters.offset_2nd_exterior, pWidth*1.5, splitGraph);
		util_computeSplitGraph_xy(o_sectionGraphs[graphId], splitGraph);

		o_trimGraphs_SlotSide[graphId] = splitGraph;
		zScalarArray scalar_slot1;
		if (funcNum >= 2)
		{
			fnField.getScalarsAsEdgeDistance(scalar_slot1, slotGraph, _printParameters.slotStartWidth, false);
			//fnField.getScalarsAsEdgeDistance(scalar_slot1, o_trimGraphs[graphId], bracingEdgeWidth, false);
		}

		zScalarArray scalar_interiorBracing;

		zScalarArray scalar_cable;
		zScalarArray scalar_cableBracing;
		zScalarArray scalar_cableBracingSlots;
		zScalarArray scalar_cableBoolean;


		if (funcNum >= 3)
		{
			
			getScalars_3dp_cable_bracing(o_sectionGraphs[graphId], o_trimGraphs_bracing[graphId], graphId % 2 == 0, scalar_cableBracingSlots,
				scalar_cableBracing, scalar_interiorBracing);
		}
		
		//fnField.smoothField(scalar_offset_inner, numSmooth);
		//fnField.smoothField(scalar_offset_outer, numSmooth);
		


		zScalarArray booleanField_0;
		if (funcNum >= 4) fnField.boolean_subtract(scalar_offset_inner, scalar_interiorBracing, booleanField_0, false);

		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(scalar_offset_outer, booleanField_0, booleanField_1, false);


		zScalarArray scalar_booleanSlot;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_1, scalar_slot1, scalar_booleanSlot, false);

		//zScalarArray booleanField_2;
		//if (funcNum >= 6) fnField.boolean_union(booleanField_1, patternField, booleanField_2, false);


		float sdfWidth = printWidth / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(scalar_offset_outer, zFieldSDF, sdfWidth);
			break;

		case 2:
			fnField.setFieldValues(scalar_slot1, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(scalar_interiorBracing, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(booleanField_0, zFieldSDF, sdfWidth);
			break;

		case 5:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, sdfWidth);
			break;

		case 7:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, printWidth / 2.0);
			break;
		}

		/*for (zItMeshScalarField f(o_field); !f.end(); f++)
		{
			cout << f.getValue() << endl;
		}*/
		
		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		int pres = 3;
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, zVector(0,0,1), pres, 0.001);
		zObjMesh o_meshTemp;
		/*fnField.getIsolineMesh(o_meshTemp, 0.0, false);
		zFnMesh fnTemp(o_meshTemp);
		printf("\n o_meshTemp : nV - nE %i - %i ", fnTemp.numVertices(), fnTemp.numEdges());*/

		/*zItMeshHalfEdge heStart, heWalk;
		for(zItMeshHalfEdge e(o_meshTemp); !e.end(); e++)
		{
			if (e.onBoundary())
			{
				heStart = e;
				break;
			}
		}


		zPointArray tempPos;
		zIntArray tempEConnect;
		tempPos.push_back(heStart.getStartVertex().getPosition());
		heWalk = heStart;
		int safetyCounter = 0;
		do
		{
			tempPos.push_back(heStart.getVertex().getPosition());
			tempEConnect.push_back(tempPos.size() - 2);
			tempEConnect.push_back(tempPos.size() - 1);
			heWalk = heWalk.getNext();
		} while (heWalk != heStart || counter++ < 10000);
		
		tempEConnect.push_back(tempPos.size() - 1);
		tempEConnect.push_back(0);

		fnIsoGraph.create(tempPos, tempEConnect, zVector(0, 0, 1));*/

		zItGraphVertexArray vArray;

		for (zItGraphVertex v(o_contourGraphs[graphId]); !v.end(); v++)
		{
			if (!v.checkValency(2))
			{
				vArray.push_back(v);
			}
		}
		if (vArray.size()>0)
		{
			printf("\n [%i] - valence 2 verts  %i \n \n ", graphId, vArray.size());

		}

		cleanContourGraph(graphId);

		zFnGraph fngraph(o_contourGraphs[graphId]);


		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);


		// transform back 

		fnIsoGraph.setTransform(t, true, true);
		transformAllGraphs_planar(graphId, false);

		/*fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		fng = zFnGraph(o_trimGraphs[graphId]);
		fng.setTransform(t, true, true);
		fng = zFnGraph(o_trimGraphs_bracing[graphId]);
		fng.setTransform(t, true, true);
		fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
		fng.setTransform(t, true, true);
		fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
		fng.setTransform(t, true, true);
		fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
		fng.setTransform(t, true, true);*/
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_BlockSDF_Planar_pentagon(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth)
	{

		if (graphId >= o_sectionGraphs.size())return;
		printf("\n fREP graphID %i | funcNum %i  ", graphId, funcNum);

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		zFnMeshScalarField fnField(o_field);


		float pWidth = printWidth;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		transformAllGraphs_planar(graphId, true);



		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));

		//Polygon and offset
		zScalarArray polyField, scalar_offset_outer, scalar_offset_inner;
		getScalars_offset(o_sectionGraphs[graphId], numSmooth, polyField, scalar_offset_outer, scalar_offset_inner);

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zObjGraph slotGraph, splitGraph;
		float graphLength = (_printParameters.offset_1st_interior + _printParameters.offset_2nd_interior + _printParameters.printWidthInterior);
		util_computeSlotGraph(planeXY, o_sectionGraphs[graphId], graphLength, graphId % 2 == 0, slotGraph);
		//splitGraph_1(planeXY, o_sectionGraphs[graphId], _printParameters.offset_2nd_exterior, pWidth * 1.5, splitGraph);
		util_computeSplitGraph_xy(o_sectionGraphs[graphId], splitGraph);
		o_trimGraphs_SlotSide[graphId] = splitGraph;
		zScalarArray scalar_slot1;
		if (funcNum >= 2)
		{
			fnField.getScalarsAsEdgeDistance(scalar_slot1, slotGraph, _printParameters.slotStartWidth, false);
			//fnField.getScalarsAsEdgeDistance(scalar_slot1, o_trimGraphs[graphId], bracingEdgeWidth, false);
		}

		zScalarArray scalar_interiorBracing;

		zScalarArray scalar_cable;
		zScalarArray scalar_cableBracing;
		zScalarArray scalar_cableBracingSlots;
		zScalarArray scalar_cableBoolean;


		if (funcNum >= 3)
		{
			getScalars_3dp_cable_bracing(o_sectionGraphs[graphId], o_trimGraphs_bracing[graphId], graphId % 2 == 0, scalar_cableBracingSlots,
				scalar_cableBracing, scalar_interiorBracing);
			
		}

		//fnField.smoothField(scalar_offset_inner, numSmooth);
		//fnField.smoothField(scalar_offset_outer, numSmooth);



		zScalarArray booleanField_0;
		if (funcNum >= 4) fnField.boolean_subtract(scalar_offset_inner, scalar_interiorBracing, booleanField_0, false);



		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(scalar_offset_outer, booleanField_0, booleanField_1, false);


		zScalarArray scalar_booleanSlot;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_1, scalar_slot1, scalar_booleanSlot, false);

		//zScalarArray booleanField_2;
		//if (funcNum >= 6) fnField.boolean_union(booleanField_1, patternField, booleanField_2, false);


		float sdfWidth = printWidth / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(scalar_offset_outer, zFieldSDF, sdfWidth);
			break;

		case 2:
			fnField.setFieldValues(scalar_slot1, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(scalar_interiorBracing, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(booleanField_0, zFieldSDF, sdfWidth);
			break;

		case 5:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, sdfWidth);
			break;

		case 7:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, printWidth / 2.0);
			break;
		}

		/*for (zItMeshScalarField f(o_field); !f.end(); f++)
		{
			cout << f.getValue() << endl;
		}*/

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		int pres = 3;
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, zVector(0, 0, 1), pres, 0.001);
		zObjMesh o_meshTemp;

		zItGraphVertexArray vArray;

		for (zItGraphVertex v(o_contourGraphs[graphId]); !v.end(); v++)
		{
			if (!v.checkValency(2))
			{
				vArray.push_back(v);
			}
		}
		if (vArray.size() > 0)
		{
			printf("\n [%i] - valence 2 verts  %i \n \n ", graphId, vArray.size());

		}


		//cleanContourGraph(o_contourGraphs[graphId]);

		zFnGraph fngraph(o_contourGraphs[graphId]);


		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);


		// transform back 

		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		transformAllGraphs_planar(graphId, false);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_BlockSDF_Planar_wall(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth)
	{

		if (graphId >= o_sectionGraphs.size())return;
		printf("\n fREP graphID %i | funcNum %i  ", graphId, funcNum);

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		zFnMeshScalarField fnField(o_field);


		float pWidth = printWidth;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		transformAllGraphs_planar(graphId, true);



		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));

		//Polygon and offset
		zScalarArray polyField, scalar_offset_outer, scalar_offset_inner;
		getScalars_offset(o_sectionGraphs[graphId], numSmooth, polyField, scalar_offset_outer, scalar_offset_inner);



		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		//create the slot for the start
		zObjGraph slotGraph, splitGraph;
		float graphLength = 0.05;
		util_computeSlotGraph(planeXY, o_sectionGraphs[graphId], graphLength, graphId % 2 == 0, slotGraph);
		util_computeSplitGraph_xy(o_sectionGraphs[graphId], splitGraph);
		o_trimGraphs_SlotSide[graphId] = splitGraph;

		zScalarArray scalar_slot1;
		if (funcNum >= 2)
		{
			fnField.getScalarsAsEdgeDistance(scalar_slot1, slotGraph, _printParameters.slotStartWidth, false);
		}

		zScalarArray scalar_interiorBracing, scalar_bracing, scalar_bracingSlots;


		if (funcNum >= 3)
		{
			getScalars_3dp_wall_bracing(o_sectionGraphs[graphId], o_trimGraphs_bracing[graphId], _printParameters.slotIterating, graphId % 2 == 0, scalar_interiorBracing, scalar_bracing, scalar_bracingSlots);
		}


		zScalarArray scalar_triangles;
		zScalarArray scalar_boolean_trianglesInner;
		zScalarArray booleanField_0;



		if (funcNum >= 4)
		{
			getScalars_3dp_wall_triangles(o_sectionGraphs[graphId], scalar_triangles);
			fnField.boolean_subtract(scalar_offset_inner, scalar_triangles, scalar_boolean_trianglesInner, false);
			if (numSmooth > 0) fnField.smoothField(scalar_boolean_trianglesInner, numSmooth); // smooth field
			fnField.boolean_subtract(scalar_boolean_trianglesInner, scalar_interiorBracing, booleanField_0, false);


		}

		zScalarArray booleanField_1;
		zScalarArray scalar_booleanSlot;

		if (funcNum >= 5)
		{
			fnField.boolean_subtract(scalar_offset_outer, booleanField_0, booleanField_1, false);
			fnField.boolean_subtract(booleanField_1, scalar_slot1, scalar_booleanSlot, false);
		}

		float sdfWidth = _printParameters.printWidthInterior / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(scalar_offset_outer, zFieldSDF, sdfWidth);
			break;

		case 2:
			fnField.setFieldValues(scalar_slot1, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(scalar_bracing, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(scalar_bracingSlots, zFieldSDF, sdfWidth);
			break;
		case 5:
			fnField.setFieldValues(booleanField_0, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 7:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, printWidth / 2.0);
			break;


		case 8:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, printWidth / 2.0);
			break;
		}

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		int pres = 3;
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, zVector(0, 0, 1), pres, 0.001);
		cleanContourGraph(graphId);
		
		zFnGraph fngraph(o_contourGraphs[graphId]);


		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);

		fnIsoGraph.setTransform(t, true, true);
		transformAllGraphs_planar(graphId, false);




	}	
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_BlockSDF_NonPlanar(int funcNum, int numSmooth, int graphId, bool alternate )
	{
		zPrintParamSDF _printParameters;

		if (graphId >= o_sectionGraphs.size())return;

		printf("\n fREP graphID %i | funcNum %i  ", graphId, funcNum);


		zFnGraph fnTmpGraph(o_sectionGraphs[graphId]);
		zPointArray positions;
		fnTmpGraph.getVertexPositions(positions);

		zPoint minBB, maxxBB;
		fnTmpGraph.getBounds(minBB, maxxBB);

		zPoint refPt = minBB;/*positions[0];*/ // update to medial point

		zObjMesh oUnrolledMesh;// = o_sectionMeshes[graphId];

		
		// Unroll
		int num_BSF = 0;

		zInt2DArray oriVertex_UnrollVertex_map;
		unordered_map<zIntPair, int, zPair_hash> oriFaceVertex_UnrollVertex;
		zItGraphVertexArray bsf_Vertices;
		zIntPairArray bsf_vertexPairs;

		zObjGraph oDualGraph;
		zObjMesh o_projectionMesh = o_sectionMeshes[graphId];
		//getPokeMesh(o_sectionMeshes[graphId], o_projectionMesh);
		creatUnrollMesh(o_projectionMesh, oUnrolledMesh, oDualGraph, oriVertex_UnrollVertex_map, oriFaceVertex_UnrollVertex, bsf_Vertices, bsf_vertexPairs);
		num_BSF = bsf_Vertices.size();
		zVector unitZ = zVector(0, 0, 1);
		zTransform newFrame;// = core.getPlaneFromOrigin_Normal(refPt, unitZ);
		unrollMesh(o_projectionMesh, oUnrolledMesh, oDualGraph, oriVertex_UnrollVertex_map, oriFaceVertex_UnrollVertex, bsf_vertexPairs, newFrame);
		mergeMesh(oUnrolledMesh);

		zObjGraph oFlatGraph;// = o_sectionGraphs[graphId];
		//setPtGraph(oFlatGraph, refPt, false, false, true);

		createBoundaryEdgeGraph(oUnrolledMesh, true, oFlatGraph);



		//sectionFrames[graphId] = newFrame;
		zTransform t = newFrame;

		// Transform
		zTransform tLocal;
		tLocal.setIdentity();

		//transform all graphs to local
		zFnGraph fng;
		/*fng = zFnGraph(oFlatGraph);
		fng.setTransform(t, true, false);
		fng.setTransform(tLocal, true, true);*/

		//make flatten version of the trimGraph
		

		//transformAllGraphs(graphId,t, true);
		zPoint refPt2(0, 0, 0);
		zObjGraph trimGraphs_bracing_flat;

		compute_TrimGraphs_BracingWall(oFlatGraph, trimGraphs_bracing_flat);
		o_trimGraphs_bracing_flat[graphId] = trimGraphs_bracing_flat;
		
		zObjGraph o_trimGraphs_slotSide_flat;
		compute_TrimGraphs_SlotSide(oFlatGraph, o_trimGraphs_slotSide_flat);





		zFnMeshScalarField fnField(o_field);




		


		zFnGraph fnGraph(oFlatGraph);
		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));
		
		//Polygon and offset
		zScalarArray polyField, scalar_offset_outer, scalar_offset_inner;
		getScalars_offset(oFlatGraph, numSmooth, polyField, scalar_offset_outer, scalar_offset_inner);


		
		
		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;
		zPoint oo(0, 0, 0);
		planeXY = core.getTransformFromOrigin_Normal(oo, zAxis);

		zObjGraph slotGraph, splitGraph;
		float graphLength = 0.05;
		util_computeSlotGraph(planeXY, oFlatGraph, graphLength, graphId % 2 == 0, slotGraph);
		//splitGraph_1(planeXY, oFlatGraph, _printParameters.offset_2nd_exterior, pWidth * 1.5, splitGraph);
		util_computeSplitGraph_xy(oFlatGraph, splitGraph);
		
		zScalarArray scalar_slot1;
		if (funcNum >= 2)
		{
			fnField.getScalarsAsEdgeDistance(scalar_slot1, slotGraph, _printParameters.slotBracingWidth , false);
		}

		zScalarArray scalar_interiorBracing;
		zScalarArray scalar_bracing;
		zScalarArray scalar_bracingSlots;


		if (funcNum >= 3)
		{
			getScalars_3dp_wall_bracing(oFlatGraph, trimGraphs_bracing_flat, _printParameters.slotIterating, graphId % 2 == 0, scalar_interiorBracing, scalar_bracing, scalar_bracingSlots);

		}


		zScalarArray scalar_triangles;
		zScalarArray scalar_boolean_trianglesInner;
		zScalarArray booleanField_0;

		if (funcNum >= 4)
		{
			getScalars_3dp_wall_triangles(oFlatGraph, scalar_triangles);
			fnField.boolean_subtract(scalar_offset_inner, scalar_triangles, scalar_boolean_trianglesInner, false);
			if (numSmooth > 0) fnField.smoothField(scalar_boolean_trianglesInner, numSmooth); // smooth field
			fnField.boolean_subtract(scalar_boolean_trianglesInner, scalar_interiorBracing, booleanField_0, false);
		}



		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(scalar_offset_outer, booleanField_0, booleanField_1, false);


		zScalarArray scalar_booleanSlot;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_1, scalar_slot1, scalar_booleanSlot, false);

		float sdfWidth = _printParameters.printWidthInterior / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(scalar_offset_inner, zFieldSDF, sdfWidth);
			break;

		case 2:
			fnField.setFieldValues(scalar_slot1, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(scalar_interiorBracing, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(scalar_triangles, zFieldSDF, sdfWidth);
			break;
		case 5:
			fnField.setFieldValues(scalar_boolean_trianglesInner, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 7:
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, sdfWidth);
			break;


		case 8:
			if (numSmooth > 0) fnField.smoothField(scalar_booleanSlot, numSmooth); // smooth field
			fnField.setFieldValues(scalar_booleanSlot, zFieldSDF, _printParameters.printWidthInterior / 2.0);
			break;
		}

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		int pres = 3;
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, zVector(0, 0, 1), pres, 0.001);
		

		zFnGraph fngraph(o_contourGraphs[graphId]);
		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);
		

		// project to  section Mesh
		zFnGraph fnContour(o_contourGraphs[graphId]);
		zPointArray contourPositions, projectedPositions;
		zIntArray faceIDs;

		zVectorArray pNorms, pNormsTemp;

		//project contour back to section mesh
		barycentericProjection_triMesh(o_contourGraphs[graphId], oUnrolledMesh, o_projectionMesh, pNorms);

		//project splitGraph back to section mesh
		barycentericProjection_triMesh(splitGraph, oUnrolledMesh, o_projectionMesh, pNormsTemp);
		o_trimGraphs_SlotSide[graphId] = splitGraph;


		//fnContour.setVertexPositions(projectedPositions);
		fnContour.getVertexPositions(contourPositions);

		if (graphId > 0)
		{
			//zFloatArray pHeights;
			//getPrintHeight(contourPositions, pNorms, o_sectionMeshes[graphId - 1], pHeights, o_contourHeightLines[graphId]);

		}

		fnContour.setEdgeColor(zBLUE);
		fnContour.setEdgeWeight(3);

		o_sectionMeshesPar[graphId] = oUnrolledMesh;
		o_contourGraphs_flatten[graphId] = oFlatGraph;
		
	}

	//EXPORT MAIN method
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_update(string pathCurrent, string dir)
	{
		string folderName = dir + "/" + to_string(blockId);

		//EXPORT blockMesh
		exportJSON_sliceMesh(pathCurrent, folderName);

		// EXPORT sections, trims, contours
		int end = isRegular? o_sectionGraphs.size() : floor(o_sectionGraphs.size() * 0.5);
		for (int i = 0; i < end; i++)
		{
			exportJSON_graphID(folderName, i, true);
			if(!isRegular) exportJSON_graphID(folderName, i, false);
		}

		return true;
	}
	//EXPORT sub method
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_sliceMesh(string pathCurrent, string folderName)
	{
		zFnMesh fn;
		json j;
		string blockPath = pathCurrent + "blockMesh_" + to_string(blockId) + ".json";
		bool fileChk = fn.json_read(blockPath, j);
		if (!fileChk) return false;
		_mkdir(folderName.c_str());
		for (const auto& entry : std::filesystem::directory_iterator(folderName)) std::filesystem::remove_all(entry.path());


		// EXPORT BlockMesh
		//string blockID_padded = coreUtils.getPaddedIndexString(blockId, 3);
		string blockMeshName = folderName + "/block" + ".json";
		ofstream myfile;
		myfile.open(blockMeshName.c_str());
		if (myfile.fail())
		{
			cout << " error in opening outPath  " << blockMeshName.c_str() << endl;
			return false;
		}
		myfile << j.dump();
		myfile.close();

		// Export Left Right Mesh
		string leftMeshName = folderName + "/block_left" + ".json";
		string rightMeshName = folderName + "/block_right" + ".json";

		zFnMesh fnMeshLeft(o_SliceMesh_Left);
		fnMeshLeft.to(leftMeshName, zJSON);

		zFnMesh fnMeshRight(o_SliceMesh_Right);
		fnMeshRight.to(rightMeshName, zJSON);
		return true;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_graphID(string folderName, int graphId, bool left)
	{
		//string blockID_padded = coreUtils.getPaddedIndexString(blockId, 3);

		int gID = left? graphId : graphId + floor(o_sectionGraphs.size() * 0.5);
		string graphID_padded = core.getPaddedIndexString(graphId, 3);
		int i = gID;
		int extraID = left? 0 : 1;

		//string extName = blockID_padded + "_" + graphID_padded + "_" + to_string(extraID) + ".json";
		string extName = graphID_padded + "_" + to_string(extraID) + ".json";

		//section graph export
		exportJSON_graphID_section(folderName, extName, gID);

		// trim graph export
		exportJSON_graphID_trims(folderName, extName, gID);
		
		// graph export
		exportJSON_graphID_contours(folderName, extName, gID);

		return true;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_graphID_trims(string folderName, string extName, int graphId)
	{
		zFnGraph fnTrimGraph;
		string outName;

		fnTrimGraph = zFnGraph(o_trimGraphs_bracing[graphId]);
		outName = folderName + "/trim_bracing_" + extName;
		fnTrimGraph.to(outName, zJSON);

		fnTrimGraph = zFnGraph(o_trimGraphs_features_hard[graphId]);
		outName = folderName + "/trim_features_hard_" + extName;
		fnTrimGraph.to(outName, zJSON);

		fnTrimGraph = zFnGraph(o_trimGraphs_features_soft[graphId]);
		outName = folderName + "/trim_features_soft_" + extName;
		fnTrimGraph.to(outName, zJSON);

		fnTrimGraph = zFnGraph(o_trimGraphs_SlotSide[graphId]);
		outName = folderName + "/trim_interiorSplit_" + extName;
		fnTrimGraph.to(outName, zJSON);

		fnTrimGraph = zFnGraph(o_trimGraphs_seamAlignment[graphId]);
		outName = folderName + "/trim_seamAlignmentt_" + extName;
		fnTrimGraph.to(outName, zJSON);

		return true;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_graphID_contours(string folderName, string extName, int graphId)
	{
		zFnGraph fnGraph(o_contourGraphs[graphId]);
		if (fnGraph.numVertices() == 0) return false;

		zIntArray vSequence;
		zItGraphVertexArray vArray;

		for (zItGraphVertex v(o_contourGraphs[graphId]); !v.end(); v++)
		{
			if (!v.checkValency(2))
			{
				vArray.push_back(v);
			}
		}

		if (vArray.size() > 0)
		{
			printf("\n contour[%i] - valence != 2 verts %i", graphId, vArray.size());
		}

		if (vArray.size() == 2)
		{

			zItGraphHalfEdge he = vArray[0].getHalfEdge();
			vSequence.push_back(vArray[0].getId());
			int safetyCounter = 0;
			do
			{
				vSequence.push_back(he.getVertex().getId());
				he = he.getNext();
			} while (he.getVertex() != vArray[1] && safetyCounter++ < 10000);

			vSequence.push_back(vArray[1].getId());
			vSequence.push_back(vArray[0].getId());

		}

		if (vArray.size() == 0)
		{
			zItGraphHalfEdge he(o_contourGraphs[graphId], 0);
			zItGraphVertex startV = he.getStartVertex();
			zItGraphHalfEdge startHe = he;
			vSequence.push_back(he.getStartVertex().getId());

			do
			{
				if (he.getVertex() == startV)
				{
					vSequence.push_back(he.getVertex().getId());
				}
				else
				{
					vSequence.push_back(he.getVertex().getId());
				}
				he = he.getNext();
			} while (he != startHe);
		}

		json contourJson;
		fnGraph.to(contourJson);
		contourJson["VertexSequence"] = vSequence;
		string outName = folderName + "/contours_" + extName;
		/*if (graphId > 0)
		{
			zPointArray projectedPositions;
			zFnGraph fnContour(o_contourGraphs[graphId]);
			fnContour.getVertexPositions(projectedPositions);
			zFloatArray pHeights;
			zVectorArray pNorms;
			getPrintHeight(projectedPositions, pNorms, o_sectionMeshes[graphId - 1], pHeights);

		}*/
		ofstream file(outName);
		file << contourJson.dump(4); // Write formatted JSON with indentation of 4 spaces
		file.close();

		


		return true;
	}
	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON_graphID_section(string folderName, string extName, int graphId)
	{
		json sectionJson;
		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		fnGraph.to(sectionJson);
		vector<zDoubleArray> frameData;
		get2DArrayFromTransform(sectionFrames[graphId], frameData);
		sectionJson["Frame"] = frameData;

		string outName = folderName + "/section_" + extName;
		ofstream file(outName);
		file << sectionJson.dump(4); // Write formatted JSON with indentation of 4 spaces
		file.close();

		if (!planarBlock)
		{
			zFnMesh fnMesh(o_sectionMeshes[graphId]);
			string meshName = folderName + "/sectionMesh_" + extName;
			fnMesh.to(meshName, zJSON);
		}

		return true;
	}

	
	
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::check_PrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset, bool runBothPlanes, bool runPlaneLeft)
	{
		printHeightDomain = _printHeightDomain;
		neopreneOffset = _neopreneOffset;

		zStringArray files;
		core.getFilesFromDirectory(files, folderDir, zJSON);
		
		for (auto& s : files)
		{
			printf("\n All File: %s ", s.c_str());
		}

		string outFileName = folderDir + "printLayerHeights_allBlocks.csv";

		ofstream myfile;
		myfile.open(outFileName.c_str());

		if (myfile.fail())
		{
			cout << " error in opening outPath  " << outFileName.c_str() << endl;
			return;
		}

		printf("\n numFiles %i ", files.size());
		zBoolArray Block_visitied;
		Block_visitied.assign(files.size(), false);

		myfile << "blockID" << ","
			<< "blockType" << ","
			<< "frameCHECKS" << ","
			<< "minHeight" << ","
			<< "maxHeight" << ","
			<< "sdfCHECK" << ","
			<< "geometryCHECK" << ","
			<< "criticalPtsMin" << ","
			<< "criticalPtsMax" << ","
			<< "planeSpacing" << ","
			<< "dihedralAngle" << endl;


		for (auto& s : files)
		{
			//printf("\n File: %s ", s.c_str());

			zStringArray split_0 = core.splitString(s, ".");
			zStringArray split_1 = core.splitString(split_0[split_0.size() - 2], "_");

			printf("\n File: %s ", s.c_str());
			int _blockID = atoi(split_1[split_1.size() - 1].c_str());

			setFromJSON(folderDir, _blockID, runBothPlanes, runPlaneLeft);
			Block_visitied[_blockID] = true;

			//if (!planarBlock) continue;


			bool frameCHECKS = false;
			bool geomCHECKS = true;
			bool sdfCHECKS = true;

			bool maxHeight = true;
			bool minHeight = true;

			int minCriticalPtsCount = INT_MAX;
			float bestPlaneSpacing = FLT_MAX;
			if (planarBlock)
			{
				for (float printPlaneSpacing = printHeightDomain.max; printPlaneSpacing >= printHeightDomain.min; printPlaneSpacing -= 0.00025)
				{


					//printf("\n printPlaneSpace %1.4f ", printPlaneSpacing);
					compute_PrintSectionFromPlaneSpacing(printPlaneSpacing, printHeightDomain, frameCHECKS, sdfCHECKS, geomCHECKS);
					if (frameCHECKS)
					{
						bestPlaneSpacing = printPlaneSpacing;
						break;
					}

					zFnPointCloud fnCritical_min(criticalMinLayer_pts);
					zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

					if (minCriticalPtsCount > (fnCritical_min.numVertices() + fnCritical_max.numVertices()))
					{
						minCriticalPtsCount = fnCritical_min.numVertices() + fnCritical_max.numVertices();
						bestPlaneSpacing = printPlaneSpacing;
					}

					if (actualPrintHeightDomain.min < printHeightDomain.min && actualPrintHeightDomain.max < printHeightDomain.max)
					{
						printf("\n minimum print height was reached and no solution was found.");
						break;
					}

				}
			}
			else
			{
				//printf("\n non-planar printPlaneSpace %1.4f ", printPlaneSpacing);
				int numBlendFrames;
				zVector norm(0, 0, 1);
				vector<zItMeshHalfEdgeArray> vLoops;
				zObjMesh oMesh_top, oMesh_bottom;
				//vector<zObjMesh> oSectionMeshes;

				computeVLoops(o_SliceMesh_Left, medialIDS, FeaturedNumStrides, norm, vLoops, oMesh_top, oMesh_bottom);

				zScalarArray scalars;
				computeGeodesicScalars(o_SliceMesh_Left, vLoops, scalars, true);

				o_sectionMeshes.clear();
				computeGeodesicContours(vLoops, scalars, 0.01, oMesh_top, oMesh_bottom, o_sectionMeshes);;
				//computeGeodesicContours(vLoops, scalars, printPlaneSpacing, oMesh_top, oMesh_bottom, o_sectionMeshes);;
				createSectionGraphs(o_sectionMeshes, o_sectionGraphs);

				zFnPointCloud fnCritical_min(criticalMinLayer_pts);
				zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

				int minHeight = FLT_MAX;
				int maxHeight = FLT_MIN;


				for (int m = 0; m < o_sectionGraphs.size(); m++)
				{
					zPointArray secPts;
					zFnGraph fnGraph(o_sectionGraphs[m]);
					fnGraph.getVertexPositions(secPts);

					zVectorArray pNorms;
					zObjGraph outPrintHeightLines;
					for (auto& p : secPts)
					{
						for (zItMeshFace f(o_sectionMeshes[m]); !f.end(); f++)
						{
							zPointArray fVerts;
							f.getVertexPositions(fVerts);
							bool ptIntTriangle = core.pointInTriangle(p, fVerts[0], fVerts[1], fVerts[2]);
							if (ptIntTriangle)
							{

								zVector n = f.getNormal();
								pNorms.push_back(n);
								break;
							}
						}
					}
					zFloatArray pHeights;
					getPrintHeight(secPts, pNorms, o_sectionMeshes[m], pHeights, outPrintHeightLines);



					//iterate throught the print width and check for all criticcal points bellow or above maximum printDomain
					
					for (int h = 0; h < pHeights.size(); h++)
					{
						if (pHeights[h] < printHeightDomain.min)fnCritical_min.addPosition(secPts[h]);
						if (pHeights[h] > printHeightDomain.max)fnCritical_max.addPosition(secPts[h]);
					
						if (pHeights[h] < minHeight) minHeight = pHeights[h];
						if (pHeights[h] > maxHeight) maxHeight = pHeights[h];
						
					}
				}
				actualPrintHeightDomain.min = minHeight;
				actualPrintHeightDomain.max = maxHeight;

			}

			

			if (!frameCHECKS)
			{
				printf("\n Layer height check FALSE!  Choose best planeXY spacing %1.4f", bestPlaneSpacing);

				compute_PrintSectionFromPlaneSpacing(bestPlaneSpacing, printHeightDomain, frameCHECKS, sdfCHECKS, geomCHECKS);
			}



			printf("\n ----------- \n BlockID %i | %s | %1.4f %1.4f \n", _blockID, (frameCHECKS) ? "True" : "False", actualPrintHeightDomain.min, actualPrintHeightDomain.max);

			int minPts, maxPts;
			zFnPointCloud fnptCloudMin, fnptCloudMax;
			fnptCloudMin = zFnPointCloud(criticalMinLayer_pts);
			fnptCloudMax = zFnPointCloud(criticalMaxLayer_pts);

			minPts = fnptCloudMin.numVertices();
			maxPts = fnptCloudMax.numVertices();


			zFloatArray ptsMin, ptsMax;

			if (fnptCloudMin.numVertices() > 0)
			{
				zPointArray pts;
				fnptCloudMin.getVertexPositions(pts);
				for (auto& p : pts)
				{
					ptsMin.push_back(p.x);
					ptsMin.push_back(p.y);
					ptsMin.push_back(p.z);
				}
			}
			if (fnptCloudMax.numVertices() > 0)
			{
				zPointArray pts;
				fnptCloudMax.getVertexPositions(pts);
				for (auto& p : pts)
				{
					ptsMax.push_back(p.x);
					ptsMax.push_back(p.y);
					ptsMax.push_back(p.z);
				}
			}

			string path = folderDir + "blockMesh_" + to_string(_blockID) + ".json";
			json j;
			core.json_read(path, j);
			//output json for critical points
			//json j;
			string planeTypes = runningType == 0 ? "BothPlanes" : runningType == 1 ? "LeftPlanes" : "RightPlanes";
			string outDir = folderDir + "/criticalPointsCheck_" + planeTypes;
			string outPath = outDir + "/outBlock_" + to_string(_blockID) + ".json";
			if (!filesystem::is_directory(outDir) || !filesystem::exists(outDir)) filesystem::create_directory(outDir);
			j["CP_Min"] = ptsMin;
			j["CP_Max"] = ptsMax;
			j["Layer_Height_Min"] = actualPrintHeightDomain.min;
			j["Layer_Height_Max"] = actualPrintHeightDomain.max;
			j["PlaneSpacing"] = bestPlaneSpacing;

			if (ptsMin.size() > 0 || ptsMax.size() > 0)
			{

				core.json_write(outPath, j);
			}


			float leftAngle = core.dihedralAngleBetweenPlanes(leftPlanes[0], leftPlanes[1]);
			float rightAngle = core.dihedralAngleBetweenPlanes(rightPlanes[0], rightPlanes[1]);

			float maxAngle = leftAngle > rightAngle ? leftAngle : rightAngle;


			myfile << _blockID << ","
				<< ((planarBlock) ? "PlanarBlock" : "NonPlanarBlock") << ","
				<< ((frameCHECKS) ? "True" : "False") << ","
				<< actualPrintHeightDomain.min * 1000 << ","
				<< actualPrintHeightDomain.max * 1000 << ","
				<< ((sdfCHECKS) ? "True" : "False") << ","
				<< ((geomCHECKS) ? "True" : "False") << ","
				<< minPts << ","
				<< maxPts << ","
				<< bestPlaneSpacing * 1000 << ","
				<< maxAngle << endl;


			string outDir2 = folderDir + "/outSections";
			if (!filesystem::is_directory(outDir2) || !filesystem::exists(outDir2)) filesystem::create_directory(outDir2);

			outDir2 += "/block_" + to_string(_blockID);
			if (!filesystem::is_directory(outDir2) || !filesystem::exists(outDir2)) filesystem::create_directory(outDir2);
			for (const auto& entry : std::filesystem::directory_iterator(outDir2)) std::filesystem::remove_all(entry.path());

			//if (planarBlock)
			if (true)
			{
				zFnGraph fnGraph;
				for (int ss = 0; ss < o_sectionGraphs.size(); ss++)
				{
					string outPath2 = outDir2 + "/outSection_" + to_string(ss) + ".json";
					fnGraph = zFnGraph(o_sectionGraphs[ss]);
					fnGraph.to(outPath2, zJSON);
				}
			}
			printf("\n Finished block %i ", _blockID);


		}

		myfile.close();

		cout << " \n outPath exported : " << outFileName.c_str() << endl;

		for (int i = 0; i < Block_visitied.size(); i++)
		{
			if (!Block_visitied[i]) printf("\n %i ", i);
		}
	}

	
	

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::compute_cable_CableSectionPoints(int graphId, zObjGraph& o_cableGraph, zPointArray& outputPoints, float threshold )
	{
		

		zFnGraph fnSection(o_sectionGraphs[graphId]);
		zPointArray graphPts;
		fnSection.getVertexPositions(graphPts);
		zPlane frame = sectionFrames[graphId];
		zVector normal(frame(2, 0), frame(2, 1), frame(2, 2));

		zPointArray intersectionPts, intersectionPts_2;
		util_intersect_graphPlane(o_cableGraph, frame, false, intersectionPts);
		zTransformationMatrix tLocal;
		zTransformationMatrix tWorld;
		tWorld.setTransform(frame);
		zTransform t = tWorld.getToTransform(tLocal);

		zItGraphHalfEdge heS(o_sectionGraphs[graphId], 0);
		zItGraphHalfEdge heW(o_sectionGraphs[graphId], 0);
		zPointArray ptsSeq;
		do
		{
			ptsSeq.push_back(heS.getStartVertex().getPosition());
			heS = heS.getNext();
		} while (heW != heS);


		if (graphId == 16)
		{
			cout << endl << "polygons: ";
			for (auto p : ptsSeq)
			{
				cout << endl << p;
			}
			cout << endl << "intersection: ";
			for (auto p : intersectionPts)
			{
				cout << endl << p;
			}

		}

		/*for(auto& p : intersectionPts)
		{
			intersectionPts_2.push_back( p * t);
		}
		for(auto&p : graphPts)
		{
			p = p * t;
		}*/

		
		
		outputPoints.clear();
		if (intersectionPts.size() > 0)
		{
			for (int i = 0; i < intersectionPts.size(); i++)
			{
				if (core.pointInPlanarPolygon(intersectionPts[i], ptsSeq, normal))
				{
					outputPoints.push_back(intersectionPts[i]);
				}
			}
		}
		else
		{
			//printf("\n section [%i] cableGraph intersection is ZERO", graphId);
		}
		//printf("\n cableGraph [%i] intersection - number of intersection %i ",graphId, outputPoints.size());

		/*if (outputPoints.size() > 0)
		{

		}*/

		//if (intersectionEvents.Count() > 0)
		//{
		//	fnField.getScalars_Circle(scalars, intersectionPoint, radius);

		//}
	}
	ZSPACE_TOOLSETS_INLINE int zTsNatpowerSDF::compute_cable_CableGraphIndexPerGraph(int graphId)
	{
		//iterate through all the graphs to get the graph the block belongs to
		int index = -1;
		int counter = 0;
		for (int i = 0; i < o_CableGraphs.size(); i++)
		{
			zPointArray intersectionPts;
			compute_cable_CableSectionPoints(graphId, o_CableGraphs[i], intersectionPts);
			if (intersectionPts.size() > 0)
			{
				index = i;
				counter++;
			}
		}
		
		if (counter != 1)
		{

			index = -1;
			if (counter == 0)
			{

				printf("\n section [%i] cableGraph intersection is ZERO, update threshold!", graphId);
				

				//try with a different threshold
				for (int i = 0; i < o_CableGraphs.size(); i++)
				{
					zPointArray intersectionPts;
					compute_cable_CableSectionPoints(graphId, o_CableGraphs[i], intersectionPts, 0.001);
					if (intersectionPts.size() > 0)
					{
						index = i;
						counter++;
					}
				}
				if (counter == 0)
				{
					index = -1;
					printf("\n section [%i] cableGraph intersection is ZERO, FAILD!", graphId);
				}
			}
			else
			{
				printf("\n section [%i] cableGraph intersection - number of intersection %i index %i", graphId, counter, index);
			}
		}
		

		//printf("\n cableGraph intersection - number of intersection %i index %i", counter, index);

		return index;


		//// create graphs
		//zFnGraph fnGraph(o_cableGraph);

		//zPoint O(frame(3, 0), frame(3, 1), frame(3, 2));
		//zVector N(frame(2, 0), frame(2, 1), frame(2, 2));

		//zVector X(frame(0, 0), frame(0, 1), frame(0, 2));
		//X.normalize();

		//zObjNurbsCurve curve;
		//zFnNurbsCurve fnnurbs(curve);
		//fnnurbs.create(o_cableGraph, 0.01, 1, false, true, 20);

		//zObjPlane oPlane;
		//zFnPlane fnPlane(oPlane);

		//fnPlane.createFromMatrix(frame);


		//ON_SimpleArray<ON_X_EVENT> intersectionEvents;
		//zFnNurbsCurve fnCurve(curve);
		//ON_NurbsCurve* onCurve = fnCurve.getRawON_Curve();
		//int intersection_result = onCurve->IntersectPlane(oPlane.on_plane.plane_equation, intersectionEvents);

		//zPoint intersectionPoint;
		//float dist = FLT_MAX;
		//for (int i = 0; i < intersectionEvents.Count(); i++)
		//{
		//	ON_3dPoint onpt = intersectionEvents[i].m_A[0];
		//	zPoint pt = zPoint(onpt.x, onpt.y, onpt.z);
		//	float d = pt.distanceTo(O);
		//	if (d < dist)
		//	{
		//		intersectionPoint = pt;
		//		dist = d;
		//	}
		//}
		////zFnMesh fn(o_SliceMesh_Left);

		//if (intersectionEvents.Count() > 0)
		//{
		//	fnField.getScalars_Circle(scalars, intersectionPoint, radius);

		//}
	}



	ZSPACE_TOOLSETS_INLINE zPoint zTsNatpowerSDF::util_getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh)
	{
		float scaleVal = core.ofMap(threshold, thresholdLow, thresholdHigh, 0.0000f, 1.0000f);
		zVector e = vertex_higher - vertex_lower;
		double edgeLen = e.length();
		e.normalize();
		return (vertex_lower + (e * edgeLen * scaleVal));
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints)
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
			zPoint pos1 = util_getContourPosition(threshold, v1, v0, s1, s0);
			contourPoints.push_back(pos1);
		}
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints, float threshold)
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
		util_isoContour(o_graph, vertexScalars, threshold, contourPoints);
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





	//---- PROTECTED UTILITY METHODS


	ZSPACE_TOOLSETS_INLINE zItMeshHalfEdge zTsNatpowerSDF::util_getStartHalfEdge(zObjMesh& o_Mesh, int startVID, int endVID)
	{

		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		zVector dir = endPoint - startPoint;
		dir.normalize();

		zItMeshVertex vStart(o_Mesh, startVID);

		//compute start half edge
		zItMeshHalfEdge heStart;

		zItMeshHalfEdgeArray cHEdges;
		vStart.getConnectedHalfEdges(cHEdges);

		float val = 10000;
		for (zItMeshHalfEdge& he : cHEdges)
		{
			zVector heVec = he.getVector();
			heVec.normalize();

			if (1 - (heVec * dir) < val)
			{
				val = 1 - (heVec * dir);
				heStart = he;
			}
		}

		return heStart;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_createGraphFromHEArray(zItGraphHalfEdgeArray& heArray, zObjGraph& outGraph)
	{
		zPointArray positions;
		zIntArray eConnects;

		
		for (int i = 0; i < heArray.size(); i++)
		{
			int indexStart = -1;
			int indexEnd = -1;
			zPoint ptSt = heArray[i].getStartVertex().getPosition();
			zPoint ptEnd = heArray[i].getVertex().getPosition();
			bool existSt = core.checkRepeatVector(ptSt, positions, indexStart);
			bool existEnd = core.checkRepeatVector(ptEnd, positions, indexEnd);

			if (!existSt)
			{
				positions.push_back(ptSt);
				indexStart = positions.size() - 1;
			}
			if (!existEnd)
			{
				positions.push_back(ptEnd);
				indexEnd = positions.size() - 1;
			}
			eConnects.push_back(indexStart);
			eConnects.push_back(indexEnd);
		}
		//for (int i = 0; i < heArray.size(); i++)
		//{
		//	positions.push_back(heArray[i].getStartVertex().getPosition());

		//	if (positions.size() > 1)
		//	{
		//		eConnects.push_back(positions.size() - 2);
		//		eConnects.push_back(positions.size() - 1);
		//	}

		//	if (i == heArray.size() - 1)
		//	{
		//		int index = -1;
		//		zPoint pt = heArray[i].getVertex().getPosition();
		//		bool exist = core.checkRepeatVector(pt, positions, index);
		//		if (!exist)
		//		{
		//			positions.push_back(heArray[i].getVertex().getPosition());

		//			eConnects.push_back(positions.size() - 2);
		//			eConnects.push_back(positions.size() - 1);
		//		}
		//		else
		//		{
		//			eConnects.push_back(positions.size() - 1);
		//			eConnects.push_back(0);
		//		}
		//	}
		//}

		zFnGraph fnGraph(outGraph);
		fnGraph.create(positions, eConnects);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_getPerpendicularVector(zPlane& plane, zVector edgeVector, zPoint midPoint, float graphLength, zObjGraph& outGraph)
	{
		zVector planeNormal(plane(2, 0), plane(2, 1), plane(2, 2));
		//zVector vector = planeNormal ^ edgeVector;
		zVector vector = edgeVector.rotateAboutAxis(planeNormal, 90.0f);
		vector.normalize();
		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(midPoint - (vector * (graphLength / 2)));
		gPts.push_back(midPoint + (vector * (graphLength)));
		gEdges.push_back(0);
		gEdges.push_back(1);
		zObjGraph oGraph;
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);
	}

	ZSPACE_TOOLSETS_INLINE zVector zTsNatpowerSDF::util_averageVectorsAtGraphVertex(zItGraphVertex& v)
	{
		zItGraphHalfEdgeArray hes;
		v.getConnectedHalfEdges(hes);

		zVector result(0, 0, 0);
		for (zItGraphHalfEdge he : hes)
		{
			zVector vec = he.getVector();
			vec.normalize();
			result += vec;
		}
		result /= hes.size();
		result.normalize();
		return result;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_combineMultipleGraphs(zObjGraphArray& inGraphs, zObjGraph& outGraph)
	{
		zPointArray positions;
		zIntArray eConnects;
		zColorArray eColors;

		zFnGraph fng;
		for (zObjGraph& g : inGraphs)
		{
			zPointArray pts;
			zIntArray es;
			zColorArray colors;
			fng = zFnGraph(g);
			fng.getVertexPositions(pts);
			fng.getEdgeData(es);
			fng.getEdgeColors(colors);

			// store the mapping of old indices to new indices for the current graph
			zIntArray indexMapping(pts.size());

			// Add points to the combined positions array and update the map
			for (int i = 0; i < pts.size(); ++i)
			{
				const zPoint& p = pts[i];
				int tempIndex = -1;
				bool exist = core.checkRepeatVector(pts[i], positions, tempIndex);
				if (!exist) positions.push_back(p);
				indexMapping[i] = exist? tempIndex : positions.size() - 1;
			}

			// Adjust edge connectivity based on new vertex indices
			for (int i = 0; i < es.size(); i += 2)
			{
				int newStartIdx = indexMapping[es[i]];
				int newEndIdx = indexMapping[es[i + 1]];
				eConnects.push_back(newStartIdx);
				eConnects.push_back(newEndIdx);
			}

			// Combine edge colors
			eColors.insert(eColors.end(), colors.begin(), colors.end());
		}

		// Create the combined graph
		fng = zFnGraph(outGraph);
		fng.create(positions, eConnects);
		fng.setEdgeColors(eColors, false);

	}

	ZSPACE_TOOLSETS_INLINE zPoint zTsNatpowerSDF::util_getGraphPointAtParameter(zObjGraph& inGraph, float normalizedPar, int& outEdgeIndex)
	{
		zFnGraph fnGraph(inGraph);
		zDoubleArray eLengths;
		double totLength = fnGraph.getEdgeLengths(eLengths);
	
		double parLength = totLength * normalizedPar;
		double currentLength = 0;
		outEdgeIndex = -1;
		zPoint point;
		for (int i = 0; i < eLengths.size(); i++)
		{
			currentLength += eLengths[i];
			if (currentLength >= parLength)
			{
				outEdgeIndex = i;
				zItGraphEdge e(inGraph, i);
				zPointArray ePts;
				e.getVertexPositions(ePts);
				//zPoint p = e.getStartVertex().getPosition();
				zPoint p = ePts[0];
				zVector v = e.getVector();
				v.normalize();
				point = p + (v * (currentLength - parLength));
			}
		}
		return point;
	}

	ZSPACE_TOOLSETS_INLINE zPoint zTsNatpowerSDF::util_getPointAtParameterHalfEdge(zItGraphHalfEdge& he, float normalizedPar)
	{
		double length = he.getLength();
		double parLength = length * normalizedPar;
		zPoint p = he.getStartVertex().getPosition();
		zVector v = he.getVector();
		v.normalize();
		p += (v * (length - parLength));
		return p;
	}


	ZSPACE_TOOLSETS_INLINE double zTsNatpowerSDF::util_normalise(double value, double min, double max)
	{
		return (value - min) / (max - min);
	}

	ZSPACE_TOOLSETS_INLINE double zTsNatpowerSDF::util_denormalise(double value, double min, double max)
	{
		return value * (max - min) + min;
	}

	
	
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_computeSlotGraph(zPlane plane, zObjGraph& inPoly, float graphLength, bool iterate, zObjGraph& outGraph)
	{
		zPrintParamSDF _printParameters;

		//poly is sectionGraph -> we don't have to specify right/left
		zFnGraph inFnGraph(inPoly);

		zVector Y(0, 1, 0);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = (blockType != zBlockType::Arch) ? _col_in_corner_st : _col_in_corner;
		zColor endColor = (blockType != zBlockType::Arch) ? _col_out_corner_st : _col_out_corner;

		if (blockType == zBlockType::Bottom)
		{
			startColor =  _col_out_corner_st ;
			endColor =  _col_in_corner_st ;
		}

		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		/*float ptOffset = edgeLength / 2.0 ;
		if (iterate) ptOffset += graphLength;*/

		float ptOffset = _printParameters.slotStart;
		if (iterate) ptOffset += _printParameters.slotIterating;

		//if (blockType != zBlockType::Arch)
		//{
		//	ptOffset = iterate ? graphLength * 2.5 : edgeLength - (graphLength * 2.5);
		//}
		//if (blockType == zBlockType::Bottom)
		//{
		//	//ptOffset = iterate ? edgeLength - (graphLength * 2.5) : edgeLength - (graphLength * 5);
		//}
		//if (isCorner)
		//{
		//	ptOffset = iterate ? (graphLength * 2.0) : edgeLength - (graphLength * 2.0);
		//}
		//if (blockType == zBlockType::Wall)
		//{
		//	ptOffset = iterate ? edgeLength - (graphLength * 2.5) : edgeLength - (graphLength * 3);

		//}
		
		startV += (edgeVector * ptOffset);

		util_getPerpendicularVector(plane, edgeVector, startV, graphLength * 2, outGraph);



	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_computeSplitGraph_plane(zPlane plane, zObjGraph& inPoly, float offset, float trim, zObjGraph& outGraph)
	{
		zFnGraph inFnGraph(inPoly);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = (blockType != zBlockType::Arch) ? _col_in_corner_st : _col_in_corner;
		zColor endColor = (blockType != zBlockType::Arch) ? _col_out_corner_st : _col_out_corner;
		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}

		
		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		//startV += (edgeVector * trim);
		//endV -= (edgeVector * trim);

		////move the graph to the offset
		//zVector n = zVector(plane(2, 0), plane(2, 1), plane(2, 2));
		//zVector dir = edgeVector.rotateAboutAxis(n, 90.0f);

		//dir.normalize();
		//dir *= offset;
		//startV += dir;
		//endV += dir;

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV);
		gPts.push_back(endV);
		gEdges.push_back(0);
		gEdges.push_back(1);
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::util_computeSplitGraph_xy(zObjGraph& inPoly, zObjGraph& outGraph)
	{
		zPrintParamSDF _printParameters;

		zFnGraph fnGraph(inPoly);

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zVector* inPositions = fnGraph.getRawVertexPositions();
		zColor* inColors = fnGraph.getRawVertexColors();

		zPoint startV, endV;

		zVector edgeVector;

		int counter;

		zColor startColor = (blockType != zBlockType::Arch) ? _col_in_corner_st : _col_in_corner;
		zColor endColor = (blockType != zBlockType::Arch) ? _col_out_corner_st : _col_out_corner;

		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				startV = v.getPosition();
				counter++;
			}
			if (v.getColor() == endColor)
			{
				endV = v.getPosition();
				counter++;
			}
			if (counter == 2)
			{
				break;
			}
		}


		edgeVector = endV - startV;
		float edgeLength = edgeVector.length();
		edgeVector.normalize();

		//Shorten the edge
		startV += (edgeVector * _printParameters.splitTrimTimming);
		endV -= (edgeVector * _printParameters.splitTrimTimming);

		//move the edge so it is in the middle between 1st and 2nd offset
		//The vector of the offset has the following:
		////1. rotate along the normal of the plane by 90 degrees and -90 degrees (vec1 and vec2)
		////2. get a vector between the center of the edge with the center of the section graph (vecC)
		////3. find the angle between the vec1/vec2 in step1 and vecC, and choose the one with the smallest angle (vec)
		////4. move using the vector in step 3

		zVector n = zVector(planeXY(2, 0), planeXY(2, 1), planeXY(2, 2));
		zVector vec1 = edgeVector.rotateAboutAxis(n, 90.0f);
		zVector vec2 = edgeVector.rotateAboutAxis(n, -90.0f);
		
		zVector vecC = fnGraph.getCenter() - ((startV+endV)/2.0);
		vecC.normalize();

		float angle1 = abs( vec1.angle(vecC));
		float angle2 = abs(vec2.angle(vecC));

		zVector vec = (angle1 < angle2) ? vec1 : vec2;
		vec.normalize();
		vec *= _printParameters.splitTrimOffset;
		startV += vec;
		endV += vec;


		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV);
		gPts.push_back(endV);
		gEdges.push_back(0);
		gEdges.push_back(1);
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);
	}


	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::util_getShortestHEsBetweenColors(zObjGraph& graph, zColor startColor, zColor endColor, zItGraphHalfEdgeArray& outHEs)
	{
		outHEs.clear();
		vector<zItGraphHalfEdgeArray> tempHEsArray;
		zFloatArray lengths;
		zFnGraph fng(graph);
		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (v.getColor() == startColor)
			{
				//walk until you reach a vertex with color _colorFeatureOuter
				zItGraphHalfEdgeArray hes;
				v.getConnectedHalfEdges(hes);
				for (auto& he : hes)
				{
					float length = 0;

					zItGraphHalfEdge heStart = he;
					zItGraphHalfEdge he = heStart;
					zItGraphHalfEdgeArray innerHE;
					int safetyCounter = 0;
					while (safetyCounter < fng.numEdges())
					{
						innerHE.push_back(he);
						length += he.getLength();
						if (he.getVertex().getColor() == endColor) break;
						he = he.getNext();
						safetyCounter++;

					}
					/*if (safetyCounter >= fng.numEdges())
					{
						printf("\n heColor was not found");
					}*/
					tempHEsArray.push_back(innerHE);
					lengths.push_back(length);
					//printf("\n slotGraph_Arch index-size %i | %i", innerHEs.size(), innerHE.size());
					//innerLengths.push_back(length);
				}
			}
		}
		if (tempHEsArray.size() == 0)
		{
			printf("\n getShotestHEsBetweenColors no tempHEsArray found! RETURN");
			return false;
		}
		int index = 0;
		int minCount = INT_MAX;
		float minLength = FLT_MAX;
		for (int i = 0; i < tempHEsArray.size(); i++)
		{
			//if (tempHEsArray[i].size() < minCount)
			if (lengths[i] < minLength)
			{
				minCount = tempHEsArray[i].size();
				minLength = lengths[i];
				index = i;
			}
		}
		outHEs = tempHEsArray[index];
		return true;
	}

	ZSPACE_TOOLSETS_INLINE int zTsNatpowerSDF::util_getGraphClosestPoint(zObjGraph& graph, zPoint& samplePoint, zPoint& outPoint, float& dist)
	{
		int index = -1;
		double minDist = DBL_MAX;
		zPoint cP;
		for (zItGraphEdge e(graph); !e.end(); e++)
		{
			zPointArray pts;
			e.getVertexPositions(pts);
			zPoint p;
			double d = core.minDist_Edge_Point(samplePoint, pts[0], pts[1],p);
			if (d < minDist)
			{
				minDist = d;
				index = e.getId();
				cP = p;
			}
		}
		dist = (float)minDist;
		outPoint = cP;
		return index;
	}

	ZSPACE_TOOLSETS_INLINE int zTsNatpowerSDF::util_getHeArrayClosestPoint(zItGraphHalfEdgeArray& hes, zPoint& samplePoint, zPoint& outPoint, float& dist)
	{
		int index = -1;
		double minDist = DBL_MAX;
		zPoint cP;
		for (int i = 0; i < hes.size(); i++)
		{
			
			zPoint tempCp;
			zPoint p0 = hes[i].getStartVertex().getPosition();
			zPoint p1 = hes[i].getVertex().getPosition();
			double d = core.minDist_Edge_Point(samplePoint, p0,p1, tempCp);
			if (d < minDist)
			{
				minDist = d;
				index = i;
				cP = tempCp;
			}
		}
		
		dist = (float)minDist;
		outPoint = cP;
		return index;
	}


	


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_cable_bracing(zObjGraph& sectionGraph, zObjGraph& bracingGraph, bool iterateChk, zScalarArray& scalar_cableBracingSlots, zScalarArray& scalar_cableBracing, zScalarArray& scalar_interiorBracing )
	{

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(sectionGraph);
		zFnMeshScalarField fnField(o_field);

		zScalarArray scalar_cable;
		
		zScalarArray scalar_cableBoolean;

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;



		float cableRadius_offset = _printParameters.cableWidth + (_printParameters.bracingEdgeWidth*2);
		float cableRadius = _printParameters.cableWidth;
		zItGraphVertex vCenter(bracingGraph, 0);
		zPoint cableCenter = vCenter.getPosition();

		fnField.getScalars_Circle(scalar_cable, cableCenter, cableRadius, 0, false);

		fnField.getScalarsAsEdgeDistance(scalar_cableBracing, bracingGraph, _printParameters.bracingEdgeWidth, false);

		//get all edges except the one with the least length. Create new graph
		//get edge of the smallest length
		//create an array of all edges except the minimum
		//get all edges

		zItGraphEdgeArray esTemp;
		double minLength = DBL_MAX;
		int minIndex = -1;

		for (zItGraphEdge e(bracingGraph); !e.end(); e++)
		{
			if (e.getLength() < minLength)
			{
				minLength = e.getLength();
				minIndex = e.getId();
			}
			esTemp.push_back(e);
		}
		if (minIndex < esTemp.size()) esTemp.erase(esTemp.begin() + minIndex); // Erase element at index
		//if (minIndex < esTemp.size()) esTemp.erase(esTemp.begin() + 2); // Erase element at index

		//get HE array for each edge in one direction. following the edge connectivity 
		zItGraphHalfEdgeArray hesTemp, hesTemp2;
		//zBoolArray isCableEdge;
		for (zItGraphEdge e : esTemp)
		{
			zItGraphVertexArray vs;
			e.getVertices(vs);
			int ind = e.getHalfEdge(0).getStartVertex().getId() == vs[0].getId() ? 1 : 0;
			hesTemp.push_back(e.getHalfEdge(ind));
			//isCableEdge.push_back(e.getHalfEdge(ind).getStartVertex().getPosition() == vCenter.getPosition());
		}


		zObjGraph cableBracingSlots;
		zObjGraphArray cableBracingSlotsArray;
		cableBracingSlotsArray.assign(hesTemp.size(), zObjGraph());
		int counter = 0;
		for (zItGraphHalfEdge& he : hesTemp)
		{
			bool isCableEdge = he.getStartVertex().getPosition() == vCenter.getPosition() || he.getVertex().getPosition() == vCenter.getPosition();
			
			float edgeStartLength = (he.getLength() - _printParameters.offset_1st_interior - _printParameters.offset_2nd_interior);
			float cableoffset = (edgeStartLength - cableRadius_offset) / 2.0;
			float bracingOffset = (edgeStartLength - _printParameters.bracingEdgeWidth) / 2.0;
			float slotOffset = isCableEdge ? cableoffset : bracingOffset;
			//if (graphId % 2 == 0) slotOffset += (bracingEdgeWidth * 2);

			if (iterateChk) slotOffset += _printParameters.slotIterating;

			zVector vec = he.getVector();
			vec.normalize();
			vec *= slotOffset;
			zPoint midP = he.getStartVertex().getPosition() + vec;



			//getPerpendicularVector(sectionFrames[graphId], vec, midP, bracingEdgeWidth * 2, cableBracingSlotsArray[counter]);
			util_getPerpendicularVector(planeXY, vec, midP, _printParameters.bracingEdgeWidth * 2, cableBracingSlotsArray[counter]);

			counter++;

		}
		util_combineMultipleGraphs(cableBracingSlotsArray, cableBracingSlots);

		fnField.getScalarsAsEdgeDistance(scalar_cableBracingSlots, cableBracingSlots, _printParameters.slotBracingWidth, false);

		fnField.boolean_union(scalar_cable, scalar_cableBracing, scalar_cableBoolean, false);

		fnField.boolean_subtract(scalar_cableBoolean, scalar_cableBracingSlots, scalar_interiorBracing, false);


	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_wall_bracing(zObjGraph& sectionGraph, zObjGraph& bracingGraph, float iterateOffset, bool iterateChk, zScalarArray & outScalar_interiorBracing, zScalarArray & outScalar_bracing, zScalarArray & outScalar_bracingSlots)
	{
		/// This method create the bracing
		/// 1. get the scalar of the bracing graph
		/// 2. get the scalar of the bracing slots
		/// 3. subtract the bracing slots from the bracing graph
		/// 4. get the scalar of the interior bracing

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(sectionGraph);
		zFnMeshScalarField fnField(o_field);


		fnField.getScalarsAsEdgeDistance(outScalar_bracing, bracingGraph, _printParameters.bracingEdgeWidth, false);
		zFnGraph fnTrimG(bracingGraph);
		int bracingCount = fnTrimG.numEdges();

		zItGraphHalfEdgeArray hesTemp;
		fnTrimG.numHalfEdges();
		for (int i = 0; i < fnTrimG.numHalfEdges(); i++)
		{
			zItGraphHalfEdge he(bracingGraph, i);
			if (he.getStartVertex().getColor() == _col_in_feature)
			{
				hesTemp.push_back(he);
			}
		}

		zPlane planeXY;
		planeXY.setIdentity();
		zVector zAxis(0, 0, 1);
		planeXY(2, 0) = 0;
		planeXY(2, 1) = 0;
		planeXY(2, 2) = 1;

		zObjGraph o_bracingSlots;
		zObjGraphArray bracingSlotsArray;
		bracingSlotsArray.assign(hesTemp.size(), zObjGraph());
		int counter = 0;
		for (zItGraphHalfEdge& he : hesTemp)
		{
			//to get the slot offset, we have the following steps
			//1. get the length of the bracing edge
			//2. get the offset of the triangle point (using the factor value)
			//3. get the offset of the 1st offset and the 2nd offset
			//4. find the middle point between step 2 and step 3
			//5. iterate
			float triangleStart = he.getLength() * _printParameters.wall_triangleOffsetFactor;
			float exteriorStart = _printParameters.offset_1st_exterior + _printParameters.offset_2nd_exterior;
			float slotOffset = triangleStart + ((he.getLength() - exteriorStart - triangleStart) / 2);

			if (iterateChk) slotOffset += _printParameters.slotIterating;
			zVector vec = he.getVector();
			vec.normalize();
			vec *= slotOffset;
			zPoint midP = he.getStartVertex().getPosition() + vec;
			float graphLength = _printParameters.bracingEdgeWidth * 4;
			util_getPerpendicularVector(planeXY, vec, midP, graphLength, bracingSlotsArray[counter]);

			counter++;

		}
		util_combineMultipleGraphs(bracingSlotsArray, o_bracingSlots);

		fnField.getScalarsAsEdgeDistance(outScalar_bracingSlots, o_bracingSlots, _printParameters.bracingEdgeWidth, false);
		fnField.boolean_subtract(outScalar_bracing, outScalar_bracingSlots, outScalar_interiorBracing, false);
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_wall_triangles(zObjGraph& sectionGraph, zScalarArray & outScalar_triangles)
	{

		zPrintParamSDF _printParameters;
		zFnGraph fnGraph(sectionGraph);
		zFnMeshScalarField fnField(o_field);

		//inner vertex are the vertex of the sectionGraph on the inner side (_col_in_feature - magenta) 
		//and outer vertex are the vertex of the sectionGraph on the outer side (_col_out_feature - orange)
		zItGraphVertexArray innerVertx, outerVertx;
		util_innerOuter(sectionGraph, true, innerVertx, outerVertx);

		zPointArray gPts;
		zIntArray eConnect;
		for (int i = 0; i < innerVertx.size(); i++)
		{
			//this is the middle point between the inner and outer vertex - this value can be changed to any value
			zVector edgeDir = outerVertx[i].getPosition() - innerVertx[i].getPosition();
			edgeDir *= _printParameters.wall_triangleOffsetFactor;
			zPoint Pt0 = innerVertx[i].getPosition() + edgeDir;

			gPts.push_back(Pt0);
			if (i < innerVertx.size() - 1)
			{
				zPoint p1 = (innerVertx[i].getPosition() + innerVertx[i + 1].getPosition()) / 2;
				zPoint p2 = (outerVertx[i].getPosition() + outerVertx[i + 1].getPosition()) / 2;
				zVector v = p2 - p1;
				v.normalize();
				zPoint Pt1 = p1 + (v * 0.001); //< this is the peak of the triangle (average) , small offset so the full polygon is not through the same top edge
				gPts.push_back(Pt1);
			}

		}
		//create the triangles sides 
		for (int i = 0; i < gPts.size() - 1; i++)
		{
			eConnect.push_back(i);
			eConnect.push_back(i + 1);
		}
		//create the top part of the triangle - this is the full of the inner edge
		for (int i = innerVertx.size() - 1; i >= 0; i--)
		{
			eConnect.push_back(gPts.size() - 1);
			//printf("\n [%i] innerVertx %i  ", i, innerVertx.size() - 1);
			gPts.push_back(innerVertx[i].getPosition());
			eConnect.push_back(gPts.size() - 1);
		}
		eConnect.push_back(gPts.size() - 1);
		eConnect.push_back(0);
		zObjGraph o_triangles;
		zFnGraph fnG(o_triangles);
		fnG.create(gPts, eConnect);
		

		fnField.getScalars_Polygon(outScalar_triangles, o_triangles, false);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_offset(zObjGraph& sectionGraph, int numSmooth, zScalarArray& outScalar_polygon, zScalarArray& outScalar_offset_outer, zScalarArray & outScalar_offset_inner)
	{
		zPrintParamSDF _printParameters;
		zFnMeshScalarField fnField(o_field);
		zFnGraph fnGraph(sectionGraph);


		// Profile polygon field
		zIntArray edgeId;
		fnField.getScalars_Polygon(outScalar_polygon, sectionGraph, edgeId, false);

		//create a map of edge offset based on the edgeId
		zFloatArray outerOffsetArray, innerOffsetArray;
		outerOffsetArray.assign(fnGraph.numEdges(), _printParameters.offset_1st_exterior);
		innerOffsetArray.assign(fnGraph.numEdges(), _printParameters.offset_1st_exterior + _printParameters.offset_2nd_exterior);


		//color the section based on the offset color
		zItGraphHalfEdgeArray hesInterior;
		util_getShortestHEsBetweenColors(sectionGraph, _col_in_corner_st, _col_in_corner, hesInterior);

		//update the offsetArray based on the interior edges
		for (zItGraphHalfEdge he : hesInterior)
		{
			outerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior;
			innerOffsetArray[he.getEdge().getId()] = _printParameters.offset_1st_interior + _printParameters.offset_2nd_interior;
		}

		outScalar_offset_outer = outScalar_polygon;
		outScalar_offset_inner = outScalar_polygon;

		//update the scalar offset based on the edgeId
		for (int sf = 0; sf < outScalar_offset_outer.size(); sf++)
		{
			outScalar_offset_outer[sf] += outerOffsetArray[edgeId[sf]];
			outScalar_offset_inner[sf] += outerOffsetArray[edgeId[sf]] + innerOffsetArray[edgeId[sf]];
		}
		//smooth fields
		fnField.smoothField(outScalar_offset_inner, numSmooth);
		fnField.smoothField(outScalar_offset_outer, numSmooth);

	}




	//-------------------------- Non planar UTILS ------------------------------

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getPokeMesh(zObjMesh& o_mesh, zObjMesh& o_TriMesh)
	{
		zFnMesh fnMesh(o_mesh);

		zColorArray vColors;
		fnMesh.getVertexColors(vColors);

		zPointArray vPositions;
		fnMesh.getVertexPositions(vPositions);

		zPointArray fCens;
		fnMesh.getCenters(zHEData::zFaceData, fCens);

		zPointArray positions;
		zIntArray pCounts, pConnects;

		positions.insert(positions.end(), vPositions.begin(), vPositions.end());
		positions.insert(positions.end(), fCens.begin(), fCens.end());

		int numOriginalVerts = vPositions.size();

		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zIntArray fVerts;
			f.getVertices(fVerts);

			int fID = f.getId();

			for (int i = 0; i < fVerts.size(); i++)
			{
				int nextID = (i + 1) % fVerts.size();

				pConnects.push_back(fVerts[i]);
				pConnects.push_back(fVerts[nextID]);
				pConnects.push_back(numOriginalVerts + fID);

				pCounts.push_back(3);
			}

			vColors.push_back(zBLACK);
		}

		zFnMesh fnPokeMesh(o_TriMesh);
		fnPokeMesh.create(positions, pCounts, pConnects);

		fnPokeMesh.setVertexColors(vColors);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getLoop(zItMeshHalfEdge& heStart, bool forward, bool corner, int vCounter, vector<zItMeshHalfEdgeArray>& v_Loops)
	{
		zItMeshHalfEdge he_U = (forward) ? heStart.getNext() : heStart.getPrev();
		if (corner) he_U = heStart;

		bool exit_1 = false;

		zItMeshHalfEdge he_V = (forward) ? he_U.getSym().getNext() : he_U.getSym().getPrev();

		zItMeshHalfEdgeArray tempV;

		bool exit_2 = false;

		for (int i = 0; i < vCounter; i++)
		{
			if (forward) tempV.push_back(he_V.getSym());
			else tempV.push_back(he_V);

			//he_V.getEdge().setColor(zBLUE);


			if (!exit_2) he_V = (forward) ? he_V.getNext().getSym().getNext() : he_V.getPrev().getSym().getPrev();
		}


		v_Loops.push_back(tempV);

		//he_U.getEdge().setColor(zMAGENTA);



	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zPointArray& fVerts, zColorArray& fVColors)
	{
		fVerts.clear();
		fVColors.clear();

		zItMeshHalfEdge he = heStart;

		do
		{
			if (forward)
			{
				fVerts.push_back(he.getVertex().getPosition());
				fVColors.push_back(he.getVertex().getColor());
				he = he.getNext();
			}
			else
			{
				fVerts.push_back(he.getStartVertex().getPosition());
				fVColors.push_back(he.getStartVertex().getColor());
				he = he.getPrev();
			}

		} while (he != heStart);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getFaceVerticesFromHalfedge(zItMeshHalfEdge& heStart, bool forward, zIntArray& fVerts)
	{
		fVerts.clear();

		zItMeshHalfEdge he = heStart;

		do
		{
			if (forward)
			{
				fVerts.push_back(he.getVertex().getId());
				he = he.getNext();
			}
			else
			{
				fVerts.push_back(he.getStartVertex().getId());
				he = he.getPrev();
			}

		} while (he != heStart);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createBoundaryEdgeGraph(zObjMesh& o_mesh, bool closeGraph, zObjGraph& o_Graph)
	{
		zPointArray positions;
		zIntArray eConnects;
		zColorArray vColors;

		zItMeshHalfEdge he;

		for (zItMeshHalfEdge tmpHE(o_mesh); !tmpHE.end(); tmpHE++)
		{
			if (tmpHE.onBoundary())
			{
				he = tmpHE;
				break;
			}
		}

		zItMeshHalfEdge startHE = he;
		positions.push_back(he.getStartVertex().getPosition());
		vColors.push_back(he.getStartVertex().getColor());

		do
		{
			positions.push_back(he.getVertex().getPosition());
			vColors.push_back(he.getVertex().getColor());

			eConnects.push_back(positions.size() - 2);
			eConnects.push_back(positions.size() - 1);

			he = he.getNext();

		} while (he != startHE);

		if (closeGraph)
		{
			eConnects.push_back(positions.size() - 1);
			eConnects.push_back(0);
		}

		zFnGraph fnGraph(o_Graph);
		fnGraph.create(positions, eConnects);

		fnGraph.setVertexColors(vColors);
		;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::colorMesh(zObjMesh& o_mesh, zFloatArray& scalars)
	{
		zFnMesh fnMesh(o_mesh);

		// color mesh
		zColor* mesh_vColors = fnMesh.getRawVertexColors();

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		zDomainFloat distanceDomain(minScalar, maxScalar);
		zDomainColor colDomain(zColor(1, 0, 0, 1), zColor(0, 1, 0, 1));

		for (int i = 0; i < fnMesh.numVertices(); i++)
		{
			mesh_vColors[i] = core.blendColor(scalars[i], distanceDomain, colDomain, zRGB);
		}

		fnMesh.computeFaceColorfromVertexColor();

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setPtGraph(zObjGraph& o_Graph, zPoint& refPt, bool setX, bool setY, bool setZ)
	{
		zFnGraph fnGraph(o_Graph);
		zPoint* positions = fnGraph.getRawVertexPositions();

		for (int i = 0; i < fnGraph.numVertices(); i++)
		{
			if (setX) positions[i].x = refPt.x;
			if (setY) positions[i].y = refPt.y;
			if (setZ) positions[i].z = refPt.z;
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setPtMesh(zObjMesh& o_Mesh, zPoint& refPt, bool setX, bool setY, bool setZ)
	{
		zFnMesh fnMesh(o_Mesh);
		zPoint* positions = fnMesh.getRawVertexPositions();

		for (int i = 0; i < fnMesh.numVertices(); i++)
		{
			if (setX) positions[i].x = refPt.x;
			if (setY) positions[i].y = refPt.y;
			if (setZ) positions[i].z = refPt.z;
		}
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getBoundaryOffset(zObjMesh& _oMesh, bool keepExistingFaces, float offset, zObjMesh& outMesh)
	{
		zFnMesh fnMesh(_oMesh);

		vector<zVector>positions;
		vector<int>polyConnects;
		vector<int>polyCounts;

		vector<zColor> vertexColors;


		vector<zVector> fCenters;
		fnMesh.getCenters(zFaceData, fCenters);

		if (keepExistingFaces)
		{
			vector<vector<int>> inVertex_newVertex;

			for (zItMeshVertex v(_oMesh); !v.end(); v++)
			{
				vector<int> temp;

				if (v.onBoundary())
				{

					temp.push_back(positions.size());


					positions.push_back(v.getPosition());



					float vertexVal = v.getColor().r;
					zColor newCol(v.getColor().r, 0, 0, 1);

					double extrudeVal = offset;


					zItMeshHalfEdge vEdge;

					vector<zItMeshHalfEdge> cEdges;
					v.getConnectedHalfEdges(cEdges);

					for (auto& he : cEdges)
					{
						if (he.onBoundary())
						{
							vEdge = he;
						}
					}

					//if (vEdge == NULL) continue;

					zItMeshVertex next = vEdge.getVertex();;
					zItMeshVertex prev = vEdge.getPrev().getSym().getVertex();

					zVector vNorm = v.getNormal();

					zVector Ori = v.getPosition();;

					zVector v1 = Ori - prev.getPosition();
					v1.normalize();

					zVector n1 = v1 ^ vNorm;
					n1.normalize();

					zVector v2 = next.getPosition() - Ori;
					v2.normalize();

					zVector n2 = v2 ^ vNorm;
					n2.normalize();

					v1 = v1 ^ v2;
					zVector v3 = (n1 + n2);

					v3 *= 0.5;
					v3.normalize();



					double cs = v3 * v2;
					double length = extrudeVal;


					zVector a1 = v2 * cs;
					zVector a2 = v3 - a1;

					double alpha = 0;
					if (a2.length() > 0) alpha = sqrt(a2.length() * a2.length());

					if (cs < 0 && a2.length() > 0) alpha *= -1;

					if (alpha > 0) length /= alpha;

					zVector offPos = Ori + (v3 * length);

					temp.push_back(positions.size());
					positions.push_back(offPos);




				}

				inVertex_newVertex.push_back(temp);

			}

			// poly connects 
			for (zItMeshHalfEdge he(_oMesh); !he.end(); he++)
			{
				if (he.onBoundary())
				{
					vector<int> eVerts;
					he.getVertices(eVerts);

					polyConnects.push_back(inVertex_newVertex[eVerts[0]][0]);
					polyConnects.push_back(inVertex_newVertex[eVerts[1]][0]);
					polyConnects.push_back(inVertex_newVertex[eVerts[1]][1]);
					polyConnects.push_back(inVertex_newVertex[eVerts[0]][1]);

					polyCounts.push_back(4);

				}
			}

		}


		if (!keepExistingFaces)
		{

			for (zItMeshVertex v(_oMesh); !v.end(); v++)
			{
				vector<int> temp;

				if (v.onBoundary())
				{

					float vertexVal = v.getColor().r;
					zColor newCol(v.getColor().r, 0, 0, 1);

					double extrudeVal = offset;


					zItMeshHalfEdge vEdge = v.getHalfEdge();

					vector<zItMeshHalfEdge> cEdges;
					v.getConnectedHalfEdges(cEdges);

					for (auto& he : cEdges)
					{
						if (he.onBoundary())
						{
							vEdge = he;
						}
					}

					//if (vEdge == NULL) continue;


					zItMeshVertex next = vEdge.getVertex();;
					zItMeshVertex prev = vEdge.getPrev().getSym().getVertex();

					zVector vNorm = v.getNormal();

					zVector Ori = v.getPosition();

					zVector v1 = Ori - prev.getPosition();
					v1.normalize();

					zVector n1 = v1 ^ vNorm;
					n1.normalize();

					zVector v2 = next.getPosition() - Ori;
					v2.normalize();

					zVector n2 = v2 ^ vNorm;
					n2.normalize();

					v1 = v1 ^ v2;
					zVector v3 = (n1 + n2);

					v3 *= 0.5;
					v3.normalize();



					double cs = v3 * v2;
					double length = extrudeVal;


					zVector a1 = v2 * cs;
					zVector a2 = v3 - a1;

					double alpha = 0;
					if (a2.length() > 0) alpha = sqrt(a2.length() * a2.length());

					if (cs < 0 && a2.length() > 0) alpha *= -1;

					if (alpha > 0) length /= alpha;

					zVector offPos = Ori + (v3 * length);


					positions.push_back(offPos);


				}

				else
				{
					positions.push_back(v.getPosition());

				}


			}

			// poly connects 
			for (zItMeshFace f(_oMesh); !f.end(); f++)
			{

				vector<int> fVerts;
				f.getVertices(fVerts);

				for (int j = 0; j < fVerts.size(); j++)
				{
					polyConnects.push_back(fVerts[j]);
				}


				polyCounts.push_back(fVerts.size());


			}

		}

		zObjMesh out;
		zFnMesh tempFn(outMesh);

		if (positions.size() > 0)
		{
			tempFn.create(positions, polyCounts, polyConnects);
		}


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::closestPointsToMesh(zPointArray& inPoints, zObjMesh oMesh, zIntArray& faceIDs, zPointArray& closestPoints, zVectorArray& printNorms)
	{

		//zObjMesh oTmpMesh;
		//getPokeMesh(oMesh, oTmpMesh);

		zFnMesh fnMesh(oMesh);

		MatrixXd V;
		MatrixXi F;

		fnMesh.getMatrices_trimesh(V, F);

		MatrixXd P(inPoints.size(), 3);
		for (int i = 0; i < inPoints.size(); i++)
		{
			P(i, 0) = inPoints[i].x;
			P(i, 1) = inPoints[i].y;
			P(i, 2) = inPoints[i].z;
		}


		VectorXd sqrD;
		VectorXi I;
		MatrixXd C;
		igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

		faceIDs.clear();
		closestPoints.clear();
		printNorms.clear();

		faceIDs.assign(inPoints.size(), int());
		closestPoints.assign(inPoints.size(), zPoint());
		printNorms.assign(inPoints.size(), zPoint());

		zVectorArray faceNormals;
		fnMesh.getFaceNormals(faceNormals);

		//cout << endl << endl;

		for (int i = 0; i < inPoints.size(); i++)
		{
			faceIDs[i] = I(i);

			closestPoints[i].x = C(i, 0);
			closestPoints[i].y = C(i, 1);
			closestPoints[i].z = C(i, 2);

			printNorms[i] = faceNormals[faceIDs[i]] * -1;
			//cout << endl << inPoints[i] << ", " << closestPoints[i] << ", " << faceIDs[i];
		}

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getPrintHeight(zPointArray& pPoints, zVectorArray& pNorms, zObjMesh& o_Mesh, zFloatArray& pHeights, zObjGraph& outPrintHeightLines)
	{
		zFnGraph fnGraph(outPrintHeightLines);
		zIntArray eConnects;
		zPointArray gPositions;
		zFnMesh fnMesh(o_Mesh);
		printf("\n pPoint count %i pNorms %i numPolygons %i", pPoints.size(), pNorms.size(), fnMesh.numPolygons());
		for (int i = 0; i < pPoints.size(); i++)
		{
			float d = 10000;

			zPoint closestPt;
			bool check = false;
			for (zItMeshFace f(o_Mesh); !f.end(); f++)
			{
				zPointArray fVPositions;
				f.getVertexPositions(fVPositions);
				zPoint cP;
				//cout << endl << fVPositions[0] << "| " << fVPositions[1] << " | " << fVPositions[2] << " || " << pNorms[i] << " || " << pPoints[i];
				bool chk = core.ray_triangleIntersection(fVPositions[0], fVPositions[1], fVPositions[2], pNorms[i], pPoints[i], cP);

				if (chk)
				{
					if (cP.distanceTo(pPoints[i]) < d)
					{
						d = cP.distanceTo(pPoints[i]);
						closestPt = cP;
						check = true;
					}
				}
				
			}

			if (check)
			{
				pHeights.push_back(d);
				gPositions.push_back(pPoints[i]);
				gPositions.push_back(closestPt);
				eConnects.push_back(gPositions.size() - 2);
				eConnects.push_back(gPositions.size() - 1);
			}
			if (!check)
			{
				printf("\n no intersection found for point %i ", i);
			}
			
		}
		fnGraph.create(gPositions, eConnects);
		printf("\n printHeightGraph v %i - e %i ", fnGraph.numVertices(), fnGraph.numEdges());
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::projectToMesh(zPointArray& pPoints, zObjMesh& o_Mesh, zPointArray& updatePts, zVectorArray& pNorm)
	{
		updatePts.clear();
		updatePts.assign(pPoints.size(), zPoint());
		zVector norm(0, 0, 1);
		pNorm.clear();
		pNorm.assign(pPoints.size(), zVector());

		for (int i = 0; i < pPoints.size(); i++)
		{
			float d = 10000;

			zPoint closestPt;
			zVector closesNormal;
			for (zItMeshFace f(o_Mesh); !f.end(); f++)
			{
				zPointArray fVPositions;
				f.getVertexPositions(fVPositions);


				zPoint cP;
				//core.ray_triangleIntersection(fVPositions[0], fVPositions[1], fVPositions[2], pNorms[i], pPoints[i], cP);
				core.ray_triangleIntersection(fVPositions[0], fVPositions[1], fVPositions[2], norm, pPoints[i], cP);

				if (cP.distanceTo(pPoints[i]) < d)
				{
					d = cP.distanceTo(pPoints[i]);
					closestPt = cP;
					closesNormal = f.getNormal();
				}
			}
			updatePts[i] = closestPt;
			pNorm[i] = closesNormal;

		}
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::UVParametrisation(zObjMesh& oTmpMesh, zObjMesh& oParamMesh)
	{
		//zObjMesh oTmpMesh;
		//getPokeMesh(oMesh, oTmpMesh);

		zFnMesh fnMesh(oTmpMesh);

		MatrixXd V;
		MatrixXi F;

		fnMesh.getMatrices_trimesh(V, F);

		// Initialize UV matrix
		Eigen::MatrixXd initial_guess_UV(V.rows(), 2);
		initial_guess_UV.setZero();

		// Boundary loop (to use as fixed points for LSCM)
		Eigen::VectorXi bnd;
		igl::boundary_loop(F, bnd);

		// Fixed positions for boundary vertices (arc length parameterization)
		Eigen::MatrixXd bc(bnd.size(), 2);
		double total_length = 0;
		Eigen::VectorXd lengths(bnd.size());
		for (int i = 0; i < bnd.size(); ++i)
		{
			lengths(i) = (V.row(bnd(i)) - V.row(bnd((i + 1) % bnd.size()))).norm();
			total_length += lengths(i);
		}
		double current_length = 0;
		for (int i = 0; i < bnd.size(); ++i)
		{
			current_length += lengths(i);
			bc(i, 0) = cos(2.0 * Z_PI * current_length / total_length);
			bc(i, 1) = sin(2.0 * Z_PI * current_length / total_length);
		}

		// Compute LSCM parameterization
		Eigen::MatrixXd UV;
		igl::lscm(V, F, bnd, bc, UV);

		// update vertex positions

		oParamMesh = oTmpMesh;
		zFnMesh fnParamMesh(oParamMesh);
		zPoint* vPositions = fnParamMesh.getRawVertexPositions();

		for (int i = 0; i < fnParamMesh.numVertices(); i++)
		{
			vPositions[i].x = UV(i, 0);
			vPositions[i].y = UV(i, 1);
			vPositions[i].z = 0;
		}

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getBaryCentricCoordinates_triangle(zPoint& pt, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& baryCoordinates)
	{
		zVector v0 = t1 - t0;
		zVector v1 = t2 - t0;
		zVector v2 = pt - t0;

		float d00 = v0 * (v0);
		float d01 = v0 * (v1);
		float d11 = v1 * (v1);
		float d20 = v2 * (v0);
		float d21 = v2 * (v1);

		float denom = d00 * d11 - d01 * d01;

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0 - v - w;

		baryCoordinates = zPoint(u, v, w);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getProjectionPoint_triangle(zPoint& baryCoordinates, zPoint& t0, zPoint& t1, zPoint& t2, zPoint& projectionPt)
	{
		projectionPt = t0 * baryCoordinates.x + t1 * baryCoordinates.y + t2 * baryCoordinates.z;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::barycentericProjection_triMesh(zObjGraph& o_graph, zObjMesh& o_inMesh, zObjMesh& o_projectionMesh, zVectorArray& outNotmals)
	{
		zFnGraph fnGraph(o_graph);

		zPointArray positions;
		fnGraph.getVertexPositions(positions);
		outNotmals.clear();
		bool check = false;
		zFnMesh fnMesh(o_inMesh);
		
		for (auto& p : positions)
		{
			for (zItMeshFace f(o_inMesh); !f.end(); f++)
			{
				
				zPointArray fVerts;
				f.getVertexPositions(fVerts);
				bool ptIntTriangle = core.pointInTriangle(p, fVerts[0], fVerts[1], fVerts[2]);
				if (ptIntTriangle)
				{
					zPoint baryCoordinates;
					getBaryCentricCoordinates_triangle(p, fVerts[0], fVerts[1], fVerts[2], baryCoordinates);


					zItMeshFace fProjection(o_projectionMesh, f.getId());

					zPointArray fVerts_projection;
					fProjection.getVertexPositions(fVerts_projection);

					zPoint projectionPt;
					getProjectionPoint_triangle(baryCoordinates, fVerts_projection[0], fVerts_projection[1], fVerts_projection[2], projectionPt);

					p = projectionPt;
					//zItMeshFace fProjection(o_projectionMesh, f.getId());
					zVector n = fProjection.getNormal();
					outNotmals.push_back(n);
					check = true;
					break;
				}
			}
			if (!check)
			{

				printf("\n barycentericProjection_triMesh check %i", check);
			}
		}

		fnGraph.setVertexPositions(positions);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::unrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zIntPairArray& bsf_vertexPairs, zTransform& outTransformStart)
	{
		zFnMesh fnMesh(o_mesh);
		zPoint* vPositions = fnMesh.getRawVertexPositions();

		zFnMesh fnMesh_unroll(o_mesh_unroll);
		zPoint* vPositions_unroll = fnMesh_unroll.getRawVertexPositions();
		//zTransform outTransformStart;
		/*for (int i = 0; i < oriVertex_UnrollVertex_map.size(); i++)
		{
			zPoint p = vPositions[i];

			for (auto vID : oriVertex_UnrollVertex_map[i])
			{
				vPositions_unroll[vID] = p;
			}
		}*/

		// unroll 
		// https://computergraphics.stackexchange.com/questions/8774/unfold-a-3d-mesh-to-a-2d-plane

		for (int i = 0; i < bsf_vertexPairs.size(); i++)
		{
			zItMeshFace f1(o_mesh, bsf_vertexPairs[i].first);
			zItMeshFace f2(o_mesh, bsf_vertexPairs[i].second);

			zItMeshFace f1_unroll(o_mesh_unroll, bsf_vertexPairs[i].first);
			zItMeshFace f2_unroll(o_mesh_unroll, bsf_vertexPairs[i].second);

			zIntPair hePair = getCommonEdge(f1, f2);

			zItMeshHalfEdge he_1(o_mesh, hePair.first);
			zItMeshHalfEdge he_2(o_mesh, hePair.second);

			zPoint A = vPositions[he_2.getStartVertex().getId()];
			zPoint B = vPositions[he_2.getVertex().getId()];

			// unroll first face
			if (i == 0)
			{
				zItMeshHalfEdge he_walker_1 = he_1;
				int f1_numV = f1.getNumVertices();

				float l_ab = he_1.getLength();

				//set the transform
				zVector y = B - A;
				zVector x = zVector(0, 0, 1) ^ y;
				zVector z = zVector(0, 0, 1);
				outTransformStart = core.getTransformFromVectors(A, x, y, z);

				zPoint a(0, 0, 0);
				zPoint b(a.x, a.y + l_ab, a.z);

				// update  postions of corresponding a & b in unroll mesh
				zIntPair hashKey_a(f1.getId(), he_1.getStartVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_a = oriFaceVertex_UnrollVertex.find(hashKey_a);
				if (got_a != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_a->second] = a;

				zIntPair hashKey_b(f1.getId(), he_1.getVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_b = oriFaceVertex_UnrollVertex.find(hashKey_b);
				if (got_b != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_b->second] = b;

				for (int j = 0; j < f1_numV; j++)
				{
					he_walker_1 = he_walker_1.getNext();
					zPoint C = vPositions[he_walker_1.getVertex().getId()];

					zVector ca = C - A;
					zVector ba = B - A;;

					float s = ((ba ^ ca).length()) / (l_ab * l_ab);
					float c = (ba * ca) / (l_ab * l_ab);

					// alternate point
					/*zPoint c1;
					c1.x = a.x + c * (b.x - a.x) - s * (b.y - a.y);
					c1.y = a.y + c * (b.y - a.y) + s * (b.x - a.x);
					c1.z = 0;*/

					zPoint c1;
					c1.x = a.x + c * (b.x - a.x) + s * (b.y - a.y);
					c1.y = a.y + c * (b.y - a.y) - s * (b.x - a.x);
					c1.z = 0;

					// update  postions of corresponding a & b in unroll mesh
					zIntPair hashKey_c(f1.getId(), he_walker_1.getVertex().getId());
					std::unordered_map<zIntPair, int>::const_iterator got_c = oriFaceVertex_UnrollVertex.find(hashKey_c);
					if (got_c != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_c->second] = c1;
				}

			}





			zItMeshHalfEdge he_walker_2 = he_2;
			int f2_numV = f2.getNumVertices();

			// get positions of the prev edge unrolled.
			zIntPair hashKey_a_prev(f1.getId(), he_2.getStartVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_a_prev = oriFaceVertex_UnrollVertex.find(hashKey_a_prev);
			zPoint a = vPositions_unroll[got_a_prev->second];

			zIntPair hashKey_b_prev(f1.getId(), he_2.getVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_b_prev = oriFaceVertex_UnrollVertex.find(hashKey_b_prev);
			zPoint b = vPositions_unroll[got_b_prev->second];

			// update  postions of corresponding a & b in unroll mesh
			zIntPair hashKey_a(f2.getId(), he_2.getStartVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_a = oriFaceVertex_UnrollVertex.find(hashKey_a);
			if (got_a != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_a->second] = a;

			zIntPair hashKey_b(f2.getId(), he_2.getVertex().getId());
			std::unordered_map<zIntPair, int>::const_iterator got_b = oriFaceVertex_UnrollVertex.find(hashKey_b);
			if (got_b != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_b->second] = b;


			for (int j = 0; j < f2_numV; j++)
			{
				he_walker_2 = he_walker_2.getNext();
				zPoint C = vPositions[he_walker_2.getVertex().getId()];

				float l_ab = he_2.getLength();

				zVector ca = C - A;
				zVector ba = B - A;;

				float s = ((ba ^ ca).length()) / (l_ab * l_ab);
				float c = (ba * ca) / (l_ab * l_ab);

				zPoint c1;
				c1.x = a.x + c * (b.x - a.x) - s * (b.y - a.y);
				c1.y = a.y + c * (b.y - a.y) + s * (b.x - a.x);
				c1.z = 0;

				// update  postions of corresponding a & b in unroll mesh
				zIntPair hashKey_c(f2.getId(), he_walker_2.getVertex().getId());
				std::unordered_map<zIntPair, int>::const_iterator got_c = oriFaceVertex_UnrollVertex.find(hashKey_c);
				if (got_c != oriFaceVertex_UnrollVertex.end()) vPositions_unroll[got_c->second] = c1;

			}


		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::creatUnrollMesh(zObjMesh& o_mesh, zObjMesh& o_mesh_unroll, zObjGraph& o_dualgraph, zInt2DArray& oriVertex_UnrollVertex_map, unordered_map<zIntPair, int, zPair_hash>& oriFaceVertex_UnrollVertex, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs)
	{
		zFnMesh fnMesh(o_mesh);

		computeDualGraph_BST(o_mesh, o_dualgraph, bsf_Vertices, bsf_vertexPairs);

		zPoint* vPositions = fnMesh.getRawVertexPositions();
		zColor* vColors = fnMesh.getRawVertexColors();
		zPointArray positions;
		zIntArray pConnects, pCounts;
		zColorArray colors;
		oriVertex_UnrollVertex_map.clear();
		oriVertex_UnrollVertex_map.assign(fnMesh.numVertices(), zIntArray());

		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zIntArray fVerts;
			f.getVertices(fVerts);

			for (auto fV : fVerts)
			{
				int numVerts = positions.size();

				pConnects.push_back(numVerts);
				oriVertex_UnrollVertex_map[fV].push_back(numVerts);

				zIntPair hashKey(f.getId(), fV);
				oriFaceVertex_UnrollVertex[hashKey] = numVerts;

				positions.push_back(vPositions[fV]);
				colors.push_back(vColors[fV]);

			}

			pCounts.push_back(fVerts.size());
		}


		zFnMesh fnMesh_unroll(o_mesh_unroll);
		fnMesh_unroll.create(positions, pCounts, pConnects);
		fnMesh_unroll.setVertexColors(colors);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeDualGraph_BST(zObjMesh& o_mesh, zObjGraph& o_graph, zItGraphVertexArray& bsf_Vertices, zIntPairArray& bsf_vertexPairs)
	{
		zFnMesh fnMesh(o_mesh);

		zIntArray inEdge_dualEdge;
		zIntArray dualEdge_inEdge;
		fnMesh.getDualGraph(o_graph, inEdge_dualEdge, dualEdge_inEdge, true, false, false);

		zFnGraph fnGraph(o_graph);
		fnGraph.setEdgeColor(zColor(1, 1, 0, 1));

		zItGraphVertex v_MaxValence;
		int maxValence = 0;;

		for (zItGraphVertex v(o_graph); !v.end(); v++)
		{
			if (v.getValence() > maxValence)
			{
				v_MaxValence = v;
				maxValence = v.getValence();
			}
		}

		maxValence += 1;

		// breadth search first sorting

		v_MaxValence.getBSF(bsf_Vertices, bsf_vertexPairs);


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::mergeMesh(zObjMesh& o_mesh)
	{
		zObjMesh oTmpMesh = o_mesh;

		zPointArray positions;
		zIntArray pCounts, pConnects;
		zColorArray colors;

		for (zItMeshFace f(o_mesh); !f.end(); f++)
		{
			zPointArray fVPositions;
			f.getVertexPositions(fVPositions);
			zItMeshVertexArray fVertices;
			f.getVertices(fVertices);
			for (auto& p : fVertices)
			{

				int id = -1;
				zPoint pPos = p.getPosition();
				if (!core.checkRepeatVector(pPos, positions, id, 2))
				{
					id = positions.size();
					positions.push_back(p.getPosition());
					colors.push_back(p.getColor());
				}

				pConnects.push_back(id);
			}

			pCounts.push_back(fVPositions.size());

		}

		zFnMesh fnMesh(o_mesh);
		fnMesh.clear();
		fnMesh.create(positions, pCounts, pConnects);
		fnMesh.setVertexColors(colors);

	}
	ZSPACE_TOOLSETS_INLINE zIntPair zTsNatpowerSDF::getCommonEdge(zItMeshFace& f1, zItMeshFace& f2)
	{
		zIntPair out;

		zItMeshHalfEdgeArray f1_HEdges;
		f1.getHalfEdges(f1_HEdges);

		zItMeshHalfEdgeArray f2_HEdges;
		f2.getHalfEdges(f2_HEdges);


		for (auto& f1HE : f1_HEdges)
		{
			for (auto& f2HE : f2_HEdges)
			{
				if (f1HE.getEdge().getId() == f2HE.getEdge().getId())
				{
					out = zIntPair(f1HE.getId(), f2HE.getId());
					break;
				}
			}
		}

		return out;
	}




	///NON-PLANAR BLOCKS

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeVLoops(zObjMesh& oMesh, zIntArray& medialIDS, zIntArray& featuredNumStrides, zVector& norm, vector<zItMeshHalfEdgeArray>& v_Loops, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom)
	{


		int stride = 0;
		for (int i = 0; i < (featuredNumStrides.size() - 1)/2; i++)
		{
			stride += featuredNumStrides[i];
		}
		featuredNumStrides[0];
		int startVID = medialIDS[0];
		int endVID = medialIDS[1];

		zItMeshVertex vStart(oMesh, startVID);
		zItMeshVertex vEnd(oMesh, endVID);

		zVector dir = vEnd.getPosition() - vStart.getPosition();

		//numFrames = ceil((vEnd.getPosition().z - vStart.getPosition().z) / spacing) + 1;
		//printf("\n numFrames %i", numFrames);

		zItMeshHalfEdgeArray hEdges_Start;
		vStart.getConnectedHalfEdges(hEdges_Start);

		float ang = FLT_MAX;
		zItMeshHalfEdge heStart;

		for (auto& he : hEdges_Start)
		{
			if (he.getVector().angle(dir) < ang)
			{
				ang = he.getVector().angle(dir);
				heStart = he;
			}
		}

		zItMeshHalfEdge he = heStart;
		norm.normalize();

		zItMeshHalfEdge he_Bottom, he_Top;
		int VCounter = 0;
		int tempCounter = 0;
		do
		{
			zVector fNorm = he.getFace().getNormal();
			fNorm.normalize();

			if (norm * fNorm > 0.98)
			{
				VCounter = tempCounter;
				he_Top = he;
			}

			if (norm * fNorm < -0.98)
			{
				he_Bottom = he;
			}

			he = he.getNext().getSym().getNext();
			tempCounter++;

		} while (he != heStart);


		printf("\n VCounter %i ", VCounter);

		zPointArray positions_top, positions_bottom;
		zIntArray pCounts_top, pCounts_bottom;
		zIntArray pConnects_top, pConnects_bottom;

		zIntArray pMap_bottom, pMap_top;

		zFnMesh fnMesh_in(oMesh);
		

		pMap_bottom.assign(fnMesh_in.numVertices(), -1);
		pMap_top.assign(fnMesh_in.numVertices(), -1);

		bool corner = true;

		for (int i = 0; i < stride; i++)
		{
			he_Top = he_Top.getNext().getNext();
			he_Bottom = he_Bottom.getNext().getNext();

			if ((i + 1) % stride != 0)
			{
				he_Bottom = he_Bottom.getSym();
				he_Top = he_Top.getSym();
			}
		}

		zColorArray vColor_top, vColors_bottom;

		zItMeshHalfEdge heWalk_Bottom = he_Bottom;
		int walkCounter = 0;
		int loopCounter;
		do
		{
			if (corner)
			{
				loopCounter = v_Loops.size();
				getLoop(heWalk_Bottom, true, corner, VCounter, v_Loops);
				corner = false;

				pMap_bottom[v_Loops[loopCounter][0].getVertex().getId()] = positions_bottom.size();
				positions_bottom.push_back(v_Loops[loopCounter][0].getVertex().getPosition());
				vColors_bottom.push_back(v_Loops[loopCounter][0].getVertex().getColor());

				pMap_top[v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getId()] = positions_top.size();
				positions_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getPosition());
				vColor_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getColor());


			}

			loopCounter = v_Loops.size();

			//zItMeshHalfEdge he = heWalk.getNext();
			getLoop(heWalk_Bottom, true, corner, VCounter, v_Loops);
			//he.getEdge().setColor(zBLUE);

			pMap_bottom[v_Loops[loopCounter][0].getVertex().getId()] = positions_bottom.size();
			positions_bottom.push_back(v_Loops[loopCounter][0].getVertex().getPosition());
			vColors_bottom.push_back(v_Loops[loopCounter][0].getVertex().getColor());

			pMap_top[v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getId()] = positions_top.size();
			positions_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getPosition());
			vColor_top.push_back(v_Loops[loopCounter][v_Loops[loopCounter].size() - 1].getStartVertex().getColor());

			heWalk_Bottom = heWalk_Bottom.getNext().getNext();

			if ((walkCounter + 1) % (2 * stride) != 0) heWalk_Bottom = heWalk_Bottom.getSym();
			else corner = true;
			walkCounter++;

		} while (heWalk_Bottom != he_Bottom);


		// create meshes

		for (int i = 0; i < stride * 2; i++)
		{

			zIntArray fVerts;
			getFaceVerticesFromHalfedge(he_Bottom, true, fVerts);
			for (auto& id : fVerts) pConnects_bottom.push_back(pMap_bottom[id]);
			pCounts_bottom.push_back(fVerts.size());


			fVerts.clear();
			getFaceVerticesFromHalfedge(he_Top, true, fVerts);
			for (auto& id : fVerts) pConnects_top.push_back(pMap_top[id]);
			pCounts_top.push_back(fVerts.size());

			he_Bottom = he_Bottom.getNext().getNext().getSym();
			he_Top = he_Top.getNext().getNext().getSym();
		}



		zFnMesh fnTop(oMesh_top);
		fnTop.create(positions_top, pCounts_top, pConnects_top);
		fnTop.setVertexColors(vColor_top);

		zFnMesh fnBottom(oMesh_bottom);
		fnBottom.create(positions_bottom, pCounts_bottom, pConnects_bottom);
		fnBottom.setVertexColors(vColors_bottom);

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeGeodesicScalars(zObjMesh& oMesh, vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, bool normalise)
	{
		zFnMesh fnMesh(oMesh);

		scalars.clear();
		scalars.assign(fnMesh.numVertices(), -1);

		float minMaxDist = 10000;
		vector<zDomainFloat> loopDomains;
		loopDomains.assign(v_Loops.size(), zDomainFloat(10000, -10000));

		for (int l = 0; l < v_Loops.size(); l++)
		{
			float length = 0;
			for (int j = 0; j < v_Loops[l].size(); j++)
			{
				if (j == 0)
				{
					scalars[v_Loops[l][j].getVertex().getId()] = length;
					loopDomains[l].min = length;
				}

				length += v_Loops[l][j].getLength();
				scalars[v_Loops[l][j].getStartVertex().getId()] = length;

				if (j == v_Loops[l].size() - 1 && length < minMaxDist) minMaxDist = length;

				if (length > loopDomains[l].max) loopDomains[l].max = length;
			}
		}

		if (normalise)
		{
			zDomainFloat outDomain(0, minMaxDist);
			for (int l = 0; l < v_Loops.size(); l++)
			{
				for (int j = 0; j < v_Loops[l].size(); j++)
				{
					scalars[v_Loops[l][j].getStartVertex().getId()] = core.ofMap(scalars[v_Loops[l][j].getStartVertex().getId()], loopDomains[l], outDomain);
				}
			}
		}

		//_colorMesh(oMesh, scalars);

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		zColor* mesh_vColors = fnMesh.getRawVertexColors();

		zDomainFloat distanceDomain(minScalar, maxScalar);
		printf("\n scalar domain %1.4f %1.4f | minMaxDist %1.4f ", minScalar, maxScalar, minMaxDist);

		zDomainColor colDomain(zColor(1, 0, 0, 1), zColor(0, 1, 0, 1));

		//for (int i = 0; i < fnMesh.numVertices(); i++)
		//{
		//	//mesh_vColors[i] = core.blendColor(scalars[i], distanceDomain, colDomain, zRGB);
		//}

		//fnMesh.computeFaceColorfromVertexColor();
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeGeodesicContours(vector<zItMeshHalfEdgeArray>& v_Loops, zScalarArray& scalars, float spacing, zObjMesh& oMesh_top, zObjMesh& oMesh_bottom, zObjMeshArray& oMeshes)
	{

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		int totalContours = ceil((maxScalar - minScalar) / spacing) ;
		float increments = (maxScalar - minScalar) / totalContours;

		printf("\n totalContours %i increments %1.2f ", totalContours, increments);
		zObjMeshArray oTmpMeshes;

		oTmpMeshes.clear();
		oTmpMeshes.assign(totalContours + 1, oMesh_bottom);

		

		for (int l = 0; l < oTmpMeshes.size(); l++)
		{
			float threshold = l * increments;
			printf("\n threshold %1.4f ", threshold);

			zFnMesh fnMesh(oTmpMeshes[l]);
			zPoint* points = fnMesh.getRawVertexPositions();

			for (int i = 0; i < v_Loops.size(); i++)
			{
				for (int j = 0; j < v_Loops[i].size(); j++)
				{
					float s0 = scalars[v_Loops[i][j].getStartVertex().getId()];
					float s1 = scalars[v_Loops[i][j].getVertex().getId()];

					bool contour = false;
					if (s0 <= threshold && s1 >= threshold)contour = true;
					if (s0 >= threshold && s1 <= threshold)contour = true;

					if (!contour) continue;

					zPoint v0 = v_Loops[i][j].getStartVertex().getPosition();
					zPoint v1 = v_Loops[i][j].getVertex().getPosition();
					

					zPoint pos1 = util_getContourPosition(threshold, v1, v0, s1, s0);
					points[i] = (pos1);
				}
			}

			
			

		}
		/*zFnMesh fnMesh_top(oMesh_top);
		oTmpMeshes.push_back(oMesh_top);*/

		// poke mesh 
		oMeshes.clear();
		oMeshes.assign(oTmpMeshes.size(), zObjMesh());

		for (int l = 0; l < oTmpMeshes.size(); l++)
		{
			getPokeMesh(oTmpMeshes[l], oMeshes[l]);
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeGeodesicContours(zObjMesh& o_mesh, zFloatArray& scalars, float spacing, zObjGraphArray& o_contourGraphs)
	{

		zScalar minScalar = core.zMin(scalars);
		zScalar maxScalar = core.zMax(scalars);

		int totalContours = ceil((maxScalar - minScalar) / spacing);


		float increments = (maxScalar - minScalar) / totalContours;

		//printf("\n totalContours %i increments %1.2f ", totalContours, increments);

		o_contourGraphs.clear();
		o_contourGraphs.assign(totalContours, zObjGraph());



		for (int i = 0; i < totalContours; i++)
		{
			// Generate the isocontour using the threshold value
			zPointArray positions;
			zIntArray edgeConnects;
			zColorArray vColors;
			int pres = 3;
			zFnMesh fnMesh(o_mesh);
			fnMesh.getIsoContour(scalars, i * increments, positions, edgeConnects, vColors, pres, pow(10, -1 * pres));

			// Create graph from the isocontour
			zFnGraph tempFn(o_contourGraphs[i]);
			tempFn.create(positions, edgeConnects);
			tempFn.setEdgeColor(zColor(1, 1, 1, 1));
			tempFn.setEdgeWeight(2);
			tempFn.setVertexColors(vColors, false);
		}



	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createSectionGraphs(zObjMeshArray& oMeshes, zObjGraphArray& o_sectionsGraphs)
	{
		o_sectionsGraphs.clear();
		o_sectionsGraphs.assign(oMeshes.size(), zObjGraph());

		int counter = 0;
		for (auto& oMesh : oMeshes)
		{
			createBoundaryEdgeGraph(oMesh, true, o_sectionsGraphs[counter]);

			zFnGraph fnGraph(o_sectionsGraphs[counter]);
			fnGraph.setEdgeColor(zGREEN);
			fnGraph.setEdgeWeight(3);

			counter++;
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::transformAllGraphs_planar(int graphId, bool toLocal)
	{
		zTransform t = sectionFrames[graphId];
		zFnGraph fng;
		if (toLocal)
		{
			// Transform
			zTransform tLocal;
			tLocal.setIdentity();

			//transform all graphs to local			
			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);

		}
		else
		{
			// transform back 

			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, true);
		}
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::transformAllGraphs(int graphId, zTransform t , bool toLocal)
	{
		zFnGraph fng;
		if (toLocal)
		{
			// Transform
			zTransform tLocal;
			tLocal.setIdentity();

			//transform all graphs to local


			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, false);
			fng.setTransform(tLocal, true, true);

		}
		else
		{
			// transform back 

			fng = zFnGraph(o_sectionGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_bracing[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_hard[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_features_soft[graphId]);
			fng.setTransform(t, true, true);
			fng = zFnGraph(o_trimGraphs_SlotSide[graphId]);
			fng.setTransform(t, true, true);
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::cleanContourGraph(int graphId)
	{
		zFnGraph fnGraph(o_contourGraphs[graphId]);
		//check for bourndary verticies
		int boundaryCount = 0;
		zItGraphVertexArray vArray;
		for (zItGraphVertex vIt(o_contourGraphs[graphId]); !vIt.end(); vIt++)
		{
			if (vIt.getValence() == 1)
			{
				boundaryCount++;
				vArray.push_back(vIt);
			}
		}

		zPointArray positions;
		zIntArray eConnects;

		fnGraph.getVertexPositions(positions);
		fnGraph.getEdgeData(eConnects);

		if (vArray.size() == 2) //two boundary verticies - we have one loop that is open - still better to create new graph with ordered positions
		{
			eConnects.push_back(vArray[1].getId());
			eConnects.push_back(vArray[0].getId());
			fnGraph.create(positions, eConnects);
			printf("\n cleaning contour graph [%i]. valence = 2 ", graphId);
			/*zItGraphHalfEdge heStart(o_contourGraphs[graphId], vArray[0].getId());
			zItGraphHalfEdge he = heStart;
			int safetyCounter = 0;
			zItGraphHalfEdgeArray loop;
			while (safetyCounter <= fnGraph.numEdges())
			{
				loop.push_back(he);
				he = he.getNext();
				if (he.getStartVertex().checkValency(1)) break;
				safetyCounter++;
			}
			createGraphFromHEArray(loop, o_contourGraphs[graphId]);*/
		}
		else if (vArray.size() == 0) //no boundary verticies - check for how many closed loops
		{
			printf("\n cleaning contour graph [%i]. valence = 0 ", graphId);


			//check for closed loops
			vector<zItGraphHalfEdgeArray> loops;

			zBoolArray visited;
			visited.assign(fnGraph.numEdges(), false);

			int edgeCounter = 0;
			zItGraphHalfEdgeArray longestLoop;
			while (edgeCounter < fnGraph.numEdges())
			{
				zItGraphHalfEdge heStart(o_contourGraphs[graphId], edgeCounter);
				zItGraphHalfEdge he = heStart;
				int safetyCounter = 0;
				zItGraphHalfEdgeArray loop;

				while (safetyCounter <= fnGraph.numEdges())
				{
					visited[he.getEdge().getId()] = true;
					loop.push_back(he);
					he = he.getNext();
					if (he == heStart) break;
					safetyCounter++;
				}
				edgeCounter += loop.size();
				if (loop.size() > longestLoop.size())
				{
					longestLoop = loop;
				}
			}
			util_createGraphFromHEArray(longestLoop, o_contourGraphs[graphId]);
		}
		else
		{
			printf("\n cleaning contour graph [%i]. other - valence = %i ", graphId, vArray.size());

			//more than two boundary verticies - we have multiple loops
			//check for closed loops
			vector<zItGraphHalfEdgeArray> loops;
			zBoolArray visited;
			visited.assign(fnGraph.numEdges(), false);
			int edgeCounter = 0;
			zItGraphHalfEdgeArray longestLoop;
			bool loopClosed = false;
			//get the loop from all boundary verticies
			for (auto& v : vArray)
			{
				zItGraphHalfEdge heStart = v.getHalfEdge();
				if (visited[heStart.getEdge().getId()])
				{
					continue;
				}
				zItGraphHalfEdge he = heStart;
				zItGraphHalfEdgeArray loop;
				int safetyCounter = 0;

				while (safetyCounter <= fnGraph.numEdges())
				{
					visited[he.getEdge().getId()] = true;
					loop.push_back(he);
					he = he.getNext();
					if (he == heStart)
					{
						loopClosed = true;
						break;
					}
					if (he.getStartVertex().checkValency(1))
					{
						loopClosed = false;
						break;
					}
					safetyCounter++;
				}
				if (loop.size() > longestLoop.size())
				{
					longestLoop = loop;
				}
			}
			util_createGraphFromHEArray(longestLoop, o_contourGraphs[graphId]);


			if (!loopClosed)
			{
				printf("\n			cleaning contour graph [%i]. valence = %i - not closed loop", graphId, vArray.size());

				fnGraph.getVertexPositions(positions);
				fnGraph.getEdgeData(eConnects);
				zItGraphVertexArray vs;
				for (zItGraphVertex vIt(o_contourGraphs[graphId]); !vIt.end(); vIt++)
				{
					if (vIt.checkValency(1))
					{
						vs.push_back(vIt);
					}
				}
				//makesure that u have two vertex
				if (vs.size() == 2)
				{
					eConnects.push_back(vs[1].getId());
					eConnects.push_back(vs[0].getId());
					fnGraph.create(positions, eConnects);
				}
				else
				{
					printf	("\n error: invalid cleaning contour graph. valence = %i ", vs.size());
					//throw std::invalid_argument(" error: invalid cleaning contour graph.");
				}
				


			}
			

		}
	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::readJSON(string path, int _blockID, bool runBothPlanes, bool runPlaneLeft)
	{
		printf("\n readJSON 0");

		json j;
		zFnMesh fnMesh(o_GuideMesh);
		bool fileChk = fnMesh.json_read(path, j);

		if (!fileChk)
		{
			throw std::invalid_argument(" error: invalid outPath. ");
			return;
		}

		fnMesh.clear();
		fnMesh.from(path, zJSON);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();


		blockId = _blockID;
		int sID = j["MedialStartEnd"][0];
		int eID = j["MedialStartEnd"][1];

		isCorner = j["IsCorner"];
		int type = j["BlockType"];
		blockType = static_cast<zBlockType>(type);

		zTransform sPlaneRight, ePlaneRight, sPlaneLeft, ePlaneLeft;

		zPoint oS, oE, xS, xE, yS, yE, zS, zE;
		zVector zVectorEndRight, zVectorEndLeft;

		oS = zPoint(j["LeftPlanes"][0][3], j["LeftPlanes"][0][7], j["LeftPlanes"][0][11]);
		oE = zPoint(j["LeftPlanes"][1][3], j["LeftPlanes"][1][7], j["LeftPlanes"][1][11]);

		xS = zPoint(j["LeftPlanes"][0][0], j["LeftPlanes"][0][4], j["LeftPlanes"][0][8]);
		xE = zPoint(j["LeftPlanes"][1][0], j["LeftPlanes"][1][4], j["LeftPlanes"][1][8]);

		yS = zPoint(j["LeftPlanes"][0][1], j["LeftPlanes"][0][5], j["LeftPlanes"][0][9]);
		yE = zPoint(j["LeftPlanes"][1][1], j["LeftPlanes"][1][5], j["LeftPlanes"][1][9]);

		zS = zPoint(j["LeftPlanes"][0][2], j["LeftPlanes"][0][6], j["LeftPlanes"][0][10]);
		zE = zPoint(j["LeftPlanes"][1][2], j["LeftPlanes"][1][6], j["LeftPlanes"][1][10]);
		zVectorEndLeft = zE;
		sPlaneLeft = core.getPlaneFromVectors(oS, xS, yS, zS);
		ePlaneLeft = core.getPlaneFromVectors(oE, xE, yE, zE);

		oS = zPoint(j["RightPlanes"][0][3], j["RightPlanes"][0][7], j["RightPlanes"][0][11]);
		oE = zPoint(j["RightPlanes"][1][3], j["RightPlanes"][1][7], j["RightPlanes"][1][11]);
		xS = zPoint(j["RightPlanes"][0][0], j["RightPlanes"][0][4], j["RightPlanes"][0][8]);
		xE = zPoint(j["RightPlanes"][1][0], j["RightPlanes"][1][4], j["RightPlanes"][1][8]);
		yS = zPoint(j["RightPlanes"][0][1], j["RightPlanes"][0][5], j["RightPlanes"][0][9]);
		yE = zPoint(j["RightPlanes"][1][1], j["RightPlanes"][1][5], j["RightPlanes"][1][9]);
		zS = zPoint(j["RightPlanes"][0][2], j["RightPlanes"][0][6], j["RightPlanes"][0][10]);
		zE = zPoint(j["RightPlanes"][1][2], j["RightPlanes"][1][6], j["RightPlanes"][1][10]);
		zVectorEndRight = zE;

		sPlaneRight = core.getPlaneFromVectors(oS, xS, yS, zS);
		ePlaneRight = core.getPlaneFromVectors(oE, xE, yE, zE);

		runningType = runBothPlanes ? 0 : runPlaneLeft ? 1 : 2;

		if (!runBothPlanes)
		{
			if (isCorner && blockType == zBlockType::Top)
			{
				bool check = abs(zVectorEndLeft.z) < abs(zVectorEndRight.z);


				check = runPlaneLeft ? check : !check;

				if (check)
				{
					sPlaneRight = sPlaneLeft;
					ePlaneRight = ePlaneLeft;
				}
				else
				{
					sPlaneLeft = sPlaneRight;
					ePlaneLeft = ePlaneRight;
				}
			}
		}


		setStartEndPlanes(sPlaneLeft, ePlaneLeft, true);
		base_world = sPlaneLeft;
		setStartEndPlanes(sPlaneRight, ePlaneRight, false);

		base_local.setIdentity();
		//medial axis
		compute_MedialGraph(o_GuideMesh, sID, eID);


		//get feature stride
		FeaturedNumStrides.clear();
		core.json_readAttribute(j, "FeaturedNumStrides", FeaturedNumStrides);
		medialIDS.clear();
		core.json_readAttribute(j, "MedialStartEnd", medialIDS);
		StartCornerVID = j["StartCornerVID"];

		// left mesh
		isRegular = !blockType == zBlockType::Bottom;
		
		if (isRegular)
		{
			compute_SliceMesh_Regular(o_GuideMesh, sID, eID, FeaturedNumStrides);
			printf("\n slice mesh regular");
		}
		else
		{
			//left Mesh
			compute_SliceMesh_Pentagon(o_GuideMesh, sID, eID, FeaturedNumStrides, true);
			printf("\n slice mesh1");

			//Right Mesh
			compute_SliceMesh_Pentagon(o_GuideMesh, sID, eID, FeaturedNumStrides, false);
			printf("\n slice mesh2");
		}
	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::get2DArrayFromTransform(zTransform& transform, vector<zDoubleArray>& arr)
	{
		//vector<zDoubleArray> arr;
		arr.clear();
		arr.assign(4, zDoubleArray());
		for (int i = 0; i < 4; i++)
		{
			arr[i].assign(4, 0.0);
			for (int j = 0; j < 4; j++)
			{
				arr[i][j] = transform(i, j);
			}
		}
	}

}

