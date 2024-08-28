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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createFieldMesh(zDomain<zPoint>& bb, int resX, int resY)
	{
		zFnMeshScalarField fnField(o_field);

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
		bool jsonCheck = coreUtils.json_read(path, j);

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

		if (planarBlock)
		{

			readJSON_planarBlock(path, _blockID, runBothPlanes, runPlaneLeft);
		}
		else
		{
			//to be implemented
		}


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

		if (numGraphs == 0)return out;

		for (auto& graph : o_trimGraphs)
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

	//---- COMPUTE METHODS 

	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::isPlanarBlock()
	{
		return planarBlock;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeFrontBackHalfEdges(zObjGraph& graph, zItGraphHalfEdgeArray& outHeInner, zItGraphHalfEdgeArray& outHEOuter)
	{
		//get HE for front and back sides
		outHeInner.clear();
		outHEOuter.clear();
		//vector<zItGraphHalfEdges> hesBack, hesFront;
		//starting from boundary vertex, walk all half edges till you reach next boundary vertex
		zItGraphHalfEdge he1, he2;
		for (zItGraphVertex v(graph); !v.end(); v++)
		{
			if (v.getColor() == _colorCornersBoundary)
			{
				zItGraphHalfEdgeArray hc;
				v.getConnectedHalfEdges(hc);
				for (int i; i < 2; i++)
				{
					if (hc[i].getVertex().getColor() == _colorCornersBoundary) he1 = hc[i].getNext();
					else he2 = hc[i];
				}
				break;
			}
		}
		bool he1IsOuter = false;
		int safetyCounter = 0;
		zItGraphHalfEdgeArray hes1, hes2;
		while (safetyCounter < 1000)
		{
			safetyCounter++;
			hes1.push_back(he1);
			he1 = he1.getNext();
			if (he1.getVertex().getColor() == _colorFeatureOuter) he1IsOuter = true;
			if (he1.getStartVertex().getColor() == _colorCornersBoundary) break;
		}
		safetyCounter = 0;
		while (safetyCounter < 1000)
		{
			safetyCounter++;
			hes2.push_back(he2);
			he2 = he2.getNext();
			if (he2.getStartVertex().getColor() == _colorCornersBoundary) break;
		}

		outHEOuter = he1IsOuter ? hes1 : hes2;
		outHeInner = !he1IsOuter ? hes1 : hes2;

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computePrintSectionFromPlaneSpacing(float printPlaneSpacing, zDomainFloat& _printHeightDomain, zDomainFloat _neopreneOffset, bool& frameCHECKS, bool& sdfCHECKS, bool& geomCHECKS)
	{
		frameCHECKS = false;
		geomCHECKS = true;
		sdfCHECKS = true;

		printf("\n printPlaneSpace %1.4f ", printPlaneSpacing);

		sectionFrames.clear();
		computePrintBlockFrames(printPlaneSpacing, neopreneOffset.min, neopreneOffset.max, true);
		if (!isRegular) computePrintBlockFrames(printPlaneSpacing, neopreneOffset.min, neopreneOffset.max, false);

		o_sectionGraphs.clear();
		o_sectionGraphs.assign(sectionFrames.size(), zObjGraph());
		computePrintBlockSections(true);
		if (!isRegular) computePrintBlockSections(false);

		frameCHECKS = checkPrintLayerHeights(sdfCHECKS, geomCHECKS);
		printf("\n");
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computePrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, zDomainFloat _neopreneOffset, bool compFrames, bool compSDF)
	{
		printHeightDomain = _printHeightDomain;
		neopreneOffset = _neopreneOffset;

		bool frameCHECKS = false;
		bool geomCHECKS = true;
		bool sdfCHECKS = true;

		int minCriticalPtsCount = INT_MAX;
		float bestPlaneSpacing = FLT_MAX;


		if (compFrames)
		{

			for (float printPlaneSpacing = printHeightDomain.max; printPlaneSpacing >= printHeightDomain.min; printPlaneSpacing -= 0.00025)
			{
				computePrintSectionFromPlaneSpacing(printPlaneSpacing, printHeightDomain, _neopreneOffset, frameCHECKS, sdfCHECKS, geomCHECKS);

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
				printf("\n Layer height check FALSE!  Choose best plane spacing %1.4f", bestPlaneSpacing);

				computePrintSectionFromPlaneSpacing(bestPlaneSpacing, printHeightDomain, _neopreneOffset, frameCHECKS, sdfCHECKS, geomCHECKS);
			}

			printf("\n frameCHECKS %s ", (frameCHECKS) ? "T" : "F");
		}
		printf("\n layerCheck %i | geomChk %i | sdCheck %i ", (frameCHECKS), geomCHECKS, sdfCHECKS);
		printf("\n sectionFrames %i | o_sectionGraphs %i", sectionFrames.size(), o_sectionGraphs.size());

		o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());
		for (int i = 1; i < o_sectionGraphs.size(); i++)
		{
			zObjGraph slotGraph0;
			float graphLength0 = 0.2;
			slotGraph_Arch(i, graphLength0, slotGraph0);
		}
		


		if (compSDF)
		{
			printf("\n \n  SDF \n \n");

			computeSDF(allSDFLayers, numSDFlayers, funcNum, numSmooth, printLayerWidth, neopreneOffset.min, raftLayerWidth);
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeMedialGraph(zObjMesh& o_Mesh, int startVID, int endVID)
	{

		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeMedial_BraceEdges(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, int braceStride)
	{
		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeSliceMesh(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides, bool left)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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

		int featureWalkNum = (FeaturedNumStrides.size() - 1) / 2;
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
				for (int k = featureStart; k < featureWalkNum; k++)
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
							bool chkRepeat = coreUtils.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								coreUtils.addToPositionMap(positionVertex, p, vID, PRECISION);
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
				for (int k = featureStart; k < featureWalkNum; k++)
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
							bool chkRepeat = coreUtils.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								coreUtils.addToPositionMap(positionVertex, p, vID, PRECISION);
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

			for (int k = featureStart; k < featureWalkNum; k++)
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
				coreUtils.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				coreUtils.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				coreUtils.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				coreUtils.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				coreUtils.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				coreUtils.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				coreUtils.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				coreUtils.vertexExists(positionVertex, p7, PRECISION, v7);

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

			for (int k = 0; k < featureWalkNum; k++)
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
				coreUtils.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				coreUtils.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				coreUtils.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				coreUtils.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				coreUtils.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				coreUtils.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				coreUtils.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				coreUtils.vertexExists(positionVertex, p7, PRECISION, v7);

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



		//walk on the slice mesh to color the edges
		//get heStart
		zFnMesh fnGuide(o_Mesh);
		zPointArray guidePts;
		fnGuide.getVertexPositions(guidePts);
		int vSIDSlice, vEIDSlice;
		vSIDSlice = coreUtils.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = coreUtils.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
		zItMeshHalfEdge heEndSlice = heStartSlice;
		int sideDivisionCount = 1;
		while (heEndSlice.getVertex().getId() != vEIDSlice)
		{
			heEndSlice.getEdge().setColor(zCYAN);

			heEndSlice = heEndSlice.getNext().getSym().getNext();
			sideDivisionCount++;
			//printf("\n sideDivisionCount %i ", sideDivisionCount);
		}
		heEndSlice.getEdge().setColor(zCYAN);

		//get start of all feature edges
		zItMeshHalfEdgeArray hesFeatureInner, hesFeatureOuter;
		//zItMeshHalfEdge heTemp = left ? heStartSlice.getSym().getNext() : heStartSlice.getPrev();
		zItMeshHalfEdge heTemp = heStartSlice;

		hesFeatureInner.push_back(heTemp);
		hesFeatureOuter.push_back(
			left ?
			heTemp.getPrev().getPrev().getSym() : heTemp.getSym().getNext().getNext());

		//hesFeatureOuter.push_back(heTemp.getPrev().getPrev().getSym());

		vector<bool> featureCheck;
		featureCheck.push_back(true);

		fnGuide.setEdgeColor(zBLACK);

		for (int k = featureStart; k < featureWalkNum; k++)
		{
			int blockStride = FeaturedNumStrides[k];
			for (int i = 0; i < blockStride; i++)
			{
				featureCheck.push_back(i == blockStride - 1);

				heTemp = left ? heTemp.getSym().getNext().getNext() : heTemp.getPrev().getPrev().getSym();

				//}
				hesFeatureInner.push_back(heTemp);

				zItMeshHalfEdge eOuter;
				if (left)
				{
					eOuter = heTemp.
						getPrev().getSym().getPrev().
						getPrev().getSym().getPrev().
						getSym();
				}
				else
				{
					eOuter = heTemp.
						getSym().
						getNext().getSym().getNext().
						getNext().getSym().getNext();
				}
				hesFeatureOuter.push_back(eOuter);
			}

		}

		//assuming bridging is always 1



		//printf("\n hesFeatureInner size %i", hesFeatureInner.size());
		//for (auto& e : hesFeatureInner)
		//{
		//	//e.getEdge().setColor(zRED);
		//}
		for (int i = 0; i < hesFeatureInner.size(); i++)
		{
			//hesFeatureInner[i].getEdge().setColor(zMAGENTA);
			//hesFeatureOuter[i].getEdge().setColor(zORANGE);
			zItMeshHalfEdge eInner = hesFeatureInner[i];
			zItMeshHalfEdge eOuter = hesFeatureOuter[i];

			//bool slotType1 = true;



			bool startChk = i == FeaturedNumStrides[0] / 2;

			bool patternChk = blockType != zBlockType::Wall ? i > FeaturedNumStrides[0] : true;
			//printf("\n patternChk %d	- checkWall %d		- i %i		- FeaturedNumStrides[0] %i", patternChk, blockType == zBlockType::Wall, i, FeaturedNumStrides[0]);

			zColor outerColor = startChk ? _colorFeatureOuter : _colorPattern;
			zColor innerColor = startChk ? _colorFeatureInner : zBLACK;

			if (i == 0)
			{
				//innerColor = _colorCornersInterior;
				//outerColor = _colorCornersInterior;
				innerColor = _colorFeatureInner;
				outerColor = _colorFeatureOuter;
			}
			else if (i == hesFeatureInner.size() - 1)
			{
				innerColor = _colorCornersBoundary;
				outerColor = _colorCornersBoundary;
			}
			else if (featureCheck[i])
			{
				innerColor = _colorFeatureInner;
				outerColor = _colorFeatureOuter;
			}
			else
			{
				innerColor = zBLACK;
				outerColor = patternChk ? _colorPattern : zBLACK;
			}

			//innerColor = featureCheck[i] ? zRED : zGREY;
			//outerColor = featureCheck[i] ? zCYAN : zGREY;

			hesFeatureInner[i].getEdge().setColor(innerColor);
			hesFeatureOuter[i].getEdge().setColor(outerColor);

			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				eInner = eInner.getNext().getSym().getNext();
				eOuter = eOuter.getNext().getSym().getNext();

				eInner.getEdge().setColor(innerColor);
				eOuter.getEdge().setColor(outerColor);

			}



		}


		//*/





		//printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeSliceMesh_Regular(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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

		const int featureWalkNum = (FeaturedNumStrides.size() - 1) / 2;
		int featureStart = 0;
		/*
		do
		{
			////right
			if (!left)
			{
				//for (int ff = FeaturedNumStrides.size()-1; ff >= 0; ff--)
				//{
					//int blockStride = FeaturedNumStrides[ff];
				zItMeshHalfEdge he_Left = he.getPrev();

				//for (int& blockStride : FeaturedNumStrides)
				for (int k = featureStart; k < featureWalkNum; k++)
				{
					int blockStride = FeaturedNumStrides[k];
					for (int i = 0; i < blockStride; i++)
					{

						zItMeshHalfEdge heTmp = he_Left;
						int pCount = 0;
						//get each face
						do
						{
							zPoint p = heTmp.getVertex().getPosition();
							int vID;
							bool chkRepeat = coreUtils.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								coreUtils.addToPositionMap(positionVertex, p, vID, PRECISION);
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

			////left
			else
			{
				zItMeshHalfEdge he_Right = he.getSym().getNext();
				//for (int& blockStride : FeaturedNumStrides)
				for (int k = featureStart; k < featureWalkNum; k++)
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
							bool chkRepeat = coreUtils.vertexExists(positionVertex, p, PRECISION, vID);

							if (!chkRepeat)
							{
								vID = positions.size();
								coreUtils.addToPositionMap(positionVertex, p, vID, PRECISION);
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

		printf("\n sliceMesh input %i %i %i ", positions.size(), pCounts.size(), pConnects.size());


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

			for (int k = featureStart; k < featureWalkNum; k++)
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
				coreUtils.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				coreUtils.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				coreUtils.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				coreUtils.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				coreUtils.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				coreUtils.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				coreUtils.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				coreUtils.vertexExists(positionVertex, p7, PRECISION, v7);

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

			for (int k = 0; k < featureWalkNum; k++)
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
				coreUtils.vertexExists(positionVertex, p0, PRECISION, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				coreUtils.vertexExists(positionVertex, p1, PRECISION, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				coreUtils.vertexExists(positionVertex, p2, PRECISION, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				coreUtils.vertexExists(positionVertex, p3, PRECISION, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				coreUtils.vertexExists(positionVertex, p4, PRECISION, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				coreUtils.vertexExists(positionVertex, p5, PRECISION, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				coreUtils.vertexExists(positionVertex, p6, PRECISION, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				coreUtils.vertexExists(positionVertex, p7, PRECISION, v7);

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
		*/



		zObjMesh* o_sliceMesh = &o_SliceMesh_Left;
		o_SliceMesh_Right.mesh.clear();

		zFnMesh fnMesh(*o_sliceMesh);


		zFnMesh fnGuideMesh(o_Mesh);
		fnGuideMesh.getVertexPositions(positions);
		fnGuideMesh.getPolygonData(pConnects, pCounts);
		fnMesh.create(positions, pCounts, pConnects);



		fnMesh.setFaceColor(zGREY);


		//walk on the slice mesh to color the edges
		//get heStart
		zFnMesh fnGuide(o_Mesh);
		zPointArray guidePts;
		fnGuide.getVertexPositions(guidePts);
		int vSIDSlice, vEIDSlice;
		vSIDSlice = coreUtils.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = coreUtils.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
		zItMeshHalfEdge heEndSlice = heStartSlice;
		int sideDivisionCount = 1;
		while (heEndSlice.getVertex().getId() != vEIDSlice)
		{
			heEndSlice.getEdge().setColor(zCYAN);

			heEndSlice = heEndSlice.getNext().getSym().getNext();
			sideDivisionCount++;
		}
		heEndSlice.getEdge().setColor(zCYAN);

		//get start of all feature edges
		zItMeshHalfEdgeArray hesFeatureInner, hesFeatureOuter;
		zItMeshHalfEdge heTemp0 = heStartSlice;
		zItMeshHalfEdge heTemp1 = heStartSlice;



		hesFeatureInner.push_back(heTemp0);
		hesFeatureInner.push_back(heTemp1);

		hesFeatureOuter.push_back(heTemp0.getPrev().getPrev().getSym());
		hesFeatureOuter.push_back(heTemp1.getSym().getNext().getNext());

		vector<bool> featureCheck;
		featureCheck.push_back(true);
		featureCheck.push_back(true);

		fnGuide.setEdgeColor(zBLACK);

		for (int k = 0; k < featureWalkNum; k++)
		{
			int blockStride = FeaturedNumStrides[k];
			//if (checkWall) blockStride--;
			for (int i = 0; i < blockStride; i++)
			{
				featureCheck.push_back(i == blockStride - 1);
				featureCheck.push_back(i == blockStride - 1);

				heTemp1 = heTemp1.getSym().getNext().getNext();
				heTemp0 = heTemp0.getPrev().getPrev().getSym();

				//}
				hesFeatureInner.push_back(heTemp0);
				hesFeatureInner.push_back(heTemp1);

				zItMeshHalfEdge eOuter0, eOuter1;
				{
					eOuter1 = heTemp1.
						getPrev().getSym().getPrev().
						getPrev().getSym().getPrev().
						getSym();
				}
				{
					eOuter0 = heTemp0.
						getSym().
						getNext().getSym().getNext().
						getNext().getSym().getNext();
				}
				hesFeatureOuter.push_back(eOuter0);
				hesFeatureOuter.push_back(eOuter1);

				/*eOuter0.getEdge().setColor(zBLUE);
				eOuter1.getEdge().setColor(zORANGE);
				heTemp0.getEdge().setColor(zGREEN);
				heTemp1.getEdge().setColor(zRED);*/
			}

		}

		//assuming bridging is always 1
		for (int m = 0; m < hesFeatureInner.size(); m++)
		{
			zItMeshHalfEdge eInner = hesFeatureInner[m];
			zItMeshHalfEdge eOuter = hesFeatureOuter[m];

			//bool slotType1 = true;

			int i = m % 2 == 0 ? m / 2 : (m - 1) / 2;

			bool startChk = i == FeaturedNumStrides[0] / 2;

			bool patternChk = blockType == zBlockType::Wall ? true : i > FeaturedNumStrides[0];

			zColor outerColor = startChk ? _colorFeatureOuter : _colorPattern;
			zColor innerColor = startChk ? _colorFeatureInner : zBLACK;
			printf("\n m %i - i %i | %d - %d - %d ", m, i, blockType == zBlockType::Wall, patternChk, i == (hesFeatureInner.size() / 2) - 1);

			if (i == 0)
			{
				innerColor = _colorCornersInterior;
				outerColor = _colorCornersInterior;
			}
			if (i == (hesFeatureInner.size() / 2) - 1)
			{
				innerColor = _colorCornersBoundary;
				outerColor = _colorCornersBoundary;
			}
			else if (startChk)
			{
				innerColor = _colorFeatureInner;
				outerColor = _colorFeatureOuter;
			}
			else
			{
				innerColor = zBLACK;
				outerColor = patternChk ? _colorPattern : zBLACK;
			}



			//innerColor = featureCheck[m] ? zRED : zGREY;
			//outerColor = featureCheck[m] ? zCYAN : zGREY;

			eInner.getEdge().setColor(innerColor);
			eOuter.getEdge().setColor(outerColor);

			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				eInner = eInner.getNext().getSym().getNext();
				eOuter = eOuter.getNext().getSym().getNext();

				eInner.getEdge().setColor(innerColor);
				eOuter.getEdge().setColor(outerColor);

			}



		}


		//*/





		//printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeSliceMesh_Top(zObjMesh& o_Mesh, int startVID, int endVID, zIntArray& FeaturedNumStrides)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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
		vSIDSlice = coreUtils.getClosest_PointCloud(guidePts[startVID], positions);
		vEIDSlice = coreUtils.getClosest_PointCloud(guidePts[endVID], positions);
		zItMeshHalfEdge heStartSlice = getStartHalfEdge(*o_sliceMesh, vSIDSlice, vEIDSlice);
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
		zColor innerColor = _colorFeatureInner;
		zColor outerColor = _colorFeatureOuter;
		bool ignoreOuterEdge = false;


		//walk on the lower part of the mesh and color corner vertices
		zItMeshHalfEdge heLowerTemp = heStartSlice.getSym().getNext();
		int cornerCounter = 0;
		while (heLowerTemp.getVertex().getId() != heStartSlice.getStartVertex().getId())
		{
			zItMeshHalfEdge heS1, heS2;

			heS1 = heLowerTemp;
			heLowerTemp = heLowerTemp.getNext().getSym().getNext();
			heS2 = heLowerTemp;
			//get the angle between the two sides heS1 and heS2
			float angle = heS1.getVector().angle(heS2.getVector());
			bool chk = angle > 20;
			if (chk)
			{
				printf("\n angle0 %1.4f", angle);

				heLowerTemp.getStartVertex().setColor(zORANGE);
				cornerCounter++;
			}
		}
		int heCornerEdgeTopId = -1;
		zItMeshHalfEdge heCornerTopStart;

		for (int k = 0; k < featureWalkNum * 2; k++)
		{
			int featureCounter = k % FeaturedNumStrides.size();
			int blockStride = FeaturedNumStrides[featureCounter];
			if (blockStride == bridgingStride)
			{
				innerCounter = 0;
				outerCounter = 0;

				//isCornerEdge = heTemp.getStartVertex().getColor() == zORANGE;
			}
			isCornerEdge = heTemp.getStartVertex().getColor() == zORANGE;


			if (isCornerEdge)
			{
				innerColor = _colorCornersInterior;
				outerColor = _colorCornersBoundary;
				heTemp.getStartVertex().setColor(isBackSide ? innerColor : outerColor);

			}
			else
			{
				innerColor = _colorFeatureInner;
				outerColor = _colorFeatureOuter;
			}

			zColor edgeColor = isBackSide ? innerColor : outerColor;
			heTemp.getEdge().setColor(edgeColor);
			zItMeshHalfEdge heRest = heTemp;
			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				heRest = heRest.getNext().getSym().getNext();
				heRest.getEdge().setColor(edgeColor);
			}
			if (isCornerEdge)
			{
				heRest.getVertex().setColor(edgeColor);
			}
			if (heCornerEdgeTopId == -1 && isCornerEdge && isBackSide)
			{
				heCornerEdgeTopId = heRest.getId();
				printf("\n heCornerEdgeTopId %i ", heCornerEdgeTopId);
				heCornerTopStart = heRest;
				//heRest.getEdge().setColor(zRED);
			}

			for (int i = 0; i < blockStride; i++)
			{
				heTemp = heTemp.getSym().getNext().getNext();
			}
			if (isBackSide)
			{
				innerCounter++;
				hesFeatureInner.push_back(heTemp);
			}
			else
			{
				outerCounter++;
				hesFeatureOuter.push_back(heTemp);
			}
			if (blockStride == bridgingStride) isBackSide = !isBackSide;
		}


		vector<zItMeshHalfEdgeArray> hesTopSides, hesStartSides;

		////get first corner
		zItMeshHalfEdge heEdgeTop = heCornerTopStart.getSym().getPrev().getSym();

		int safetyCounter = 0;
		isBackSide = true;
		int tempCounter = 0;
		for (int i = 0; i < cornerCounter; i++)
		{
			zItMeshHalfEdgeArray hesSide, heStart;
			zColor startColor = heEdgeTop.getStartVertex().getColor();
			printf("\n startColor %1.4f | %1.4f | %1.4f ", startColor.r, startColor.g, startColor.b);
			zColor c0;
			while (true)
			{
				hesSide.push_back(heEdgeTop);

				heEdgeTop = heEdgeTop.getVertex().checkValency(3) ?
					heEdgeTop.getNext() : heEdgeTop.getNext().getSym().getNext();

				c0 = heEdgeTop.getStartVertex().getColor();
				if (c0 == _colorCornersInterior || c0 == _colorCornersBoundary)
				{
					if (!(c0 == startColor))
					{
						for (auto& h : hesSide)
						{
							h.getEdge().setColor(zBLACK);
							heStart.push_back(h);
						}
						hesStartSides.push_back(heStart);
					}
					else
					{
						for (auto& h : hesSide)
						{
							h.getEdge().setColor(startColor);
						}
					}
					hesTopSides.push_back(hesSide);

					break;
				}
				safetyCounter++;
				if (safetyCounter > fnMesh.numEdges())
				{
					printf("\n End 00 was not reached!");
					break;
				}
			}
			printf("\n hesTopSides %i | %i", hesTopSides.size(), hesSide.size());
		}
		int indexStart = -1;
		int minCount = INT_MAX;
		for (int i = 0; i < hesStartSides.size(); i++)
		{
			//get the smallest length of the side ones and color the faces of that one
			if (hesStartSides[i].size() < minCount)
			{
				indexStart = i;
				minCount = hesStartSides[i].size();
			}
		}

		for (auto& h : hesStartSides[indexStart])
		{
			zItMeshHalfEdge heTemp = h.getSym();
			zItMeshHalfEdge heRest = heTemp;
			for (int j = 0; j < sideDivisionCount - 1; j++)
			{
				heRest.getFace().setColor(zCYAN);
				heRest = heRest.getNext().getNext().getSym();
				//heRest.getEdge().setColor(edgeColor);
			}

			/*zItMeshFaceArray fs;
			h.getFaces(fs);
			for (auto& f : fs) f.setColor(zCYAN);*/
		}

		printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());
	}



	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computePrintBlockFrames(float printPlaneSpacing, float neopreneOffset_start, float neopreneOffset_end, bool leftBlock)
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

				float dStart = coreUtils.minDist_Point_Plane(pOnCurve, startPlanePoint, startPlaneNormal);

				if (abs(dStart) >= neopreneOffset_start)  right = true;


			}

			if (!left)
			{
				zPoint startPlanePoint = zVector(leftPlanes[0](3, 0), leftPlanes[0](3, 1), leftPlanes[0](3, 2));
				zPoint startPlaneNormal = zVector(leftPlanes[0](2, 0), leftPlanes[0](2, 1), leftPlanes[0](2, 2));

				float dStart = coreUtils.minDist_Point_Plane(pOnCurve, startPlanePoint, startPlaneNormal);
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

				float dEnd = coreUtils.minDist_Point_Plane(pOnCurve, endPlanePoint, endPlaneNormal);

				if (abs(dEnd) >= neopreneOffset_end)  right = true;

			}

			if (!left)
			{
				zPoint endPlanePoint = zVector(leftPlanes[1](3, 0), leftPlanes[1](3, 1), leftPlanes[1](3, 2));
				zPoint endPlaneNormal = zVector(leftPlanes[1](2, 0), leftPlanes[1](2, 1), leftPlanes[1](2, 2));

				float dEnd = coreUtils.minDist_Point_Plane(pOnCurve, endPlanePoint, endPlaneNormal);

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
			if (weight >= weights[l] && weight <= weights[l + 1])   mult = coreUtils.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
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
		zTransform pFrame = coreUtils.getTransformFromOrigin_Normal(O, tempZ, zVector(1, 1, 0));
		sectionFrames.push_back(pFrame);

		pOnCurve = O;
		dStart = coreUtils.minDist_Point_Plane(pOnCurve, startOrig, startNorm);


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
				if (weight >= weights[l] && weight <= weights[l + 1])   mult = coreUtils.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
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
			pFrame = coreUtils.getTransformFromOrigin_Normal(O, tempZ, zVector(1, 1, 0));
			sectionFrames.push_back(pFrame);

			pOnCurve = O;

			if (j == numLayers - 1)
			{
				float dEnd = coreUtils.minDist_Point_Plane(pOnCurve, endOrig, endNorm);

				if (leftBlock) printf(" \n left sD %1.4f eD %1.4f ", dStart, dEnd);
				else printf(" \n right sD %1.4f eD %1.4f ", dStart, dEnd);
			}
		}




		//supsitute frames in a new calculation way
		//zPlane startPlane = coreUtils.getPlaneFromOrigin_Normal(startOrig, startNorm);
		zPointArray origins;
		for (int i = 0; i < numLayers; i++)
		{
			origins.push_back(zVector(sectionFrames[i](3, 0), sectionFrames[i](3, 1), sectionFrames[i](3, 2)));
		}
		vector<zPlane> leftFrames, rightFrames;
		coreUtils.interpolatePlanes_slerp(leftPlanes[0], leftPlanes[1], numLayers, origins, leftFrames);
		if (!isRegular) coreUtils.interpolatePlanes_slerp(rightPlanes[0], rightPlanes[1], numLayers, origins, rightFrames);
		sectionFrames.clear();
		sectionFrames.insert(sectionFrames.end(), leftFrames.begin(), leftFrames.end());
		if (!isRegular) sectionFrames.insert(sectionFrames.end(), rightFrames.begin(), rightFrames.end());





	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computePrintBlockSections(bool left)
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
				float minDist_Plane = coreUtils.minDist_Point_Plane(P, O, N);
				scalars.push_back(minDist_Plane);
			}

			zPointArray positions;
			zIntArray edgeConnects;
			zColorArray vColors;
			int pres = 3;
			fn_sliceMesh.getIsoContour(scalars, 0.0, positions, edgeConnects, vColors, pres, pow(10, -1 * pres));

			// create graphs
			zFnGraph tempFn(o_sectionGraphs[i]);
			tempFn.create(positions, edgeConnects);;

			//tempFn.setEdgeColor(zGREY);
			tempFn.setVertexColors(vColors, false);

		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeSDF(bool allSDFLayers, int& numSDFlayers, int funcNum, int numSmooth, float printWidth, float neopreneOffset, float raftWidth)
	{

		o_contourGraphs.clear();
		o_contourGraphs.assign(o_sectionGraphs.size(), zObjGraph());

		o_trimGraphs.clear();
		o_trimGraphs.assign(o_sectionGraphs.size(), zObjGraph());

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
			//int kStart = (!blockPlanarCheck && !rightPlaneExists) ? 1 : 0;
			//int kEnd = (!blockPlanarCheck && !leftPlaneExists) ? 1 : 2;


			int kStart = (blockPlanarCheck) ? 0 : 1;
			int kEnd = (blockPlanarCheck) ? 2 : 1;

			printf("\n SDF kStart-kEnd %i - %i", kStart, kEnd);

			for (int k = kStart; k < kEnd; k++)
			{
				printf("\n SDF k %i", k);

				int i = (k == 0) ? j : j + end;

				if (!blockPlanarCheck) i = j;

				if (!blockPlanarCheck)
				{
					if (i == r0)
					{
						printf("\n SDF i %i", i);

						//raft 
					}
					else
					{
						//To be implemented
						printf("\n SDF i %i", i);

						//computeBlockSDF_NonPlanar(funcNum, numSmooth, i, (j % 2 == 0), printWidth, neopreneOffset, false, 0, raftWidth);

					}
				}
				else
				{
					if (i == r0)
					{
						printf("\n SDF i %i", i);

						//raft 
					}
					else if (i == r2)
					{
						printf("\n SDF i %i", i);

						//raft 
					}
					else
					{
						printf("\n SDF else");

						if (blockType == zBlockType::Wall)
						{
							printf("\n SDF Wall");

							computeBlockSDF_Planar_wall(funcNum, numSmooth, i, (j % 2 == 0), printWidth, neopreneOffset, false, 0, raftWidth);

						}
						else
						{
							printf("\n SDF bracing");
							computeBlockSDF_Planar_bracing(funcNum, numSmooth, i, (j % 2 == 0), printWidth, neopreneOffset, false, 0, raftWidth);
						}
					}
				}

			}
		}



	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computePrintBlockTrimGraphs(zObjGraph& inPolyObj, zObjGraph& o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE)
	{
		zFnGraph inFnGraph(inPolyObj);
		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPointArray gPositions;
		zIntArray gEdgeCOnnects;

		for (auto& he : bottomHE)
		{
			if (inColors[he.getStartVertex().getId()] == zBLUE) gPositions.push_back(inPositions[he.getStartVertex().getId()]);
			if (inColors[he.getVertex().getId()] == zMAGENTA) gPositions.push_back(inPositions[he.getVertex().getId()]);
			if (inColors[he.getVertex().getId()] == zORANGE) gPositions.push_back(inPositions[he.getVertex().getId()]);


		}

		for (auto& he : topHE)
		{
			if (inColors[he.getStartVertex().getId()] == zBLUE) gPositions.push_back(inPositions[he.getStartVertex().getId()]);
			if (inColors[he.getVertex().getId()] == zMAGENTA) gPositions.push_back(inPositions[he.getVertex().getId()]);
			if (inColors[he.getVertex().getId()] == zORANGE) gPositions.push_back(inPositions[he.getVertex().getId()]);
		}

		int end = floor(gPositions.size() * 0.5);
		//zDoubleArray eWeights;
		for (int i = 0; i < end; i++)
		{
			gEdgeCOnnects.push_back(i);
			gEdgeCOnnects.push_back(i + end);

			//(i == end - 1) ? eWeights.push_back(5) : eWeights.push_back(1);

		}

		zFnGraph fnOutGraph(o_outGraph);
		fnOutGraph.create(gPositions, gEdgeCOnnects);

		//fnOutGraph.setEdgeWeights(eWeights);

		//printf("\n trimGraph:  %i %i ", fnOutGraph.numVertices(), fnOutGraph.numEdges());

		//for (auto& p : gPositions) cout << endl << p;

	}

	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::checkPrintLayerHeights(bool& checkSDF, bool& checkGeometry)
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

				for (zItGraphVertex v(o_sectionGraphs[i]); !v.end(); v++)
				{
					zPoint p = v.getPosition();

					if (v.getColor() == zORANGE) orangeCounter++;
					if (v.getColor() == zMAGENTA) magentaCounter++;

					if (v.getColor() == _colorFeatureInner) innerSideCounter++;
					if (v.getColor() == _colorFeatureOuter) outerSideCounter++;


					zPoint p1 = p + norm * 1.0;

					zPoint intPt;
					bool check = coreUtils.line_PlaneIntersection(p, p1, prevNorm, prevOrigin, intPt);

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

		int precision = 3;
		bool checkLeftPlanar = checkInterfacePoints(true, pow(10, -1 * precision));
		bool checkRightPlanar = isRegular ? true : checkInterfacePoints(false, pow(10, -1 * precision));

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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::checkPrintLayerHeights_Folder(string folderDir, zDomainFloat& _printHeightDomain, zDomainFloat& _neopreneOffset, bool runBothPlanes, bool runPlaneLeft)
	{
		printHeightDomain = _printHeightDomain;
		neopreneOffset = _neopreneOffset;

		zStringArray files;
		coreUtils.getFilesFromDirectory(files, folderDir, zJSON);

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
			zStringArray split_0 = coreUtils.splitString(s, ".");
			zStringArray split_1 = coreUtils.splitString(split_0[split_0.size() - 2], "_");

			printf("\n File: %s ", s.c_str());
			int _blockID = atoi(split_1[split_1.size() - 1].c_str());

			setFromJSON(folderDir, _blockID, runBothPlanes, runPlaneLeft);
			Block_visitied[_blockID] = true;

			if (!planarBlock) continue;


			bool frameCHECKS = false;
			bool geomCHECKS = true;
			bool sdfCHECKS = true;

			bool maxHeight = true;
			bool minHeight = true;

			int minCriticalPtsCount = INT_MAX;
			float bestPlaneSpacing = FLT_MAX;

			for (float printPlaneSpacing = printHeightDomain.max; printPlaneSpacing >= printHeightDomain.min; printPlaneSpacing -= 0.00025)
			{
				if (planarBlock)
				{
					//printf("\n printPlaneSpace %1.4f ", printPlaneSpacing);
					computePrintSectionFromPlaneSpacing(printPlaneSpacing, printHeightDomain, _neopreneOffset, frameCHECKS, sdfCHECKS, geomCHECKS);
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
				else
				{

				}

			}

			if (!frameCHECKS)
			{
				printf("\n Layer height check FALSE!  Choose best plane spacing %1.4f", bestPlaneSpacing);

				computePrintSectionFromPlaneSpacing(bestPlaneSpacing, printHeightDomain, _neopreneOffset, frameCHECKS, sdfCHECKS, geomCHECKS);
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
			coreUtils.json_read(path, j);
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

				coreUtils.json_write(outPath, j);
			}


			float leftAngle = coreUtils.dihedralAngleBetweenPlanes(leftPlanes[0], leftPlanes[1]);
			float rightAngle = coreUtils.dihedralAngleBetweenPlanes(rightPlanes[0], rightPlanes[1]);

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

			printf("\n Finished block %i ", _blockID);





		}

		myfile.close();

		cout << " \n outPath exported : " << outFileName.c_str() << endl;

		for (int i = 0; i < Block_visitied.size(); i++)
		{
			if (!Block_visitied[i]) printf("\n %i ", i);
		}
	}

	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::checkInterfacePoints(bool left, float distTolerance)
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
			float d = coreUtils.minDist_Point_Plane(vPositions[vId], sPoint, sNorm);
			//printf("\n start %1.4f ", d);
			if (d > distTolerance) chkStart = false;
			if (!chkStart) break;
		}
		for (auto& vId : interfaceVerts_end)
		{
			float d = coreUtils.minDist_Point_Plane(vPositions[vId], ePoint, eNorm);
			//printf("\n end %1.4f ", d);
			if (d > distTolerance) chkEnd = false;
			if (!chkEnd) break;
		}
		bool out = (!chkStart || !chkEnd) ? false : true;
		if (!out) printf("\n checkGeom Fail! %s | start-End %s - %s ", left ? "Left" : "Right", chkStart ? "True" : "false", chkEnd ? "True" : "false");

		return out;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeBlockSDF_Planar_wall(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth)
	{
		if (graphId >= o_sectionGraphs.size())return;


		printf("\n fREP graphID %i  funcNum %i ", graphId, funcNum);

		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		float pWidth = (addRaft) ? printWidth : printWidth;
		float inc3dWidth = (addRaft) ? 0.025 : 0.025;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		fnGraph.setTransform(t, true, false);

		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));




		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_inner = 1.5 * pWidth;

		// Transform
		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// TOP BOTTOM HALF EDGES
		zItGraphHalfEdgeArray topHE, bottomHE;
		float topLength, bottomLength;


		slotGraph_1(sectionFrames[graphId], o_sectionGraphs[graphId], bottomLength, o_trimGraphs[graphId]);

		if (false)
		{
			polyTopBottomEdges(o_sectionGraphs[graphId], topHE, bottomHE, topLength, bottomLength);
			//trim graph
			computePrintBlockTrimGraphs(o_sectionGraphs[graphId], o_trimGraphs[graphId], topHE, bottomHE);
		}
		zFnGraph fnTrimGraph(o_trimGraphs[graphId]);

		// Profile polygon field
		zScalarArray polyField;
		if (funcNum >= 0) fnField.getScalars_Polygon(polyField, o_sectionGraphs[graphId], false);

		zScalarArray polyField_offset_outer, polyField_offset_inner;
		if (funcNum >= 1)
		{
			polyField_offset_outer = polyField;;
			for (auto& s : polyField_offset_outer) s += offset_outer;

			polyField_offset_inner = polyField;;
			for (auto& s : polyField_offset_inner) s += offset_inner;
		}

		zScalarArray slotField;
		float slotRadius = 0.02;
		zObjGraph slotGraph;
		if (funcNum >= 2) fnField.getScalarsAsEdgeDistance(slotField, o_trimGraphs[graphId], 0.20 * pWidth, false);

		zScalarArray patternField, patternField_1, patternField_2;
		float patternOffset = 0.025;
		float patternOffset2 = 0.015;
		float patternMove = offset_outer;
		if (funcNum >= 3)
		{
			getScalars_3dp_pattern_2(patternField_1, o_sectionGraphs[graphId], patternOffset, patternMove, true, graphId % 2 == 0);
			getScalars_3dp_pattern_2(patternField_2, o_sectionGraphs[graphId], patternOffset2, patternMove, true, graphId % 2 != 0);

			fnField.boolean_union(patternField_1, patternField_2, patternField, false);
		}


		zScalarArray booleanField_0;
		if (funcNum >= 4) fnField.boolean_subtract(polyField_offset_outer, polyField_offset_inner, booleanField_0, false);

		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_0, slotField, booleanField_1, false);

		zScalarArray booleanField_2;
		if (funcNum >= 6) fnField.boolean_union(booleanField_1, patternField, booleanField_2, false);


		float sdfWidth = printWidth / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(polyField_offset_outer, zFieldSDF, sdfWidth);
			break;


		case 2:

			fnField.setFieldValues(slotField, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(patternField, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(booleanField_0, zFieldSDF, sdfWidth);
			break;

		case 5:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(booleanField_2, zFieldSDF, sdfWidth);
			break;

		case 7:
			if (numSmooth > 0) fnField.smoothField(booleanField_1, numSmooth); // smooth field
			fnField.setFieldValues(booleanField_1, zFieldSDF, printWidth / 2.0);
			break;
		}

		/*for (zItMeshScalarField f(o_field); !f.end(); f++)
		{
			cout << f.getValue() << endl;
		}*/
		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, PRECISION, distanceTolerance);

		zFnGraph fngraph(o_contourGraphs[graphId]);


		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);


		// transform back 

		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		fnTrimGraph.setTransform(t, true, true);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeBlockSDF_Planar_bracing(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth)
	{
		printf("\n 0 fREP graphID %i  o_sectionGraphs.size() %i", graphId, o_sectionGraphs.size());

		if (graphId >= o_sectionGraphs.size())return;


		printf("\n fREP graphID %i  funcNum %i ", graphId, funcNum);

		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		float pWidth = (addRaft) ? printWidth : printWidth;
		float inc3dWidth = (addRaft) ? 0.025 : 0.025;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		fnGraph.setTransform(t, true, false);

		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));




		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_outer_2 = 0.25 * pWidth;
		float offset_inner = 1.5 * pWidth;

		// Transform
		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// TOP BOTTOM HALF EDGES
		zItGraphHalfEdgeArray topHE, bottomHE;

		zObjGraph slotGraph0;
		float graphLength0 = 0.2;
		slotGraph_Arch(graphId, graphLength0, slotGraph0);

		if (false)
		{
			//polyTopBottomEdges(o_sectionGraphs[graphId], topHE, bottomHE, topLength, bottomLength);
			//trim graph
			//computePrintBlockTrimGraphs(o_sectionGraphs[graphId], o_trimGraphs[graphId], topHE, bottomHE);
		}
		zFnGraph fnTrimGraph(o_trimGraphs[graphId]);

		// Profile polygon field
		zScalarArray polyField;
		zIntArray edgeId;
		if (funcNum >= 0) fnField.getScalars_Polygon(polyField, o_sectionGraphs[graphId], edgeId, false);

		zScalarArray polyField_offset_outer, polyField_offset_inner;

		//create a map of edge offset based on the edgeId
		zFloatArray outerOffsetArray;
		outerOffsetArray.assign(fnGraph.numEdges(), offset_outer);
		for (zItGraphEdge e(o_sectionGraphs[graphId]); !e.end(); e++)
		{
			zItGraphVertexArray vs;
			e.getVertices(vs);
			if (vs[0].getColor() == _colorCornersInterior && vs[1].getColor() == _colorCornersInterior)
			{
				outerOffsetArray[e.getId()] = offset_outer_2;
			}
		}

		if (funcNum >= 1)
		{
			polyField_offset_outer = polyField;;
			for (int sf = 0; sf < polyField_offset_outer.size(); sf++)
			{
				polyField_offset_outer[sf] += outerOffsetArray[edgeId[sf]];

			}
			/*for (auto& s : polyField_offset_outer)
			{
				s += offset_outer;
			}*/

			polyField_offset_inner = polyField;;
			for (auto& s : polyField_offset_inner) s += offset_inner;
		}


		zObjGraph slotGraph;
		float graphLength = 2 / pWidth;
		slotGraph_1(sectionFrames[graphId], o_sectionGraphs[graphId], graphLength, slotGraph);

		zScalarArray slotField;
		float slotRadius = 0.02;
		if (funcNum >= 2) fnField.getScalarsAsEdgeDistance(slotField, slotGraph, 0.20 * pWidth, false);

		zScalarArray patternField, patternField_1, patternField_2;
		float patternOffset = 0.025;
		float patternOffset2 = 0.015;
		float patternMove = offset_outer;
		/*if (funcNum >= 3)
		{
			getScalars_3dp_pattern_2(patternField_1, o_sectionGraphs[graphId], patternOffset,	patternMove, true, graphId % 2 == 0);
			getScalars_3dp_pattern_2(patternField_2, o_sectionGraphs[graphId], patternOffset2,	patternMove, true, graphId % 2 != 0);

			fnField.boolean_union(patternField_1, patternField_2, patternField, false);
		}*/


		//computeCableSectionPoints(graphId, )



		zScalarArray booleanField_0;
		if (funcNum >= 4) fnField.boolean_subtract(polyField_offset_outer, polyField_offset_inner, booleanField_0, false);

		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_subtract(booleanField_0, slotField, booleanField_1, false);

		zScalarArray booleanField_2;
		if (funcNum >= 6) fnField.boolean_union(booleanField_1, patternField, booleanField_2, false);


		float sdfWidth = printWidth / 2.0;
		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField, zFieldSDF, sdfWidth);
			break;

		case 1:
			fnField.setFieldValues(polyField_offset_outer, zFieldSDF, sdfWidth);
			break;

		case 2:

			fnField.setFieldValues(slotField, zFieldSDF, sdfWidth);
			break;

		case 3:
			fnField.setFieldValues(patternField, zFieldSDF, sdfWidth);
			break;

		case 4:
			fnField.setFieldValues(booleanField_0, zFieldSDF, sdfWidth);
			break;

		case 5:
			fnField.setFieldValues(booleanField_1, zFieldSDF, sdfWidth);
			break;

		case 6:
			fnField.setFieldValues(booleanField_2, zFieldSDF, sdfWidth);
			break;

		case 7:
			if (numSmooth > 0) fnField.smoothField(booleanField_1, numSmooth); // smooth field
			fnField.setFieldValues(booleanField_1, zFieldSDF, printWidth / 2.0);
			break;
		}

		/*for (zItMeshScalarField f(o_field); !f.end(); f++)
		{
			cout << f.getValue() << endl;
		}*/
		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0, PRECISION, distanceTolerance);

		zFnGraph fngraph(o_contourGraphs[graphId]);


		printf("\n o_contourGraphs[%i] : nV - nE %i - %i ", graphId, fngraph.numVertices(), fngraph.numEdges());
		fnIsoGraph.setEdgeWeight(2);


		// transform back 

		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		fnTrimGraph.setTransform(t, true, true);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeBlockSDF_NonPlanar(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth)
	{
		if (graphId >= o_sectionGraphs.size())return;


		printf("\n fREP graphID %i  funcNum %i ", graphId, funcNum);

		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		float pWidth = (addRaft) ? printWidth : printWidth;
		float inc3dWidth = (addRaft) ? 0.025 : 0.025;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		fnGraph.setTransform(t, true, false);

		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));

		// TOP BOTTOM HALF EDGES
		zItGraphHalfEdgeArray topHE, bottomHE;
		float topLength, bottomLength;
		if (funcNum >= 2) polyTopBottomEdges(o_sectionGraphs[graphId], topHE, bottomHE, topLength, bottomLength);


		// Transform
		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_inner = 1.5 * pWidth;

		//trim graph
		zFnGraph fnTrimGraph;
		if (funcNum >= 2)
		{
			computePrintBlockTrimGraphs(o_sectionGraphs[graphId], o_trimGraphs[graphId], topHE, bottomHE);
			fnTrimGraph = zFnGraph(o_trimGraphs[graphId]);
		}



		// Profile polygon field
		zScalarArray polyField;
		if (funcNum >= 0) fnField.getScalars_Polygon(polyField, o_sectionGraphs[graphId], false);

		zScalarArray polyField_offset_outer;
		if (funcNum >= 1)
		{
			polyField_offset_outer = polyField;;
			for (auto& s : polyField_offset_outer) s += offset_outer;
		}

		zScalarArray braceField;
		if (funcNum >= 2) getScalars_3dp_brace(braceField, o_trimGraphs[graphId], 0, 0.23 * pWidth, alternate);

		zScalarArray trimField;
		if (funcNum >= 3) getScalars_3dp_trim(trimField, o_trimGraphs[graphId], 0.20 * pWidth, alternate);

		zScalarArray booleanField_0;
		if (funcNum >= 4) fnField.boolean_subtract(polyField_offset_outer, braceField, booleanField_0, false);

		zScalarArray booleanField_1;
		if (funcNum >= 5) fnField.boolean_union(booleanField_0, trimField, booleanField_1, false);

		// RESULT FIELDS
		switch (funcNum)
		{
		case 0:
			fnField.setFieldValues(polyField);
			break;

		case 1:
			fnField.setFieldValues(polyField_offset_outer);
			break;

		case 2:
			fnField.setFieldValues(braceField);
			break;

		case 3:
			fnField.setFieldValues(trimField);
			break;

		case 4:
			fnField.setFieldValues(booleanField_0);
			break;

		case 5:
			if (numSmooth > 0) fnField.smoothField(booleanField_1, numSmooth);
			fnField.setFieldValues(booleanField_1);
			break;
		}



		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0);

		fnIsoGraph.setEdgeWeight(2);

		// transform back 
		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		if (funcNum >= 2) fnTrimGraph.setTransform(t, true, true);
	}

	ZSPACE_TOOLSETS_INLINE bool zTsNatpowerSDF::exportJSON(string pathCurrent, string dir, string filename, float printLyerWidth, float raftLayerWidth)
	{
		zFnMesh fn;
		json j;
		string blockPath = pathCurrent + "blockMesh_" + to_string(blockId) + ".json";

		bool fileChk = fn.json_read(blockPath, j);

		if (!fileChk) return false;

		string folderName = dir + "/" + to_string(blockId);
		_mkdir(folderName.c_str());
		for (const auto& entry : std::filesystem::directory_iterator(folderName)) std::filesystem::remove_all(entry.path());

		// EXPORT	BlockMesh

		string blockID_padded = coreUtils.getPaddedIndexString(blockId, 3);

		string blockMeshName = folderName + "/block_" + blockID_padded + ".json";

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
		string leftMeshName = folderName + "/block_left_" + blockID_padded + ".json";
		string rightMeshName = folderName + "/block_right_" + blockID_padded + ".json";

		zFnMesh fnMeshLeft(o_SliceMesh_Left);
		fnMeshLeft.to(leftMeshName, zJSON);

		zFnMesh fnMeshRight(o_SliceMesh_Right);
		fnMeshRight.to(rightMeshName, zJSON);

		// EXPORT graph

		float minLayerHeight = 10;
		float maxLayerHeight = 0;

		int r0 = 0;
		int r1 = floor(o_contourGraphs.size() * 0.5) - 1;
		int r2 = floor(o_contourGraphs.size() * 0.5);
		int r3 = (o_contourGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		//bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;
		bool deckBlock = true;

		int end = (deckBlock) ? floor(o_contourGraphs.size() * 0.5) : o_contourGraphs.size();

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

				string graphID_padded = coreUtils.getPaddedIndexString(j, 3);

				// trim graph export
				zFnGraph fnTrimGraph(o_trimGraphs[i]);
				if (fnTrimGraph.numVertices() == 0) continue;

				zFnGraph fnIsoGraph(o_contourGraphs[i]);
				if (fnIsoGraph.numVertices() == 0) continue;

				string outName1 = folderName;
				outName1 += "/";
				outName1 += "trim";
				outName1 += "_";
				outName1 += blockID_padded + "_" + graphID_padded + "_" + to_string(k) + ".json";

				fnTrimGraph.to(outName1, zJSON);

				// graph export

				string outName = folderName;
				outName += "/";
				outName += filename;
				outName += "_";
				outName += blockID_padded + "_" + graphID_padded + "_" + to_string(k) + ".json";
				fnIsoGraph.to(outName, zJSON);

				// read existing data in the json 
				json jSON;

				ifstream in_myfile;
				in_myfile.open(outName.c_str());

				int lineCnt = 0;

				if (in_myfile.fail())
				{
					cout << " error in opening outPath  " << outName.c_str() << endl;
				}

				in_myfile >> jSON;
				in_myfile.close();

				// CREATE JSON FILE
				zUtilsJsonHE graphJSON;


				graphJSON.vertexAttributes.clear();
				graphJSON.vertexAttributes = (jSON["VertexAttributes"].get<vector<vector<double>>>());



				// attribute export
				zVector norm(sectionFrames[i](2, 0), sectionFrames[i](2, 1), sectionFrames[i](2, 2));
				norm *= -1;
				norm.normalize();

				zVector prevNorm(sectionFrames[i - 1](2, 0), sectionFrames[i - 1](2, 1), sectionFrames[i - 1](2, 2));
				zVector prevOrigin(sectionFrames[i - 1](3, 0), sectionFrames[i - 1](3, 1), sectionFrames[i - 1](3, 2));

				for (zItGraphVertex v(o_contourGraphs[i]); !v.end(); v++)
				{
					vector<double> v_attrib = graphJSON.vertexAttributes[v.getId()];


					zPoint p = v.getPosition();

					zPoint p1 = p + norm * 1.0;

					zPoint intPt;
					bool check = coreUtils.line_PlaneIntersection(p, p1, prevNorm, prevOrigin, intPt);

					float layerHeight = intPt.distanceTo(p);
					maxLayerHeight = (layerHeight > maxLayerHeight) ? layerHeight : maxLayerHeight;
					minLayerHeight = (layerHeight < minLayerHeight) ? layerHeight : minLayerHeight;


					v_attrib.push_back(norm.x);
					v_attrib.push_back(norm.y);
					v_attrib.push_back(norm.z);

					v_attrib.push_back(printLyerWidth);
					v_attrib.push_back(layerHeight);

					graphJSON.vertexAttributes[v.getId()] = v_attrib;

				}

				// vertex sequence
				zIntArray vSequence;
				zItGraphVertexArray vArray;

				for (zItGraphVertex v(o_contourGraphs[i]); !v.end(); v++)
				{
					if (!v.checkValency(2))
					{
						vArray.push_back(v);
					}
				}

				printf("\n %s - valence 2 verts  %i", graphID_padded.c_str(), vArray.size());

				if (vArray.size() == 2)
				{
					zItGraphHalfEdge he = vArray[0].getHalfEdge();
					vSequence.push_back(vArray[0].getId());

					do
					{

						vSequence.push_back(he.getVertex().getId());

						he = he.getNext();
					} while (he.getVertex() != vArray[1]);

					vSequence.push_back(vArray[1].getId());
					vSequence.push_back(vArray[0].getId());

				}

				if (vArray.size() == 0)
				{
					zItGraphHalfEdge he(o_contourGraphs[i], 0);

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

				// Json file 

				jSON["VertexAttributes"] = graphJSON.vertexAttributes;
				jSON["VertexSequence"] = vSequence;

				// EXPORT	
				ofstream myfile;
				myfile.open(outName.c_str());

				if (myfile.fail())
				{
					cout << " error in opening outPath  " << outName.c_str() << endl;
					return false;
				}

				//myfile.precision(16);
				myfile << jSON.dump();
				myfile.close();


				// Print Length
				for (zItGraphEdge e(o_contourGraphs[i]); !e.end(); e++)
				{
					printLength += e.getLength();
				}


			}
		}

		printf("\n block| %1.4f %1.4f| %1.1f ", minLayerHeight, maxLayerHeight, printLength);



		return true;
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeClosestPointToGradientMesh(zPointArray& inPoints, zIntArray& faceIDs, zPointArray& closestPoints)
	{
		//MatrixXd P(inPoints.size(), 3);
		//for (int i = 0; i < inPoints.size(); i++)
		//{
		//	P(i, 0) = inPoints[i].x;
		//	P(i, 1) = inPoints[i].y;
		//	P(i, 2) = inPoints[i].z;
		//}
		//

		//VectorXd sqrD;
		//VectorXi I;
		//MatrixXd C;
		//igl::point_mesh_squared_distance(P, gradientTriMesh_V, gradientTriMesh_FTris, sqrD, I, C);

		//faceIDs.clear();
		//closestPoints.clear();

		//faceIDs.assign(inPoints.size(), int());
		//closestPoints.assign(inPoints.size(), zPoint());
		//
		////cout << endl << endl;

		//for (int i = 0; i< inPoints.size(); i++)
		//{
		//	faceIDs[i] = I(i);

		//	closestPoints[i].x = C(i, 0);
		//	closestPoints[i].y = C(i, 1);
		//	closestPoints[i].z = C(i, 2);		

		//	//cout << endl << inPoints[i] << ", " << closestPoints[i] << ", " << faceIDs[i];
		//}

		//


	}

	ZSPACE_TOOLSETS_INLINE float zTsNatpowerSDF::computeWeightedGradientValue(int& faceID, zPoint& closestPt)
	{
		zItMeshFace f(o_gradientTriMesh, faceID);

		zItMeshVertexArray fVerts;
		f.getVertices(fVerts);

		zPointArray fVPositions;
		zColorArray fVColors;

		for (auto& v : fVerts)
		{
			fVPositions.push_back(v.getPosition());
			fVColors.push_back(v.getColor());
		}

		zDoubleArray weights;
		coreUtils.getDistanceWeights(closestPt, fVPositions, 2.0, weights);

		//printf("\n weights %i ", weights.size());

		float out = 0;
		double w = 0;
		for (int i = 0; i < weights.size(); i++)
		{
			out += fVColors[i].r * weights[i];
			w += weights[i];
		}

		out = (out == 0.0) ? 0.0 : out / w;
		out = (out > 1.0) ? 1.0 : out;

		return out;
	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeCableSectionPoints(int graphId, zObjGraph& o_cableGraph, zPointArray& outputPoints)
	{
		zFnGraph fnSection(o_sectionGraphs[graphId]);
		zPointArray graphPts;
		fnSection.getVertexPositions(graphPts);
		zPlane frame = sectionFrames[graphId];
		zVector normal(frame(2, 0), frame(2, 1), frame(2, 2));

		zPointArray intersectionPts;
		intersect_graphPlane(o_cableGraph, frame, false, intersectionPts);

		outputPoints.clear();
		if (intersectionPts.size() > 0)
		{
			for (auto& p : intersectionPts)
			{
				if (coreUtils.pointInPlanarPolygon(p, graphPts, normal))
				{
					outputPoints.push_back(p);
				}
			}
		}
		if (outputPoints.size() > 0)
		{
			printf("\n cableGraph intersection - number of intersection %i ", outputPoints.size());

		}

		//if (intersectionEvents.Count() > 0)
		//{
		//	fnField.getScalars_Circle(scalars, intersectionPoint, radius);

		//}
	}
	ZSPACE_TOOLSETS_INLINE int zTsNatpowerSDF::getCableGraphIndexPerGraph(int graphId)
	{
		//iterate through all the graphs to get the graph the block belongs to
		int index = -1;
		int counter = 0;
		for (int i = 0; i < o_CableGraphs.size(); i++)
		{
			zPointArray intersectionPts;
			computeCableSectionPoints(graphId, o_CableGraphs[i], intersectionPts);
			if (intersectionPts.size() > 0)
			{
				index = i;
				counter++;
			}
		}

		if (counter != 1) index = -1;

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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::computeGraphEdgesForSlot_Arch()
	{
		//iterate through all the section graphs, and save the one u want
	}


	ZSPACE_TOOLSETS_INLINE zPoint zTsNatpowerSDF::getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh)
	{
		float scaleVal = coreUtils.ofMap(threshold, thresholdLow, thresholdHigh, 0.0000f, 1.0000f);
		zVector e = vertex_higher - vertex_lower;
		double edgeLen = e.length();
		e.normalize();
		return (vertex_lower + (e * edgeLen * scaleVal));
	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints)
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
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints)
	{
		outPoints.clear();
		zScalarArray vertexScalars;
		zPoint O(inPlane(3, 0), inPlane(3, 1), inPlane(3, 2));
		zVector N(inPlane(2, 0), inPlane(2, 1), inPlane(2, 2));
		for (zItGraphVertex v(o_graph); !v.end(); v++)
		{
			zPoint P = v.getPosition();
			float minDist_Plane = coreUtils.minDist_Point_Plane(P, O, N);
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





	//---- PROTECTED UTILITY METHODS


	ZSPACE_TOOLSETS_INLINE zItMeshHalfEdge zTsNatpowerSDF::getStartHalfEdge(zObjMesh& o_Mesh, int startVID, int endVID)
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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::createGraphFromHEArray(zItGraphHalfEdgeArray& heArray, zObjGraph& outGraph)
	{
		zPointArray positions;
		zIntArray eConnects;

		for (int i = 0; i < heArray.size(); i++)
		{
			positions.push_back(heArray[i].getStartVertex().getPosition());

			if (positions.size() > 1)
			{
				eConnects.push_back(positions.size() - 2);
				eConnects.push_back(positions.size() - 1);
			}

			if (i == heArray.size() - 1)
			{
				positions.push_back(heArray[i].getVertex().getPosition());

				eConnects.push_back(positions.size() - 2);
				eConnects.push_back(positions.size() - 1);
			}
		}

		zFnGraph fnGraph(outGraph);
		fnGraph.create(positions, eConnects);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getPerpendicularVector(zPlane& plane, zVector edgeVector, zPoint midPoint, float graphLength, zObjGraph& outGraph)
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

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength)
	{
		zFnGraph inFnGraph(inPoly);
		inFnGraph.setEdgeColor(zGREY);

		bottomHE.clear();
		topHE.clear();

		zVector Y(0, 1, 0);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zIntArray startVerts;
		for (int i = 0; i < inFnGraph.numVertices(); i++)
		{
			if (inColors[i] == zBLUE) startVerts.push_back(i);;
		}

		// Flip if needed, to have lower point in 0 index
		if (inPositions[startVerts[1]].z < inPositions[startVerts[0]].z)
		{
			int temp = startVerts[0];

			startVerts[0] = startVerts[1];
			startVerts[1] = temp;
		}

		zItGraphHalfEdge tHe;
		zItGraphHalfEdge bHe;

		for (int i = 0; i < 2; i++)
		{
			zItGraphVertex v(inPoly, startVerts[i]);

			zItGraphHalfEdgeArray cHEdges;
			v.getConnectedHalfEdges(cHEdges);

			if (i == 0)
			{
				//bHe = (cHEdges[0].getVertex().getId() == startVerts[1]) ? cHEdges[1] : cHEdges[0];

				zVector he = cHEdges[0].getVector();
				he.normalize();
				bHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
			}

			if (i == 1)
			{
				//tHe = (cHEdges[0].getVertex().getId() == startVerts[0]) ? cHEdges[1] : cHEdges[0];

				zVector he = cHEdges[0].getVector();
				he.normalize();
				tHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
			}

		}

		bool exit = false;
		while (!exit)
		{
			bottomHE.push_back(bHe);
			bottomLength += bHe.getLength();

			bHe.getEdge().setColor(zColor(0, 0, 1, 1));

			zItGraphVertex v = bHe.getVertex();

			if (v.getColor() == zORANGE) exit = true;

			bHe = bHe.getNext();
		}

		exit = false;
		while (!exit)
		{
			topHE.push_back(tHe);
			topLength += tHe.getLength();

			tHe.getEdge().setColor(zColor(1, 0, 0, 1));
			zItGraphVertex v = tHe.getVertex();
			if (v.getColor() == zORANGE) exit = true;

			tHe = tHe.getNext();
		}


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::slotGraph_0(zObjGraph& inPoly, zObjGraph& innerHE, float& innerLength)
	{
		printf("\n poly 0");
		//poly is sectionGraph -> we don't have to specify right/left
		zFnGraph inFnGraph(inPoly);
		inFnGraph.setEdgeColor(zGREY);


		zVector Y(0, 1, 0);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zItGraphVertex startV, endV;
		for (zItGraphVertex v(inPoly); !v.end(); v++)
		{
			if (v.getColor() == _colorFeatureInner)
			{

				startV = v;
				//printf("\n  start  found %1.4f, %1.4, %1.4f", startV.getPosition().x, startV.getPosition().y, startV.getPosition().z);

			}
			if (v.getColor() == _colorFeatureOuter)
			{
				endV = v;
				//printf("\n  end  found %1.4f, %1.4, %1.4f", endV.getPosition().x, endV.getPosition().y, endV.getPosition().z);

			}
			/*if (startV.getRawPosition() && endV.getRawPosition())
			{
				printf("\n both start and end found");
				break;
			}*/
		}

		zVector vector = endV.getPosition() - startV.getPosition();
		float length = vector.length() / 2;
		vector.normalize();

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV.getPosition() - (vector * (length / 2)));
		gPts.push_back(startV.getPosition() + (vector * (length)));
		gEdges.push_back(0);
		gEdges.push_back(1);

		zObjGraph oGraph;
		zFnGraph fnG(innerHE);
		fnG.create(gPts, gEdges);
		innerLength = length;

		//// Flip if needed, to have lower point in 0 index
		//if (inPositions[startVerts[1]].z < inPositions[startVerts[0]].z)
		//{
		//	int temp = startVerts[0];

		//	startVerts[0] = startVerts[1];
		//	startVerts[1] = temp;
		//}

		//zItGraphHalfEdge oHe;
		//zItGraphHalfEdge iHe;
		//printf("\n poly startVerts %i", startVerts.size());
		//if (startVerts.size() <2)
		//{
		//	printf("\n poly startVerts FAIL %i", startVerts.size());
		//	throw;
		//}
		//

		//for (int i = 0; i < 2; i++)
		//{
		//	zItGraphVertex v(inPoly, startVerts[i]);

		//	zItGraphHalfEdgeArray cHEdges;
		//	v.getConnectedHalfEdges(cHEdges);
		//	printf("\n poly cHEdges %i", cHEdges.size());

		//	if (i == 0)
		//	{
		//		//bHe = (cHEdges[0].getVertex().getId() == startVerts[1]) ? cHEdges[1] : cHEdges[0];

		//		zVector he = cHEdges[0].getVector();
		//		he.normalize();
		//		iHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
		//	}

		//	if (i == 1)
		//	{
		//		//tHe = (cHEdges[0].getVertex().getId() == startVerts[0]) ? cHEdges[1] : cHEdges[0];

		//		zVector he = cHEdges[0].getVector();
		//		he.normalize();
		//		oHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
		//	}

		//}

		//bool exit = false;
		//while (!exit)
		//{
		//	innerHE.push_back(iHe);
		//	innerLength += iHe.getLength();

		//	iHe.getEdge().setColor(zColor(0, 0, 1, 1));

		//	zItGraphVertex v = iHe.getVertex();

		//	if (v.getColor() == zORANGE) exit = true;

		//	iHe = iHe.getNext();
		//}
		//printf("\n poly 5");

		//exit = false;
		//while (!exit)
		//{
		//	outerHE.push_back(oHe);
		//	outerLength += oHe.getLength();

		//	oHe.getEdge().setColor(zColor(1, 0, 0, 1));
		//	zItGraphVertex v = oHe.getVertex();
		//	if (v.getColor() == zORANGE) exit = true;

		//	oHe = oHe.getNext();
		//}


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::slotGraph_1(zPlane plane, zObjGraph& inPoly, float graphLength, zObjGraph& outGraph)
	{
		//poly is sectionGraph -> we don't have to specify right/left
		zFnGraph inFnGraph(inPoly);
		inFnGraph.setEdgeColor(zGREY);


		zVector Y(0, 1, 0);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPoint startV, endV;

		int startEdgeId = INT_MAX;
		int endEdgeId = INT_MAX;
		zVector edgeVector;

		for (zItGraphEdge e(inPoly); !e.end(); e++)
		{
			zItGraphVertexArray vs;
			e.getVertices(vs);
			if (vs[0].getColor() == _colorCornersBoundary && vs[1].getColor() == _colorCornersBoundary)
			{
				if (startEdgeId == INT_MAX)
				{
					startV = e.getCenter();
					startEdgeId = e.getId();
					edgeVector = e.getVector();
					printf("\n start slot found");
				}
				else if (blockType == zBlockType::Wall)
				{
					endV = e.getCenter();
					endEdgeId = e.getId();
					printf("\n end slot found");

					break;
				}

			}
			if (!blockType == zBlockType::Wall)
			{
				if (vs[0].getColor() == _colorCornersInterior && vs[1].getColor() == _colorCornersInterior)
				{
					endV = e.getCenter();
				}
			}


		}


		zVector planeNormal(plane(2, 0), plane(2, 1), plane(2, 2));
		//printf("\n \n slotGraph planeNormal %1.4f - %1.4f - %1.4f", planeNormal.x, planeNormal.y, planeNormal.z);
		//printf("\n slotGraph edgeVector %1.4f - %1.4f - %1.4f ", edgeVector.x, edgeVector.y, edgeVector.z);

		zVector vector = planeNormal ^ edgeVector;
		//printf("\n slotGraph vector %1.4f - %1.4f - %1.4f \n", vector.x, vector.y, vector.z);

		//float length = 0.02;
		vector.normalize();

		zPointArray gPts;
		zIntArray gEdges;
		gPts.push_back(startV - (vector * (graphLength / 2)));
		gPts.push_back(startV + (vector * (graphLength)));
		gEdges.push_back(0);
		gEdges.push_back(1);

		zObjGraph oGraph;
		zFnGraph fnG(outGraph);
		fnG.create(gPts, gEdges);

		//printf("\n \n slotGraph gPts %1.4f - %1.4f - %1.4f", gPts[0].x, gPts[0].y, gPts[0].z);
		//printf("\n slotGraph gPts %1.4f - %1.4f - %1.4f \n", gPts[1].x, gPts[1].y, gPts[1].z);

		//// Flip if needed, to have lower point in 0 index
		//if (inPositions[startVerts[1]].z < inPositions[startVerts[0]].z)
		//{
		//	int temp = startVerts[0];

		//	startVerts[0] = startVerts[1];
		//	startVerts[1] = temp;
		//}

		//zItGraphHalfEdge oHe;
		//zItGraphHalfEdge iHe;
		//printf("\n poly startVerts %i", startVerts.size());
		//if (startVerts.size() <2)
		//{
		//	printf("\n poly startVerts FAIL %i", startVerts.size());
		//	throw;
		//}
		//

		//for (int i = 0; i < 2; i++)
		//{
		//	zItGraphVertex v(inPoly, startVerts[i]);

		//	zItGraphHalfEdgeArray cHEdges;
		//	v.getConnectedHalfEdges(cHEdges);
		//	printf("\n poly cHEdges %i", cHEdges.size());

		//	if (i == 0)
		//	{
		//		//bHe = (cHEdges[0].getVertex().getId() == startVerts[1]) ? cHEdges[1] : cHEdges[0];

		//		zVector he = cHEdges[0].getVector();
		//		he.normalize();
		//		iHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
		//	}

		//	if (i == 1)
		//	{
		//		//tHe = (cHEdges[0].getVertex().getId() == startVerts[0]) ? cHEdges[1] : cHEdges[0];

		//		zVector he = cHEdges[0].getVector();
		//		he.normalize();
		//		oHe = (abs(he * Y) > 0.9) ? cHEdges[0] : cHEdges[1];
		//	}

		//}

		//bool exit = false;
		//while (!exit)
		//{
		//	innerHE.push_back(iHe);
		//	innerLength += iHe.getLength();

		//	iHe.getEdge().setColor(zColor(0, 0, 1, 1));

		//	zItGraphVertex v = iHe.getVertex();

		//	if (v.getColor() == zORANGE) exit = true;

		//	iHe = iHe.getNext();
		//}
		//printf("\n poly 5");

		//exit = false;
		//while (!exit)
		//{
		//	outerHE.push_back(oHe);
		//	outerLength += oHe.getLength();

		//	oHe.getEdge().setColor(zColor(1, 0, 0, 1));
		//	zItGraphVertex v = oHe.getVertex();
		//	if (v.getColor() == zORANGE) exit = true;

		//	oHe = oHe.getNext();
		//}


	}


	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::slotGraph_Arch(int graphId, float graphLength, zObjGraph& outGraph)
	{
		//poly is sectionGraph -> we don't have to specify right/left
		zFnGraph fngraph(o_sectionGraphs[graphId]);

		//get graph edge that is between two colors. Then choose the one you want. Still not sure how
		vector<zItGraphHalfEdgeArray> innerHEs;
		vector<int> counts;
		int vertexFoundCounter = 0;
		zItGraphEdge startSlotE;
		bool found = false;
		cout << endl;
		//printf("\n (o_sectionGraphs[%i]) ne = %i ", graphId, fngraph.numEdges());

		//for (zItGraphHalfEdge he (o_sectionGraphs[graphId]); !he.end(); he++)
		//{
		//	//get the he where the start is inetior and end is boundary
		//	if (he.getStartVertex().getColor() == _colorCornersInterior && he.getVertex().getColor() == _colorCornersBoundary )
		//	{
		//		startSlotHE = he;
		//		printf("\n slot graph he found");
		//		found = true;

		//		break;
		//	}
		//}
		for (zItGraphEdge e(o_sectionGraphs[graphId]); !e.end(); e++)
		{
			//get the he where the start is inetior and end is boundary
			zItGraphVertexArray verts;
			e.getVertices(verts);
			if (verts[0].getColor() == _colorCornersInterior && verts[1].getColor() == _colorCornersBoundary)
			{
				startSlotE = e;
				//printf("\n slot graph he found");
				found = true;

				break;
			}
		}
		if (!found)
		{
			printf("\n slot graph, no edge was found! counter %i | %i | graphId = %i THROW!", fngraph.numEdges(), graphId);
			//throw;
			return;
		}

		//for (zItGraphVertex v(o_sectionGraphs[graphId]); !v.end(); v++)
		//{
		//	if (v.getColor() == _colorCornersInterior)
		//	{
		//		vertexFoundCounter++;
		//		//walk until you reach a vertex with color _colorFeatureOuter
		//		zItGraphHalfEdgeArray hes;
		//		v.getConnectedHalfEdges(hes);
		//		for (auto& he : hes)
		//		{
		//			zItGraphHalfEdge heStart = he;
		//			zItGraphHalfEdge he = heStart;
		//			zItGraphHalfEdgeArray innerHE;
		//			while (true)
		//			{
		//				innerHE.push_back(he);
		//				if (he.getVertex().getColor() == _colorCornersBoundary) break;
		//				he = he.getNext();
		//			}
		//			innerHEs.push_back(innerHE);
		//			counts.push_back(innerHE.size());
		//			printf("\n slotGraph_Arch index-size %i | %i", innerHEs.size(), innerHE.size());
		//			//innerLengths.push_back(length);
		//		}
		//	}
		//}
		////get the array with the smallest number of vertices (there might be another check needed, but this is for now)
		//if (innerHEs.size() == 0)
		//{
		//	printf("\n slotGraph_Arch no innerHEs found! THROW");
		//	throw;
		//}
		//int index = 0;
		//int minCount = INT_MAX;
		//for (int i = 0; i < innerHEs.size(); i++)
		//{
		//	if (counts[i] < minCount)
		//	{
		//		//minCount = innerHEs[i].size();
		//		minCount = counts[i];
		//		index = i;
		//	}
		//}
		//printf("\n slotGraph_Arch index %i | %i", index, innerHEs.size());
		//zItGraphHalfEdgeArray heEdge = innerHEs[index];
		/*for (auto& he : heEdge)
		{
			he.getEdge().setColor(zBLUE);
		}*/

		zObjGraph edgeGraph;
		zFnGraph fng(edgeGraph);
		/*zPointArray pts;
		pts.push_back(startSlotE.getStartVertex().getPosition());
		pts.push_back(startSlotE.getVertex().getPosition());*/
		//= { startSlotHE.getStartVertex().getPosition(), startSlotHE.getVertex().getPosition() };
		zPointArray pts;
		startSlotE.getVertexPositions(pts);

		zIntArray connect = { 0, 1 };
		fng.create(pts, connect);
		//createGraphFromHEArray(heEdge, edgeGraph);

		zPoint mid = fng.getCenter();
		zItGraphVertex v0(edgeGraph, 0);
		zItGraphVertex v1(edgeGraph, fng.numVertices() - 1);
		zVector vec = v1.getPosition() - v0.getPosition();

		getPerpendicularVector(sectionFrames[graphId], vec, mid, graphLength, outGraph);

		o_trimGraphs[graphId] = outGraph;

		//zFnGraph fnOutGraph(outGraph);
		//fnOutGraph.setEdgeColor(zBLUE);


	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_slot(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnTrimGraph(o_trimGraph);
		zPoint* trimPositions = fnTrimGraph.getRawVertexPositions();

		zPointArray tmpPositions;
		zPointArray gPositions;

		int end = floor(fnTrimGraph.numVertices() * 0.5);

		for (int i = 1; i < end - 1; i++)
		{

			zPoint p0 = (trimPositions[i] + trimPositions[i + end]) * 0.5;
			gPositions.push_back(p0);
		}

		fnField.getScalarsAsVertexDistance(scalars, gPositions, offset, false);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_pattern(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, bool alternate, bool alternateChk)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnSectionGraph(o_sectionGraph);
		zPoint* trimPositions = fnSectionGraph.getRawVertexPositions();

		zPointArray tmpPositions;
		zPointArray gPositions;

		bool chk = alternateChk;
		for (zItGraphVertex v(o_sectionGraph); !v.end(); v++)
		{
			//if (v.getColor() == _colorPattern || v.getColor() == _colorStartOuter)
			if (v.getColor() == _colorPattern)
			{
				if (alternate)
				{
					if (chk)
					{
						gPositions.push_back(v.getPosition());
					}
					chk = !chk;
				}
				else
				{
					gPositions.push_back(v.getPosition());
				}
			}

		}

		fnField.getScalarsAsVertexDistance(scalars, gPositions, offset, false);


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_pattern_2(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, float offset2, bool alternate, bool alternateChk)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnSectionGraph(o_sectionGraph);
		zPoint* trimPositions = fnSectionGraph.getRawVertexPositions();

		zPointArray tmpPositions;
		zPointArray gPositions;

		bool chk = alternateChk;
		for (zItGraphVertex v(o_sectionGraph); !v.end(); v++)
		{
			bool colorChk = blockType == zBlockType::Wall ? v.getColor() == _colorPattern || v.getColor() == _colorFeatureOuter : v.getColor() == _colorPattern;
			//if (v.getColor() == _colorPattern || v.getColor() == _colorStartOuter)
			if (v.getColor() == _colorPattern)
			{
				if (alternate)
				{
					if (chk)
					{
						//get vector and normal
						zItGraphHalfEdgeArray hes;
						v.getConnectedHalfEdges(hes);

						zVector v0 = hes[0].getVector();
						zVector v1 = hes[1].getVector();

						v0.normalize();
						v1.normalize();

						zVector n = (v0 + v1) / 2.0;
						n.normalize();
						n *= offset2;

						zPoint p = v.getPosition() - n;
						gPositions.push_back(p);
					}
					chk = !chk;
				}
				else
				{
					gPositions.push_back(v.getPosition());
				}
			}

		}

		fnField.getScalarsAsVertexDistance(scalars, gPositions, offset, false);


	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_pattern_3(zScalarArray& scalars, zObjGraph& o_sectionGraph, float offset, float offset2, bool alternate, bool alternateChk)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnSectionGraph(o_sectionGraph);
		zPoint* trimPositions = fnSectionGraph.getRawVertexPositions();

		zPointArray tmpPositions;
		zPointArray gPositions;

		bool chk = alternateChk;

		for (zItGraphVertex v(o_sectionGraph); !v.end(); v++)
		{
			bool colorChk = blockType == zBlockType::Wall ? v.getColor() == _colorPattern || v.getColor() == _colorFeatureOuter : v.getColor() == _colorPattern;
			//if (v.getColor() == _colorPattern || v.getColor() == _colorStartOuter)
			if (v.getColor() == _colorPattern)
			{
				if (alternate)
				{
					if (chk)
					{
						//get vector and normal
						zItGraphHalfEdgeArray hes;
						v.getConnectedHalfEdges(hes);

						zVector v0 = hes[0].getVector();
						zVector v1 = hes[1].getVector();

						v0.normalize();
						v1.normalize();

						zVector n = (v0 + v1) / 2.0;
						n.normalize();
						n *= offset2;

						zPoint p = v.getPosition() - n;
						gPositions.push_back(p);
					}
					chk = !chk;
				}
				else
				{
					gPositions.push_back(v.getPosition());
				}
			}

		}

		fnField.getScalarsAsVertexDistance(scalars, gPositions, offset, false);


	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_patternMesh(zScalarArray& scalars, zObjGraph& o_sectionGraph, zIntArray& faceIDs, zPointArray& closestPoints, zDomainFloat& offsetDomain)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnSectionGraph(o_sectionGraph);
		zPoint* trimPositions = fnSectionGraph.getRawVertexPositions();

		zPointArray tmpPositions;
		zPointArray gPositions;

		for (zItGraphVertex v(o_sectionGraph); !v.end(); v++)
		{
			if (v.getColor() == zMAGENTA)
			{
				gPositions.push_back(v.getPosition());
			}

		}

		//zScalarArray scalars_temp, booleanField;

		//printf("\n gPositions %i faceIDs %i ", gPositions.size(), faceIDs.size());


		zFnMesh fnGradMesh(o_gradientTriMesh);
		zColorArray gradColors;
		fnGradMesh.getFaceColors(gradColors);
		zDomainFloat inDomain(0.0, 1.0);

		zColor* vColors = fnGradMesh.getRawVertexColors();

		//printf("\n numM %i ", gPositions.size());

		//for (int i = 0; i < gPositions.size(); i++)
		//{
		//	float inVal = gradColors[faceIDs[i]].r;
		//	float offset = coreUtils.ofMap(inVal, inDomain, offsetDomain);

		//	printf("\n %i | %i | %1.2f | %1.2f ",i, faceIDs[i], gradColors[faceIDs[i]].r, offset);

		//	/*if (i == 0)
		//	{
		//		fnField.getScalars_Circle(scalars, gPositions[i], offset, 0, false);
		//	}
		//	else
		//	{*/
		//		scalars_temp.clear();
		//		booleanField.clear();
		//		fnField.getScalars_Circle(scalars_temp, gPositions[i], offset, 0, false);

		//		fnField.boolean_union(booleanField, scalars, scalars_temp, false);

		//		//scalars.clear();
		//		scalars = booleanField;
		//	//}
		//}

		scalars.clear();
		zPoint* fieldMeshPositions = fnField.fnMesh.getRawVertexPositions();


		zFloatArray offsets;
		for (int j = 0; j < faceIDs.size(); j++)
		{
			//weight from face
			float offset = computeWeightedGradientValue(faceIDs[j], closestPoints[j]);
			offsets.push_back(offset);

			// closest point of the face
			//zItMeshFace f(o_gradientTriMesh, faceIDs[j]);
			//
			//zPointArray fVPositions;
			//f.getVertexPositions(fVPositions);

			//zIntArray fVerts;
			//f.getVertices(fVerts);


			//int vId = coreUtils.getClosest_PointCloud(closestPoints[j], fVPositions);
			////printf("\n %i %1.2f ", fVerts[vId], vColors[fVerts[vId]].r);


			//offsets.push_back(vColors[fVerts[vId]].r);
		}


		for (int i = 0; i < fnField.fnMesh.numVertices(); i++)
		{
			double d = 0.0;
			double tempDist = 10000;

			for (int j = 0; j < gPositions.size(); j++)
			{
				double r = fieldMeshPositions[i].distanceTo(gPositions[j]);

				//float inVal = gradColors[faceIDs[j]].r;

				float inVal = offsets[j];

				float offset = coreUtils.ofMap(inVal, inDomain, offsetDomain);
				//printf("\n %1.2 %1.2f ", offset, offsets[j]);

				r = r - offset;

				if (r < tempDist)
				{
					d = r;
					tempDist = r;
				}

			}

			scalars.push_back(d);
		}



	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_brace(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset, bool alternate)
	{
		zFnMeshScalarField fnField(o_field);


		zPointArray gPositions;
		zIntArray gEdgeCOnnects;

		zFnGraph fnTrimGraph(o_trimGraph);
		zPoint* trimPositions = fnTrimGraph.getRawVertexPositions();
		for (int i = 2; i < fnTrimGraph.numVertices() - 2; i++)
		{
			gPositions.push_back(trimPositions[i]);
		}

		int end = floor(fnTrimGraph.numVertices() * 0.5);

		for (int i = 1; i < end - 1; i++)
		{
			zPoint p0 = trimPositions[i];
			zPoint p1 = trimPositions[i + end];

			if (outer_printWidth != 0)
			{
				if (alternate)
				{
					zVector newE = p1 - p0;
					float newELen = newE.length();
					newE.normalize();

					float dist = outer_printWidth + offset;
					newELen -= dist;

					p1 = p0 + newE * newELen;


				}
				else
				{
					zVector newE = p0 - p1;
					float newELen = newE.length();
					newE.normalize();

					float dist = outer_printWidth + offset;
					newELen -= dist;

					p0 = p1 + newE * newELen;
				}
			}

			gPositions.push_back(p0);
			gPositions.push_back(p1);

			gEdgeCOnnects.push_back(gPositions.size() - 2);
			gEdgeCOnnects.push_back(gPositions.size() - 1);


		}


		zObjGraph o_tempGraph;
		zFnGraph fnTempGraph(o_tempGraph);
		fnTempGraph.create(gPositions, gEdgeCOnnects);

		fnField.getScalarsAsEdgeDistance(scalars, o_tempGraph, offset, false);

	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_brace2(zScalarArray& scalars, zObjGraph& poly, float outer_printWidth, float offset, bool alternate)
	{
		zFnMeshScalarField fnField(o_field);

		zFnGraph fnGraph(poly);



		//pick points on the inner as the tip of the triangle



	}
	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_trim(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate)
	{
		zFnMeshScalarField fnField(o_field);
		zIntArray tEdgeConnects;
		zPointArray tPositions;

		zFnGraph fnTrimGraph(o_trimGraph);
		zPoint* trimPositions = fnTrimGraph.getRawVertexPositions();

		int end = floor(fnTrimGraph.numVertices() * 0.5);
		zVector upVec(0, 0, 1);

		//printf("\n end %i ", end);

		int offsetMult = 2.15;

		for (int i = 1; i < end - 1; i++)
		{
			zPoint tmp0 = trimPositions[i];
			zPoint tmp1 = trimPositions[i + end];

			zVector e = tmp1 - tmp0;
			float eLen = e.length();
			e.normalize();

			zPoint midPt = tmp0 + e * eLen * 0.5;

			zPoint p0 = (alternate) ? midPt + e * (offset * offsetMult) : midPt - e * (offset * offsetMult);
			zVector dir = upVec ^ e;
			dir.normalize();

			tPositions.push_back(p0 + dir * offset * 1);
			tPositions.push_back(p0 - dir * offset * 1);

			tEdgeConnects.push_back(tPositions.size() - 2);
			tEdgeConnects.push_back(tPositions.size() - 1);

		}
		zObjGraph tempGraph;
		zFnGraph fnTempGraph(tempGraph);
		fnTempGraph.create(tPositions, tEdgeConnects);

		fnField.getScalarsAsEdgeDistance(scalars, tempGraph, offset, false);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_trimStart(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate)
	{
		zFnMeshScalarField fnField(o_field);
		zIntArray tEdgeConnects;
		zPointArray tPositions;

		zFnGraph fnTrimGraph(o_trimGraph);
		zPoint* trimPositions = fnTrimGraph.getRawVertexPositions();

		int end = floor(fnTrimGraph.numVertices() * 0.5);
		zVector upVec(0, 0, 1);



		//printf("\n end %i ", end);

		int offsetMult = 2.15;

		for (int i = 1; i < end - 1; i++)
		{
			zPoint tmp0 = trimPositions[i];
			zPoint tmp1 = trimPositions[i + end];

			zVector e = tmp1 - tmp0;
			float eLen = e.length();
			e.normalize();

			zPoint midPt = tmp0 + e * eLen * 0.5;

			zPoint p0 = (alternate) ? midPt + e * (offset * offsetMult) : midPt - e * (offset * offsetMult);
			zVector dir = upVec ^ e;
			dir.normalize();

			tPositions.push_back(p0 + dir * offset * 1);
			tPositions.push_back(p0 - dir * offset * 1);

			tEdgeConnects.push_back(tPositions.size() - 2);
			tEdgeConnects.push_back(tPositions.size() - 1);

		}
		zObjGraph tempGraph;
		zFnGraph fnTempGraph(tempGraph);
		fnTempGraph.create(tPositions, tEdgeConnects);

		fnField.getScalarsAsEdgeDistance(scalars, tempGraph, offset, false);
	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::getScalars_3dp_cable(zScalarArray& scalars, zPoint cablePoint, float radius)
	{
		zFnMeshScalarField fnField(o_field);
		scalars.clear();
		fnField.getScalars_Circle(scalars, cablePoint, radius);
	}




	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::readJSON_planarBlock(string path, int _blockID, bool runBothPlanes, bool runPlaneLeft)
	{
		printf("\n readJSON_planarBlock 0");

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
		//int blockStride = j["BlockStride"];



		isCorner = j["IsCorner"];
		//isFront = j["IsFront"];
		int type = j["BlockType"];
		blockType = static_cast<zBlockType>(type);




		//zFloatArray arrayL, eArrayR;
		zTransform sPlaneRight, ePlaneRight, sPlaneLeft, ePlaneLeft;

		//vector<zFloatArray> arrayL = j["LeftPlane"];
		//vector<zFloatArray> arrayR = j["RightPlane"];

		zPoint oS, oE, xS, xE, yS, yE, zS, zE;
		zVector zVectorEndRight, zVectorEndLeft;

		//zFloatArray sArrayL = j["leftPlanes"][0];
		//zFloatArray eArrayL = j["leftPlanes"][1];

		//sPlane = coreUtils.getPlaneFromArray(sArrayL);
		//ePlane = coreUtils.getPlaneFromArray(eArrayL);

		oS = zPoint(j["LeftPlanes"][0][3], j["LeftPlanes"][0][7], j["LeftPlanes"][0][11]);
		oE = zPoint(j["LeftPlanes"][1][3], j["LeftPlanes"][1][7], j["LeftPlanes"][1][11]);

		xS = zPoint(j["LeftPlanes"][0][0], j["LeftPlanes"][0][4], j["LeftPlanes"][0][8]);
		xE = zPoint(j["LeftPlanes"][1][0], j["LeftPlanes"][1][4], j["LeftPlanes"][1][8]);

		yS = zPoint(j["LeftPlanes"][0][1], j["LeftPlanes"][0][5], j["LeftPlanes"][0][9]);
		yE = zPoint(j["LeftPlanes"][1][1], j["LeftPlanes"][1][5], j["LeftPlanes"][1][9]);

		zS = zPoint(j["LeftPlanes"][0][2], j["LeftPlanes"][0][6], j["LeftPlanes"][0][10]);
		zE = zPoint(j["LeftPlanes"][1][2], j["LeftPlanes"][1][6], j["LeftPlanes"][1][10]);
		zVectorEndLeft = zE;
		sPlaneLeft = coreUtils.getPlaneFromVectors(oS, xS, yS, zS);
		ePlaneLeft = coreUtils.getPlaneFromVectors(oE, xE, yE, zE);



		//zFloatArray sArrayR = j["rightPlanes"][0];
		//zFloatArray eArrayR = j["rightPlanes"][1];
		//sPlane = coreUtils.getPlaneFromArray(sArrayR);
		//ePlane = coreUtils.getPlaneFromArray(eArrayR);

		oS = zPoint(j["RightPlanes"][0][3], j["RightPlanes"][0][7], j["RightPlanes"][0][11]);
		oE = zPoint(j["RightPlanes"][1][3], j["RightPlanes"][1][7], j["RightPlanes"][1][11]);
		xS = zPoint(j["RightPlanes"][0][0], j["RightPlanes"][0][4], j["RightPlanes"][0][8]);
		xE = zPoint(j["RightPlanes"][1][0], j["RightPlanes"][1][4], j["RightPlanes"][1][8]);
		yS = zPoint(j["RightPlanes"][0][1], j["RightPlanes"][0][5], j["RightPlanes"][0][9]);
		yE = zPoint(j["RightPlanes"][1][1], j["RightPlanes"][1][5], j["RightPlanes"][1][9]);
		zS = zPoint(j["RightPlanes"][0][2], j["RightPlanes"][0][6], j["RightPlanes"][0][10]);
		zE = zPoint(j["RightPlanes"][1][2], j["RightPlanes"][1][6], j["RightPlanes"][1][10]);
		zVectorEndRight = zE;

		sPlaneRight = coreUtils.getPlaneFromVectors(oS, xS, yS, zS);
		ePlaneRight = coreUtils.getPlaneFromVectors(oE, xE, yE, zE);

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
		computeMedialGraph(o_GuideMesh, sID, eID);


		//get feature stride
		zIntArray FeaturedNumStrides, FeaturedNumStrides2;
		coreUtils.json_readAttribute(j, "FeaturedNumStrides", FeaturedNumStrides);
		//coreUtils.json_readAttribute(j, "FeaturedNumStrides_2", FeaturedNumStrides2);


		//checkWall = blockType == zBlockType::Wall;


		/*
		zPoint startPoint = tmpPositions[sID];
		zPoint endPoint = tmpPositions[eID];

		zVector fNormStart(j["LeftPlanes"][0][0], j["LeftPlanes"][0][1], j["LeftPlanes"][0][2]);
		//fNormStart *= -1;
		zVector baseVector(1, 1, 0);
		zTransform sPlane = coreUtils.getTransformFromOrigin_Normal(startPoint, fNormStart, baseVector);

		zVector fNormEnd(j["LeftPlanes"][1][0], j["LeftPlanes"][1][1], j["LeftPlanes"][1][2]);
		zTransform ePlane = coreUtils.getTransformFromOrigin_Normal(endPoint, fNormEnd, baseVector);

		cout << "\n left sNorm " << fNormStart;
		cout << "\n left eNorm " << fNormEnd;

		setStartEndPlanes(sPlane, ePlane, true);

		base_world = sPlane;


		fNormStart = zVector(j["RightPlanes"][0][0], j["RightPlanes"][0][1], j["RightPlanes"][0][2]);
		//fNormStart *= -1;
		sPlane = coreUtils.getTransformFromOrigin_Normal(startPoint, fNormStart);

		 fNormEnd = zVector(j["RightPlanes"][1][0], j["RightPlanes"][1][1], j["RightPlanes"][1][2]);
		 ePlane = coreUtils.getTransformFromOrigin_Normal(endPoint, fNormEnd); ;

		cout << "\n right sNorm " << fNormStart;
		cout << "\n right eNorm " << fNormEnd;

		setStartEndPlanes(sPlane, ePlane, false);




		base_local.setIdentity();
		//medial axis
		computeMedialGraph(o_GuideMesh, sID, eID);

		*/

		//left mesh
		//if (leftPlaneExists)
		{
			// left mesh
			if (planarBlock)
			{
				//isRegular = blockType == zBlockType::Wall || blockType == zBlockType::Arch || (blockType == zBlockType::Top && !isCorner);
				isRegular = blockType == zBlockType::Wall || blockType == zBlockType::Arch || (blockType == zBlockType::Top);

				if (blockType != zBlockType::Top)
				{
					printf("\n slice mesh corner");
					StartCornerVID = j["StartCornerVID"];
					computeSliceMesh_Top(o_GuideMesh, sID, eID, FeaturedNumStrides);
				}
				else if (isRegular)
				{
					computeSliceMesh_Regular(o_GuideMesh, sID, eID, FeaturedNumStrides);
					printf("\n slice mesh regular");


				}
				else
				{

					//left Mesh
					computeSliceMesh(o_GuideMesh, sID, eID, FeaturedNumStrides, true);
					printf("\n slice mesh1");

					//Right Mesh
					computeSliceMesh(o_GuideMesh, sID, eID, FeaturedNumStrides, false);
					printf("\n slice mesh2");
				}


				//computeMedial_BraceEdges(o_SliceMesh_Left, 3, 0, blockStride, braceStride);
			}
			else
			{
				/*zFnMesh fnMesh_Left(o_SliceMesh_Left);
				fnMesh_Left.clear();
				fnMesh_Left.from(path, zJSON);

				computeMedial_BraceEdges(o_SliceMesh_Left, sID, eID, blockStride, braceStride);*/
			}

		}

		//right mesh
		//if (rightPlaneExists)
		{
			// right mesh
			if (planarBlock)
			{
				//computeSliceMesh(o_GuideMesh, sID, eID, blockStride, false);
				//computeMedial_BraceEdges(o_SliceMesh_Right, 0, 1, blockStride, braceStride);
			}
			else
			{
				/*zFnMesh fnMesh_Right(o_SliceMesh_Right);
				fnMesh_Right.clear();
				fnMesh_Right.from(path, zJSON);

				computeMedial_BraceEdges(o_SliceMesh_Right, sID, eID, blockStride, braceStride);*/
			}
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsNatpowerSDF::setCableGraph(string folderDir)
	{
		//read all files in directory
		zStringArray files;
		coreUtils.getFilesFromDirectory(files, folderDir, zJSON);
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


}


