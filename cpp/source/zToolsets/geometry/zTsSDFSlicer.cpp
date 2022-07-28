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


#include<headers/zToolsets/geometry/zTsSDFSlicer.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSDFSlicer::zTsSDFSlicer()
	{
		red = zColor(1, 0, 0, 1);
		yellow = zColor(1, 1, 0, 1); 
		green = zColor(0, 1, 0, 1);
		cyan = zColor(0, 1, 1, 1); 
		blue = zColor(0, 0, 1, 1); 
		magenta = zColor(1, 0, 1, 1);

		grey = zColor(0.5, 0.5, 0.5, 1);

		orange = zColor(1, 0.5, 0, 1);

		printHeightDomain = zDomainFloat(0.006, 0.012);
	
	}


	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSDFSlicer::~zTsSDFSlicer() {}

	//---- CREATE METHODS

	ZSPACE_INLINE void zTsSDFSlicer::createFieldMesh(zDomain<zPoint>& bb, int resX, int resY)
	{
		zFnMeshScalarField fnField(o_field);

		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);


		zDomainColor dCol(red, green);
		fnField.setFieldColorDomain(dCol);
	}

	//--- SET METHODS 

	ZSPACE_INLINE void zTsSDFSlicer::setFromJSON(string path, int blockStride, int braceStride)
	{
		json j;
		bool fileChk = readJSON(path, j);

		if (!fileChk) return;


		zFnMesh fnMesh(o_GuideMesh);
		fnMesh.from(path, zJSON);


		zPoint* tmpPositions = fnMesh.getRawVertexPositions();

		string id = j["BlockID"];
		//blockId = j["BlockAttributes"][0];

		int sID = j["MedialStartEnd"][0];
		int eID = j["MedialStartEnd"][1];

		zPoint startPoint = tmpPositions[sID];
		zPoint endPoint = tmpPositions[eID];

		computeMedialGraph(o_GuideMesh, sID, eID);
		

		//left plane
		if (j["LeftPlanes"][0].is_null())
		{
			//do nothing  left plane doesnt exist
			leftPlaneExists = false;
		}
		else
		{
			leftPlaneExists = true;
			
			zVector fNormStart(j["LeftPlanes"][0][0], j["LeftPlanes"][0][1], j["LeftPlanes"][0][2]);
			fNormStart *= -1;
			zTransform sPlane = setTransformFromOrigin_Normal(startPoint, fNormStart);				
			
			zVector fNormEnd(j["LeftPlanes"][1][0], j["LeftPlanes"][1][1], j["LeftPlanes"][1][2]);
			zTransform ePlane = setTransformFromOrigin_Normal(endPoint, fNormEnd);			

			cout << "\n left sNorm " << fNormStart;
			cout << "\n left eNorm " << fNormEnd;

			setStartEndPlanes(sPlane, ePlane, true);
			
		}

		//right planes
		if (j["RightPlanes"][0].is_null())
		{
			//do nothing  right plane doesnt exist
			rightPlaneExists = false;

		}
		else
		{
			rightPlaneExists = true;		

			zVector fNormStart(j["RightPlanes"][0][0], j["RightPlanes"][0][1], j["RightPlanes"][0][2]);
			fNormStart *= -1;
			zTransform sPlane = setTransformFromOrigin_Normal(startPoint, fNormStart);

			zVector fNormEnd(j["RightPlanes"][1][0], j["RightPlanes"][1][1], j["RightPlanes"][1][2]);
			zTransform ePlane = setTransformFromOrigin_Normal(endPoint, fNormEnd); ;

			cout << "\n right sNorm " << fNormStart;
			cout << "\n right eNorm " << fNormEnd;

			setStartEndPlanes(sPlane, ePlane, false);

		}

		//left mesh
		if (leftPlaneExists)
		{
			// left mesh
			computeSliceMesh(o_GuideMesh, sID, eID, blockStride, true);
			computeMedial_BraceEdges(o_SliceMesh_Left, 3, 0, blockStride, braceStride);
		}

		//right mesh
		if (rightPlaneExists)
		{
			//right mesh
			computeSliceMesh(o_GuideMesh, sID, eID, blockStride, false);
			computeMedial_BraceEdges(o_SliceMesh_Right, 0, 1, blockStride, braceStride);
		}

	}

	ZSPACE_INLINE void zTsSDFSlicer::setFromJSON(string dir, int _blockID)
	{
		string fileDeck = dir + "deck_" + to_string(_blockID) + ".json";
		bool checkDeck = fileExists(fileDeck);
		if (checkDeck) deckBlock = true;
		if (checkDeck) printf("\n %s exists", fileDeck.c_str());
		else printf("\n %s doesnt exists", fileDeck.c_str());

		string fileBalustrade = dir + "balustrade_" + to_string(_blockID) + ".json";
		bool checkBalustrade = fileExists(fileBalustrade);
		if (checkBalustrade) deckBlock = false;
		if (checkBalustrade) printf("\n %s exists", fileBalustrade.c_str());
		else printf("\n %s doesnt exists", fileBalustrade.c_str());


		json j;
		string path = (checkDeck) ? fileDeck : fileBalustrade;
		bool fileChk = readJSON(path, j);

		if (!fileChk) return;



		zFnMesh fnMesh(o_GuideMesh);
		fnMesh.clear();
		fnMesh.from(path, zJSON);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();


		blockId = _blockID;
		int sID = j["MedialStartEnd"][0];
		int eID = j["MedialStartEnd"][1];
		int blockStride = j["BlockStride"];
		int braceStride = j["BraceStride"];

		numMagentaLoops = blockStride - 1;
		
		zPoint startPoint = tmpPositions[sID];
		zPoint endPoint = tmpPositions[eID];


		//left plane
		if (j["LeftPlanes"][0].is_null()) 
		{
			//do nothing  left plane doesn't exist
			leftPlaneExists = false;
		}
		else
		{
			leftPlaneExists = true;

			zVector fNormStart(j["LeftPlanes"][0][0], j["LeftPlanes"][0][1], j["LeftPlanes"][0][2]);
			fNormStart *= -1;
			zTransform sPlane = setTransformFromOrigin_Normal(startPoint, fNormStart);

			zVector fNormEnd(j["LeftPlanes"][1][0], j["LeftPlanes"][1][1], j["LeftPlanes"][1][2]);
			zTransform ePlane = setTransformFromOrigin_Normal(endPoint, fNormEnd);

			cout << "\n left sNorm " << fNormStart;
			cout << "\n left eNorm " << fNormEnd;

			setStartEndPlanes(sPlane, ePlane, true);

			base_world = sPlane;

		}

		//right planes
		if (j["RightPlanes"][0].is_null()) {
			//do nothing  right plane doesn't exist.
			rightPlaneExists = false;
		}
		else
		{
			rightPlaneExists = true;

			zVector fNormStart(j["RightPlanes"][0][0], j["RightPlanes"][0][1], j["RightPlanes"][0][2]);
			fNormStart *= -1;
			zTransform sPlane = setTransformFromOrigin_Normal(startPoint, fNormStart);

			zVector fNormEnd(j["RightPlanes"][1][0], j["RightPlanes"][1][1], j["RightPlanes"][1][2]);
			zTransform ePlane = setTransformFromOrigin_Normal(endPoint, fNormEnd); ;

			cout << "\n right sNorm " << fNormStart;
			cout << "\n right eNorm " << fNormEnd;

			setStartEndPlanes(sPlane, ePlane, false);

			if (!leftPlaneExists) base_world = sPlane;

		}

		base_local.setIdentity();

		//medial axis
		computeMedialGraph(o_GuideMesh, sID, eID);


		//left mesh
		if (leftPlaneExists)
		{
			// left mesh
			if (deckBlock)
			{
				computeSliceMesh(o_GuideMesh, sID, eID, blockStride, true);
				computeMedial_BraceEdges(o_SliceMesh_Left, 3, 0, blockStride, braceStride);
			}
			else
			{
				zFnMesh fnMesh_Left(o_SliceMesh_Left);
				fnMesh_Left.clear();
				fnMesh_Left.from(path, zJSON);

				computeMedial_BraceEdges(o_SliceMesh_Left, sID, eID, blockStride, braceStride);
			}
			
		}

		//right mesh
		if (rightPlaneExists)
		{			
			// right mesh
			if (deckBlock)
			{
				computeSliceMesh(o_GuideMesh, sID, eID, blockStride, false);
				computeMedial_BraceEdges(o_SliceMesh_Right, 0, 1, blockStride, braceStride);
			}
			else
			{
				zFnMesh fnMesh_Right(o_SliceMesh_Right);
				fnMesh_Right.clear();
				fnMesh_Right.from(path, zJSON);

				computeMedial_BraceEdges(o_SliceMesh_Right, sID, eID, blockStride, braceStride);
			}
		}


	}

	ZSPACE_INLINE void zTsSDFSlicer::setSliceMesh(zObjMesh& _o_SliceMesh, bool left)
	{
		(left) ? o_SliceMesh_Left = _o_SliceMesh : o_SliceMesh_Right = _o_SliceMesh;
	}

	ZSPACE_INLINE void zTsSDFSlicer::setMedialGraph(zObjGraph& _o_MedialGraph)
	{
		o_MedialGraph = _o_MedialGraph;
	}

	ZSPACE_INLINE void zTsSDFSlicer::setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane, bool left)
	{
		(left) ? leftPlanes[0] = _sPlane : rightPlanes[0] = _sPlane;
		(left) ? leftPlanes[1] = _ePlane : rightPlanes[1] = _ePlane;

		(left) ? leftPlaneExists = true : rightPlaneExists = true;
	}

	ZSPACE_INLINE void zTsSDFSlicer::setTransforms(bool toLocal)
	{
		if (toLocal)
		{
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

		}
		else
		{
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

	//---- GET METHODS

	ZSPACE_INLINE vector<zTransform> zTsSDFSlicer::getBlockFrames()
	{	
		return sectionFrames;
	}

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockSectionGraphs( int& numGraphs)
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

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockRaftGraphs(int& numGraphs)
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

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockContourGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_contourGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_contourGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockTrimGraphs(int& numGraphs)
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

	ZSPACE_INLINE zObjPointCloud* zTsSDFSlicer::getRawCriticalPoints(bool minHeight)
	{
		return (minHeight) ? &criticalMinLayer_pts : &criticalMaxLayer_pts;
	}

	ZSPACE_INLINE zObjMeshScalarField* zTsSDFSlicer::getRawFieldMesh()
	{
		return &o_field;
	}

	ZSPACE_INLINE zObjGraph* zTsSDFSlicer::getRawMedialGraph()
	{
		return &o_MedialGraph;
	}

	ZSPACE_INLINE zObjMesh* zTsSDFSlicer::getRawLeftMesh()
	{
		return &o_SliceMesh_Left;
	}

	ZSPACE_INLINE zObjMesh* zTsSDFSlicer::getRawRightMesh()
	{
		return &o_SliceMesh_Right;
	}

	ZSPACE_INLINE zObjMesh* zTsSDFSlicer::getRawGuideMesh()
	{
		return &o_GuideMesh;
	}

	//---- COMPUTE METHODS

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlocks(zDomainFloat& _printHeightDomain, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFlayers, int funcNum , int numSmooth, zDomainFloat neopreneOffset ,  bool compFrames , bool compSDF)
	{
		printHeightDomain = _printHeightDomain;

		bool frameCHECKS = false;

		if (compFrames)
		{			
			for (float printPlaneSpacing = printHeightDomain.min; printPlaneSpacing <= printHeightDomain.max; printPlaneSpacing += 0.0005)
			{
				printf("\n printPlaneSpace %1.4f ", printPlaneSpacing);

				sectionFrames.clear();
				if (leftPlaneExists) computePrintBlockFrames(printPlaneSpacing, neopreneOffset.min, neopreneOffset.max, true);
				if (rightPlaneExists) computePrintBlockFrames(printPlaneSpacing, neopreneOffset.min, neopreneOffset.max, false);

				o_sectionGraphs.clear();
				o_sectionGraphs.assign(sectionFrames.size(), zObjGraph());
				if (leftPlaneExists) computePrintBlockSections(true);
				if (rightPlaneExists) computePrintBlockSections(false);				

				frameCHECKS = checkPrintLayerHeights();
				printf("\n");

				if (frameCHECKS) break;				
			}

			printf("\n frameCHECKS %s ", (frameCHECKS) ? "T" : "F");
		}

		if (compSDF)
		{
			computeSDF(allSDFLayers, numSDFlayers, funcNum, numSmooth, printLayerWidth, neopreneOffset.min, raftLayerWidth);
		}

	}

	ZSPACE_INLINE bool zTsSDFSlicer::onDeckBlock()
	{
		return deckBlock;
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeMedialGraph(zObjMesh& o_Mesh, int startVID, int endVID)
	{
		zFnMesh fnMesh(o_Mesh);

		zPoint* tmpPositions = fnMesh.getRawVertexPositions();
		zPoint startPoint = tmpPositions[startVID];
		zPoint endPoint = tmpPositions[endVID];

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

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
		fnMedial.setEdgeColor(green, false);
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeMedial_BraceEdges(zObjMesh& o_Mesh,int startVID, int endVID,  int blockStride, int braceStride)
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
					(flipHE) ? he_U.getPrev().getEdge().setColor(magenta) : he_U.getNext().getEdge().setColor(magenta);
				}
				if (i == blockStride - 1)
				{
					(flipHE) ? he_U.getNext().getEdge().setColor(orange) : he_U.getPrev().getEdge().setColor(orange);
				}

				he_U = (flipHE) ? he_U.getNext().getSym().getNext() : he_U.getPrev().getSym().getPrev();
			}
			

			//spine
			he.getEdge().setColor(blue);
			he = (flipHE) ? he.getPrev().getSym().getPrev(): he.getNext().getSym().getNext();

			if (he == heStart) exit = true;

		} while (!exit);


	}

	ZSPACE_INLINE void zTsSDFSlicer::computeSliceMesh(zObjMesh& o_Mesh, int startVID, int endVID, int blockStride, bool left)
	{
		unordered_map<string, int> positionVertex;
		zPointArray positions;
		zIntArray pCounts;
		zIntArray pConnects;

		//compute start half edge
		zItMeshHalfEdge heStart = getStartHalfEdge(o_Mesh, startVID, endVID);

		// walk along spine
		bool exit = false;
		zItMeshHalfEdge he = heStart;

		do
		{
			//right
			if (!left)
			{
				zItMeshHalfEdge he_Left = he.getPrev();
								
				for (int i = 0; i < blockStride; i++)
				{
					zItMeshHalfEdge heTmp = he_Left;
					int pCount = 0;

					do
					{
						zPoint p = heTmp.getVertex().getPosition();
						int vID;
						bool chkRepeat = vertexExistsinPositionMap(positionVertex, p, vID);
						
						if (!chkRepeat)
						{
							vID = positions.size();
							addVertexToPositionMap(positionVertex, p, vID);
							positions.push_back(p);							
						}

						pConnects.push_back(vID);
						pCount++;

						heTmp = heTmp.getNext();

					} while (heTmp != he_Left);

					pCounts.push_back(pCount);

					he_Left = he_Left.getPrev().getSym().getPrev();
				}
			}			

			//left
			else
			{
				zItMeshHalfEdge he_Right = he.getSym().getNext();

				for (int i = 0; i < blockStride; i++)
				{
					zItMeshHalfEdge heTmp = he_Right;
					int pCount = 0;

					do
					{
						zPoint p = heTmp.getStartVertex().getPosition();
						int vID;
						bool chkRepeat = vertexExistsinPositionMap(positionVertex, p, vID);

						if (!chkRepeat)
						{
							vID = positions.size();
							addVertexToPositionMap(positionVertex, p, vID);
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
			

			//spine walk
			he = he.getNext().getSym().getNext();

			if (he == heStart) exit = true;



		} while (!exit);

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

			for (int i = 0; i < blockStride; i++)
			{
				heTop_corner = heTop_corner.getPrev().getPrev().getSym();
				heBottom_corner = heBottom_corner.getPrev().getPrev().getSym();
			}


			exit = false;
			do
			{
				if (heTop.getVertex().getId() == endVID) exit = true;

				zPoint p0 = heTop.getVertex().getPosition();
				int v0;
				vertexExistsinPositionMap(positionVertex, p0, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				vertexExistsinPositionMap(positionVertex, p1, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				vertexExistsinPositionMap(positionVertex, p2, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				vertexExistsinPositionMap(positionVertex, p3, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				vertexExistsinPositionMap(positionVertex, p4, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				vertexExistsinPositionMap(positionVertex, p5, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				vertexExistsinPositionMap(positionVertex, p6, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				vertexExistsinPositionMap(positionVertex, p7, v7);
				
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

			for (int i = 0; i < blockStride; i++)
			{
				heTop_corner = heTop_corner.getPrev().getPrev().getSym();
				heBottom_corner = heBottom_corner.getPrev().getPrev().getSym();
			}

			exit = false;
			do
			{
				if (heTop.getStartVertex().getId() == endVID) exit = true;

				zPoint p0 = heTop.getVertex().getPosition();
				int v0;
				vertexExistsinPositionMap(positionVertex, p0, v0);

				zPoint p1 = heTop.getStartVertex().getPosition();
				int v1;
				vertexExistsinPositionMap(positionVertex, p1, v1);

				zPoint p2 = heBottom.getVertex().getPosition();
				int v2;
				vertexExistsinPositionMap(positionVertex, p2, v2);

				zPoint p3 = heBottom.getStartVertex().getPosition();
				int v3;
				vertexExistsinPositionMap(positionVertex, p3, v3);

				pConnects.push_back(v0);
				pConnects.push_back(v1);
				pConnects.push_back(v2);
				pConnects.push_back(v3);

				pCounts.push_back(4);

				// corner
				zPoint p4 = heTop_corner.getStartVertex().getPosition();
				int v4;
				vertexExistsinPositionMap(positionVertex, p4, v4);

				zPoint p5 = heTop_corner.getVertex().getPosition();
				int v5;
				vertexExistsinPositionMap(positionVertex, p5, v5);

				zPoint p6 = heBottom_corner.getStartVertex().getPosition();
				int v6;
				vertexExistsinPositionMap(positionVertex, p6, v6);

				zPoint p7 = heBottom_corner.getVertex().getPosition();
				int v7;
				vertexExistsinPositionMap(positionVertex, p7, v7);

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

		(left) ? fnMesh.setFaceColor(green) : fnMesh.setFaceColor(magenta);

		printf("\n sliceMesh %i %i %i ", fnMesh.numVertices(), fnMesh.numEdges(), fnMesh.numPolygons());

	}

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlockFrames(float printPlaneSpacing, float neopreneOffset_start, float neopreneOffset_end, bool leftBlock)
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
	
		zVector startNorm = (leftBlock)? zVector(leftPlanes[0](2, 0), leftPlanes[0](2, 1), leftPlanes[0](2, 2)) : zVector(rightPlanes[0](2, 0), rightPlanes[0](2, 1), rightPlanes[0](2, 2));
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

		bool left = (leftPlaneExists) ? false : true;
		bool right = (rightPlaneExists) ? false : true;

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
		left = (leftPlaneExists) ? false : true;
		right = (rightPlaneExists) ? false : true;

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


		// Start point
		
		zPoint O = start;
		zVector Z = startNorm;

		zVector tempZ = Z;
		tempZ.normalize();

		zVector X;
		zVector Y(0,1,0);

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
		zTransform pFrame = setTransformFromOrigin_Normal(O, tempZ, zVector(0, 0, 1));
		sectionFrames.push_back(pFrame);

		pOnCurve = O;
		dStart = coreUtils.minDist_Point_Plane(pOnCurve, startOrig, startNorm);


		// in between points
		walkHe = startHe;
		for (int j = 0; j < numLayers; j++)
		{
			zPoint prevPoint = pOnCurve;

			zPoint eEndPoint =  walkHe.getVertex().getPosition();

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

			tempZ = Z;
			tempZ.normalize();

			Y = zVector(0, 1, 0);

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
			pFrame = setTransformFromOrigin_Normal(O, tempZ, zVector(0, 0, 1));
			sectionFrames.push_back(pFrame);

			pOnCurve = O;

			if (j == numLayers - 1)
			{
				float dEnd = coreUtils.minDist_Point_Plane(pOnCurve, endOrig, endNorm);

				if(leftBlock) printf(" \n left sD %1.4f eD %1.4f ", dStart, dEnd);
				else printf(" \n right sD %1.4f eD %1.4f ", dStart, dEnd);
			}
		}

	}

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlockSections(bool left)
	{
		

		zScalarArray scalars;

		int start = 0;
		int end = sectionFrames.size();

		if (leftPlaneExists && rightPlaneExists)
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
			fn_sliceMesh.getIsoContour(scalars, 0.0, positions, edgeConnects, vColors);

			// create graphs
			zFnGraph tempFn(o_sectionGraphs[i]);
			tempFn.create(positions, edgeConnects);;
		
			tempFn.setEdgeColor(grey);
			tempFn.setVertexColors(vColors, false);
		
		}
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeSDF(bool allSDFLayers, int& numSDFlayers,int funcNum, int numSmooth, float printWidth, float neopreneOffset, float raftWidth)
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

		bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;

		int end = (!deckBlock) ? o_sectionGraphs.size() : floor(o_sectionGraphs.size() * 0.5);
		numSDFlayers = (numSDFlayers > end) ? end : numSDFlayers;
		numSDFlayers = (allSDFLayers) ? end : numSDFlayers;

		for (int j = 1; j < numSDFlayers; j++)
		{
			int kStart = (!deckBlock && !rightPlaneExists) ? 1 : 0;
			int kEnd = (!deckBlock && !leftPlaneExists) ? 1 : 2;


			for (int k = kStart; k < kEnd; k++)
			{
				int i = (k == 0) ? j : j + end;

				if (!deckBlock) i = j;

				if (!deckBlock)
				{
					if (i == r0)
					{
						//raft 
					}
					else
					{
						computeBlockSDF_Balustrade(funcNum, numSmooth, i, (j % 2 == 0), printWidth, neopreneOffset, false, 0, raftWidth);

					}
				}
				else
				{
					if (i == r0) 
					{ 
						//raft 
					}
					else if (i == r2)
					{
						//raft 
					}
					else
					{
						computeBlockSDF_Deck(funcNum, numSmooth, i, (j % 2 == 0), printWidth, neopreneOffset, false, 0, raftWidth);
					}
				}

			}
		}

		

	}

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlockTrimGraphs(zObjGraph& inPolyObj, zObjGraph& o_outGraph, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE)
	{
		zFnGraph inFnGraph(inPolyObj);
		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zPointArray gPositions;
		zIntArray gEdgeCOnnects;

		for (auto& he : bottomHE)
		{
			if (inColors[he.getStartVertex().getId()] == orange) gPositions.push_back(inPositions[he.getStartVertex().getId()]);
			if (inColors[he.getVertex().getId()] == magenta) gPositions.push_back(inPositions[he.getVertex().getId()]);
			if (inColors[he.getVertex().getId()] == blue) gPositions.push_back(inPositions[he.getVertex().getId()]);

		
		}

		for (auto& he : topHE)
		{
			if (inColors[he.getStartVertex().getId()] == orange) gPositions.push_back(inPositions[he.getStartVertex().getId()]);
			if (inColors[he.getVertex().getId()] == magenta) gPositions.push_back(inPositions[he.getVertex().getId()]);
			if (inColors[he.getVertex().getId()] == blue) gPositions.push_back(inPositions[he.getVertex().getId()]);
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

	ZSPACE_INLINE bool zTsSDFSlicer::checkPrintLayerHeights()
	{
		float minLayerHeight = 10;
		float maxLayerHeight = 0;
		int minHeightGraphID = -1;
		bool checkOranges = true;
		bool checkMagentas = true;

		zFnPointCloud fnCritical_min(criticalMinLayer_pts);
		zFnPointCloud fnCritical_max(criticalMaxLayer_pts);

		fnCritical_min.clear();
		fnCritical_max.clear();

		int r0 = 0;
		int r1 = floor(o_sectionGraphs.size() * 0.5) - 1;
		int r2 = floor(o_sectionGraphs.size() * 0.5);
		int r3 = (o_sectionGraphs.size()) - 1;
		printf("\n r: %i  %i %i %i ", r0, r1, r2, r3);

		bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;

		int end = (leftPlaneExists && rightPlaneExists) ? floor(o_sectionGraphs.size() * 0.5) : o_sectionGraphs.size();

		float printLength = 0;

		// PRINT LAYERS
		for (int j = 0; j < end; j++)
		{
			int kStart = (!deckBlock && !rightPlaneExists) ? 1 : 0;
			int kEnd = (!deckBlock && !leftPlaneExists) ? 1 : 2;

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

				for (zItGraphVertex v(o_sectionGraphs[i]); !v.end(); v++)
				{
					zPoint p = v.getPosition();

					if (v.getColor() == orange) orangeCounter++;
					if (v.getColor() == magenta) magentaCounter++;

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
					checkOranges = false;	

					
				}

				if (magentaCounter != numMagentaLoops * 2) checkMagentas = false;	

				for (zItGraphEdge e(o_sectionGraphs[i]); !e.end(); e++)
				{
					printLength += e.getLength();
				}


			}
		}

		fnCritical_min.setVertexColor(magenta);		
		fnCritical_max.setVertexColor(magenta);

		printf("\n block| %1.4f %1.4f| %1.1f | chkOrange %s | chkMagenta %s | < min ht pts %i  | > max ht pts %i ", minLayerHeight, maxLayerHeight, printLength, (checkOranges) ? "T" : "F", (checkMagentas) ? "T" : "F", fnCritical_min.numVertices(), fnCritical_max.numVertices());

		bool out = true;
		if (!checkOranges || !checkMagentas) out = false;
		if (out)
		{
			if (minLayerHeight < printHeightDomain.min) out = false;
			if (minLayerHeight > printHeightDomain.max) out = false;

			if (maxLayerHeight < printHeightDomain.min) out = false;
			if (maxLayerHeight > printHeightDomain.max) out = false;
		}

		return out;
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeBlockSDF_Deck(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset,  bool addRaft, int raftId, float raftWidth)
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
		polyTopBottomEdges(o_sectionGraphs[graphId], topHE, bottomHE, topLength, bottomLength);


		// Transform
		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_inner = 1.5 * pWidth;

		//trim graph
		computePrintBlockTrimGraphs(o_sectionGraphs[graphId], o_trimGraphs[graphId], topHE, bottomHE);
		zFnGraph fnTrimGraph(o_trimGraphs[graphId]);
		

		// Profile polygon field
		zScalarArray polyField;
		if (funcNum >= 0) fnField.getScalars_Polygon(polyField, o_sectionGraphs[graphId], false);

		zScalarArray polyField_offset_outer;
		if (funcNum >= 1)
		{
			polyField_offset_outer = polyField;;
			for (auto& s : polyField_offset_outer) s += offset_outer;
		}

		//zScalarArray braceField;
		//if (funcNum >= 2) getScalars_3dp_brace(braceField, o_trimGraphs[graphId], pWidth, 0.25 * pWidth, alternate);

		//zScalarArray slotField;
		//float slotRadius = 0.02;
		//if (funcNum >= 3) getScalars_3dp_slot(slotField, o_trimGraphs[graphId], slotRadius + 0.25 * pWidth); 

		//zScalarArray booleanField_0;
		//if (funcNum >= 4) fnField.boolean_subtract(polyField_offset_outer, braceField, booleanField_0, false);

		//zScalarArray booleanField_1;
		//if (funcNum >= 5) fnField.boolean_subtract(booleanField_0, slotField, booleanField_1, false);

		// No SLOT, same as Balustrade
		zScalarArray braceField;
		if (funcNum >= 2) getScalars_3dp_brace(braceField, o_trimGraphs[graphId], 0, 0.25 * pWidth, alternate);

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
			if(numSmooth > 0) fnField.smoothField(booleanField_1, numSmooth); // smooth field
			fnField.setFieldValues(booleanField_1);
			break;
		}
		

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0);

		fnIsoGraph.setEdgeWeight(2);

		// transform back 
		
		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		fnTrimGraph.setTransform(t, true, true);
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeBlockSDF_Balustrade(int funcNum, int numSmooth, int graphId, bool alternate, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth)
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
		polyTopBottomEdges(o_sectionGraphs[graphId], topHE, bottomHE, topLength, bottomLength);


		// Transform
		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_inner = 1.5 * pWidth;

		//trim graph
		computePrintBlockTrimGraphs(o_sectionGraphs[graphId], o_trimGraphs[graphId], topHE, bottomHE);
		zFnGraph fnTrimGraph(o_trimGraphs[graphId]);
		

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
		if (funcNum >= 2) getScalars_3dp_brace(braceField, o_trimGraphs[graphId], 0, 0.25 * pWidth, alternate);

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
			if(numSmooth > 0) fnField.smoothField(booleanField_1, numSmooth);
			fnField.setFieldValues(booleanField_1);
			break;
		}


		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0);

		fnIsoGraph.setEdgeWeight(2);

		// transform back 
		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
		fnTrimGraph.setTransform(t, true, true);
	}


	ZSPACE_INLINE bool zTsSDFSlicer::exportJSON(string pathCurrent, string dir, string filename, float printLyerWidth, float raftLayerWidth)
	{

		json j;
		bool fileChk = readJSON(pathCurrent, j);

		if (!fileChk) return false;

		string folderName = dir + "/" + to_string(blockId);
		_mkdir(folderName.c_str());
		for (const auto& entry : std::filesystem::directory_iterator(folderName)) std::filesystem::remove_all(entry.path());

		// EXPORT	BlockMesh
	
		string blockID_padded = coreUtils.getPaddedIndexString(blockId,3);

		string blockMeshName = folderName + "/block_" + blockID_padded + ".json";

		ofstream myfile;
		myfile.open(blockMeshName.c_str());

		if (myfile.fail())
		{
			cout << " error in opening file  " << blockMeshName.c_str() << endl;
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

		bool deckBlock = (leftPlaneExists && rightPlaneExists) ? true : false;

		int end = (leftPlaneExists && rightPlaneExists) ? floor(o_contourGraphs.size() * 0.5) : o_contourGraphs.size();

		float printLength = 0;

		// PRINT LAYERS
		for (int j = 0; j < end; j++)
		{
			int kStart = (!deckBlock && !rightPlaneExists) ? 1 : 0;
			int kEnd = (!deckBlock && !leftPlaneExists) ? 1 : 2;

			for (int k = kStart; k < kEnd; k++)
			{
				int i = (k == 0) ? j : j + end;

				if (!deckBlock) i = j;

				if (i == r0) continue;

				if (deckBlock && i == r2) continue;

				string graphID_padded = coreUtils.getPaddedIndexString(j, 3);

				// trim graph export
				string outName1 = folderName;
				outName1 += "/";
				outName1 += "trim";
				outName1 += "_";
				outName1 += blockID_padded + "_" + graphID_padded + "_" + to_string(k) + ".json";

				zFnGraph fnTrimGraph(o_trimGraphs[i]);
				fnTrimGraph.to(outName1, zJSON);

				// graph export

				string outName = folderName;
				outName += "/";
				outName += filename;
				outName += "_";
				outName += blockID_padded + "_" + graphID_padded + "_" + to_string(k) + ".json";

				zFnGraph fnIsoGraph(o_contourGraphs[i]);
				fnIsoGraph.to(outName, zJSON);

				// read existing data in the json 
				json jSON;

				ifstream in_myfile;
				in_myfile.open(outName.c_str());

				int lineCnt = 0;

				if (in_myfile.fail())
				{
					cout << " error in opening file  " << outName.c_str() << endl;
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
					cout << " error in opening file  " << outName.c_str() << endl;
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

	
	//---- UTILITY METHODS


	ZSPACE_INLINE zTransform zTsSDFSlicer::setTransformFromVectors(zPoint& O, zVector& X, zVector& Y, zVector& Z)
	{
		zTransform out;

		out(0, 0) = X.x; out(0, 1) = X.y; out(0, 2) = X.z; out(0, 3) = 1;
		out(1, 0) = Y.x; out(1, 1) = Y.y; out(1, 2) = Y.z; out(1, 3) = 1;
		out(2, 0) = Z.x; out(2, 1) = Z.y; out(2, 2) = Z.z; out(2, 3) = 1;
		out(3, 0) = O.x; out(3, 1) = O.y; out(3, 2) = O.z; out(3, 3) = 1;


		return out;
	}

	zTransform zTsSDFSlicer::setTransformFromOrigin_Normal(zPoint& O, zVector& Z, zVector Basis)
	{
		zVector X = Basis ^ Z;		

		zVector Y = Z ^ X;
		Y.normalize();

		X = Y ^ Z;
		X.normalize();

		return setTransformFromVectors(O,X,Y,Z);
	}

	ZSPACE_INLINE zItMeshHalfEdge zTsSDFSlicer::getStartHalfEdge(zObjMesh& o_Mesh, int startVID, int endVID)
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
			if (1 - (heVec * dir) < val)
			{
				val = 1 - (heVec * dir);
				heStart = he;
			}
		}

		return heStart;
	}

	ZSPACE_INLINE void zTsSDFSlicer::polyTopBottomEdges(zObjGraph& inPoly, zItGraphHalfEdgeArray& topHE, zItGraphHalfEdgeArray& bottomHE, float& topLength, float& bottomLength)
	{
		zFnGraph inFnGraph(inPoly);
		inFnGraph.setEdgeColor(grey);

		bottomHE.clear();
		topHE.clear();

		zVector Y(0, 1, 0);

		zVector* inPositions = inFnGraph.getRawVertexPositions();
		zColor* inColors = inFnGraph.getRawVertexColors();

		zIntArray startVerts;
		for (int i = 0; i < inFnGraph.numVertices(); i++)
		{
			if (inColors[i] == orange) startVerts.push_back(i);;
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

			if (v.getColor() == blue) exit = true;		

			bHe = bHe.getNext();
		}

		exit = false;
		while (!exit)
		{
			topHE.push_back(tHe);
			topLength += tHe.getLength();

			tHe.getEdge().setColor(zColor(1, 0, 0, 1));
			zItGraphVertex v = tHe.getVertex();
			if (v.getColor() == blue) exit = true;		

			tHe = tHe.getNext();
		}


	}

	ZSPACE_INLINE void zTsSDFSlicer::getScalars_3dp_slot(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset)
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

	ZSPACE_INLINE void zTsSDFSlicer::getScalars_3dp_brace(zScalarArray& scalars, zObjGraph& o_trimGraph, float outer_printWidth, float offset, bool alternate)
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

	ZSPACE_INLINE void zTsSDFSlicer::getScalars_3dp_trim(zScalarArray& scalars, zObjGraph& o_trimGraph, float offset, bool alternate)
	{
		zFnMeshScalarField fnField(o_field);
		zIntArray tEdgeConnects;
		zPointArray tPositions;

		zFnGraph fnTrimGraph(o_trimGraph);
		zPoint* trimPositions = fnTrimGraph.getRawVertexPositions();

		int end = floor(fnTrimGraph.numVertices() * 0.5);
		zVector upVec(0, 0, 1);

		//printf("\n end %i ", end);

		int offsetMult = 2;

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

	ZSPACE_INLINE void zTsSDFSlicer::addVertexToPositionMap(unordered_map<string, int>& positionVertex, zPoint& pos, int index, int precisionfactor)
	{
		double factor = pow(10, precisionfactor);
		double x = std::round(pos.x * factor) / factor;
		double y = std::round(pos.y * factor) / factor;
		double z = std::round(pos.z * factor) / factor;

		string hashKey = (to_string(x) + "," + to_string(y) + "," + to_string(z));
		positionVertex[hashKey] = index;
	}

	ZSPACE_INLINE bool zTsSDFSlicer::vertexExistsinPositionMap(unordered_map<string, int>& positionVertex, zPoint& pos, int& outVertexId, int precisionfactor)
	{
		bool out = false;;
		outVertexId = -1;

		double factor = pow(10, precisionfactor);
		double x = std::round(pos.x * factor) / factor;
		double y = std::round(pos.y * factor) / factor;
		double z = std::round(pos.z * factor) / factor;

		string hashKey = (to_string(x) + "," + to_string(y) + "," + to_string(z));
		std::unordered_map<std::string, int>::const_iterator got = positionVertex.find(hashKey);


		if (got != positionVertex.end())
		{
			out = true;
			outVertexId = got->second;
		}


		return out;
	}

	ZSPACE_INLINE bool zTsSDFSlicer::readJSON(string path, json& j)
	{		
		j.clear();
		ifstream in_myfile;
		in_myfile.open(path.c_str());

		int lineCnt = 0;

		if (in_myfile.fail())
		{
			cout << " error in opening file  " << path.c_str() << endl;
			return false;
		}

		in_myfile >> j;
		in_myfile.close();

		return true;
	}

	ZSPACE_INLINE bool zTsSDFSlicer::fileExists(string& path)
	{
		ifstream f(path.c_str());
		return f.good();

	}

	

}