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


#include<headers/zToolsets/0_externalMethods/zExtSDFSlicer.h>


namespace zSpace
{
	ZSPACE_TOOLSETS_INLINE void ext_sdf_SetFromJSON(zExtTsSDFSlicer& extSlicer, char* path, int blockID)
	{
		std::string pathSt(path);
		//extSlicer = zExtTsSDFSlicer();
		extSlicer.slicer = new zTsSDFSlicer();
		cout << pathSt << endl;
		extSlicer.slicer->setFromJSON(pathSt, blockID);
		extSlicer.updateFields();

		printf("\n ext_sdf_SetFromJSON: BlockID %i is %s - ", extSlicer.blockID, (extSlicer.deckBlock) ? "DeckBlock" : "Balustrade");
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_createFieldMesh(zExtTsSDFSlicer& extSlicer, float* pointDomainMin, float* pointDomainMax, int resX, int resY, zExtMesh& extFieldMesh)
	{
		if (!extSlicer.slicer)
		{
			printf("\n slicer is null - FieldMesh");
			return 0;
		}
		//zTsSDFSlicer* slicer = extSlicer.slicer;
		zDomain<zPoint> bb = zDomain<zPoint>(zPoint(pointDomainMin[0], pointDomainMin[1], pointDomainMin[2]), zPoint(pointDomainMax[0], pointDomainMax[1], pointDomainMax[2]));
		extSlicer.slicer->createFieldMesh(bb, resX, resY);
		extSlicer.updateFields();

		extFieldMesh.mesh = extSlicer.slicer->getRawFieldMesh();
		extFieldMesh.updateFields();


		return 1;
	}

	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& rightMesh, zExtMesh& leftMesh)
	{
		if (!extSlicer.slicer)
		{
			return 0;
		}
		rightMesh.mesh = extSlicer.slicer->getRawRightMesh();
		rightMesh.updateFields();
		leftMesh.mesh = extSlicer.slicer->getRawLeftMesh();
		leftMesh.updateFields();

		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockSectionGraphs(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet)
	{
		if (!extSlicer.slicer)
		{
			return 0;
		}
		graphSet.graphSet = new zObjGraphPointerArray();
		int numGraphs = 0;
		*graphSet.graphSet = extSlicer.slicer->getBlockSectionGraphs(numGraphs);
		graphSet.graphsCount = graphSet.graphSet->size();

		return 1;
	}


	ZSPACE_TOOLSETS_INLINE void ext_sdf_SetFromJSON3(zTsSDFSlicer*& slicer, char* path, int blockID, bool& isDeckBlock)
	{
		std::string pathSt(path);
		if (!slicer)
		{
			slicer = new zTsSDFSlicer();
		}
		cout << pathSt << endl;
		slicer->setFromJSON(pathSt, blockID);
		isDeckBlock = slicer->onDeckBlock();

		printf("\n ext_sdf_SetFromJSON: BlockID %i is %s - ", slicer->blockId, (isDeckBlock)? "DeckBlock" : "Balustrade");
	}
	ZSPACE_TOOLSETS_INLINE void ext_sdf_SetFromJSON2(zTsSDFSlicer*& slicer, char* path, int blockStride, int braceStride)
	{
		//printf("\n ext_sdf_JSON: %i", 0);

		std::string pathSt(path);
		if (!slicer)
		{
			slicer = new zTsSDFSlicer();
		}
	
		cout << pathSt << endl;
		slicer->setFromJSON(pathSt, blockStride, braceStride);
		

	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_checkLayerHeight2(zTsSDFSlicer*& slicer, bool& check)
	{
		if (!slicer)
		{
			printf("\n slicer is null - CheckLayerHeights");
			return 0;
		}
		bool checkSDF; bool checkGeo;
		check = slicer->checkPrintLayerHeights(checkSDF, checkGeo);
	}

	ZSPACE_TOOLSETS_INLINE void ext_sdf_createFieldMesh2(zTsSDFSlicer*& slicer, float* pointDomainMin, float* pointDomainMax, int resX, int resY)
	{
		if (!slicer)
		{
			printf("\n slicer is null - FieldMesh");
			slicer = new zTsSDFSlicer();
		}
		zDomain<zPoint> bb = zDomain<zPoint>(zPoint(pointDomainMin[0], pointDomainMin[1], pointDomainMin[2]), zPoint(pointDomainMax[0], pointDomainMax[1], pointDomainMax[2]));
		slicer->createFieldMesh(bb, resX, resY);
		//printf("\n ext_sdf_createFieldMesh: BlockID %i", slicer->blockId);
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_computePrintBlocks2(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF)
	{
		//printf("\n ext_sdf_computePrintBlocks %i", slicer);
		//printf("\n o_sectionGraphs %i", slicer->o_sectionGraphs.size());
		if (!slicer)
		{
			return 0;
		}
		zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);

		zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
		float printPlaneSpace;
		/*for ( printPlaneSpace = printHeightDomain.min; printPlaneSpace <= printHeightDomain.max; printPlaneSpace += printPlaneStep)
		{
			printf("\n printPlaneSpace %1.4f ", printPlaneSpace);

			slicer->computePrintBlocks(printHeightDomain, printLayerWidth, raftLayerWidth, allSDFLayers, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, true, false);

			bool layerHeight;

			bool checkSDF; bool checkGeo;
			bool frameCHECKS = slicer->checkPrintLayerHeights(checkSDF, checkGeo);

			if (frameCHECKS) break;

			printf("\n");
		}*/

		slicer->computePrintBlocks(printHeightDomain, printLayerWidth, raftLayerWidth, allSDFLayers, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, compFrames, compSDF);



		//printf("\n compute 1");
		//slicer->computePrintBlocks(printPlaneSpacing, printLayerWidth, raftLayerWidth, allSDFLayer, computeGraphCount, funcNum,  neopreneOffset, compFrames, compSDF);
		//slicer->checkPrintLayerHeights();
		//slicer->computePrintBlocks(printPlaneSpacing, printLayerWidth, raftLayerWidth, allSDFLayer, computeGraphCount, funcNum, neopreneOffset, true, compSDF);

		//printf("\n compute 3");
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_computePrintBlocksPar2(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF)
	{

		if (!slicer)
		{
			return 0;
		}
		zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);

		zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
		float printPlaneSpace;
		

		slicer->computePrintBlocksPar(printHeightDomain, printLayerWidth, raftLayerWidth, allSDFLayers, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, compFrames, compSDF);



		return 1;
	}

	ZSPACE_TOOLSETS_INLINE void ext_sdf_computeSDFLayer2(zTsSDFSlicer*& slicer, int SDFLayerNumber, float printLayerWidth, float raftLayerWidth,  int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames)
	{
		
		if (compFrames)
		{
			zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);
			zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
			int numSDFLayers = 5;
			slicer->computePrintBlocks(printHeightDomain, printLayerWidth, raftLayerWidth, false, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, true, false);
		}

		slicer->computeSDFSingleLayer(SDFLayerNumber, SDFFunc_Num, SDFFunc_NumSmooth, printLayerWidth, neopreneOffsetMin, raftLayerWidth);
	}
	ZSPACE_TOOLSETS_INLINE void ext_sdf_duplicateSlicer2(zTsSDFSlicer* inSlicer, zTsSDFSlicer*& outSlicer)
	{
		//if (!inSlicer)  inSlicer = new zTsSDFSlicer();
		outSlicer = inSlicer->copy();

	}

	//Get Methods
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawMesh2(zTsSDFSlicer* slicer, zObjMesh*& objMesh, bool rightMesh)
	{
		if (!slicer)
		{
			return 0;
		}
		if (rightMesh) objMesh = slicer->getRawRightMesh();
		else objMesh = slicer->getRawLeftMesh();

		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawFieldMesh2(zTsSDFSlicer* slicer, zObjMesh*& objMesh)
	{
		if (!slicer)
		{
			return 0;
		}
		objMesh = slicer->getRawFieldMesh();
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockSectionGraphs2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount)
	{
		if (!slicer)
		{
			return 0;
		}
		graph = new zObjGraphPointerArray();
		*graph = slicer->getBlockSectionGraphs(numGraphs);
		outGraphsCount = graph->size();

		

		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockContourGraphs2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount)
	{
		if (!slicer)
		{
			return 0;
		}
		graph = new zObjGraphPointerArray();
		*graph = slicer->getBlockContourGraphs(numGraphs);
		outGraphsCount = graph->size();
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawMedialGraph2(zTsSDFSlicer* slicer, zObjGraph*& graph)
	{
		if (!slicer)
		{
			return 0;
		}
		graph = slicer->getRawMedialGraph();
		
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockFrame2s(zTsSDFSlicer* slicer, vector<zTransform>*& planes, int& count)
	{
		if (!slicer)
		{
			return 0;
		}
		printf("\n getBlockFrames - 0");
		planes = new vector<zTransform>;
		*planes = slicer->getBlockFrames();
		count = planes->size();
		printf("\n planes count  %i", count);
		return 1;
	}

	ZSPACE_TOOLSETS_INLINE int ext_sdf_getTrimGraph2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graphs, int& Count)
	{
		if (!slicer)
		{
			return 0;
		}
		graphs = new zObjGraphPointerArray();
		*graphs = slicer->getBlockTrimGraphs(Count);
		//Count = slicer->o_trimGraphs.size();
		return 1;

	}

	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getPlanesData2(vector<zTransform>* graph, float* matrix)
	//{
	//	for (int i = 0; i < graph->size(); i++)
	//	{
	//		//outPlanes[i * 4 + 0] = 
	//		zTransform frame = graph->at(i);

	//		zTransform frameTranspose = frame.transpose();

	//		float* m = frameTranspose.data();
	//		
	//		//float* m = graph->at(i).transpose().data(); //doesn't work
	//		for (int j = 0; j < 16; j++)
	//		{
	//			matrix[i * 16 + j] = m[j];
	//		}

	//		
	//	}
	//}





	////Graph Data
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getGraphsSetFromPointersVector(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray)
	//{
	//	//printf("\n C++ EXT (getGraphSet) num of graphs %i", graphs->size());
	//	for (int i = 0; i < graphs->size(); i++)
	//	{
	//		outGraphArray[i] = graphs->at(i);
	//	}
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getGraphsSetFromVector(zObjGraphArray* graphs, zObjGraph** outGraphArray)
	//{
	//	for (int i = 0; i < graphs->size(); i++)
	//	{
	//		

	//		outGraphArray[i] = &graphs->at(i);

	//	}
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getGraphCounts(zObjGraph* graph, int& outvCount, int& outeCount)
	//{
	//	zFnGraph g(*graph);
	//	
	//	outeCount = g.numEdges();
	//	outvCount = g.numVertices();
	//	
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getGraphData(zObjGraph* graph, float* vPositions, float* vColors, int* ePairs, float* eColors)
	//{
	//	zFnGraph g(*graph);
	//	int vCount = g.numVertices();
	//	int eCount = g.numEdges();


	//	zPointArray inVerticies;
	//	g.getVertexPositions(inVerticies);

	//	zColorArray inVColors;
	//	g.getVertexColors(inVColors);

	//	zColorArray ineColors;
	//	g.getEdgeColors(ineColors);

	//	zIntArray inEdges;
	//	g.getEdgeData(inEdges);
	//	//printf("\n graph edges size %i ", vCount);
	//	for (int i = 0; i < vCount; i++)
	//	{
	//		vPositions[i * 3 + 0] = inVerticies[i].x;
	//		vPositions[i * 3 + 1] = inVerticies[i].y;
	//		vPositions[i * 3 + 2] = inVerticies[i].z;

	//		vColors[i * 4 + 0] = inVColors[i].r;
	//		vColors[i * 4 + 1] = inVColors[i].g;
	//		vColors[i * 4 + 2] = inVColors[i].b;
	//		vColors[i * 4 + 3] = inVColors[i].a;
	//	}
	//	//printf("\n graph edges size %i ", g.numEdges());
	//	//printf("\n graph edgeVerticies size %i ", inEdges.size());
	//	//printf("\n graph vCount %i ", vCount);
	//	//printf("\n graph edgeColors size %i ", inVerticies.size());
	//	//printf("\n graph edgeColors size %i ", inVColors.size());

	//	//printf("\n graph eCount %i ", eCount);
	//	//printf("\n graph inEdges size %i ", inEdges.size());
	//	//printf("\n graph ineColors size %i ", ineColors.size());
	//	for (int64_t i = 0; i < eCount; i++)
	//	{

	//		ePairs[i * 2 + 0] = inEdges[i * 2 + 0];
	//		ePairs[i * 2 + 1] = inEdges[i * 2 + 1];

	//		eColors[i * 4 + 0] = ineColors[i].r;
	//		eColors[i * 4 + 1] = ineColors[i].g;
	//		eColors[i * 4 + 2] = ineColors[i].b;
	//		eColors[i * 4 + 3] = ineColors[i].a;

	//	}


	//	//printf("\n C++ EDGES: %i ", eCount);
	//	//for (int i = 0; i < eCount; i++)
	//	//{
	//	//	//printf("\n edge pair %i - %i", ePairs[i * 2], ePairs[i * 2 + 1]);
	//	//}



	//}

	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getGraphSequence(zObjGraph* graph, int* outSequence)//the size of the area = the number of vertices + 1
	//{
	//	zIntArray seq;
	//	seq = zTsSDFSlicer::getGraphSequence(*graph);
	//	//printf("\n sequence: \n");
	//	for (int i = 0; i < seq.size(); i++)
	//	{
	//		outSequence[i] = seq[i];
	//		//printf("%i , ", outSequence[i]);
	//	}
	//}
	//
	////Mesh Data
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getMeshCounts(zObjMesh* objMesh, int& out_vCount, int& out_fCount)
	//{
	//	zFnMesh fn(*objMesh);
	//	out_vCount = fn.numVertices();
	//	out_fCount = fn.numPolygons();
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getMeshPosition(zObjMesh* objMesh, float* outVPostions, float* outVColors)
	//{
	//	zFnMesh fn(*objMesh);
	//	zPoint* pts = fn.getRawVertexPositions();
	//	zColor* colors = fn.getRawVertexColors();
	//	for (int i = 0; i < fn.numVertices(); i++)
	//	{
	//		outVPostions[i * 3 + 0] = pts[i].x;
	//		outVPostions[i * 3 + 1] = pts[i].y;
	//		outVPostions[i * 3 + 2] = pts[i].z;

	//		outVColors[i * 4 + 0] = colors[i].r;
	//		outVColors[i * 4 + 1] = colors[i].g;
	//		outVColors[i * 4 + 2] = colors[i].b;
	//		outVColors[i * 4 + 3] = colors[i].a;
	//	}
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getMeshFaceCount(zObjMesh* objMesh, int* outfCounts)
	//{
	//	zFnMesh fn(*objMesh);
	//	zIntArray pCounts;
	//	zIntArray pConnects;
	//	fn.getPolygonData(pConnects, pCounts);
	//	for (int i = 0; i < pCounts.size(); i++)
	//	{
	//		outfCounts[i] = pCounts[i];
	//	}
	//}
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_getMeshFaceConnect(zObjMesh* objMesh, int* outfConnects)
	//{
	//	zFnMesh fn(*objMesh);
	//	zIntArray pCounts;
	//	zIntArray pConnects;
	//	fn.getPolygonData(pConnects, pCounts);
	//	for (int i = 0; i < pConnects.size(); i++)
	//	{
	//		outfConnects[i] = pConnects[i];
	//	}
	//}

	//
	//ZSPACE_TOOLSETS_INLINE void ext_sdf_ExportJSON(zTsSDFSlicer* slicer, char* fileCurrent, int fileCount, char* fileNew, int newCount,  char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth)
	//{

	//	string file = "";
	//	string exportDir="";
	//	string exportName="";
	//	for (int i = 0; i < fileCount; i++) { file += fileCurrent[i]; }
	//	for (int i = 0; i < newCount; i++) { exportDir += fileNew[i]; }
	//	for (int i = 0; i < nameCount; i++) { exportName += fileName[i]; }


	//	cout << "\n Export_2 entryDir" << endl;
	//	cout << file << endl;
	//	cout << "\n Export_2 exportDir" << endl;
	//	cout << exportDir << endl;
	//	cout << "\n Export_2 exportName" << endl;
	//	cout << exportName << endl;

	//	if (!slicer)
	//	{
	//		printf("\n slicer is null - export");
	//	}

	//	slicer->exportJSON(file, exportDir, exportName, printLayerWidth, raftLayerWidth);
	//}

	//ZSPACE_TOOLSETS_INLINE void ext_sdf_CheckFolder(char* folderDirectoryChar, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax)
	//{
	//	std::string folder(folderDirectoryChar);
	//	zTsSDFSlicer slicer = zTsSDFSlicer();
	//	
	//	zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);
	//	zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
	//	slicer.checkPrintLayerHeights_Folder(folder, printHeightDomain, neopreneOffset); 
	//}

	//ZSPACE_TOOLSETS_INLINE void ext_sdf_DrawMesh(zObjMesh* objMesh)
	//{
	//	printf("\n c++ before draw method");

	//	glColor3f(0, 0, 0);
	//	glPointSize(5);

	//	glBegin(GL_POINTS);
	//	glVertex3f(1, 1, 1);
	//	glEnd();

	//	//objMesh->draw();
	//	printf("\n c++ after draw method");

	//}

	//Export Data

	

	

}