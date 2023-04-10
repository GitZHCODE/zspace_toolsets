//// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
//// data analysis & visualization framework.
////
//// Copyright (C) 2019 ZSPACE 
//// 
//// This Source Code Form is subject to the terms of the MIT License 
//// If a copy of the MIT License was not distributed with this file, You can 
//// obtain one at https://opensource.org/licenses/MIT.
////
//// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
////
//
//
//#include<headers/zToolsets/externalMethods/zExtSDFSlicer.h>
//
//
//namespace zSpace
//{
//	ZSPACE_TOOLSETS_INLINE  zExtTsSDFSlicer::zExtTsSDFSlicer(zTsSDFSlicer* s)
//	{
//		slicer = s;
//		updateAttributes();
//	}
//	ZSPACE_TOOLSETS_INLINE void zExtTsSDFSlicer::updateAttributes()
//	{
//		blockID = slicer->blockId;
//		deckBlock = slicer->onDeckBlock()? 1:0;
//		leftPlaneExist = slicer->leftPlaneExists ? 1 : 0;
//		rightPlaneExist = slicer->rightPlaneExists ? 1 : 0;
//		sectionFrameCount = slicer->getBlockFrames().size();
//		slicer->getBlockSectionGraphs(sectionGraphCount);
//		slicer->getBlockContourGraphs(contourGraphCount);
//		slicer->getBlockTrimGraphs(trimGraphCount);
//
//		printf("\n \n SLICER ATTRIBUTES C++ - {update}");
//		printf("\n blockID				%i", blockID);
//		printf("\n deck Block			%i", deckBlock );
//		printf("\n left PlaneExist		%i", leftPlaneExist );
//		printf("\n right PlaneExist		%i", rightPlaneExist );
//		printf("\n section Frame Count	%i", sectionFrameCount	);
//		printf("\n section Graph Count	%i", sectionGraphCount	);
//		printf("\n contour Graph Count	%i", contourGraphCount	);
//		printf("\n trim Graph Count		%i", trimGraphCount		);
//
//	}
//
//
//	//SET
//	ZSPACE_TOOLSETS_INLINE void ext_sdf_SetFromJSON(zExtTsSDFSlicer& extSlicer, char* path, int blockID)
//	{
//		std::string pathSt(path);
//		//extSlicer = zExtTsSDFSlicer();
//		cout << pathSt << endl;
//		//zTsSDFSlicer s;// = new zTsSDFSlicer();
//		//s.setFromJSON(pathSt, blockID);
//
//		extSlicer.slicer = new zTsSDFSlicer();
//		extSlicer.slicer->setFromJSON(pathSt, blockID);
//		extSlicer.updateAttributes();
//
//		/*printf("\n \n SLICER ATTRIBUTES C++ - {SET}");
//		printf("\n blockID					%i", extSlicer.blockID);
//		printf("\n deck Block				%s", extSlicer.deckBlock ? "T" : "F");
//		printf("\n left PlaneExist			%s", extSlicer.leftPlaneExist ? "T" : "F");
//		printf("\n right PlaneExist			%s", extSlicer.rightPlaneExist ? "T" : "F");
//		printf("\n section Frame Count		%i", extSlicer.sectionFrameCount);
//		printf("\n section Graph Count		%i", extSlicer.sectionGraphCount);
//		printf("\n contour Graph Count		%i", extSlicer.contourGraphCount);
//		printf("\n trim Graph Count			%i", extSlicer.trimGraphCount);*/
//
//		printf("\n ext_sdf_SetFromJSON: BlockID %i is %s - ", extSlicer.blockID, (extSlicer.deckBlock) ? "DeckBlock" : "Balustrade");
//	}
//	
//	//COMPUTE
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_createFieldMesh(zExtTsSDFSlicer& extSlicer, float* pointDomainMin, float* pointDomainMax, int resX, int resY, zExtMesh& extFieldMesh)
//	{
//		if (!extSlicer.slicer)
//		{
//			printf("\n slicer is null - FieldMesh");
//			return 0;
//		}
//		//zTsSDFSlicer* slicer = extSlicer.slicer;
//		zDomain<zPoint> bb = zDomain<zPoint>(zPoint(pointDomainMin[0], pointDomainMin[1], pointDomainMin[2]), zPoint(pointDomainMax[0], pointDomainMax[1], pointDomainMax[2]));
//		extSlicer.slicer->createFieldMesh(bb, resX, resY);
//		extSlicer.updateAttributes();
//
//		extFieldMesh.mesh = extSlicer.slicer->getRawFieldMesh();
//		extFieldMesh.updateAttributes();
//
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_computePrintBlocks(zExtTsSDFSlicer& extSlicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);
//		zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
//		
//		float printPlaneSpace;
//
//		extSlicer.slicer->computePrintBlocks(printHeightDomain, printLayerWidth, raftLayerWidth, allSDFLayers, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, compFrames, compSDF);
//		
//		extSlicer.updateAttributes();
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_computeSDFSingleLayer(zExtTsSDFSlicer& extSlicer, int SDFLayerNumber, float printLayerWidth, float raftLayerWidth, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames)
//	{
//		if (!extSlicer.slicer)
//		{
//			printf("\n slicer is null - computeSDFLayer");
//			return 0;
//		}
//		if (compFrames)
//		{
//			zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);
//			zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
//			int numSDFLayers = 5;
//			extSlicer.slicer->computePrintBlocks(printHeightDomain, printLayerWidth, raftLayerWidth, false, numSDFLayers, SDFFunc_Num, SDFFunc_NumSmooth, neopreneOffset, true, false);
//		}
//
//		extSlicer.slicer->computeSDFSingleLayer(SDFLayerNumber, SDFFunc_Num, SDFFunc_NumSmooth, printLayerWidth, neopreneOffsetMin, raftLayerWidth);
//		extSlicer.updateAttributes();
//
//		return 1;
//	}
//
//	//CHECK
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_checkLayerHeight(zExtTsSDFSlicer& extSlicer, bool& check)
//	{
//		if (!extSlicer.slicer)
//		{
//			printf("\n slicer is null - CheckLayerHeights");
//			return 0;
//		}
//		bool checkSDF; bool checkGeo;
//		check = extSlicer.slicer->checkPrintLayerHeights(checkSDF, checkGeo);
//		extSlicer.updateAttributes();
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE void ext_sdf_CheckFolder(char* folderDirectoryChar, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax)
//	{
//		std::string folder(folderDirectoryChar);
//		zTsSDFSlicer slicer = zTsSDFSlicer();
//
//		zDomainFloat neopreneOffset = zDomainFloat(neopreneOffsetMin, neopreneOffsetMax);
//		zDomainFloat printHeightDomain(heightDomainMin, heightDomainMax);
//		slicer.checkPrintLayerHeights_Folder(folder, printHeightDomain, neopreneOffset);
//	}
//
//	//GET
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& rightMesh, zExtMesh& leftMesh)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		rightMesh.mesh = extSlicer.slicer->getRawRightMesh();
//		rightMesh.updateAttributes();
//		leftMesh.mesh = extSlicer.slicer->getRawLeftMesh();
//		leftMesh.updateAttributes();
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawFieldMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& fieldMesh)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		fieldMesh.mesh = extSlicer.slicer->getRawFieldMesh();
//		fieldMesh.updateAttributes();
//		
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getRawMedialGraph(zExtTsSDFSlicer& extSlicer, zExtGraph& graph)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		int numGraphs = 0;
//		graph.graph = extSlicer.slicer->getRawMedialGraph();
//		graph.updateAttributes();
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockFrames(zExtTsSDFSlicer& extSlicer, zExtPlane* frames)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		vector<zTransform> planes;
//		planes = extSlicer.slicer->getBlockFrames();
//
//		for (int i = 0; i < planes.size(); i++)
//		{
//			frames[i] = zExtPlane(planes.at(i));
//		}
//
//		return 1;
//	}
//	
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockSectionGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//
//		zObjGraphPointerArray graphs = extSlicer.slicer->getBlockSectionGraphs(numGraphs);
//		printf("\n extSlicer.sectionGraphCount : numGraphs - %i : %i", extSlicer.sectionGraphCount, numGraphs);
//
//		if (extSlicer.sectionGraphCount != numGraphs)
//		{
//			printf("\n extSlicer.sectionGraphCount != numGraphs - %i != %i", extSlicer.sectionGraphCount, numGraphs);
//			return 0;
//		}
//		for (int i = 0; i < graphs.size(); i++)
//		{
//			graphSet[i] = zExtGraph(graphs[i]);
//		}
//
//
//		
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockContourGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//		zObjGraphPointerArray graphs = extSlicer.slicer->getBlockContourGraphs(numGraphs);
//		if (extSlicer.contourGraphCount != numGraphs)
//		{
//			printf("extSlicer.contourGraphCount != numGraphs - %i != %i", extSlicer.contourGraphCount, numGraphs);
//			return 0;
//		}
//		for (int i = 0; i < graphs.size(); i++)
//		{
//			graphSet[i] = zExtGraph(graphs[i]);
//		}
//
//
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getTrimGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//		zObjGraphPointerArray graphs = extSlicer.slicer->getBlockTrimGraphs(numGraphs);
//		if (extSlicer.trimGraphCount != numGraphs)
//		{
//			printf("extSlicer.trimGraphCount != numGraphs - %i != %i", extSlicer.trimGraphCount, numGraphs);
//			return 0;
//		}
//		for (int i = 0; i < graphs.size(); i++)
//		{
//			graphSet[i] = zExtGraph(graphs[i]);
//		}
//
//
//
//		return 1;
//	}
//	
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockSectionGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//		*graphSet.graphSet = extSlicer.slicer->getBlockSectionGraphs(numGraphs);
//		graphSet.updateAttributes();
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getBlockContourGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//		*graphSet.graphSet = extSlicer.slicer->getBlockContourGraphs(numGraphs);
//		graphSet.updateAttributes();
//
//		return 1;
//	}
//	ZSPACE_TOOLSETS_INLINE int ext_sdf_getTrimGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet)
//	{
//		if (!extSlicer.slicer)
//		{
//			return 0;
//		}
//		//graphSet.graphSet = new zObjGraphPointerArray();
//		int numGraphs = 0;
//		*graphSet.graphSet = extSlicer.slicer->getBlockTrimGraphs(numGraphs);
//		graphSet.updateAttributes();
//
//		return 1;
//	}
//	
//	//EXPORT
//	ZSPACE_TOOLSETS_INLINE void ext_sdf_ExportJSON(zExtTsSDFSlicer extSlicer, char* fileCurrent, int fileCount, char* fileNew, int newCount, char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth)
//	{
//
//		string file = "";
//		string exportDir = "";
//		string exportName = "";
//		for (int i = 0; i < fileCount; i++) { file += fileCurrent[i]; }
//		for (int i = 0; i < newCount; i++) { exportDir += fileNew[i]; }
//		for (int i = 0; i < nameCount; i++) { exportName += fileName[i]; }
//
//
//		cout << "\n Export_2 entryDir" << endl;
//		cout << file << endl;
//		cout << "\n Export_2 exportDir" << endl;
//		cout << exportDir << endl;
//		cout << "\n Export_2 exportName" << endl;
//		cout << exportName << endl;
//
//		if (!extSlicer.slicer)
//		{
//			printf("\n slicer is null - export");
//		}
//
//		extSlicer.slicer->exportJSON(file, exportDir, exportName, printLayerWidth, raftLayerWidth);
//	}
//
//
//
//}