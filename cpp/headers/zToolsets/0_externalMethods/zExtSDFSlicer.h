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

#ifndef ZSPACE_EXT_TS_GEOMETRY_SDFSLICER_H
#define ZSPACE_EXT_TS_GEOMETRY_SDFSLICER_H



#pragma once

#include <headers/zToolsets/geometry/zTsSDFSlicer.h>
#include <headers/zToolsets/0_externalMethods/zExtMeshUtility.h>
#include <headers/zToolsets/0_externalMethods/zExtGraph.h>

namespace zSpace
{
	struct zExtTsSDFSlicer
	{
		zTsSDFSlicer* slicer;
		int blockID;
		bool deckBlock;
		bool leftPlaneExist;
		bool rightPlaneExist;
		int sectionGraphCount;
		int contourGraphCount;
		int trimGraphCount;

		void updateFields()
		{
			blockID = slicer->blockId;
			deckBlock = slicer->onDeckBlock();
			leftPlaneExist = slicer->leftPlaneExists;
			rightPlaneExist = slicer->rightPlaneExists;
			sectionGraphCount = slicer->o_sectionGraphs.size();
			contourGraphCount = slicer->o_contourGraphs.size();
			trimGraphCount = slicer->o_trimGraphs.size();
		}

		/*zExtMesh rightMesh;
		zExtMesh leftMesh;
		zExtMesh fieldMesh;
		zExtGraphSet sectionGraphsSet;
		zExtGraphSet contourGraphsSet;
		zExtGraphSet trimGraphsSet;*/
	};
	
	ZSPACE_TOOLSETS_EXT
	{
		ZSPACE_TOOLSETS void ext_sdf_SetFromJSON(zExtTsSDFSlicer& extSlicer, char* path, int blockID);
		ZSPACE_TOOLSETS int ext_sdf_createFieldMesh(zExtTsSDFSlicer& extSlicer, float* domainMin, float* domainMax, int resX, int resY, zExtMesh& extFieldMesh);
		ZSPACE_TOOLSETS int ext_sdf_getRawMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& rightMesh, zExtMesh& leftMesh);
		ZSPACE_TOOLSETS int ext_sdf_getBlockSectionGraphs(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet);




		ZSPACE_TOOLSETS void ext_sdf_SetFromJSON3(zTsSDFSlicer*& slicer, char* path, int blockID, bool& isDeckBlock);
		ZSPACE_TOOLSETS void ext_sdf_SetFromJSON2(zTsSDFSlicer*& slicer, char* path, int blockStride, int braceStride);

		ZSPACE_TOOLSETS void ext_sdf_createFieldMesh2(zTsSDFSlicer*& slicer, float* domainMin, float* domainMax, int resX, int resY);
		ZSPACE_TOOLSETS int ext_sdf_checkLayerHeight2(zTsSDFSlicer*& slicer, bool& check);


		ZSPACE_TOOLSETS int ext_sdf_computePrintBlocks2(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF);
		ZSPACE_TOOLSETS int ext_sdf_computePrintBlocksPar2(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF);

		ZSPACE_TOOLSETS void ext_sdf_computeSDFLayer2(zTsSDFSlicer*& slicer, int SDFLayerNumber, float printLayerWidth, float raftLayerWidth, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames);
		ZSPACE_TOOLSETS void ext_sdf_duplicateSlicer2(zTsSDFSlicer* inSlicer, zTsSDFSlicer*& outSlicer);
		//Get Methods
		ZSPACE_TOOLSETS int ext_sdf_getRawMesh2(zTsSDFSlicer* slicer, zObjMesh*& objMesh, bool rightMesh);
		ZSPACE_TOOLSETS int ext_sdf_getRawFieldMesh2(zTsSDFSlicer* slicer, zObjMesh*& objMesh);
		ZSPACE_TOOLSETS int ext_sdf_getBlockSectionGraphs2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int ext_sdf_getBlockContourGraphs2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int ext_sdf_getRawMedialGraph2(zTsSDFSlicer* slicer, zObjGraph*& graph);
		ZSPACE_TOOLSETS int ext_sdf_getBlockFrames2(zTsSDFSlicer* slicer, vector<zTransform>*& planes, int& count);
		ZSPACE_TOOLSETS int ext_sdf_getTrimGraph2(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graphs, int& Count);

		////Plane Data
		////ZSPACE_TOOLSETS void ext_sdf_getPlanesData(vector<zTransform>* graph, float* outOrigin, float* outNormal, float* outXAxis, float* outYAxis);
		//ZSPACE_TOOLSETS void ext_sdf_getPlanesData2(vector<zTransform>* graph, float* matrix);


		////Graph Data
		//ZSPACE_TOOLSETS void ext_sdf_getGraphsSetFromPointersVector(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray);
		//ZSPACE_TOOLSETS void ext_sdf_getGraphsSetFromVector(zObjGraphArray* graphs, zObjGraph** outGraphArray);

		//ZSPACE_TOOLSETS void ext_sdf_getGraphCounts(zObjGraph* graph, int& outvCount, int& outeCount);
		//ZSPACE_TOOLSETS void ext_sdf_getGraphData(zObjGraph* graph, float* outVPostions, float* outvColors, int* outePair, float* outeColors);
		//ZSPACE_TOOLSETS void ext_sdf_getGraphSequence(zObjGraph* graph, int* outSequence);

		////Mesh Data
		//ZSPACE_TOOLSETS void ext_sdf_getMeshCounts(zObjMesh* objMesh, int& out_vCount, int& out_fCount);
		//ZSPACE_TOOLSETS void ext_sdf_getMeshPosition(zObjMesh* objMesh, float* outVPostions, float* outVColors);
		//ZSPACE_TOOLSETS void ext_sdf_getMeshFaceCount(zObjMesh* objMesh, int* outfCounts);
		//ZSPACE_TOOLSETS void ext_sdf_getMeshFaceConnect(zObjMesh* objMesh, int* outfConnects);

		////Export JSON
		//ZSPACE_TOOLSETS void ext_sdf_ExportJSON(zTsSDFSlicer* slicer, char* fileCurrent, int fileCount, char* fileNew, int newCount, char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth);


		////Check


		////Iterate
		//ZSPACE_TOOLSETS void ext_sdf_CheckFolder(char* folderDirectoryChar, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax);

		////Draw Mesh
		//ZSPACE_TOOLSETS void ext_sdf_DrawMesh(zObjMesh* objMesh);



	}

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/0_externalMethods/zExtSDFSlicer.cpp>
#endif

#endif