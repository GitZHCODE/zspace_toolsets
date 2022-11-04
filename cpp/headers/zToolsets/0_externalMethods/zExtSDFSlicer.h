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


namespace zSpace
{
	ZSPACE_TOOLSETS_EXT{
ZSPACE_TOOLSETS void EXT_SetFromJSON(zTsSDFSlicer*& slicer, char* path, int blockID, bool& isDeckBlock, bool& rightExist, bool& leftExist);
		ZSPACE_TOOLSETS void EXT_SetFromJSON2(zTsSDFSlicer*& slicer, char* path, int blockStride, int braceStride);

		ZSPACE_TOOLSETS void EXT_createFieldMesh(zTsSDFSlicer*& slicer, float* domainMin, float* domainMax, int resX, int resY);
		ZSPACE_TOOLSETS int EXT_checkLayerHeight(zTsSDFSlicer*& slicer, bool& check);


		ZSPACE_TOOLSETS int EXT_computePrintBlocks(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF);
		ZSPACE_TOOLSETS int EXT_computePrintBlocksPar(zTsSDFSlicer*& slicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF);

		ZSPACE_TOOLSETS void EXT_computeSDFLayer(zTsSDFSlicer*& slicer, int SDFLayerNumber, float printLayerWidth, float raftLayerWidth, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames);
		ZSPACE_TOOLSETS void EXT_duplicateSlicer(zTsSDFSlicer* inSlicer, zTsSDFSlicer*& outSlicer);
		//Get Methods
		ZSPACE_TOOLSETS int EXT_getRawMesh(zTsSDFSlicer* slicer, zObjMesh*& objMesh, bool rightMesh);
		ZSPACE_TOOLSETS int EXT_getRawFieldMesh(zTsSDFSlicer* slicer, zObjMesh*& objMesh);
		ZSPACE_TOOLSETS int EXT_getBlockSectionGraphs(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int EXT_getBlockContourGraphs(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graph, int numGraphs, int& outGraphsCount);
		ZSPACE_TOOLSETS int EXT_getRawMedialGraph(zTsSDFSlicer* slicer, zObjGraph*& graph);
		ZSPACE_TOOLSETS int EXT_getBlockFrames(zTsSDFSlicer* slicer, vector<zTransform>*& frames, vector<zTransform>*& rightPlanes, vector<zTransform>*& leftPlanes, int& countFrames, int& countRightPlanes, int& countLefPlanes);
		ZSPACE_TOOLSETS int EXT_getTrimGraph(zTsSDFSlicer* slicer, zObjGraphPointerArray*& graphs, int& Count);

		//Plane Data
		//ZSPACE_TOOLSETS void EXT_getPlanesData(vector<zTransform>* graph, float* outOrigin, float* outNormal, float* outXAxis, float* outYAxis);
		ZSPACE_TOOLSETS void EXT_getPlanesData(vector<zTransform>* graph, float* matrix);


		//Graph Data
		ZSPACE_TOOLSETS void EXT_getGraphsSetFromPointersVector(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray);
		ZSPACE_TOOLSETS void EXT_getGraphsSetFromVector(zObjGraphArray* graphs, zObjGraph** outGraphArray);

		ZSPACE_TOOLSETS void EXT_getGraphCounts(zObjGraph* graph, int& outvCount, int& outeCount);
		ZSPACE_TOOLSETS void EXT_getGraphData(zObjGraph* graph, float* outVPostions, float* outvColors, int* outePair, float* outeColors);
		ZSPACE_TOOLSETS void EXT_getGraphSequence(zObjGraph* graph, int* outSequence);

		//Mesh Data
		ZSPACE_TOOLSETS void EXT_getMeshCounts(zObjMesh* objMesh, int& out_vCount, int& out_fCount);
		ZSPACE_TOOLSETS void EXT_getMeshPosition(zObjMesh* objMesh, float* outVPostions, float* outVColors);
		ZSPACE_TOOLSETS void EXT_getMeshFaceCount(zObjMesh* objMesh, int* outfCounts);
		ZSPACE_TOOLSETS void EXT_getMeshFaceConnect(zObjMesh* objMesh, int* outfConnects);

		//Export JSON
		ZSPACE_TOOLSETS void EXT_ExportJSON(zTsSDFSlicer* slicer, char* fileCurrent, int fileCount, char* fileNew, int newCount, char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth);


		//Check


		//Iterate
		ZSPACE_TOOLSETS void EXT_CheckFolder(char* folderDirectoryChar, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax);

		//Draw Mesh
		ZSPACE_TOOLSETS void EXT_DrawMesh(zObjMesh* objMesh);


		
	}

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/0_externalMethods/zExtSDFSlicer.cpp>
#endif

#endif