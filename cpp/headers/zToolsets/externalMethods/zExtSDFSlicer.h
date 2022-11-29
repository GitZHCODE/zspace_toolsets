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

#include <headers/zToolsets/externalMethods/zExtMesh.h>
#include <headers/zToolsets/externalMethods/zExtGraph.h>
#include <headers/zToolsets/externalMethods/zExtPlane.h>

namespace zSpace
{
	struct zExtTsSDFSlicer
	{
		zTsSDFSlicer* slicer;
		int blockID;
		int deckBlock;
		int leftPlaneExist;
		int rightPlaneExist;
		int sectionFrameCount;
		int sectionGraphCount;
		int contourGraphCount;
		int trimGraphCount;
		
		zExtTsSDFSlicer(zTsSDFSlicer* slicer);
		void updateAttributes();

		/*zExtMesh rightMesh;
		zExtMesh leftMesh;
		zExtMesh fieldMesh;
		zExtGraphSet sectionGraphsSet;
		zExtGraphSet contourGraphsSet;
		zExtGraphSet trimGraphsSet;*/
	};
	
	ZSPACE_TOOLSETS_EXT
	{
		//SET
		ZSPACE_TOOLSETS void ext_sdf_SetFromJSON(zExtTsSDFSlicer& extSlicer, char* path, int blockID);
		
		//COMPUTE
		ZSPACE_TOOLSETS int ext_sdf_createFieldMesh(zExtTsSDFSlicer& extSlicer, float* pointDomainMin, float* pointDomainMax, int resX, int resY, zExtMesh& extFieldMesh);
		ZSPACE_TOOLSETS int ext_sdf_computePrintBlocks(zExtTsSDFSlicer& extSlicer, float printLayerWidth, float raftLayerWidth, bool allSDFLayers, int& numSDFLayers, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames, bool compSDF);
		ZSPACE_TOOLSETS int ext_sdf_computeSDFSingleLayer(zExtTsSDFSlicer& extSlicer, int SDFLayerNumber, float printLayerWidth, float raftLayerWidth, int SDFFunc_Num, int SDFFunc_NumSmooth, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax, float printPlaneStep, bool compFrames);

		//CHECK
		ZSPACE_TOOLSETS int ext_sdf_checkLayerHeight(zExtTsSDFSlicer& extSlicer, bool& check);
		ZSPACE_TOOLSETS void ext_sdf_CheckFolder(char* folderDirectoryChar, float neopreneOffsetMin, float neopreneOffsetMax, float heightDomainMin, float heightDomainMax);

		//GET
		ZSPACE_TOOLSETS int ext_sdf_getRawMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& rightMesh, zExtMesh& leftMesh);
		ZSPACE_TOOLSETS int ext_sdf_getRawFieldMesh(zExtTsSDFSlicer& extSlicer, zExtMesh& fieldMesh);
		ZSPACE_TOOLSETS int ext_sdf_getRawMedialGraph(zExtTsSDFSlicer& extSlicer, zExtGraph& graph);
		ZSPACE_TOOLSETS int ext_sdf_getBlockFrames(zExtTsSDFSlicer& extSlicer, zExtPlane* frames);

		ZSPACE_TOOLSETS int ext_sdf_getBlockSectionGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet);
		ZSPACE_TOOLSETS int ext_sdf_getBlockContourGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet);
		ZSPACE_TOOLSETS int ext_sdf_getTrimGraphs(zExtTsSDFSlicer& extSlicer, zExtGraph* graphSet);

		ZSPACE_TOOLSETS int ext_sdf_getBlockSectionGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet);
		ZSPACE_TOOLSETS int ext_sdf_getBlockContourGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet);
		ZSPACE_TOOLSETS int ext_sdf_getTrimGraphSet(zExtTsSDFSlicer& extSlicer, zExtGraphSet& graphSet);

		//EXPORT
		ZSPACE_TOOLSETS void ext_sdf_ExportJSON(zExtTsSDFSlicer extSlicer, char* fileCurrent, int fileCount, char* fileNew, int newCount, char* fileName, int nameCount, float printLayerWidth, float raftLayerWidth);
	}

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/externalMethods/zExtSDFSlicer.cpp>
#endif

#endif