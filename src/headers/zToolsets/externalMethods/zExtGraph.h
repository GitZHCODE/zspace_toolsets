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

#ifndef ZSPACE_EXT_TS_GRAPH_H
#define ZSPACE_EXT_TS_GRAPH_H



#pragma once

//#include "headers/zToolsets/geometry/zTsSDFSlicer.h"

#include "headers/base/zSpace_Toolsets.h"

#include <headers/zCore/base/zExtern.h>
#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include<execution>

using namespace std;

namespace zSpace
{
	struct zExtGraph
	{
		zObjGraph* graph;
		int vCount;
		int eCount;
		zExtGraph(zObjGraph* g);
		void updateAttributes();
		zIntArray getGraphSequence();
	};
	struct zExtGraphSet
	{
		zObjGraphPointerArray* graphSet;
		int graphsCount;
		void updateAttributes();
	};
	
	ZSPACE_TOOLSETS_EXT
	{

		ZSPACE_TOOLSETS void ext_graphUtil_getGraphData(zExtGraph extGraph, float* vPositions, float* vColors, int* ePairs, float* eColors);
		ZSPACE_TOOLSETS void ext_graphUtil_getGraphSequence(zExtGraph extGraph, int* outSequence);
		ZSPACE_TOOLSETS void ext_graphUtil_getGraphsSet(zExtGraphSet graphSet, zExtGraph* outGraphArray);

		////Graph Data
		//ZSPACE_TOOLSETS void ext_graphUtil_getGraphsSetFromPointersVector2(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray);
		//ZSPACE_TOOLSETS void ext_graphUtil_getGraphsSetFromVector2(zObjGraphArray* graphs, zObjGraph** outGraphArray);
		//ZSPACE_TOOLSETS void ext_graphUtil_getGraphCounts2(zObjGraph* graph, int& outvCount, int& outeCount);
		//ZSPACE_TOOLSETS void ext_graphUtil_getGraphData2(zObjGraph* graph, float* outVPostions, float* outvColors, int* outePair, float* outeColors);

		
	}

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/externalMethods/zExtGraph.cpp>
#endif

#endif