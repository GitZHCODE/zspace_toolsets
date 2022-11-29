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

#ifndef ZSPACE_EXT_TS_DATA_KMEANS_H
#define ZSPACE_EXT_TS_DATA_KMEANS_H



#pragma once

#include <headers/zToolsets/data/zTsKMeans.h>


namespace zSpace
{
	ZSPACE_TOOLSETS_EXT{
		ZSPACE_TOOLSETS void ext_KMean_computeKmeansManualInput(zTsKMeans*& kmean, double* _data, int dataCount, int strideCount, int& numCluster, int& numIterations, double* initMeans, int* outClusterID);
		ZSPACE_TOOLSETS void ext_KMean_computeKmeans(zTsKMeans*& kmean, double* _data, int dataCount, int strideCount, int& numCluster, int& numIterations, int initMethod, int seed1, int seed2, int* outClusterID);
		ZSPACE_TOOLSETS void ext_KMean_getMeans(zTsKMeans*& kmean, double* outMeans);
		ZSPACE_TOOLSETS int ext_KMean_findOptimalK(double* _data, int dataCount, int strideCount, int& numIterations, int initMethod, int seed1, int seed2, int min, int max, int increment, int optimalCountMethod, int* outK, double* outScore);
		
	}

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/0_externalMethods/zExtKMeans.cpp>
#endif

#endif