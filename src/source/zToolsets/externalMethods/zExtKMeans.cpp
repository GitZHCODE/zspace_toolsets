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


#include<headers/zToolsets/externalMethods/zExtKMeans.h>

namespace zSpace
{
	ZSPACE_TOOLSETS_INLINE void ext_KMean_computeKmeansManualInput(zTsKMeans*& kmean, double* _data, int dataCount, int strideCount, int& numCluster, int& numIterations, double* initMeans, int* outClusterID)
	{
		//create matrix from 1d data
		MatrixXf data(dataCount, strideCount);
		MatrixXf manualMeans(numCluster, strideCount);
		for (int i = 0; i < dataCount; i++)
		{
			for (int j = 0; j < strideCount; j++)
			{
				data(i, j) = _data[i * strideCount + j];
			}
		}
		for (int i = 0; i < numCluster; i++)
		{
			for (int j = 0; j < strideCount; j++)
			{
				manualMeans(i, j) = initMeans[i * strideCount + j];
			}
		}
		kmean = new zTsKMeans(data, numCluster, numIterations);
		int actualNumCluster = numCluster;
		numIterations = kmean->getKMeansClusters(actualNumCluster, manualMeans);
		numCluster = actualNumCluster;
		zIntArray clusterIDs = kmean->clusterIDS;
		for (int i = 0; i < dataCount; i++)
		{
			outClusterID[i] = clusterIDs[i];
		}
	}

	ZSPACE_TOOLSETS_INLINE void ext_KMean_computeKmeans(zTsKMeans*& kmean, double* _data, int dataCount, int strideCount, int& numCluster, int& numIterations, int initMethod, int seed1, int seed2, int* outClusterID)
	{
		//create matrix from 1d data
		MatrixXf data(dataCount, strideCount);
		for (int i = 0; i < dataCount; i++)
		{
			for (int j = 0; j < strideCount; j++)
			{
				data(i, j) = _data[i * strideCount + j];
			}
		}
		kmean = new zTsKMeans(data, numCluster, numIterations);
		int actualNumCluster = numCluster;
		if (initMethod > 1) initMethod = 1;
		zTsKMeans::initialisationMethod init = static_cast<zTsKMeans::initialisationMethod>(initMethod);
		numIterations = kmean->getKMeansClusters(actualNumCluster, init, seed1, seed2);
		numCluster = actualNumCluster;
		zIntArray clusterIDs = kmean->clusterIDS;
		for (int i = 0; i < dataCount; i++)
		{
			outClusterID[i] = clusterIDs[i];
		}
	}

	ZSPACE_TOOLSETS_INLINE void ext_KMean_getMeans(zTsKMeans*& kmean, double* outMeans)
	{
		if (!kmean)
		{
			return;
		}
		MatrixXf means =  kmean->means;
		int stride = means.cols();
		for (int i = 0; i < means.rows(); i++)
		{
			for (int j = 0; j < stride; j++)
			{
				outMeans[i * means.cols() + j] = means(i, j);
			}
		}
	}

	ZSPACE_TOOLSETS_INLINE int ext_KMean_findOptimalK(double* _data, int dataCount, int strideCount, int& numIterations,  int initMethod,  int seed1, int seed2, int min, int max, int increment, int optimalCountMethod, int* outK, double* outScore)
	{
		printf("\n ext_KMean_findOptimalK - dataCount = %i strideCount %i ", dataCount, strideCount);
		if (!_data)
		{
			printf("\n ext_KMean_findOptimalK - _data null ");
			return 0;
		}
		printf("\n min:%i max:%i optimalCountMethod:%i ", min, max, optimalCountMethod);

		MatrixXf data(dataCount, strideCount);
		for (int i = 0; i < dataCount; i++)
		{
			for (int j = 0; j < strideCount; j++)
			{
				//printf("\n _data[%i]	%f ", i * strideCount + j, _data[i * strideCount + j]);

				data(i, j) = _data[i * strideCount + j];
				//printf("	--- data(%i, %i)	%f ",i, j, data(i, j));

			}
		}	
		zTsKMeans kmean;
		kmean = zTsKMeans(data, min, numIterations);
		zTsKMeans::initialisationMethod init = static_cast<zTsKMeans::initialisationMethod>(initMethod);
		if (initMethod == 0) init = zTsKMeans::initialisationMethod::random;
		else if (initMethod == 1) init = zTsKMeans::initialisationMethod::kmeansPlusPlus;
		vector<pair<int, float>> KScorePair;
		if (optimalCountMethod > 2) optimalCountMethod == 2;
		int optimalK = min;
		/*if(optimalCountMethod == 0) optimalK = kmean.findOptimalK_Elbow(init, true, min, max, seed1, seed2, KScorePair);
		if(optimalCountMethod == 1) optimalK = kmean.findOptimalK_Elbow(init, false, min, max, seed1, seed2, KScorePair);
		if(optimalCountMethod == 2) optimalK = kmean.findOptimalK_Silhouette(init, min, max, seed1, seed2, KScorePair);*/
		
		switch (optimalCountMethod)
		{
		case 0:
			optimalK = kmean.findOptimalK_Elbow(init, true, min, max, increment, seed1, seed2, KScorePair); 
			break;
		case 1:
			optimalK = kmean.findOptimalK_Elbow(init, false, min, max, increment, seed1, seed2, KScorePair); 
			break;
		case 2:
			optimalK = kmean.findOptimalK_Silhouette(init, min, max, increment, seed1, seed2, KScorePair); 
			break;
		}
		int count = max - min +1;
		

		//outK = new int[count];
		//outScore = new double[count];
		for (int i = 0; i < KScorePair.size(); i++)
		{
			outK[i] = KScorePair[i].first;

			outScore[i] = KScorePair[i].second;
		}
		

		return optimalK;
	}

}