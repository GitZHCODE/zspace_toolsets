// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Heba Eiz <heba.eiz@zaha-hadid.com>
//

#ifndef ZSPACE_TS_DATA_KMEANS_H
#define ZSPACE_TS_DATA_KMEANS_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zCore/base/zMatrix.h>

namespace zSpace
{
	
	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsData
	*	\brief toolsets for data related utilities.
	*  @{
	*/

	/** \addtogroup zClustering
	*	\brief tool sets of clustering algorithms.
	*  @{
	*/

	/*! \class zTsKMeans
	*	\brief A tool set for doing K-Means clustering.
	*	\details Based on https://www.geeksforgeeks.org/k-means-clustering-introduction/
	*
	*	\since version 0.0.2
	*/

	/** @}*/

	/** @}*/

	/** @}*/

	

	class ZSPACE_TOOLSETS zTsKMeans
	{
		

	public:

		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------
		/*! \brief core utilities object */
		zUtilsCore coreUtils;

		/*!<\brief number of clusters.*/
		int numClusters;

		/*!<\brief number of maximum iterations.*/
		int numIterations;

		/*!<\brief minimum length of stream.*/
		double *minLength;

		/*!<\brief maximum length of stream.*/
		double *maxLength;

		/*!<\brief The sum of squared distances of samples to their closest cluster center.*/
		//double inertia;

		/*!<\brief The average of the squared distances between clusters' center.*/
		//double distortion;

		/*!<\brief mean or average data.*/
		MatrixXf means;

		/*!	\brief matrix data*/
		MatrixXf dataPoints;	
		
		/*!	\brief 2 dimensional container of cluster items*/
		vector<vector<int>> clusters;

		/*!	\brief  container of cluster id for each item*/
		vector<int> clusterIDS;
		
		/*!	\brief  container of tolerance for each dimension - if the difference between the item and the mean in that dimensions exceeds the tolerance in that dimension, create new cluster*/
		vector<float> tolerances;
	
		enum initialisationMethod
		{
			random = 0,
			kmeansPlusPlus = 1,
			manual = 2
		};


		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.2
		*/
		zTsKMeans();

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_data			- input matrix data.
		*	\since version 0.0.2
		*/
		zTsKMeans(MatrixXf &_dataPoints);
		
		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_dataPoints				- input matrix data.
		*	\param		[in]	_numClusters		- input snumber of clusters.
		*	\param		[in]	_numIterations		- input number of iterations.
		*	\since version 0.0.2
		*/
		zTsKMeans(MatrixXf &_dataPoints, int &_numClusters, int &_numIterations);

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_dataPoints				- input matrix data.
		*	\param		[in]	_numClusters		- input snumber of clusters.
		*	\param		[in]	_numIterations		- input number of iterations.
		*	\param		[in]	dimsTolerances		- the tolerance for each dimension.
		*	\since version 0.0.2
		*/
		zTsKMeans(MatrixXf& _dataPoints, int& _numClusters, int& _numIterations, zFloatArray& dimsTolerances);

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsKMeans();
		
		//--------------------------
		//----  SET METHODS
		//--------------------------

		/*! \brief This method sets the number of clusters.
		*
		*	\param		[in]	_numClusters		- input number of clusters.
		*	\since version 0.0.2
		*/
		void setNumClusters(int &_numClusters);

		/*! \brief This method sets the number of iterations.
		*
		*	\param		[in]	_numIterations		- input number of iterations.
		*	\since version 0.0.2
		*/
		void setNumIterations(int &_numIterations);

		/*! \brief This method sets the number of iterations.
		*
		*	\param		[in]	_numIterations		- input number of iterations.
		*	\since version 0.0.2
		*/
		void setTolerances(vector<float> tolerances);

		//--------------------------
		//---- CLUSTERING METHODS
		//--------------------------
			   
		/*! \brief This method computes classify the input data into input number of clusters using the K-Means Algorithm.
		*
		*	\param	[out]	actualNumClusters		- actual number of clusters after removing clusters of size 0.
		*	\return			int						- number of iterations the algorithm ran.
		*/
		int getKMeansClusters(int& actualNumClusters, initialisationMethod initMethod, int seed1, int seed2, float tolerance);
		int getKMeansClusters(int& actualNumClusters, MatrixXf manualInitMeans, float tolerance = FLT_MAX);

		int getKMeansClustersWithTolerance(int& actualNumClusters, initialisationMethod initMethod, int maxClusters, int seed1, int seed2 = 1);
		int getKMeansClustersWithTolerance_2nd(int& actualNumClusters, int maxClusters, int seed1, int seed2 = 1);

		int runKMeansClusters(int maxClusters, int& actualNumClusters, zFloatArray tolerance());


		/*
		reference from https://www.geeksforgeeks.org/elbow-method-for-optimal-value-of-k-in-kmeans/
		*/

		int findOptimalK_Elbow(initialisationMethod initMethod, bool distortionMethod, int min, int max, int increment, int seed1, int seed2, vector<pair<int, float>>& KScorePair, float tolerance);
		int findOptimalK_Silhouette(initialisationMethod initMethod, int min, int max, int increment, int seed1, int seed2, vector<pair<int, float>>& KScorePair, float tolerance);

		/*
		! \brief This method computes Distortion Score of the input. It is calculated as the average of the squared distances from the cluster centers of the respective clusters. Typically, the Euclidean distance metric is used.
		* reference from https://www.geeksforgeeks.org/elbow-method-for-optimal-value-of-k-in-kmeans/
		*	\param	[out]	actualNumClusters		- actual number of clusters after removing clusters of size 0.
		*	\return			int						- number of iterations the algorithm ran. 
		*/
		double calculateDistortion();

		double calculateIntertia();
		double calculateSilhouette();
		vector<float> matrixToVector(MatrixXf m, int num, bool row);
		int probabilitySelection(vector<float> list, int seed);

		void createEvaluationGraph(vector<float, float> pairs, string xTitle, string yTitle);
		
		
		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------
	protected:
			
		vector<int> mapUniqueClusterIDs(const std::vector<int>& v);

		/*! \brief This method initialises the means based on the minimum and maximum value in the data points.
		*
		*	\param	[out]	minVal			- input minimum value in the data.
		*	\param	[out]	maxVal			- input maximum value in the data.
		* 	\param	seed					- seed value for initialisation.

		*	\return			int				- index of cluster.
		*/
		MatrixXf intialiseMeansRandom(float &minVal, float &maxVal, int seed);

		MatrixXf intialiseMeansManual(float& minVal, float& maxVal, MatrixXf manualInput);

		/*! \brief This method initialises the means based on the minimum and maximum value in the data points.
		*
		*	\param	[out]	minVal			- input minimum value in the data.
		*	\param	[out]	maxVal			- input maximum value in the data.
		*	\return			int				- index of cluster.
		*/
		MatrixXf intialiseMeansPlusPlus(float& minVal, float& maxVal, int seed1, int seed2);

		/*! \brief This method computes the cluster index based on the least euclidean distance between input data point and mean values.
		*
		*	\param	[in]	data			- input row matrix of data.
		*	\return			int				- index of cluster.
		*/
		int getClusterIndex(MatrixXf& data, MatrixXf& means, float tolerance);

		/*! \brief This method computes the cluster index based on the least euclidean distance between input data point and mean values considering the tolerance value of each dimension.
		*
		*	\param	[in]	data			- input row matrix of data.
		*	\param	[in]	mean			- input row matrix of means.
		*	\return			int				- index of cluster.
		*/
		int getClusterIndexWithTolerance(MatrixXf& data, MatrixXf& means, bool checkTolerance = true);
		int getClusterIndex(MatrixXf& data, MatrixXf& means);

		/*! \brief This method updates the mean value of the cluster based on the input data point and cluster size.
		*
		*	\param	[in]	data			- input row matrix of datapoints.
		*	\param	[in]	mean			- input row matrix of means.
		*	\param	[in]	clusterSize		- current cluster size.
		*/
		void updateMean(MatrixXf &data, MatrixXf &mean, int clusterSize);
		

		

		
	};

	//int FindOptimalClusterCount(int min, int max, int seed1, int seed2);





}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/data/zTsKMeans.cpp>
#endif

#endif