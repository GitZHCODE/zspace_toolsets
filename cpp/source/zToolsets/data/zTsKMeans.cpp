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


#include<headers/zToolsets/data/zTsKMeans.h>

namespace zSpace
{


	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsKMeans::zTsKMeans()
	{
		numClusters = 2;
		numIterations = 100;
	}

	ZSPACE_TOOLSETS_INLINE zTsKMeans::zTsKMeans(MatrixXf& _dataPoints)
	{
		dataPoints = _dataPoints;


		numClusters = 2;
		numIterations = 100;
	}

	ZSPACE_TOOLSETS_INLINE zTsKMeans::zTsKMeans(MatrixXf& _dataPoints, int& _numClusters, int& _numIterations)
	{
		dataPoints = _dataPoints;

		numClusters = _numClusters;
		numIterations = _numIterations;
	}
	ZSPACE_TOOLSETS_INLINE zTsKMeans::zTsKMeans(MatrixXf& _dataPoints, int& _numClusters, int& _numIterations, zFloatArray& dimsTolerances)
	{
		dataPoints = _dataPoints;
		numClusters = _numClusters;
		numIterations = _numIterations;
		tolerances = dimsTolerances;
	}

	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsKMeans::~zTsKMeans() {}

	//----  SET METHODS

	ZSPACE_TOOLSETS_INLINE void zTsKMeans::setNumClusters(int& _numClusters)
	{
		numClusters = _numClusters;
	}
	ZSPACE_TOOLSETS_INLINE void zTsKMeans::setNumIterations(int& _numIterations)
	{
		numIterations = _numIterations;
	}

	ZSPACE_TOOLSETS_INLINE void zTsKMeans::setTolerances(vector<float> dimsTolerances)
	{
		tolerances = dimsTolerances;

	}

	//---- CLUSTERING METHODS
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getKMeansClusters(int& actualNumClusters, MatrixXf manualInitMeans, float tolerance)
	{
		numClusters = manualInitMeans.rows();
		int numRows = dataPoints.rows();
		int numCols = dataPoints.cols();

		// get min max value of the datapoints
		int minIndex;

		float minVal = dataPoints.minCoeff();

		int maxIndex;
		float maxVal = dataPoints.maxCoeff();

		// Initialise means
		MatrixXf tempMeans = manualInitMeans;

		// Initialise container to store cluster index per data point
		//vector<int> clusterIDS;
		for (int i = 0; i < numRows; i++)
		{
			clusterIDS.push_back(0);
		}


		// Initialise container to store item indicies per cluster
		vector<vector<int>> tempClusters;

		for (int i = 0; i < numClusters; i++)
		{
			vector<int> temp;
			tempClusters.push_back(temp);
		}


		// compute means

		int numIters = 0;
		bool exit = false;
		for (int i = 0; i < numIterations; i++)
		{
			numIters = i;

			exit = true;
			//exit = false;

			for (int j = 0; j < numRows; j++)
			{
				MatrixXf data = dataPoints.row(j);

				int clusterID = getClusterIndex(data, tempMeans, tolerance);


				// check if data point changed cluster
				if (clusterID != clusterIDS[j]) exit = false;

				clusterIDS[j] = clusterID;



				// update clusters
				tempClusters[clusterID].push_back(j);


			}

			// update mean			
			for (int j = 0; j < numClusters; j++)
			{
				MatrixXf  mean = tempMeans.row(j);

				for (int l = 0; l < tempClusters[j].size(); l++)
				{
					int dataId = tempClusters[j][l];
					MatrixXf data = dataPoints.row(dataId);

					updateMean(data, mean, l + 1);
				}

				tempMeans.row(j) = mean;
			}



			if (exit) break;
			else
			{
				// clear clusters
				for (int j = 0; j < numClusters; j++)
				{
					tempClusters[j].clear();
				}
			}


		}


		// remove cluster with zero elements
		actualNumClusters = numClusters;
		clusters.clear();

		for (int i = 0; i < numClusters; i++)
		{
			if (tempClusters[i].size() != 0)
			{
				clusters.push_back(tempClusters[i]);
			}
			else actualNumClusters--;

		}

		means = MatrixXf(actualNumClusters, tempMeans.cols());
		int id = 0;
		for (int i = 0; i < tempMeans.rows(); i++)
		{
			if (tempClusters[i].size() != 0)
			{
				for (int j = 0; j < tempMeans.cols(); j++)
				{
					means(id, j) = tempMeans(i, j);
				}

				id++;
			}
		}
		//calculateDistortion();
		//calculateIntertia();
		return numIters;

	}


	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getKMeansClusters(int& actualNumClusters, initialisationMethod initMethod, int seed1, int seed2, float tolerance)
	{
		int numRows = dataPoints.rows();
		int numCols = dataPoints.cols();

		// get min max value of the datapoints
		int minIndex;

		float minVal = dataPoints.minCoeff();

		int maxIndex;
		float maxVal = dataPoints.maxCoeff();
		clusterIDS.clear();

		// Initialise means
		MatrixXf tempMeans;
		switch (initMethod)
		{
		case initialisationMethod::random:
			tempMeans = intialiseMeansRandom(minVal, maxVal, seed1);
			break;

		case initialisationMethod::kmeansPlusPlus:
			tempMeans = intialiseMeansPlusPlus(minVal, maxVal, seed1, seed2);
			break;
		}


		// Initialise container to store cluster index per data point
		//vector<int> clusterIDS;
		for (int i = 0; i < numRows; i++)
		{
			clusterIDS.push_back(0);
		}


		// Initialise container to store item indicies per cluster
		vector<vector<int>> tempClusters;
		for (int i = 0; i < numClusters; i++)
		{
			vector<int> temp;
			tempClusters.push_back(temp);
		}


		// compute means

		int numIters = 0;
		bool exit = false;
		for (int i = 0; i < numIterations; i++)
		{
			numIters = i;

			exit = true;
			//exit = false;


			for (int j = 0; j < numRows; j++)
			{
				MatrixXf data = dataPoints.row(j);

				int clusterID = getClusterIndex(data, tempMeans, tolerance);
				//cout << "numRows:" << numRows << endl;
				//cout << "tempMeans.rows():" << tempMeans.rows() << endl;

				if (clusterID == -1)
				{
					//MatrixXf newMatrix(tempMeans.rows() + 1, dataPoints.cols());
					////newMatrix.block(0, 0, tempMeans.rows(), tempMeans.cols()) = tempMeans;

					//for (int l = 0; l < tempMeans.rows(); l++)
					//{
					//	newMatrix.row(l) = tempMeans.row(l);
					//}
					//newMatrix.row(tempMeans.rows()) = data;

					//tempMeans = MatrixXf(newMatrix.rows(), dataPoints.cols());
					//tempMeans = newMatrix;

					//clusterID = getClusterIndex(data, tempMeans, tolerance);

					//numClusters++;

					tempMeans.conservativeResize(tempMeans.rows() + 1, tempMeans.cols());
					tempMeans.row(tempMeans.rows() - 1) = data;

					clusterID = tempMeans.rows() - 1;/*getClusterIndex(data, tempMeans, tolerance);*/


					tempClusters.push_back(vector<int>());

					numClusters++;
				}

				// check if data point changed cluster
				if (clusterID != clusterIDS[j]) exit = false;

				clusterIDS[j] = clusterID;



				// update clusters
				//cout << "clusterID:" << clusterID << endl;

				tempClusters[clusterID].push_back(j);


			}

			// update mean			
			for (int j = 0; j < numClusters; j++)
			{
				MatrixXf  mean = tempMeans.row(j);

				for (int l = 0; l < tempClusters[j].size(); l++)
				{
					int dataId = tempClusters[j][l];
					//cout << "dataId:" << dataId << endl;

					MatrixXf data = dataPoints.row(dataId);

					updateMean(data, mean, l + 1);
				}

				tempMeans.row(j) = mean;
			}



			if (exit) break;
			else
			{
				// clear clusters
				for (int j = 0; j < numClusters; j++)
				{
					tempClusters[j].clear();
				}
			}


		}


		// remove cluster with zero elements
		actualNumClusters = numClusters;
		clusters.clear();

		for (int i = 0; i < numClusters; i++)
		{
			if (tempClusters[i].size() != 0)
			{
				clusters.push_back(tempClusters[i]);
			}
			else actualNumClusters--;

		}

		means = MatrixXf(actualNumClusters, tempMeans.cols());
		int id = 0;
		for (int i = 0; i < tempMeans.rows(); i++)
		{
			if (tempClusters[i].size() != 0)
			{
				for (int j = 0; j < tempMeans.cols(); j++)
				{
					means(id, j) = tempMeans(i, j);
				}

				id++;
			}
		}

		clusterIDS = mapUniqueClusterIDs(clusterIDS);

		//calculateDistortion();
		//calculateIntertia();
		return numIters;

	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getKMeansClustersWithTolerance(int& actualNumClusters, initialisationMethod initMethod, int maxClusters, int seed1, int seed2)
	{
		int numRows = dataPoints.rows();
		int numCols = dataPoints.cols();

		// get min max value of the datapoints
		float minVal = dataPoints.minCoeff();
		float maxVal = dataPoints.maxCoeff();
		clusterIDS.clear();

		// Initialise means
		MatrixXf tempMeans;
		switch (initMethod)
		{
		case initialisationMethod::random:
			tempMeans = intialiseMeansRandom(minVal, maxVal, seed1);
			break;

		case initialisationMethod::kmeansPlusPlus:
			tempMeans = intialiseMeansPlusPlus(minVal, maxVal, seed1, seed2);
			break;
		}


		// Initialize container to store cluster index per data point
		for (int i = 0; i < numRows; i++)
		{
			clusterIDS.push_back(0);
		}


		// Initialize container to store item indicies per cluster
		vector<vector<int>> tempClusters;
		for (int i = 0; i < numClusters; i++)
		{
			vector<int> temp;
			tempClusters.push_back(temp);
		}


		// compute means

		int numIters = 0;
		bool exit = false;

		int maxClustersAdded = 0;
		for (int i = 0; i < numIterations; i++)
		{
			numIters = i;

			exit = true;
			//exit = false;


			for (int j = 0; j < numRows; j++)
			{
				MatrixXf data = dataPoints.row(j);
				int clusterID = numClusters < maxClusters ? getClusterIndexWithTolerance(data, tempMeans) : getClusterIndexWithTolerance(data, tempMeans, false);
				if (clusterID == -1)
				{

					printf("\n create new cluster %i", numClusters + 1);
					tempMeans.conservativeResize(tempMeans.rows() + 1, tempMeans.cols());
					tempMeans.row(tempMeans.rows() - 1) = data;

					clusterID = tempMeans.rows() - 1;/*getClusterIndex(data, tempMeans, tolerance);*/


					tempClusters.push_back(vector<int>());

					numClusters++;
				}



				// check if data point changed cluster
				if (clusterID != clusterIDS[j]) exit = false;

				clusterIDS[j] = clusterID;



				// update clusters
				//cout << "clusterID:" << clusterID << endl;

				tempClusters[clusterID].push_back(j);


			}

			// update mean			
			for (int j = 0; j < numClusters; j++)
			{
				MatrixXf  mean = tempMeans.row(j);

				for (int l = 0; l < tempClusters[j].size(); l++)
				{
					int dataId = tempClusters[j][l];
					//cout << "dataId:" << dataId << endl;

					MatrixXf data = dataPoints.row(dataId);

					updateMean(data, mean, l + 1);
				}

				tempMeans.row(j) = mean;
			}



			if (exit) break;
			else
			{
				// clear clusters
				for (int j = 0; j < numClusters; j++)
				{
					tempClusters[j].clear();
				}
			}


		}


		// remove cluster with zero elements
		actualNumClusters = numClusters;
		clusters.clear();

		for (int i = 0; i < numClusters; i++)
		{
			if (tempClusters[i].size() != 0)
			{
				clusters.push_back(tempClusters[i]);
			}
			else actualNumClusters--;

		}

		means = MatrixXf(actualNumClusters, tempMeans.cols());
		int id = 0;
		for (int i = 0; i < tempMeans.rows(); i++)
		{
			if (tempClusters[i].size() != 0)
			{
				for (int j = 0; j < tempMeans.cols(); j++)
				{
					means(id, j) = tempMeans(i, j);
				}

				id++;
			}
		}

		clusterIDS = mapUniqueClusterIDs(clusterIDS);

		//calculateDistortion();
		//calculateIntertia();
		return numIters;
	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getKMeansClustersWithTolerance(int& actualNumClusters, initialisationMethod initMethod, int maxClusters, MatrixXf manualInitMeans, int seed1, int seed2)
	{
		int numRows = dataPoints.rows();
		int numCols = dataPoints.cols();

		// get min max value of the datapoints
		float minVal = dataPoints.minCoeff();
		float maxVal = dataPoints.maxCoeff();
		clusterIDS.clear();

		// Initialise means
		MatrixXf tempMeans;
		switch (initMethod)
		{
		case initialisationMethod::random:
			tempMeans = intialiseMeansRandom(minVal, maxVal, seed1);
			break;

		case initialisationMethod::kmeansPlusPlus:
			tempMeans = intialiseMeansPlusPlus(minVal, maxVal, seed1, seed2);
			break;

		case initialisationMethod::manual:
			tempMeans = manualInitMeans;
				break;
		}

		// Initialize container to store cluster index per data point
		for (int i = 0; i < numRows; i++)
		{
			clusterIDS.push_back(0);
		}


		// Initialize container to store item indicies per cluster
		vector<vector<int>> tempClusters;
		for (int i = 0; i < numClusters; i++)
		{
			vector<int> temp;
			tempClusters.push_back(temp);
		}


		// compute means

		int numIters = 0;
		bool exit = false;

		int maxClustersAdded = 0;
		for (int i = 0; i < numIterations; i++)
		{
			numIters = i;

			exit = true;
			//exit = false;


			for (int j = 0; j < numRows; j++)
			{
				MatrixXf data = dataPoints.row(j);
				int clusterID = numClusters < maxClusters ? getClusterIndexWithTolerance(data, tempMeans) : getClusterIndexWithTolerance(data, tempMeans, false);
				if (clusterID == -1)
				{

					printf("\n create new cluster %i", numClusters + 1);
					tempMeans.conservativeResize(tempMeans.rows() + 1, tempMeans.cols());
					tempMeans.row(tempMeans.rows() - 1) = data;

					clusterID = tempMeans.rows() - 1;/*getClusterIndex(data, tempMeans, tolerance);*/


					tempClusters.push_back(vector<int>());

					numClusters++;
				}



				// check if data point changed cluster
				if (clusterID != clusterIDS[j]) exit = false;

				clusterIDS[j] = clusterID;



				// update clusters
				//cout << "clusterID:" << clusterID << endl;

				tempClusters[clusterID].push_back(j);


			}

			// update mean			
			for (int j = 0; j < numClusters; j++)
			{
				MatrixXf  mean = tempMeans.row(j);

				for (int l = 0; l < tempClusters[j].size(); l++)
				{
					int dataId = tempClusters[j][l];
					//cout << "dataId:" << dataId << endl;

					MatrixXf data = dataPoints.row(dataId);

					updateMean(data, mean, l + 1);
				}

				tempMeans.row(j) = mean;
			}



			if (exit) break;
			else
			{
				// clear clusters
				for (int j = 0; j < numClusters; j++)
				{
					tempClusters[j].clear();
				}
			}


		}


		// remove cluster with zero elements
		actualNumClusters = numClusters;
		clusters.clear();

		for (int i = 0; i < numClusters; i++)
		{
			if (tempClusters[i].size() != 0)
			{
				clusters.push_back(tempClusters[i]);
			}
			else actualNumClusters--;

		}

		means = MatrixXf(actualNumClusters, tempMeans.cols());
		int id = 0;
		for (int i = 0; i < tempMeans.rows(); i++)
		{
			if (tempClusters[i].size() != 0)
			{
				for (int j = 0; j < tempMeans.cols(); j++)
				{
					means(id, j) = tempMeans(i, j);
				}

				id++;
			}
		}

		clusterIDS = mapUniqueClusterIDs(clusterIDS);

		//calculateDistortion();
		//calculateIntertia();
		return numIters;



	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getKMeansClustersWithTolerance_2nd(int& actualNumClusters, int maxClusters, int seed1, int seed2)
	{
		int numRows = dataPoints.rows();
		int numCols = dataPoints.cols();

		// get min max value of the datapoints
		float minVal = dataPoints.minCoeff();
		float maxVal = dataPoints.maxCoeff();
		clusterIDS.clear();
		// Initialise means
		MatrixXf tempMeans = means;



		// Initialize container to store cluster index per data point
		for (int i = 0; i < numRows; i++)
		{
			clusterIDS.push_back(0);
		}


		// Initialize container to store item indicies per cluster
		vector<vector<int>> tempClusters;
		for (int i = 0; i < numClusters; i++)
		{
			vector<int> temp;
			tempClusters.push_back(temp);
		}


		// compute means

		int numIters = 0;
		bool exit = false;

		int maxClustersAdded = 0;
		for (int i = 0; i < numIterations; i++)
		{
			numIters = i;

			exit = true;
			//exit = false;


			for (int j = 0; j < numRows; j++)
			{
				MatrixXf data = dataPoints.row(j);
				int clusterID = numClusters < maxClusters ? getClusterIndexWithTolerance(data, tempMeans) : getClusterIndexWithTolerance(data, tempMeans, false);
				if (clusterID == -1)
				{

					printf("\n create new cluster %i", numClusters + 1);
					tempMeans.conservativeResize(tempMeans.rows() + 1, tempMeans.cols());
					tempMeans.row(tempMeans.rows() - 1) = data;

					clusterID = tempMeans.rows() - 1;/*getClusterIndex(data, tempMeans, tolerance);*/


					tempClusters.push_back(vector<int>());

					numClusters++;
				}



				// check if data point changed cluster
				if (clusterID != clusterIDS[j]) exit = false;

				clusterIDS[j] = clusterID;



				// update clusters
				//cout << "clusterID:" << clusterID << endl;

				tempClusters[clusterID].push_back(j);


			}

			// update mean			
			for (int j = 0; j < numClusters; j++)
			{
				MatrixXf  mean = tempMeans.row(j);

				for (int l = 0; l < tempClusters[j].size(); l++)
				{
					int dataId = tempClusters[j][l];
					//cout << "dataId:" << dataId << endl;

					MatrixXf data = dataPoints.row(dataId);

					updateMean(data, mean, l + 1);
				}

				tempMeans.row(j) = mean;
			}



			if (exit) break;
			else
			{
				// clear clusters
				for (int j = 0; j < numClusters; j++)
				{
					tempClusters[j].clear();
				}
			}


		}


		// remove cluster with zero elements
		actualNumClusters = numClusters;
		clusters.clear();

		for (int i = 0; i < numClusters; i++)
		{
			if (tempClusters[i].size() != 0)
			{
				clusters.push_back(tempClusters[i]);
			}
			else actualNumClusters--;

		}

		means = MatrixXf(actualNumClusters, tempMeans.cols());
		int id = 0;
		for (int i = 0; i < tempMeans.rows(); i++)
		{
			if (tempClusters[i].size() != 0)
			{
				for (int j = 0; j < tempMeans.cols(); j++)
				{
					means(id, j) = tempMeans(i, j);
				}

				id++;
			}
		}
		printf("\n means.rows(): %i", means.rows());
		clusterIDS = mapUniqueClusterIDs(clusterIDS);

		//calculateDistortion();
		//calculateIntertia();
		return numIters;
	}

	int zTsKMeans::runKMeansClusters(int maxClusters, int& actualNumClusters, zFloatArray tolerance())
	{
		return 0;
	}

	//---- PROTECTED METHODS

	ZSPACE_TOOLSETS_INLINE vector<int> zTsKMeans::mapUniqueClusterIDs(const std::vector<int>& v)
	{
		std::map<int, std::size_t> m;
		for (auto e : v)
		{
			m[e];
		}
		int index = 0;
		for (auto& [key, value] : m)
		{
			value = index++;
		}
		std::vector<int> res;
		for (auto e : v)
		{
			res.push_back(m[e]);
		}
		return res;
	}
	ZSPACE_TOOLSETS_INLINE MatrixXf zTsKMeans::intialiseMeansRandom(float& minVal, float& maxVal, int seed)
	{
		printf("\n init random 0");
		MatrixXf out(numClusters, dataPoints.cols());
		printf("\n init random 1");

		// to generate different random number every time the program runs
		srand(seed);

		vector<double> randNumbers;

		for (int i = 0; i < out.rows() * out.cols(); i++)
		{
			double v = coreUtils.randomNumber_double(minVal, maxVal);
			randNumbers.push_back(v);
		}
		printf("\n init random 3");

		int id = 0;

		for (int i = 0; i < out.rows(); i++)
		{

			for (int j = 0; j < out.cols(); j++)
			{
				out(i, j) = randNumbers[id];
				id++;
			}
		}
		printf("\n init random 4");

		return out;
	}
	ZSPACE_TOOLSETS_INLINE MatrixXf zTsKMeans::intialiseMeansManual(float& minVal, float& maxVal, MatrixXf manualInitMeans)
	{
		numClusters = manualInitMeans.rows();
		return manualInitMeans;
	}


	ZSPACE_TOOLSETS_INLINE MatrixXf zTsKMeans::intialiseMeansPlusPlus(float& minVal, float& maxVal, int seed1, int seed2)
	{
		MatrixXf out(numClusters, dataPoints.cols());
		MatrixXf current(1, dataPoints.cols());
		//Randomly select the first centroid as one data item
		srand(seed1);
		int index = coreUtils.randomNumber(0, dataPoints.rows() - 1);
		/*for (int i = 0; i < dataPoints.cols(); i++)
		{
			out[0, i] = dataPoints[index, i];
			current[0, i] = dataPoints[index, i];
		}*/
		out.row(0) = dataPoints.row(index); // is this correct to access row?

		//create a list to contain all the distances to the centroid we found so far
		MatrixXf pDist(dataPoints.rows(), numClusters);

		for (int k = 1; k < numClusters; k++)
		{
			//find the distance with the previous cluster index (to avoid calculating the same distance with the same centroids multiple times to find the smallest)
			for (int r = 0; r < dataPoints.rows(); r++) //for each data point
			{
				MatrixXf p(1, dataPoints.cols());
				p.row(0) = dataPoints.row(r);
				MatrixXf c(1, dataPoints.cols());
				c.row(0) = out.row(k - 1);
				float d = coreUtils.getEuclideanDistance(p, c);
				pDist(r, k - 1) = d;
			}

			//create a list with the smallest distance for each datapoint
			vector<float> pSmallestDist;
			vector<int> pSmallestIndex;

			for (int r = 0; r < dataPoints.rows(); r++) //for each data point
			{
				float distance = FLT_MAX;
				int pointIndex = 0;
				for (int ki = 0; ki < k; ki++)
				{
					if (pDist(r, ki) < distance)
					{
						distance = pDist(r, ki);
						pointIndex = r;
					}
				}
				pSmallestDist.push_back(distance);
				pSmallestIndex.push_back(pointIndex);

			}

			int newMeanIndex = pSmallestIndex[probabilitySelection(pSmallestDist, seed2 + k)];
			out.row(k) = dataPoints.row(newMeanIndex);

			//for (int r = 0; r < dataPoints.rows(); r++) //for each data point
			//{
			//	MatrixXf p(1, dataPoints.cols());
			//	p.row(0) = dataPoints.row(r);
			//	for (int ki = 0; ki < k; ki++)
			//	{
			//		MatrixXf c(1, dataPoints.cols());
			//		c.row(0) = out.row(ki);
			//		float d = coreUtils.getEuclideanDistance(p, c);
			//		if (d > distance)
			//		{
			//			distance = d;
			//			pointIndex = r;
			//		}
			//	}
			//}
			//out.row(k) = dataPoints.row(pointIndex);

		}
		return out;
	}

	ZSPACE_TOOLSETS_INLINE double zTsKMeans::calculateDistortion()
	{
		double sum = 0;
		for (int r = 0; r < means.rows(); r++)
		{
			for (int j = r + 1; j < means.rows(); j++)
			{
				MatrixXf p(1, means.cols());
				p.row(0) = means.row(r);
				MatrixXf c(1, means.cols());
				c.row(0) = means.row(j);
				float d = coreUtils.getEuclideanDistance(p, c);
				sum += d * d;
			}
		}
		return sum / means.rows();
	}
	ZSPACE_TOOLSETS_INLINE double zTsKMeans::calculateIntertia()
	{
		double sum = 0;
		for (int r = 0; r < means.rows(); r++) //for each cluster mean 
		{
			MatrixXf p(1, means.cols());
			p.row(0) = means.row(r);
			for (int index : clusters[r]) //for each sample point in the cluster 
			{
				MatrixXf c(1, dataPoints.cols());
				c.row(0) = dataPoints.row(index);
				float d = coreUtils.getEuclideanDistance(p, c);
				sum += d * d;
			}
		}
		return sum;
	}

	ZSPACE_TOOLSETS_INLINE double zTsKMeans::calculateSilhouette()
	{
		//The silhouette value measures how similar a point is to its own cluster (cohesion) compared to other clusters (separation). Ref. https://medium.com/analytics-vidhya/how-to-determine-the-optimal-k-for-k-means-708505d204eb and https://www.geeksforgeeks.org/silhouette-algorithm-to-determine-the-optimal-value-of-k/
		//The maximum 
		//average silhouette of k
		//double avgSilhouette = 0;
		double SilhouetteScore = 0;

		//distance between all data points with each others
		MatrixXf pointsDistances(dataPoints.rows(), dataPoints.rows());
		pointsDistances.fill(std::numeric_limits<float>::max());

		for (int i = 0; i < dataPoints.rows(); i++)
		{

			MatrixXf pi(1, dataPoints.cols());
			pi.row(0) = dataPoints.row(i); // point at i
			for (int j = 0; j < dataPoints.rows(); j++)
			{

				if (i != j && pointsDistances(i, j) == std::numeric_limits<float>::max())
				{

					MatrixXf pj(1, dataPoints.cols());
					pj.row(0) = dataPoints.row(j); // point at j
					float d = coreUtils.getEuclideanDistance(pi, pj);
					pointsDistances(i, j) = d;
					pointsDistances(j, i) = d;
				}
				if (i == j)
				{

					pointsDistances(i, j) = 0;
					pointsDistances(j, i) = 0;
				}
			}
		}
		/*printf("\n dataPoints.rows() %i", dataPoints.rows());
		printf("\n pointsDistances.rows() %i", pointsDistances.rows());
		printf("\n pointsDistances.cols() %i", pointsDistances.cols());
		printf("\n s-5");*/
		printf("\n i | ai | bi | si ");

		for (int i = 0; i < dataPoints.rows(); i++)
		{
			//cluster of point i 
			int id = clusterIDS[i];
			double si = 0;
			printf("\n id %i", id);
			//printf("\n clusters.size() %i", clusters.size());
			//printf("\n clusters[id].size() %i", clusters[id].size());

			if (clusters[id].size() > 1)
			{
				double ai = 0; //measure of similarity of the point i to its own cluster. It is measured as the average distance of i from other points in the cluster. 
				for (int j = 0; j < clusters[id].size(); j++)
				{
					int k = clusters[id][j];
					ai += pointsDistances(i, k);
				}
				ai = ai / (clusters[id].size() - 1);

				vector<pair<float, int>> distClusterPair;
				for (int m = 0; m < pointsDistances.cols(); m++)
				{
					distClusterPair.push_back(make_pair(pointsDistances(i, m), clusterIDS[m])); //pair.first = distance to m, pair.second = which clusterId m belongs to
				}

				std::sort(distClusterPair.begin(), distClusterPair.end());

				int datapointId = 0;

				while (distClusterPair[datapointId].second == id && datapointId < distClusterPair.size())
				{
					datapointId++;
				}

				printf("\n datapointId %i", datapointId);
				int closestClusterId = clusterIDS[datapointId];
				/*while (clusters[closestClusterId].size() == 00 && datapointId < distClusterPair.size())
				{
					datapointId++;
				}*/
				//closestClusterId = clusterIDS[datapointId];
				printf("\n closestClusterId %i", closestClusterId);
				printf("\n clusters[closestClusterId].size() %i", clusters[closestClusterId].size());

				//the closest cluster is closestClusterId
				double bi = 0; //measure of the average dissimilarity to the closest cluster which is not it’s cluster 
				for (int j = 0; j < clusters[closestClusterId].size(); j++)
				{
					printf("\n j %i clusters[closestClusterId][j] %i", j, clusters[closestClusterId][j]);

					int k = clusters[closestClusterId][j];
					bi += pointsDistances(i, k);
				}
				bi = bi / (clusters[closestClusterId].size());

				si = (bi - ai) / std::max(ai, bi);

				printf("\n %i | %f | %f | %f ", i, ai, bi, si);


			}
			SilhouetteScore += si;
		}
		SilhouetteScore = SilhouetteScore / dataPoints.rows();

		printf("\n avg %f ", SilhouetteScore);


		return SilhouetteScore;
	}

	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getClusterIndex(MatrixXf& data, MatrixXf& means, float tolerance)
	{
		double minDist = 10000000;
		int out = -1;

		for (int i = 0; i < means.rows(); i++)
		{
			MatrixXf mean = means.row(i);

			double dist = coreUtils.getEuclideanDistance(data, mean);

			if (dist < minDist && dist < tolerance)
			{
				minDist = dist;
				out = i;
			}
		}

		return out;
	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getClusterIndexWithTolerance(MatrixXf& data, MatrixXf& means, bool checkTolerance)
	{
		double minDist = DBL_MAX;
		int out = -1;

		for (int i = 0; i < means.rows(); i++)
		{
			MatrixXf mean = means.row(i);


			//find the distance in each dimension separately (the values should NOT be normalized)
			bool allDimsInTolerance = true;
			if (checkTolerance)
			{
				for (int j = 0; j < mean.cols(); j++)
				{
					float distInDim = std::abs(mean(0, j) - data(0, j));
					if (distInDim > tolerances[j])
					{
						allDimsInTolerance = false;
						//printf("\n out of tolerance dimension %i : %.3f | %.3f ", j, distInDim, tolerances[j]);
						//printf("\n mean %.3f |  %.3f ", mean(0, j), data(0, j));
						break;
					}
				}
			}


			double dist = coreUtils.getEuclideanDistance(data, mean);
			if (dist < minDist && allDimsInTolerance)
			{
				minDist = dist;
				out = i;
			}
		}
		//printf("\n clusterID %i ", out);
		return out;
	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::getClusterIndex(MatrixXf& data, MatrixXf& means)
	{
		double minDist = DBL_MAX;
		int out = -1;

		for (int i = 0; i < means.rows(); i++)
		{
			MatrixXf mean = means.row(i);


			//find the distance in each dimension separately (the values should NOT be normalized)
			bool allDimsInTolerance = true;
			if (tolerances.size() == mean.cols())
			{
				for (int j = 0; j < mean.cols(); j++)
				{
					float distInDim = abs(mean(0, j) - data(0, j));
					if (distInDim > tolerances[j])
					{
						allDimsInTolerance = false;
						//printf("\n out of tolerance dimension %i : %.3f | %.3f ", j, distInDim, tolerances[j]);
						//printf("\n mean %.3f |  %.3f ", mean(0, j), data(0, j));
						break;
					}
				}
			}


			double dist = coreUtils.getEuclideanDistance(data, mean);
			if (dist < minDist && allDimsInTolerance)
			{
				minDist = dist;
				out = i;
			}
		}
		//printf("\n clusterID %i ", out);
		return out;
	}

	ZSPACE_TOOLSETS_INLINE void zTsKMeans::updateMean(MatrixXf& data, MatrixXf& mean, int clusterSize)
	{
		for (int i = 0; i < mean.cols(); i++)
		{
			mean(0, i) = (mean(0, i) * (clusterSize - 1) + data(0, i)) / clusterSize;
		}
	}

	ZSPACE_TOOLSETS_INLINE int zTsKMeans::probabilitySelection(vector<float> list, int seed)
	{
		int n = list.size();
		vector<pair<float, int> > pair;
		for (int i = 0; i < n; i++) pair.push_back(make_pair(list[i], i));

		std::sort(pair.rbegin(), pair.rend());

		float sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += pair[i].first;
		}

		float cumulativeProb = 0;
		srand(seed);
		double random = coreUtils.randomNumber_double(0, 1);

		for (int i = 0; i < n; i++)
		{
			cumulativeProb += (pair[i].first / sum);
			if (cumulativeProb > random) return pair[i].second;
		}

		return pair[n - 1].second;
	}

	//ZSPACE_TOOLSETS_INLINE void zTsKMeans::createEvaluationGraph(vector<float, float> pairs, string xTitle, string yTitle)
	//{
	//	//draw axis 
	//}

	ZSPACE_TOOLSETS_INLINE int zTsKMeans::findOptimalK_Elbow(initialisationMethod initMethod, bool distortionMethod, int min, int max, int increment, int seed1, int seed2, vector<pair<int, float>>& KScorePair, float tolerance)
	{
		KScorePair.clear();
		//vector<int> ks;
		vector<float> slopes;
		//vector<float> scores;
		vector<float> slopeDifference;
		int counter = 0;
		int optimalCount = min;
		float deviation = 0;
		for (int k = min; k <= max; k += increment)
		{
			int numCluster = k;
			setNumClusters(k);
			getKMeansClusters(numCluster, initMethod, seed1, seed2, tolerance);
			float score = distortionMethod ? calculateDistortion() : calculateIntertia();

			KScorePair.push_back(make_pair(k, score));

			if (k - increment >= min && k + increment <= max)
			{
				int prvs = counter - 1; int next = counter + 1;
				float s0 = KScorePair[prvs].second - KScorePair[counter].second;
				float s2 = KScorePair[counter].second - KScorePair[next].second;

				slopeDifference.push_back(abs(s2 - s0));
			}
			counter++;

		}

		for (int i = 1; i < slopeDifference.size(); i++)
		{
			if (slopeDifference[i] > deviation)
			{
				deviation = slopeDifference[i];
				optimalCount = min + i;
			}
		}

		return optimalCount;
	}
	ZSPACE_TOOLSETS_INLINE int zTsKMeans::findOptimalK_Silhouette(initialisationMethod initMethod, int min, int max, int increment, int seed1, int seed2, vector<pair<int, float>>& KScorePair, float tolerance)
	{
		KScorePair.clear();
		float* mesh;
		//vector<int> ks;
		vector<float> slopes;
		//vector<float> scores;
		vector<float> slopeDifference;
		int counter = 0;
		int optimalCount = min;
		float deviation = 0;

		vector<pair<float, int>> scoreKpair;
		//vector<pair<int, float>> KScorePair;

		for (int k = min; k <= max; k += increment)
		{
			printf("\n k %i", k);

			int numCluster = k;
			setNumClusters(k);

			getKMeansClusters(numCluster, initMethod, seed1, seed2, tolerance);

			float score = calculateSilhouette();

			scoreKpair.push_back(make_pair(score, k));

			KScorePair.push_back(make_pair(k, score));


		}

		std::sort(scoreKpair.begin(), scoreKpair.end());

		optimalCount = scoreKpair[scoreKpair.size() - 1].second;
		/*printf("\n optimalCount = %i ", optimalCount);

		for (int i = 0; i < KScorePair.size(); i++)
		{
			printf("\n KScorePair[%i].first	= %i ", i, KScorePair[i].first);
			printf("\n KScorePair[%i].second	= %f", i, KScorePair[i].second);
		}
		for (int i = 0; i < scoreKpair.size(); i++)
		{
			printf("\n scoreKpair[%i].first	= %f ", i, scoreKpair[i].first);
			printf("\n scoreKpair[%i].second	= %i", i, scoreKpair[i].second);
		}*/


		return optimalCount;
	}

	ZSPACE_TOOLSETS_INLINE vector<float> zTsKMeans::matrixToVector(MatrixXf m, int num, bool row)
	{
		vector<float> v;
		if (row)
		{
			if (m.rows() > num)
			{

				for (int i = 0; i < m.cols(); i++)
				{
					v.push_back(m(num, i));
				}
			}
		}
		else
		{
			if (m.cols() <= num)
			{
				for (int i = 0; i < m.rows(); i++)
				{
					v.push_back(m(i, num));
				}
			}
		}
		return v;
	}



	//int findOptimalClusterCount(MatrixXf& _dataPoints, int& _numIterations, bool plusPlusMethod, bool distortionMethod, int min, int max, int seed1, int seed2)
	//{
	//	vector<int> ks;
	//	vector<float> slopes;
	//	vector<float> scores;
	//	vector<float> slopeDifference;
	//	int counter = 0;
	//	int optimalCount = min;
	//	float deviation = 0;
	//	for (int k = min; k <= max; k++)
	//	{
	//		zTsKMeans km;
	//		km = zTsKMeans(_dataPoints, k, _numIterations);
	//		km.getKMeansClusters(k, plusPlusMethod, seed1, seed2);
	//		
	//		float score = distortionMethod ? km.ge : km.inertia;
	//		scores.push_back(score);

	//		if (k-1 >= min && k+1 <= max)
	//		{
	//			float s0 = scores[counter-1] - scores[counter];
	//			float s2 = scores[counter] - scores[counter + 1];

	//			slopeDifference.push_back(abs(s2 - s0));
	//		}
	//		counter++;

	//	}
	//	for (int i = 1; i < slopeDifference.size(); i++)
	//	{
	//		if (slopeDifference[i] > deviation)
	//		{
	//			deviation = slopeDifference[i];
	//			optimalCount = min + i;
	//		}
	//	}


	//	//find the k where the slope decrease rabidly (the biggest difference in slope)
	//	
	//	//for (int i = 0; i < scores.size()-1; i++)
	//	//{
	//	//	 slopes.push_back(scores[i] - scores[i + 1]);

	//	//	
	//	//}
	//	//for (int i = 1; i < slopes.size(); i++)
	//	//{
	//	//	float diff = abs(slopes[i-1] - slopes[i]);
	//	//	if (diff > deviation)
	//	//	{
	//	//		deviation = diff;
	//	//		optimalCount = min + i;
	//	//	}
	//	//}

	//	return optimalCount;
	//}


}