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


#include "headers/zToolsets/geometry/zTsSpatialGraph.h"

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSpatialGraph::zTsSpatialGraph()
	{
		
	
	}


	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSpatialGraph::~zTsSpatialGraph() {}

	//---- CREATE METHODS

	//--- SET METHODS 

	ZSPACE_INLINE void zTsSpatialGraph::setFromFile(string path, zFileTpye fileType)
	{
		zFnGraph fnGraph(o_SpatialGraph);
		fnGraph.from(path, fileType);
	}

	//---- GET METHODS

	ZSPACE_INLINE zObjGraph* zTsSpatialGraph::getRawSpatialGraph()
	{
		return &o_SpatialGraph;
	}

	ZSPACE_INLINE void zTsSpatialGraph::computeConvexHulls(float edgeFactor0, float edgeFactor1, int profileVerts, float profileEdgeLen)
	{
		zFnGraph fnGraph(o_SpatialGraph);

		o_ConvexHulls.clear();
		o_ConvexHulls.assign(fnGraph.numVertices(), zObjMesh());

		zPointArray profile = computeProfilePoints(profileVerts, profileEdgeLen);

		for (zItGraphVertex v(o_SpatialGraph); !v.end(); v++)
		{
			zItGraphHalfEdgeArray cHEdges;
			v.getConnectedHalfEdges(cHEdges);

			for (auto& he : cHEdges)
			{
				zVector heVec = he.getVector();
				

			}

		}

	}




	//---- COMPUTE METHODS

	
	//---- PROTECTED UTILITY METHODS

	ZSPACE_INLINE zPointArray zTsSpatialGraph::computeProfilePoints(int numVerts, float edgeLen)
	{
		zPointArray out;
		
		double halfElen = edgeLen / 2;
		double theta = 0;

		for (int i = 0; i < numVerts; i++)
		{
			zPoint pos;
			pos.x = (halfElen * cos(theta + Z_HALF_PI));
			pos.y = (halfElen * sin(theta + Z_HALF_PI));
			pos.z = 0;

			out.push_back(pos);

			theta += (Z_TWO_PI / numVerts);
		}

		return out;
	}


}