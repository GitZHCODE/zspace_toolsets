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


#include<headers/zToolsets/0_externalMethods/zExtGraph.h>


namespace zSpace
{
	ZSPACE_TOOLSETS_INLINE void zExtGraph::updateFields()
	{
		zFnGraph g(*graph);
		eCount = g.numEdges();
		vCount = g.numVertices();
	}
	ZSPACE_TOOLSETS_INLINE void zExtGraphSet::updateFields()
	{
		graphsCount = graphSet->size();
	}

	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraph(zExtGraph extGraph, float* vPositions, float* vColors, int* ePairs, float* eColors)
	{
		zObjGraph* graph;
		graph = extGraph.graph;
		zFnGraph g(*graph);
		int vCount = g.numVertices();
		int eCount = g.numEdges();


		zPointArray inVerticies;
		g.getVertexPositions(inVerticies);

		zColorArray inVColors;
		g.getVertexColors(inVColors);

		zColorArray ineColors;
		g.getEdgeColors(ineColors);

		zIntArray inEdges;
		g.getEdgeData(inEdges);
		for (int i = 0; i < vCount; i++)
		{
			vPositions[i * 3 + 0] = inVerticies[i].x;
			vPositions[i * 3 + 1] = inVerticies[i].y;
			vPositions[i * 3 + 2] = inVerticies[i].z;

			vColors[i * 4 + 0] = inVColors[i].r;
			vColors[i * 4 + 1] = inVColors[i].g;
			vColors[i * 4 + 2] = inVColors[i].b;
			vColors[i * 4 + 3] = inVColors[i].a;
		}
		for (int64_t i = 0; i < eCount; i++)
		{

			ePairs[i * 2 + 0] = inEdges[i * 2 + 0];
			ePairs[i * 2 + 1] = inEdges[i * 2 + 1];

			eColors[i * 4 + 0] = ineColors[i].r;
			eColors[i * 4 + 1] = ineColors[i].g;
			eColors[i * 4 + 2] = ineColors[i].b;
			eColors[i * 4 + 3] = ineColors[i].a;

		}
	}
	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSet(zExtGraphSet graphSet, zExtGraph* outGraphArray)
	{

		for (int i = 0; i < graphSet.graphsCount; i++)
		{
			outGraphArray[i].graph = graphSet.graphSet->at(i);
			outGraphArray[i].updateFields();
		}
		/*zObjGraphArray* graphs;
		graphs = graphSet.graphSet;
		for (int i = 0; i < graphs->size(); i++)
		{
			outGraphArray[i] = graphSet.graphSet->at(i);
		}*/
	}






	ZSPACE_TOOLSETS_INLINE void ext_sdf_getPlanesData(vector<zTransform>* graph, float* matrix)
	{
		for (int i = 0; i < graph->size(); i++)
		{
			//outPlanes[i * 4 + 0] = 
			zTransform frame = graph->at(i);

			zTransform frameTranspose = frame.transpose();

			float* m = frameTranspose.data();
			
			//float* m = graph->at(i).transpose().data(); //doesn't work
			for (int j = 0; j < 16; j++)
			{
				matrix[i * 16 + j] = m[j];
			}
		}
	}

	//Graph Data
	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSetFromPointersVector2(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray)
	{
		for (int i = 0; i < graphs->size(); i++)
		{
			outGraphArray[i] = graphs->at(i);
		}
	}
	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSetFromVector2(zObjGraphArray* graphs, zObjGraph** outGraphArray)
	{
		for (int i = 0; i < graphs->size(); i++)
		{
			outGraphArray[i] = &graphs->at(i);
		}
	}
	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphCounts2(zObjGraph* graph, int& outvCount, int& outeCount)
	{
		zFnGraph g(*graph);
		outeCount = g.numEdges();
		outvCount = g.numVertices();
	}
	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphData2(zObjGraph* graph, float* vPositions, float* vColors, int* ePairs, float* eColors)
	{
		zFnGraph g(*graph);
		int vCount = g.numVertices();
		int eCount = g.numEdges();


		zPointArray inVerticies;
		g.getVertexPositions(inVerticies);

		zColorArray inVColors;
		g.getVertexColors(inVColors);

		zColorArray ineColors;
		g.getEdgeColors(ineColors);

		zIntArray inEdges;
		g.getEdgeData(inEdges);
		for (int i = 0; i < vCount; i++)
		{
			vPositions[i * 3 + 0] = inVerticies[i].x;
			vPositions[i * 3 + 1] = inVerticies[i].y;
			vPositions[i * 3 + 2] = inVerticies[i].z;

			vColors[i * 4 + 0] = inVColors[i].r;
			vColors[i * 4 + 1] = inVColors[i].g;
			vColors[i * 4 + 2] = inVColors[i].b;
			vColors[i * 4 + 3] = inVColors[i].a;
		}
		for (int64_t i = 0; i < eCount; i++)
		{

			ePairs[i * 2 + 0] = inEdges[i * 2 + 0];
			ePairs[i * 2 + 1] = inEdges[i * 2 + 1];

			eColors[i * 4 + 0] = ineColors[i].r;
			eColors[i * 4 + 1] = ineColors[i].g;
			eColors[i * 4 + 2] = ineColors[i].b;
			eColors[i * 4 + 3] = ineColors[i].a;

		}
	}

}