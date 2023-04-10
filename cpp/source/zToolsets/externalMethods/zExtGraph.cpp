//// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
//// data analysis & visualization framework.
////
//// Copyright (C) 2019 ZSPACE 
//// 
//// This Source Code Form is subject to the terms of the MIT License 
//// If a copy of the MIT License was not distributed with this file, You can 
//// obtain one at https://opensource.org/licenses/MIT.
////
//// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
////
//
//
//#include<headers/zToolsets/externalMethods/zExtGraph.h>
//
//
//namespace zSpace
//{
//	ZSPACE_TOOLSETS_INLINE zExtGraph::zExtGraph(zObjGraph* g)
//	{
//		graph = g;
//		updateAttributes();
//	}
//	ZSPACE_TOOLSETS_INLINE void zExtGraph::updateAttributes()
//	{
//		zFnGraph g(*graph);
//		eCount = g.numEdges();
//		vCount = g.numVertices();
//	}
//	ZSPACE_TOOLSETS_INLINE zIntArray zExtGraph::getGraphSequence()
//	{
//		zIntArray sequence;
//		zItGraphVertexArray vArray;
//
//		zFnGraph fnGraph(*graph);
//		if (fnGraph.numVertices() == 0) return sequence;
//
//		for (zItGraphVertex v(*graph); !v.end(); v++)
//		{
//			if (!v.checkValency(2))
//			{
//				vArray.push_back(v);
//			}
//		}
//
//		printf("\n vArray.size() - %i", vArray.size());
//
//		if (vArray.size() == 2)
//		{
//			zItGraphHalfEdge he = vArray[0].getHalfEdge();
//			sequence.push_back(vArray[0].getId());
//			do
//			{
//				sequence.push_back(he.getVertex().getId());
//				he = he.getNext();
//			} while (he.getVertex() != vArray[1]);
//
//			sequence.push_back(vArray[1].getId());
//			sequence.push_back(vArray[0].getId());
//		}
//		if (vArray.size() == 0)
//		{
//
//			zItGraphHalfEdge he(*graph, 0);
//
//			zItGraphVertex startV = he.getStartVertex();
//
//			zItGraphHalfEdge startHe = he;
//
//			sequence.push_back(he.getStartVertex().getId());
//
//			do
//			{
//
//				if (he.getVertex() == startV)
//				{
//
//					sequence.push_back(he.getVertex().getId());
//				}
//				else
//				{
//
//					sequence.push_back(he.getVertex().getId());
//				}
//
//				he = he.getNext();
//
//			} while (he != startHe);
//		}
//
//		printf("\n num of vertices : num of sequence ---- %i : %i \n ", fnGraph.numVertices(), sequence.size());
//
//		return sequence;
//	}
//
//
//	ZSPACE_TOOLSETS_INLINE void zExtGraphSet::updateAttributes()
//	{
//		graphsCount = graphSet->size();
//	}
//
//
//	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphData(zExtGraph extGraph, float* vPositions, float* vColors, int* ePairs, float* eColors)
//	{
//		zObjGraph* graph;
//		graph = extGraph.graph;
//		zFnGraph g(*graph);
//		int vCount = g.numVertices();
//		int eCount = g.numEdges();
//
//
//		zPointArray inVerticies;
//		g.getVertexPositions(inVerticies);
//
//		zColorArray inVColors;
//		g.getVertexColors(inVColors);
//
//		zColorArray ineColors;
//		g.getEdgeColors(ineColors);
//
//		zIntArray inEdges;
//		g.getEdgeData(inEdges);
//		for (int i = 0; i < vCount; i++)
//		{
//			vPositions[i * 3 + 0] = inVerticies[i].x;
//			vPositions[i * 3 + 1] = inVerticies[i].y;
//			vPositions[i * 3 + 2] = inVerticies[i].z;
//
//			vColors[i * 4 + 0] = inVColors[i].r;
//			vColors[i * 4 + 1] = inVColors[i].g;
//			vColors[i * 4 + 2] = inVColors[i].b;
//			vColors[i * 4 + 3] = inVColors[i].a;
//		}
//		for (int64_t i = 0; i < eCount; i++)
//		{
//
//			ePairs[i * 2 + 0] = inEdges[i * 2 + 0];
//			ePairs[i * 2 + 1] = inEdges[i * 2 + 1];
//
//			eColors[i * 4 + 0] = ineColors[i].r;
//			eColors[i * 4 + 1] = ineColors[i].g;
//			eColors[i * 4 + 2] = ineColors[i].b;
//			eColors[i * 4 + 3] = ineColors[i].a;
//
//		}
//	}
//	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphSequence(zExtGraph extGraph, int* outSequence)//the size of the area = the number of vertices + 1
//	{
//		zIntArray seq;
//		seq = extGraph.getGraphSequence();
//		for (int i = 0; i < seq.size(); i++)
//		{
//			outSequence[i] = seq[i];
//		}
//	}
//	ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSet(zExtGraphSet graphSet, zExtGraph* outGraphArray)
//	{
//
//		for (int i = 0; i < graphSet.graphsCount; i++)
//		{
//			outGraphArray[i].graph = graphSet.graphSet->at(i);
//			outGraphArray[i].updateAttributes();
//		}
//		/*zObjGraphArray* graphs;
//		graphs = graphSet.graphSet;
//		for (int i = 0; i < graphs->size(); i++)
//		{
//			outGraphArray[i] = graphSet.graphSet->at(i);
//		}*/
//	}
//
//
//	////Graph Data
//	//ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSetFromPointersVector2(zObjGraphPointerArray* graphs, zObjGraph** outGraphArray)
//	//{
//	//	for (int i = 0; i < graphs->size(); i++)
//	//	{
//	//		outGraphArray[i] = graphs->at(i);
//	//	}
//	//}
//	//ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphsSetFromVector2(zObjGraphArray* graphs, zObjGraph** outGraphArray)
//	//{
//	//	for (int i = 0; i < graphs->size(); i++)
//	//	{
//	//		outGraphArray[i] = &graphs->at(i);
//	//	}
//	//}
//	//ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphCounts2(zObjGraph* graph, int& outvCount, int& outeCount)
//	//{
//	//	zFnGraph g(*graph);
//	//	outeCount = g.numEdges();
//	//	outvCount = g.numVertices();
//	//}
//	//ZSPACE_TOOLSETS_INLINE void ext_graphUtil_getGraphData2(zObjGraph* graph, float* vPositions, float* vColors, int* ePairs, float* eColors)
//	//{
//	//	zFnGraph g(*graph);
//	//	int vCount = g.numVertices();
//	//	int eCount = g.numEdges();
//
//
//	//	zPointArray inVerticies;
//	//	g.getVertexPositions(inVerticies);
//
//	//	zColorArray inVColors;
//	//	g.getVertexColors(inVColors);
//
//	//	zColorArray ineColors;
//	//	g.getEdgeColors(ineColors);
//
//	//	zIntArray inEdges;
//	//	g.getEdgeData(inEdges);
//	//	for (int i = 0; i < vCount; i++)
//	//	{
//	//		vPositions[i * 3 + 0] = inVerticies[i].x;
//	//		vPositions[i * 3 + 1] = inVerticies[i].y;
//	//		vPositions[i * 3 + 2] = inVerticies[i].z;
//
//	//		vColors[i * 4 + 0] = inVColors[i].r;
//	//		vColors[i * 4 + 1] = inVColors[i].g;
//	//		vColors[i * 4 + 2] = inVColors[i].b;
//	//		vColors[i * 4 + 3] = inVColors[i].a;
//	//	}
//	//	for (int64_t i = 0; i < eCount; i++)
//	//	{
//
//	//		ePairs[i * 2 + 0] = inEdges[i * 2 + 0];
//	//		ePairs[i * 2 + 1] = inEdges[i * 2 + 1];
//
//	//		eColors[i * 4 + 0] = ineColors[i].r;
//	//		eColors[i * 4 + 1] = ineColors[i].g;
//	//		eColors[i * 4 + 2] = ineColors[i].b;
//	//		eColors[i * 4 + 3] = ineColors[i].a;
//
//	//	}
//	//}
//	//
//}