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
//#ifndef ZSPACE_EXT_TS_MESH_UTILITY_H
//#define ZSPACE_EXT_TS_MESH_UTILITY_H
//
//
//
//#pragma once
//#include <headers/base/zSpace_Toolsets.h>
//
//#include <headers/zCore/base/zExtern.h>
//#include <headers/zInterface/functionsets/zFnMesh.h>
//#include <headers/zInterface/functionsets/zFnGraph.h>
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <iostream>
//#include <sstream>
//
//#include<execution>
//
//using namespace std;
//
//
//namespace zSpace
//{
//	struct zExtMesh
//	{
//		zObjMesh* mesh;
//		int vCount;
//		int fCount;
//
//		zExtMesh(zObjMesh* m);
//		void updateAttributes();
//	};
//
//
//	ZSPACE_TOOLSETS_EXT
//	{
//		//from array //from file //from graph
//		ZSPACE_TOOLSETS void ext_meshUtil_createMeshOBJ(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zExtMesh& out_mesh);
//		ZSPACE_TOOLSETS int ext_meshUtil_getMeshPosition(zExtMesh objMesh, float* outVPostions, float* outVColors);
//		ZSPACE_TOOLSETS int ext_meshUtil_getMeshFaceCount(zExtMesh objMesh, int* outfCounts);
//		ZSPACE_TOOLSETS int ext_meshUtil_getMeshFaceConnect(zExtMesh objMesh, int* outfConnects);
//
//
//
//
//		//ZSPACE_TOOLSETS void ext_meshUtil_createMeshOBJ2(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zObjMesh& outMesh);
//
//		////Mesh Data
//		//ZSPACE_TOOLSETS int ext_meshUtil_getMeshCounts2(zObjMesh* objMesh, int& out_vCount, int& out_fCount);
//		//ZSPACE_TOOLSETS int ext_meshUtil_getMeshPosition2(zObjMesh* objMesh, float* outVPostions, float* outVColors);
//		//ZSPACE_TOOLSETS int ext_meshUtil_getMeshFaceCount2(zObjMesh* objMesh, int* outfCounts);
//		//ZSPACE_TOOLSETS int ext_meshUtil_getMeshFaceConnect2(zObjMesh* objMesh, int* outfConnects);
//
//
//		//Draw Mesh
//		//ZSPACE_TOOLSETS void ext_sdf_DrawMesh2(zObjMesh* objMesh);
//
//
//		
//	}
//
//}
//
//
//
//
//#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
//// All defined OK so do nothing
//#else
//#include<source/zToolsets/externalMethods/zExtMesh.cpp>
//#endif
//
//#endif