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


#include "headers/zToolsets/externalMethods/zExtMeshUtility.h"


namespace zSpace
{
	
	ZSPACE_TOOLSETS_INLINE void zExtMesh::updateFields()
	{
		zFnMesh fn(*mesh);
		vCount = fn.numVertices();
		fCount = fn.numPolygons();
	}
	


	ZSPACE_TOOLSETS_INLINE void ext_meshUtil_createMeshOBJ(double* _vertexPositions, int* _polyCounts, int* _polyConnects, int numVerts, int numFaces, zExtMesh& out_mesh)
	{
		if (!_vertexPositions || !_polyCounts || !_polyConnects) throw std::invalid_argument(" error: mesh container is empty.");

		zPointArray vPos;
		zIntArray pConnects;
		zIntArray pCounts;

		for (int i = 0; i < numVerts; i++)
		{
			zVector v;
			v = zVector(_vertexPositions[i * 3 + 0], _vertexPositions[i * 3 + 1], _vertexPositions[i * 3 + 2]);
			vPos.push_back(v);
		}
		int polyconnectsCurrentIndex = 0;
		for (int i = 0; i < numFaces; i++)
		{
			int num_faceVerts = _polyCounts[i];
			pCounts.push_back(_polyCounts[i]);

			for (int j = 0; j < num_faceVerts; j++)
			{
				pConnects.push_back(_polyConnects[polyconnectsCurrentIndex + j]);
			}

			polyconnectsCurrentIndex += num_faceVerts;
		}
		zFnMesh fnMesh(*out_mesh.mesh);
		fnMesh.create(vPos, pCounts, pConnects);
	}
	ZSPACE_TOOLSETS_INLINE int ext_meshUtil_getMeshPosition(zExtMesh objMesh, float* outVPostions, float* outVColors)
	{
		if (!objMesh.mesh)
		{
			"/n meshPointer is null";
			return 0;
		}
		zFnMesh fn(*objMesh.mesh);
		zPoint* pts = fn.getRawVertexPositions();
		zColor* colors = fn.getRawVertexColors();
		for (int i = 0; i < fn.numVertices(); i++)
		{
			outVPostions[i * 3 + 0] = pts[i].x;
			outVPostions[i * 3 + 1] = pts[i].y;
			outVPostions[i * 3 + 2] = pts[i].z;

			outVColors[i * 4 + 0] = colors[i].r;
			outVColors[i * 4 + 1] = colors[i].g;
			outVColors[i * 4 + 2] = colors[i].b;
			outVColors[i * 4 + 3] = colors[i].a;
		}
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_meshUtil_getMeshFaceCount(zExtMesh objMesh, int* outfCounts)
	{
		if (!objMesh.mesh)
		{
			"/n meshPointer is null";
			return 0;
		}
		zFnMesh fn(*objMesh.mesh);
		zIntArray pCounts;
		zIntArray pConnects;
		fn.getPolygonData(pConnects, pCounts);
		for (int i = 0; i < pCounts.size(); i++)
		{
			outfCounts[i] = pCounts[i];
		}
		return 1;
	}
	ZSPACE_TOOLSETS_INLINE int ext_meshUtil_getMeshFaceConnect(zExtMesh objMesh, int* outfConnects)
	{
		if (!objMesh.mesh)
		{
			"/n meshPointer is null";
			return 0;
		}
		zFnMesh fn(*objMesh.mesh);
		zIntArray pCounts;
		zIntArray pConnects;
		fn.getPolygonData(pConnects, pCounts);
		for (int i = 0; i < pConnects.size(); i++)
		{
			outfConnects[i] = pConnects[i];
		}
		return 1;
	}


}