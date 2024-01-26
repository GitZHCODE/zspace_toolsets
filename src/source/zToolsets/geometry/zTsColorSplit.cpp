// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Taizhong Chen <taizhong.chen@zaha-hadid.com>
//


#include <zToolsets/geometry/zTsColorSplit.h>

namespace zSpace
{

	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsColorSplit::zTsColorSplit() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsColorSplit::~zTsColorSplit() {}

	//---- SET METHODS
	ZSPACE_INLINE void zTsColorSplit::setInMesh(zObjMesh& _meshObj)
	{
		oMesh = &_meshObj;
		fnMesh = zFnMesh(_meshObj);
	}

	//---- GET METHODS

	ZSPACE_INLINE zObjMeshPointerArray zTsColorSplit::getRawSplitMesh(int& _numMesh)
	{
		zObjMeshPointerArray out;
		numMesh = 0;

		numMesh = splits.size();
		_numMesh = numMesh;

		if (numMesh == 0)return out;

		for (auto& m : splits)
		{
			out.push_back(&m);
		}

		return out;
	}

	//---- EXPORT METHODS

	ZSPACE_INLINE void zTsColorSplit::exportTo(string& _pth, zFileTpye _type)
	{
		for (int i = 0;i<splits.size();i++)
		{
			zFnMesh fnMesh(splits[i]);
			string path = _pth;
			path += "m_";
			path += to_string(i);
			if (_type == zOBJ) path += ".obj";
			if (_type == zJSON) path += ".json";
			fnMesh.to(path, _type);
		}
	}

	//---- COMPUTE METHODS

	ZSPACE_INLINE void zTsColorSplit::compute()
	{
		fnMesh.getFaceColors(fCols);
		zIntArray inToDual;
		zIntArray dualToIn;

		fnMesh.getDualGraph(dualGraph, inToDual, dualToIn, true, false, false);
		const int numF = fnMesh.numPolygons();

		//initialise data
		fIDs.assign(numF, -1);

		//check splits and add face id to the list
		bool tested = false;
		int vNow = 0;
		int meshCount = 0;
		do
		{
			tested = checkExisted(vNow, fIDs);
			if (!tested)
			{
				zIntArray recorder;
				int meshID = meshCount;

				zItGraphVertex v(dualGraph, vNow);
				zIntArray connectedV;

				bool isDead;
				int id;
				do
				{
					id = v.getId();
					recorder.push_back(id);

					isDead = true;
					v.getConnectedVertices(connectedV);

					for (int i = 0; i < connectedV.size(); i++) //check all neighours
					{
						int cv = connectedV[i];
						if (fCols[cv] == fCols[id] && !isRecorded(cv, recorder)) //find neighbour has same color and not walked on
						{
							if (fIDs[cv] != -1) //if that neighbour belongs to a recorded mesh
							{
								meshID = fIDs[cv];
							}
							else //same color, not walked on, not belongs to a recorded mesh, go to that neighbour
							{
								isDead = false;
								v = zItGraphVertex(dualGraph, cv);
							}
						}
					}
				} while (isDead != true);
				for (auto& o : recorder) fIDs[o] = meshID; //push temp recorder to the mesh id list
				if (meshCount == meshID)meshCount++;
			}
			else vNow++;
		} while (vNow != numF);

		resultMeshID = meshContainer(fIDs); //convert mesh id list to a 2d container: meshID + faceID

		generateMesh();
	}

	//---- PRIVATE METHODS

	ZSPACE_INLINE bool zTsColorSplit::isRecorded(int _id, zIntArray& _recorder)
	{
		bool isRecorded = false;
		for (auto &id : _recorder)
			if (_id == id)
				isRecorded = true;
		goto end;
	end:
		return  isRecorded;
	}

	ZSPACE_INLINE bool zTsColorSplit::checkExisted(int _id, zIntArray& _faceID)
	{
		bool found = false;
		for (int i = 0; i < _faceID.size(); i++)
			if (_faceID[_id] != -1)
			{
				found = true;
				goto end;
			}

	end:
		return found;
	}

	ZSPACE_INLINE vector<vector<int>> zTsColorSplit::meshContainer(zIntArray faceID)
	{
		set<int> set(begin(faceID), end(faceID));

		vector<vector<int>> result;
		for (auto it = set.begin(); it != set.end(); it++)
		{
			vector<int> faces;
			for (int j = 0; j < faceID.size(); j++)
				if (faceID[j] == *it)
					faces.push_back(j);

			result.push_back(faces);
		}
		return result;
	}

	ZSPACE_INLINE void zTsColorSplit::generateMesh()
	{
		//make split mesh and export
		splits.clear();
		int counter = 0;
		zPointArray vertices;
		fnMesh.getVertexPositions(vertices, false);

		splits.assign(resultMeshID.size(), zObjMesh());
		numMesh = resultMeshID.size();

		for (auto meshID = resultMeshID.begin(); meshID != resultMeshID.end(); meshID++)
		{
			zFnMesh fm(splits[counter]);

			zPointArray pVertices;
			zIntArray pCounts;
			zIntArray pConnects;
			zColor color;


			vector<int> meshFaces = *meshID;

			for (auto faceID = meshFaces.begin(); faceID != meshFaces.end(); faceID++)
			{
				zItMeshFace f(*oMesh, (int)*faceID);
				zPointArray fVertices;
				f.getVertexPositions(fVertices);
				color = f.getColor();

				for (auto& v : fVertices)
				{
					int vID;
					bool check = coreUtils.checkRepeatVector(v, pVertices, vID);
					
					if (!check)
					{
						vID = pVertices.size();
						pVertices.push_back(v);
					}
					pConnects.push_back(vID);
				}
				pCounts.push_back(f.getNumVertices());
			}

			fm.create(pVertices,pCounts,pConnects);
			fm.setVertexColor(color, true);

			cout << endl;
			cout << endl;
			cout << "---"<<"MESH_" << counter << "---" << endl;
			cout << "V_" << fm.numVertices() << ",";
			cout << "E_" << fm.numEdges() << ",";
			cout << "F_" << fm.numPolygons() << ",";

			//cout << endl;
			//cout << "pVertices" << endl;
			//for (auto& o : pVertices) cout << o << ",";
			//cout << endl;
			//cout << "pCounts" << endl;
			//for (auto& o : pCounts) cout << o << ",";
			//cout << endl;
			//cout << "pConnects" << endl;
			//for (auto& o : pConnects) cout << o << ",";

			counter++;
		}
	}

	//old method
	/*ZSPACE_INLINE zIntArray zTsColorSplit::transformID(zIntArray & v)
	{
		map<int, size_t> m;
		for (auto e : v) m[e];

		int index = 0;
		for (auto& [key, value] : m) value = index++;

		zIntArray res;
		for (auto e : v) res.push_back(m[e]);
			
		return res;
	}*/

}