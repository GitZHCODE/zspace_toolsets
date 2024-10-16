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


#include "zToolsets/configurator/zTsConfigurator.h"

namespace zSpace
{

	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsConfigurator::zTsConfigurator() {}

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsConfigurator::~zTsConfigurator() {}

	//---- SET METHODS

	ZSPACE_INLINE void zTsConfigurator::setDirectory(string& _path)
	{
		mainDir = _path;
	}

	ZSPACE_INLINE void zTsConfigurator::setBaseMeshFromFile(string& _path)
	{
		zFnMesh fn(baseMesh);
		fn.from(_path, zJSON, true);
	}

	//---- GET METHODS
	ZSPACE_INLINE zObjMesh* zTsConfigurator::getRawBaseMesh()
	{
		return &baseMesh;
	}

	ZSPACE_INLINE zObjGraph* zTsConfigurator::getRawBaseGraph()
	{
		return &baseGraph;
	}

	ZSPACE_INLINE vector<zGameObj> zTsConfigurator::getGameObjs()
	{
		return gameObjs;
	}

	//---- EXPORT METHODS


	//---- COMPUTE METHODS
	ZSPACE_INLINE void zTsConfigurator::initialise()
	{
		map_id_programme.clear();
		
		loadColorProgrammeMap();
		computeBaseGraph();
	}


	ZSPACE_INLINE void zTsConfigurator::compute()
	{
		computeGameObjs();

		string file = mainDir + "/config.json";
		exportTo(file);
	}

	ZSPACE_INLINE void zTsConfigurator::loadConfig()
	{
		string file_config = mainDir + "/config.json";

		json j;
		bool chk = coreUtils.json_read(file_config, j);

		if (chk)
		{
			json allObjectsJson = j;

			gameObjs.clear();

			for (const auto& item : allObjectsJson) 
			{
				zGameObj gameObject;
				gameObject.id = item["id"].get<int>();
				gameObject.programme = item["programme"].get<string>();
				gameObject.type = item["type"].get<zType>();

				vector<float> tData = item["transform"].get<vector<float>>();
				zTransform t;
				for (int i = 0; i < 16; ++i)
					t(i % 4, i / 4) = tData[i];

				gameObject.transform = t;

				gameObjs.push_back(gameObject);
			}
		}
	}

	ZSPACE_INLINE void zTsConfigurator::loadMesh()
	{
		for (auto& obj : gameObjs)
		{
			string file_mesh = mainDir + "/Assets/";
			file_mesh = file_mesh + obj.programme + "/m_" + to_string(obj.type) + ".json";

			zFnMesh fn(obj.oMesh);
			fn.from(file_mesh, zJSON);
			fn.setTransform(obj.transform);
		}
	}



	//---- PRIVATE METHODS

	ZSPACE_INLINE void zTsConfigurator::loadColorProgrammeMap()
	{
		json j;
		string file = "/color_programme.json";
		bool chk = coreUtils.json_read(mainDir + file, j);

		if (chk)
		{
			const json& colors = j["Colors"].get<vector<vector<float>>>();
			const json& programmes = j["Programmes"].get<vector<string>>();

			fCols.assign(colors.size(),zColor());

			for (int i = 0; i < colors.size(); i++)
			{
				size_t st;
				fCols[i] = zColor((colors[i][0]), colors[i][1], colors[i][2], 1);
				map_id_programme[i] = programmes[i];
			}
		}
	}

	ZSPACE_INLINE void zTsConfigurator::computeBaseGraph()
	{
		int numSplitMeshes;
		spliter.setInMesh(baseMesh);
		spliter.compute();
		zObjMeshPointerArray splits = spliter.getRawSplitMesh(numSplitMeshes);

		cout << "numSplitMeshes:" << numSplitMeshes << endl;
		zPointArray positions;
		zIntArray edgeConnects;
		zColorArray cols;
		int v_offset = 0;

		for (int i = 0; i < numSplitMeshes; i++)
		{
			zObjGraph tempGraph;
			zIntArray garbageA, garbageB;
			zFnMesh fnMesh(*splits[i]);
			fnMesh.getDualGraph(tempGraph, garbageA, garbageB, true, false, false);

			zFnGraph fnTempGraph(tempGraph);

			zItGraphVertex v(tempGraph);
			for (; !v.end(); v++)
			{
				zItMeshFace f(*splits[i], v.getId());
				v.setColor(f.getColor());
			}

			//add to base graph
			zPointArray positions_temp;
			zIntArray edgeConnects_temp;
			zColorArray cols_temp;
			fnTempGraph.getVertexPositions(positions_temp);
			fnTempGraph.getEdgeData(edgeConnects_temp);
			fnTempGraph.getVertexColors(cols_temp);

			for (auto& id : edgeConnects_temp)
				id += v_offset;

			positions.insert(positions.end(), positions_temp.begin(), positions_temp.end());
			edgeConnects.insert(edgeConnects.end(), edgeConnects_temp.begin(), edgeConnects_temp.end());
			cols.insert(cols.end(), cols_temp.begin(), cols_temp.end());

			v_offset += positions_temp.size();
		}

		zFnGraph fnGraph(baseGraph);
		fnGraph.create(positions, edgeConnects);
		fnGraph.setVertexColors(cols);

		numGameObjs = fnGraph.numVertices();

		//zIntArray garbageA, garbageB;
		//zFnMesh fnMesh(baseMesh);
		//fnMesh.getDualGraph(baseGraph, garbageA, garbageB, true, false, false);

		//zFnGraph fnGraph(baseGraph);
		//numGameObjs = fnGraph.numVertices();

		//zItGraphVertex v(baseGraph);
		//for (; !v.end(); v++)
		//{
		//	zItMeshFace f(baseMesh, v.getId());
		//	v.setColor(f.getColor());
		//}
	}

	ZSPACE_INLINE void zTsConfigurator::computeGameObjs()
	{
		//assign empty game objs
		gameObjs.assign(numGameObjs, zGameObj());


		for (int i = 0; i < numGameObjs; i++)
		{
			zGameObj* obj = &gameObjs[i];
			zItGraphVertex v(baseGraph, i);
			zColor col = v.getColor();
			zPoint pos = v.getPosition();

			//compute cell type
			zType objType;
			zVector alignVector;
			zVector axisX(1, 0, 0);
			zVector axisZ(0, 0, 1);
			checkType(v, objType, alignVector);

			float angle = axisZ.dihedralAngle(axisX, alignVector);
			angle *= DEG_TO_RAD;
			angle *= -1;

			//set transform
			zTransform tranform;
			tranform.setIdentity();
			tranform(0, 0) = cos(angle);	tranform(0, 1) = -sin(angle);	tranform(0, 2) = 0;
			tranform(1, 0) = sin(angle);	tranform(1, 1) = cos(angle);	tranform(1, 2) = 0;
			tranform(2, 0) = 0;				tranform(2, 1) = 0;				tranform(2, 2) = 1;
			tranform(3, 0) = pos.x;			tranform(3, 1) = pos.y;			tranform(3, 2) = pos.z;

			//assign to obj
			obj->id = i;
			obj->programme = colorToProgramme(col);
			obj->type = objType;
			obj->transform = tranform;

			cout << endl;
			//cout << "id:" << obj->id << endl;
			//cout << "color:" << col.r << "," << col.g << "," << col.b << endl;
			////cout << "colorToProgramme:" << colorToProgramme(col) << endl;
			//cout << "programme:" << obj->programme << endl;
			//cout << "type:" << obj->type << endl;
			//cout << "transform:" << obj->transform << endl;

		}
	}

	ZSPACE_INLINE string zTsConfigurator::colorToProgramme(zColor& _col)
	{


		for (int i = 0; i < fCols.size(); i++)
		{
			if (_col == fCols[i])
				return map_id_programme[i];
		}
			
		return "NOT FOUND";
			/*auto it = map_id_programme.find(id);
		if (it != map_id_programme.end())
		{
			return it->second;
		}*/
	}

	ZSPACE_INLINE void zTsConfigurator::checkType(zItGraphVertex& _v, zType& _type, zVector& _alignVector)
	{
		zType type = LINE;

		zVector alignVector;
		zVector axisX(1, 0, 0);
		zVector axisZ(0, 0, 1);

		zItGraphHalfEdgeArray cHes;
		_v.getConnectedHalfEdges(cHes);

		if (_v.checkValency(1))
		{
			type = END;
			alignVector = cHes[0].getVector();
		}
		else if (_v.checkValency(2))
		{
			zVector vec_a = cHes[0].getVector();
			zVector vec_b = cHes[1].getVector();

			if (vec_a.angle(vec_b) > 90)
			{
				type = LINE;
				alignVector = vec_a;
			}
			else
			{
				type = CORNER;
				zVector temp = vec_a ^ vec_b;
				alignVector = (temp.angle(axisZ) < 1) ? vec_a : vec_b;
			}
		}
		else if (_v.checkValency(3))
		{
			type = TRI;

			zVector vec_a = cHes[0].getVector();
			zVector vec_b = cHes[1].getVector();
			zVector vec_c = cHes[2].getVector();

			float tol = 0.01f;
			if (abs((vec_a*vec_b) + (vec_a * vec_c)) < tol)
				alignVector = vec_a;
			else if (abs((vec_b * vec_a) + (vec_b * vec_c)) < tol)
				alignVector = vec_b;
			else if (abs((vec_c * vec_a) + (vec_c * vec_b)) < tol)
				alignVector = vec_c;

			//if (vec_a.angle(vec_b) <= 91 && vec_a.angle(vec_c) <= 91)
			//	alignVector = vec_a;
			//else if (vec_b.angle(vec_c) <= 91)
			//	alignVector = vec_b;
			//else
			//	alignVector = vec_c;
		}
		else if (_v.checkValency(4))
		{
			type = CROSS;
			alignVector = cHes[0].getVector();
		}

		_alignVector = alignVector;
		_type = type;
	}

	ZSPACE_INLINE void zTsConfigurator::exportTo(string& _pth, zFileTpye _type)
	{
		json allObjectsJson = json::array();
		for (auto& obj : gameObjs)
		{
			json j;
			j["id"] = obj.id;
			j["type"] = obj.type;
			j["programme"] = obj.programme;
			to_json(j["transform"], obj.transform);

			allObjectsJson.push_back(j);
		}

		ofstream file(_pth, ios::out);
		if (file.is_open())
		{
			file << allObjectsJson.dump();
			file.close();
		}
	}

	ZSPACE_INLINE void zTsConfigurator::to_json(json& j, const zTransform& t)
	{
		for (int i = 0; i < t.size(); i++)
			j.push_back(t.data()[i]);
	}

	//---- DRAW METHODS

	ZSPACE_INLINE void zTsConfigurator::draw()
	{

	}


}