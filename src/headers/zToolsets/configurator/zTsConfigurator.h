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

#ifndef ZSPACE_TS_CONFIGURATOR_H
#define ZSPACE_TS_CONFIGURATOR_H

#pragma once
#include <zInterface/functionsets/zFnMesh.h>
#include <zApp/include/zFnSets.h>
#include <zToolsets/geometry/zTsColorSplit.h>


namespace zSpace
{
	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsConfigurator
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsConfigurator
	*	\brief Configurator.
	*	\method based on Breadth-first search (BFS) method
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/


	//enum zProgramme { office = 0, resi = 1, pub = 2 };
	enum zType { END = 0, CORNER = 1, LINE = 2, TRI = 3, CROSS = 4};

	//struct zColorHash {
	//	size_t operator()(const zColor& color) const {
	//		return std::hash<int>()(color.r) ^ std::hash<int>()(color.g) ^ std::hash<int>()(color.b);
	//	}
	//};

	class ZSPACE_TOOLSETS zGameObj
	{
	public:
		int id;
		string programme;
		zType type;
		zTransform transform;
		zObjMesh oMesh;

		// Constructor
		//zGameObj();

		//zGameObj(int id, string programme, zType type, zTransform transform)
		//	: id(id), programme(programme), type(type), transform(transform) {}
	};

	class ZSPACE_TOOLSETS zTsConfigurator
	{
	protected:		

		string mainDir;

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;

		/*!	\brief pointer to input Object  */
		zObjMesh baseMesh;

		/*!	\brief container to dual graph of input Object  */
		zObjGraph baseGraph;

		/*!	\brief container to result face indices under each mesh obj  */
		vector<zGameObj> gameObjs;

		/*!	\brief numver of result mesh splits */
		int numGameObjs;

		vector<zColor> fCols;

		unordered_map <int, string> map_id_programme;

		zTsColorSplit spliter;

	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------

		/*!	\brief input mesh function set  */
		zFnMesh fnMesh;

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.2
		*/
		zTsConfigurator();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsConfigurator();

		//--------------------------
		//---- SET METHODS
		//--------------------------

				/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_meshObj			- input mesh object.
		*	\since version 0.0.2
		*/
		void setDirectory(string& _path);

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_meshObj			- input mesh object.
		*	\since version 0.0.2
		*/
		void setBaseMeshFromFile(string& _path);

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the splits of an input mesh based on the face colors.
		*
		*	\since version 0.0.2
		*/
		void initialise();


		/*! \brief This method computes the splits of an input mesh based on the face colors.
		*
		*	\since version 0.0.2
		*/
		void compute();

		void loadConfig();

		void loadMesh();

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_meshObj			- input mesh object.
		*	\since version 0.0.2
		*/
		zObjMesh* getRawBaseMesh();

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_meshObj			- input mesh object.
		*	\since version 0.0.2
		*/
		zObjGraph* getRawBaseGraph();

		vector<zGameObj> getGameObjs();
		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method returns the split meshes.
		*
		*	\since version 0.0.2
		*/
		zObjMeshPointerArray getRawSplitMesh(int& _numMesh);

		//--------------------------
		//---- DRAW METHODS
		//--------------------------
		void draw();

	private:

		//--------------------------
		//---- PRIVATE UTILITY METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void computeBaseGraph();

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void computeGameObjs();

		string colorToProgramme(zColor& _col);

		void loadColorProgrammeMap();


		void checkType(zItGraphVertex& _v, zType& _type, zVector& _alignVector);

		/*! \brief This method exports color split meshes to the input directory.
		*
		*	\param		[in]	_pth	- export directory.
		*	\param		[in]	_type	- export file type.
		*	\since version 0.0.2
		*/
		void exportTo(string& _pth, zFileType _type = zJSON);

		void to_json(json& j, const zTransform& t);

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		bool isRecorded(int _id, zIntArray& _recorder);

	};
}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/configurator/zTsConfigurator.cpp>
#endif

#endif