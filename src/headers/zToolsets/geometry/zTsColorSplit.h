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

#ifndef ZSPACE_TS_GEOMETRY_COLOR_SPLIT_H
#define ZSPACE_TS_GEOMETRY_COLOR_SPLIT_H

#pragma once
#include <zInterface/functionsets/zFnMesh.h>
#include <zApp/include/zFnSets.h>


namespace zSpace
{
	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsColorSplit
	*	\brief A mesh split tool set class depends on mesh face colors.
	*	\method based on Breadth-first search (BFS) method
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/


	class ZSPACE_TOOLS zTsColorSplit 
	{
	protected:		

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;

		/*!	\brief pointer to input Object  */
		zObjMesh *oMesh;

		/*!	\brief container to dual graph of input Object  */
		zObjGraph dualGraph;

		/*!	\brief container to face colors  */
		zColorArray fCols;

		/*!	\brief container to face indices  */
		zIntArray fIDs;

		/*!	\brief container to result face indices under each mesh obj  */
		vector<vector<int>> resultMeshID;

		/*!	\brief container to result mesh obj */
		zObjMeshArray splits;

		/*!	\brief numver of result mesh splits */
		int numMesh;

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
		zTsColorSplit();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsColorSplit();

		//--------------------------
		//---- SET METHODS
		//--------------------------

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_meshObj			- input mesh object.
		*	\since version 0.0.2
		*/
		void setInMesh(zObjMesh& _meshObj);

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the splits of an input mesh based on the face colors.
		*
		*	\since version 0.0.2
		*/
		void compute();

		/*! \brief This method exports color split meshes to the input directory.
		*
		*	\param		[in]	_pth	- export directory.
		*	\param		[in]	_type	- export file type.
		*	\since version 0.0.2
		*/
		void exportTo(string& _pth, zFileTpye _type);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method returns the split meshes.
		*
		*	\since version 0.0.2
		*/
		zObjMeshPointerArray getRawSplitMesh(int& _numMesh);


	private:

		//--------------------------
		//---- PRIVATE UTILITY METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		bool isRecorded(int _id, zIntArray& _recorder);

		bool checkExisted(int _id, zIntArray& _faceID);

		vector<vector<int>> meshContainer(zIntArray faceID);

		void generateMesh();

		zIntArray transformID(zIntArray& v);

	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsColorSplit.cpp>
#endif

#endif