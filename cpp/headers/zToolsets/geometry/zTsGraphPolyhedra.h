// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>, Leo Bieling <leo.bieling@zaha-hadid.com>
//

#ifndef ZSPACE_TS_GEOMETRY_GRAPH_POLYHEDRA_H
#define ZSPACE_TS_GEOMETRY_GRAPH_POLYHEDRA_H

#pragma once
#include <headers/base/zSpace_Toolsets.h>
#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>
#include <headers/zInterface/model/zModel.h>

namespace zSpace
{
	struct zPair_hash {
		template <class T1, class T2>
		size_t operator()(const pair<T1, T2>& p) const
		{
			auto hash1 = hash<T1>{}(p.first);
			auto hash2 = hash<T2>{}(p.second);

			if (hash1 != hash2) {
				return hash1 ^ hash2;
			}

			// If hash1 == hash2, their XOR is zero.
			return hash1;
		}
	};

	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsGraphPolyhedra
	*	\brief A function set to convert graph data to polyhedra.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLSETS zTsGraphPolyhedra
	{
	protected:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief pointer to graph Object  */
		zObjGraph o_formGraph;

		/*!	\brief container of  form particle objects  */
		zObjParticleArray o_formParticles;

		/*!	\brief container of form particle function set  */
		vector<zFnParticle> fnFormParticles;

		/*!	\brief DISCRIPTION  */
		zUtilsCore coreUtils;

		/*!	\brief DISCRIPTION  */
		zObjMeshArray o_convexHullMeshes;

		/*!	\brief DISCRIPTION  */
		zObjMeshArray o_forceMeshes;

		/*!	\brief container of  force particle objects  */
		vector<zObjParticleArray> o_forceParticles;

		/*!	\brief container of force particle function set  */
		vector<vector<zFnParticle>> fnForceParticles;

		/*!	\brief DISCRIPTION  */

		zColorArray colors;
		unordered_map<int, zIntPair> formHalfedge_forceCellFace;
		unordered_map<zIntPair, int, zPair_hash> forceCellFace_formHalfedge;

		//vector<zIntPairArray> c_graphHalfEdge_dualCellFace;
		
		zVectorArray  form_targetHalfEdges;

		vector<zVectorArray> force_targetNormals;
		vector<zDoubleArray> force_targetAreas;

		/*!	\brief DISCRIPTION  */
		vector<pair<zPoint, zPoint>> dualConnectivityLines;

		/*!	\brief DISCRIPTION  */
		zIntArray nodeId;

	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------

		/*!	\brief form function set  */
		zDomainDouble deviations[2];


		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsGraphPolyhedra();

	
		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsGraphPolyhedra();

		//--------------------------
		//---- SET METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void setFormGraphFromFile(string& _path, zFileTpye _type, bool _staticGeom = false);

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void setFormGraphFromMesh(zObjMesh& _inMeshObj, zVector& _verticalForce);

		//--------------------------
		//---- GET METHODS
		//--------------------------

		/*! \brief This method gets pointer to the internal form graph object.
		*
		*	\return				zObjGraph*					- pointer to internal graph object.
		*	\since version 0.0.4
		*/
		zObjGraph* getRawFormGraph();

		/*! \brief This method gets the block section graphs
		*
		*	\param		[out]	numMeshes				- output number of meshes.
		*	\return				zObjMeshPointerArray	- pointer conatiner of meshes if they exist.
		*	\since version 0.0.2
		*/
		zObjMeshPointerArray getRawForceMeshes(int& numMeshes);


		//--------------------------
		//---- CREATE METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void createForceMeshes();
	
		//--------------------------
		//---- UPDATE METHODS
		//--------------------------

		bool equilibrium(bool& compTargets, float formWeight, float areaScale,  double dT, zIntergrationType type, float minMax_formEdge = 0.1, float minMax_forceEdge = 0.1, int numIterations = 1000, double angleTolerance = EPS, double areaTolerance = EPS, bool printInfo = false);




	private:
		//--------------------------
		//---- PRIVATE CREATE METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void createForceMesh(zItGraphVertex& _graphVertex, double _tol);

		//--------------------------
		//---- PRIVATE COMPUTE / UPDATE METHODS
		//--------------------------

		void computeTargets(float formWeight, float areaScale);

		void updateFormDiagram(float minmax_Edge, float dT, zIntergrationType type, int numIterations = 1);

		void updateForceDiagram(float minmax_Edge, float dT, zIntergrationType type, int numIterations = 1);

		bool checkDeviations(zDomainDouble& angleDeviations, double& angleTolerance, zDomainDouble& areaDeviations, double& areaTolerance, bool& printInfo);

		//--------------------------
		//---- PRIVATE UTILITY METHODS
		//--------------------------

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void cleanConvexHull(zItGraphVertex& _vIt);

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void colorDualFaceConnectivity();

		/*! \brief DISCRIPTION
		*
		*	\since version 0.0.4
		*/
		void colorGraphEdges();

	};
}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsGraphPolyhedra.cpp>
#endif

#endif