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

#ifndef ZSPACE_TS_GEOMETRY_MESHPARAMETERIZATION_H
#define ZSPACE_TS_GEOMETRY_MESHPARAMETERIZATION_H



#pragma once

#include <headers/base/zSpace_Toolsets.h>

#include <headers/zCore/base/zExtern.h>

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zInterface/functionsets/zFnMeshField.h>


#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>

#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/arap.h>

#include <igl/heat_geodesics.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/avg_edge_length.h>
#include <igl/isolines_map.h>

#include <igl/point_mesh_squared_distance.h>

#include <igl/copyleft/comiso/nrosy.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include<execution>

//#include<atlsafe.h>
//#include <comdef.h>
//#include<oaidl.h>

using namespace std;


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

	
	/*! \class zTsMeshParam
	*	\brief A toolset for Mesh Parameterization.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLSETS zTsMeshParam
	{
	//protected:
	public:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;			

		/*!	\brief input mesh object  */
		zObjMesh o_TriMesh;

		MatrixXd triMesh_V;
		MatrixXi triMesh_FTris;	

		/*!	\brief input mesh object  */
		zObjMesh o_paramTriMesh;
			

	public:

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zTsMeshParam();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zTsMeshParam();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------




		//--------------------------
		//--- SET METHODS 
		//--------------------------

		void setFromFile(string path, zFileTpye fType);

		//--------------------------
		//--- GET METHODS 
		//--------------------------
		
		/*! \brief This method gets pointer to the internal input triangle mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawInMesh();

		/*! \brief This method gets pointer to the internal parameterised triangle mesh object.
		*
		*	\return				zObjMesh*					- pointer to internal mesh object.
		*	\since version 0.0.4
		*/
		zObjMesh* getRawParamMesh();

		//--------------------------
		//--- COMPUTE METHODS 
		//--------------------------

		/*! \brief This method computes the harmonic parametrization of the input triangle mesh.
		*
		*	\since version 0.0.4
		*/
		void computeParam_Harmonic();

		/*! \brief This method computes the as rigid as possible (ARAP) parametrization of the input triangle mesh.
		*
		*	\since version 0.0.4
		*/
		void computeParam_ARAP();

		/*! \brief This method computes the least square conformal map (LSCM) of the input triangle mesh.
		*
		*	\since version 0.0.4
		*/
		void computeParam_LSCM();

		/*! \brief This method computes the least square conformal map (LSCM) of the input triangle mesh.
		*
		*	\since version 0.0.4
		*/
		void compute_NRosy();

		void computeGeodesics_Heat(zIntArray &_vertexIDs, zFloatArray &dScalars);

	};

}





#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsMeshParameterization.cpp>
#endif

#endif