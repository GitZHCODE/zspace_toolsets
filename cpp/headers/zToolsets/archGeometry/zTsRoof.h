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

#ifndef ZSPACE_TS_GEOMETRY_ROOFANALYSIS_H
#define ZSPACE_TS_GEOMETRY_ROOFANALYSIS_H


#pragma once

//#include <headers/zCore/base/zExtern.h>
#include <headers/base/zSpace_Toolsets.h>


#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <headers/zToolsets/archGeometry/zTsPaneling.h>
#include <headers/zApp/include/zInterOp.h>

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/gaussian_curvature.h>
#include <igl/read_triangle_mesh.h>

#include <headers/zInterOp/objects/zObjNurbsCurve.h>
#include <headers/zInterOp/functionSets/zFnNurbsCurve.h>

//#include <rhinoSdkStdafxPreamble.h>
//#include <rhinoSdk.h>

using namespace std;

namespace zSpace
{
	///** \addtogroup zToolsets
	//*	\brief Collection of toolsets for applications.
	//*  @{
	//*/

	///** \addtogroup zTsArchGeometry
	//*	\brief tool sets for geometry related utilities.
	//*  @{
	//*/

	///*! \class zTsRoof
	//*	\brief A toolset for panels analysis on freeform mesh.
	//*	\since version 0.0.4
	//*/

	///** @}*/

	///** @}*/


	class ZSPACE_TOOLS zTsRoof
	{
	protected:
		/*!	\brief coreUtils utilities Object  */
		zUtilsCore coreUtils;
		
		//<zPanel> panels;
		//zObjGraph dualGraph;
		//zIntArray cols;
	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------
		zObjMesh* o_guideMesh; //why is this zObjMesh* o_guideMesh;
		zObjGraph* grid;
		zObjNurbsCurve* nurbsGuideCurve;
		/*!	\brief form function set  */
		zFnMesh fnMesh;
		vector<zPanel> panels;
		zFourColorsMesh c_mesh;

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------
		/*! \brief Default constructor.
		*
		*	\since version 0.0.2
		*/
		zTsRoof();
		//--------------------------
		//---- DESTRUCTOR
		//--------------------------
		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsRoof();

		//--------------------------
        //--- SET METHODS 
        //--------------------------
		void setGuideMesh(zObjMesh& o_mesh);
		void setGuideGrid(zObjGraph& grid);
		void setGuideCurves(zObjNurbsCurve& nurbsGuide);

		//--------------------------
		//--- GET METHODS 
		//--------------------------
		
		zObjMesh* getRawGuideMesh();
		
		int getNumPanels();

		vector<zPanelType> getPanelTypes();

		zPanel* getRawPanels(int &numPanels);

		//zObjMesh* getRawColorMesh();

		//--------------------------
		//---- CREATE METHODS
		//--------------------------
		void colorGuideMesh();

		void colorGuideMeshFour();

		void createPanels();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the curvature type.
		*
		*	\param		[in]	x			- input x.
		*	\since version 0.0.4
		*/
		void computePanelType(float planarityDeviation);

		//--------------------------
		//---- EXPORT METHODS
		//--------------------------
		void exportPanelTypes(string folderDir);

		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:
		void createPanelMesh(zItMeshVertex& v, int panelID );
	};
}

#if defined(ZSPACE_STATIC_LIBRARY)  || defined(ZSPACE_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/archGeometry/zTsRoof.cpp>
#endif

#endif