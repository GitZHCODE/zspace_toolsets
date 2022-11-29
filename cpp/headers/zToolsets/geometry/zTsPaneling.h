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

#ifndef ZSPACE_TS_GEOMETRY_PANELANALYSIS_H
#define ZSPACE_TS_GEOMETRY_PANELANALYSIS_H


#pragma once

#include <headers/zCore/base/zExtern.h>

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


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

	/*! \class zTsPanelAnalysis
	*	\brief A toolset for panels analysis on freeform mesh.
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/
	enum zPanelType { zPlanar = 1000, zConcaveEllipsoid, zConcaveCylinder, zHyperboloid, zConvexCylinder, zConvexEllipsoid, zCustom };

	class ZSPACE_TOOLS zPanel
	{
	private:

		/*!	\brief core utilities object  */
		zUtilsCore coreUtils;
		zObjMesh o_panelMesh;
		zObjMesh o_panelMesh_tri;
		zCurvatureArray curvature;
		zColor red, yellow, green, cyan, blue, magenta, grey, orange;
		zDomainDouble tolerance;

	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------

		zFnMesh fnMesh;
		double k1, k2;
		zPanelType panelType;
		string filename;

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.4
		*/
		zPanel();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zPanel();

		//--------------------------
		//--- SET METHODS 
		//--------------------------

		//--------------------------
		//--- GET METHODS 
		//--------------------------

		/*! \brief This method gets the raw panel quad mesh.
		*
		*	\since version 0.0.4
		*/
		zObjMesh* getRawPanelMesh();

		/*! \brief This method gets the raw panel triangulated mesh.
		*
		*	\since version 0.0.4
		*/
		zObjMesh* getRawPanelMesh_Tri();

		/*! \brief This method gets panel Curvature.
		*
		*	\since version 0.0.4
		*/
		zCurvatureArray  getCurvature();

		zPanelType getType();

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes panel curvature.
		*
		*	\param		[in]	x			- input x.
		*	\since version 0.0.4
		*/
		void computeCurvature();

		/*! \brief This method computes panel curvature type.
	*
	*	\param		[in]	x			- input x.
	*	\since version 0.0.4
	*/
		void  computeType();

		/*! \brief This method returns the curvature type of a panel.
		*
		*	\param		[in]	fCenters			- input fCenters
		*	\since version 0.0.4
		*/
		bool isPanelPlanar(float tol = EPS);

		bool isPanelConcaveEllipsoid();

		bool isPanelConcaveCylinder();

		bool isPanelHyperboloid();

		bool isPanelConvexCylinder();

		bool isPanelConvexEllipsoid();

		//--------------------------
		//---- UTILITY METHODS
		//--------------------------
		void color_Type();

		//--------------------------
		//---- PROTECTED UTILITY METHODS
		//--------------------------
	protected:
	};

	class zFourColorsMesh
	{
	protected:
			
	public:

		zObjMesh* o_colorMesh;
		zObjGraph dualGraph;
		zIntArray cols;
		zColor red, yellow, green, cyan, blue, magenta, grey, orange;

	//--------------------------
	//---- CONSTRUCTOR
	//--------------------------

	/*! \brief Default constructor.
	*
	*	\since version 0.0.4
	*/
		zFourColorsMesh();

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.4
		*/
		~zFourColorsMesh();

		//--------------------------
		//--- SET METHODS 
		//--------------------------
		void setMesh(zObjMesh& o_mesh);

		//--------------------------
		//--- COMPUTE METHODS 
		//--------------------------

		/*! \brief This method gets the raw panel quad mesh.
		*
		*	\since version 0.0.4
		*/
		void computeColors();

		void setMeshColor();

	protected:

		bool checkExisted(int _id);
		
		int findColor(int _id);
		
	};

	class ZSPACE_TOOLS zTsPaneling
	{
	protected:
		/*!	\brief coreUtils utilities Object  */
		zUtilsCore coreUtils;
		zObjMesh* o_guideMesh; //why is this zObjMesh* o_guideMesh;
		//<zPanel> panels;
		zObjGraph dualGraph;
		zIntArray cols;
	public:
		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------
		zColor red, yellow, green, cyan, blue, magenta, grey, orange;
		
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
		zTsPaneling();
		//--------------------------
		//---- DESTRUCTOR
		//--------------------------
		/*! \brief Default destructor.
		*
		*	\since version 0.0.2
		*/
		~zTsPaneling();

		//--------------------------
        //--- SET METHODS 
        //--------------------------
		void setGuideMesh(zObjMesh& o_mesh);

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
		void computePanelType();

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
#include<source/zToolsets/geometry/zTsPaneling.cpp>
#endif

#endif