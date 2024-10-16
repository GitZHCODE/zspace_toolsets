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

#ifndef ZSPACE_TS_GEOMETRY_NATPOWER_H
#define ZSPACE_TS_GEOMETRY_NATPOWER_H

#pragma once
#include "base/zSpace_Toolsets.h"

#include <zInterface/functionsets/zFnMesh.h>
#include <zApp/include/zFnSets.h>
#include <zInterOp/functionSets/zFnPlane.h>

//#include <zInterOp/functionSets/zFnNurbsCurve.h>

#include <cmath>


namespace zSpace
{

	struct ZSPACE_TOOLSETS zPrintBlock
	{
		int blockID;
		int medialVID;
		int stride;

		float minLayerHeight = 0;
		float maxLayerHeight = 0;
		float totalLength = 0;

		zObjNurbsCurve medial;
		zObjPlane startPlane;
		zObjPlane endPlane;
		zObjMesh blockMesh;
		zObjGraphArray trimGraphs;
		zObjGraphArray contourGraph;

	private:
		void interpolatePlanes(zObjPlane startPlane, zObjPlane endPlane, float t, zObjPlane& outPlane);


	};

	enum zMedialType
	{
		//columnUp, 
		//column,
		Bottom, //<< lower part of the bay (column and arch)
		Top //<< top part of the bay (top arch)
	};
	struct zHalfMedial
	{


	
		zObjGraph o_medial;
		zObjGraphArray o_features;
		vector<zPlane> splitPlanes;

		zItMeshHalfEdgeArray meshHEStartEnd;


		zObjNurbsCurve o_medial_nurbs;
		zObjNurbsCurveArray o_features_nurbs;
		zObjPlaneArray o_splitPlanes;


		zHalfMedial* nextMedial = nullptr;
		zHalfMedial* prevMedial = nullptr;
		zMedialType type;

		zIntArray blockIDs;


		zHalfMedial* symMedial = nullptr;

		//zHalfMedial();
	};

	


	/** \addtogroup zToolsets
	*	\brief Collection of toolsets for applications.
	*  @{
	*/

	/** \addtogroup zTsGeometry
	*	\brief tool sets for geometry related utilities.
	*  @{
	*/

	/*! \class zTsNatPower
	*	\brief .
	*	\method based on 
	*	\since version 0.0.4
	*/

	/** @}*/

	/** @}*/


	class ZSPACE_TOOLSETS zTsNatPower
	{
	protected:
		/*!	\brief core utilities Object  */
		zUtilsCore core;
	public:

		vector<zTransform> interfacePlanes;
		vector<zObjPlaneArray> o_interfacePlanes;

		vector<zPrintBlock> printBlocks;


		//zObjMeshArray blocksMeshes;

		zObjMesh* o_inMesh_inner;
		zObjMesh* o_inMesh_outer;

		//zItMeshHalfEdgeArray itMedialHE;
		//zItMeshHalfEdgeArray medial_inner, medial_outer;
		zIntArray medialMap_inner, medialMap_outer;


		vector<zHalfMedial> medials_inner, medials_outer;



		//vector <vector<zItMeshHalfEdgeArray>> feature1_inner, feature1_outer;
		zIntArray feature1Map_inner, feature1Map_outer;

		//vector<zObjNurbsCurveArray> feature_inner,  feature_outer;





		int startVertexId_inner, startVertexId_outer;
		int meshStride;
		int baysCount = 7;
		
	public:
		//zTsNatPower();
		/*~zTsNatPower();*/
		void setMeshFromPath(string& innerMeshPath, string& outerMeshPath);
	//private:
		void computeMedials();
		void computeInterfacePlanes(float interfaceAngle, float dividLength);
		void getFeatureNurbs(vector<zObjNurbsCurveArray> & outInnerFeature, vector<zObjNurbsCurveArray> & outOuterFeature);
		//void getMedialGraph();

		void getMesh(zObjMesh& innerMesh, zObjMesh& outerMesh);

	private:

		/*
		* PRIVATE METHODS
		*/

		void readMeshJson(string& meshPath, zObjMesh* outMesh, int& outStartVertexInd, int& outStride);
		void computeAllMedialGraphs(zObjMesh& o_mesh, int startVID, int numFeatureLoops, int numStride, vector<int>& _He_Medial_Map, vector<zHalfMedial>& _medials);
		void computeInterfacePlanesOnMedial(zHalfMedial& medialInner, zHalfMedial& medialOuter, float interfaceAngle, float dividLength);
		void computePlaneIntersectionOnMedial(zHalfMedial& medialInner, zHalfMedial& medialOuter, zObjPlaneArray& inPlanes, zObjGraphArray& outGraph);
		void createNurbsFromMedial(zHalfMedial& medial, zObjNurbsCurveArray& outCurve);

		/*
		* GENERIC METHODS
		*/

		void createNurbsFromHEArray(zItMeshHalfEdgeArray& heArray, zObjNurbsCurve& outCurve);
		void intersect_graphPlane(zObjGraph& o_graph, zPlane& inPlane, bool closestPoint, zPointArray& outPoints);
		zPoint getContourPosition(float& threshold, zVector& vertex_lower, zVector& vertex_higher, float& thresholdLow, float& thresholdHigh);
		void isoContour(zObjGraph& o_graph, zScalarArray& vertexScalars, float threshold, zPointArray& contourPoints);
		void interpolatePlanes(zObjPlane startPlane, zObjPlane endPlane, float t, zObjPlane& outPlane);
		void interpolatePlanes(zObjPlane startPlane, zObjPlane endPlane, zFloatArray& params, zObjPlaneArray& outPlanes, zPointArray newOrigins = zPointArray() );
		void createGraphFromHalfEdgeLoop(zObjGraph& o_graph, zItMeshHalfEdgeArray& heLoop, zColor eColor = zBLACK);

		float getGraphCumulativeLengths(zObjGraph& graph, zFloatArray& cumulativeLengths, zFloatArray& cumulativeLengthsReparms);

		void getPointsAtGraphParams(zObjGraph& graph, zFloatArray& params, zPointArray& outPoints, zIntArray& outEdgeIndex);
		void getNurbsMidPointFromGraphVertex(zObjNurbsCurve& crv, zObjGraph& graph, double& outPar, zPoint& outPoint);


	public:
		void curvePlaneIntersection(zObjNurbsCurve& o_curve, zObjPlane& o_plane, zPointArray& intersectionPoint, zFloatArray& crvParameters);
		
	};
}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsNatPower.cpp>
#endif

#endif