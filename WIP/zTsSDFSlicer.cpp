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


#include<headers/zToolsets/geometry/zTsSDFSlicer.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_INLINE zTsSDFSlicer::zTsSDFSlicer()
	{
		red = zColor(1, 0, 0, 1);
		yellow = zColor(1, 1, 0, 1); 
		green = zColor(0, 1, 0, 1);
		cyan = zColor(0, 1, 1, 1); 
		blue = zColor(0, 0, 1, 1); 
		magenta = zColor(1, 0, 1, 1);

		grey = zColor(0.5, 0.5, 0.5, 1);

		orange = zColor(1, 0.5, 0, 1);
	
	}
	

	//---- DESTRUCTOR

	ZSPACE_INLINE zTsSDFSlicer::~zTsSDFSlicer() {}

	//---- CREATE METHODS


	ZSPACE_INLINE void zTsSDFSlicer::createFieldMesh(zDomain<zPoint>& bb, int resX, int resY)
	{
		zFnMeshScalarField fnField(o_field);

		fnField.create(bb.min, bb.max, resX, resY, 1, true, false);


		zDomainColor dCol(red, green);
		fnField.setFieldColorDomain(dCol);
	}

	//--- SET METHODS 

	ZSPACE_INLINE void zTsSDFSlicer::setSliceMesh(zObjMesh& _o_SliceMesh)
	{
		o_SliceMesh = &_o_SliceMesh;
	}

	ZSPACE_INLINE void zTsSDFSlicer::setGuideGraph(zObjGraph& _o_GuideGraph)
	{
		o_GuideGraph = &_o_GuideGraph;
	}

	ZSPACE_INLINE void zTsSDFSlicer::setStartEndPlanes(zTransform& _sPlane, zTransform& _ePlane)
	{
		startPlane = _sPlane;
		endPlane = _ePlane;
	}

	//---- GET METHODS

	ZSPACE_INLINE vector<zTransform> zTsSDFSlicer::getBlockFrames()
	{	
		return sectionFrames;
	}

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockSectionGraphs( int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;				

		numGraphs = o_sectionGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_sectionGraphs)
		{
			out.push_back(&graph);
		}				

		return out;
	}

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockRaftGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_raftGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_raftGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_INLINE zObjGraphPointerArray zTsSDFSlicer::getBlockContourGraphs(int& numGraphs)
	{
		zObjGraphPointerArray out;
		numGraphs = 0;

		numGraphs = o_contourGraphs.size();

		if (numGraphs == 0)return out;

		for (auto& graph : o_contourGraphs)
		{
			out.push_back(&graph);
		}

		return out;
	}

	ZSPACE_INLINE zObjMeshScalarField* zTsSDFSlicer::getRawFieldMesh()
	{
		return &o_field;
	}

	//---- COMPUTE METHODS

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlocks(float printPlaneSpacing, float printLayerWidth, float raftLayerWidth, zDomainFloat neopreneOffset, bool compFrames, bool compSDF)
	{
		if (compFrames)
		{
			computePrintBlockFrames(printPlaneSpacing, neopreneOffset.min, neopreneOffset.max);
			computePrintBlockSections();		
		}

		if (compSDF)
		{
			computeSDF(printLayerWidth, neopreneOffset.min, raftLayerWidth);
		}

	}

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlockFrames(float printPlaneSpacing, float neopreneOffset_start, float neopreneOffset_end)
	{
		// getLength of guide graph
		zDoubleArray eLens;
		zFnGraph fnGraph(*o_GuideGraph);
		float totalLength = fnGraph.getEdgeLengths(eLens);

		zFloatArray weights = { 0.0, 0.5, 1.0 };
		zFloatArray multVals = { 0.0, 0.5,  1.0 };

		float len = totalLength - (neopreneOffset_start + neopreneOffset_end);

		int numLayers = floor(len / printPlaneSpacing);

		float equalisedPlaneSpacing = len / numLayers;

		printf("\n %i %1.2f %1.2f ", numLayers, len, equalisedPlaneSpacing);
	
		zVector startNorm(startPlane(2, 0), startPlane(2, 1), startPlane(2, 2));
		zVector endNorm(endPlane(2, 0), endPlane(2, 1), endPlane(2, 2));
	
		zPoint startOrig(startPlane(3, 0), startPlane(3, 1), startPlane(3, 2));
		zPoint endOrig(endPlane(3, 0), endPlane(3, 1), endPlane(3, 2));

		zItGraphVertex v(*o_GuideGraph, 0);
		zItGraphHalfEdge startHe = v.getHalfEdge();

		

		// Start point
		
		zPoint O = startOrig;
		zVector Z = startNorm;

		zVector tempZ = Z;
		tempZ.normalize();

		zVector X;
		zVector Y(0,1,0);

		float weight = 0;

		float mult;

		for (int l = 0; l < weights.size() - 1; l++)
		{
			if (weight >= weights[l] && weight <= weights[l + 1])   mult = coreUtils.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
		}

		tempZ.x = (startNorm.x * (1 - mult)) + (endNorm.x * mult);
		tempZ.y = (startNorm.y * (1 - mult)) + (endNorm.y * mult);
		tempZ.z = (startNorm.z * (1 - mult)) + (endNorm.z * mult);
		tempZ.normalize();

		X = Y ^ tempZ;
		X.normalize();

		Y = tempZ ^ X;
		Y.normalize();

		zTransform pFrame = setTransformFromVectors(O, X, Y, tempZ);
		sectionFrames.push_back(pFrame);

		zPoint pOnCurve = O;

		// in between points
		zItGraphHalfEdge walkHe = startHe;
		for (int j = 0; j < numLayers; j++)
		{
			zPoint prevPoint = pOnCurve;

			zPoint eEndPoint =  walkHe.getVertex().getPosition();

			float distance_increment = equalisedPlaneSpacing;

			while (pOnCurve.distanceTo(eEndPoint) < distance_increment)
			{
				distance_increment = distance_increment - pOnCurve.distanceTo(eEndPoint);
				pOnCurve = eEndPoint;

				walkHe = (walkHe.onBoundary()) ? walkHe.getNext() : walkHe.getNext().getSym().getNext();
				eEndPoint = walkHe.getVertex().getPosition();
			}

			zVector he_vec = walkHe.getVector();
			he_vec.normalize();

			//O
			O = pOnCurve + he_vec * distance_increment;
			
			if (j == numLayers - 1)
			{
				O = endOrig;
				Z = endNorm;
			}

			tempZ = Z;
			tempZ.normalize();

			Y = zVector(0, 1, 0);

			weight = (float)(j + 1) / numLayers;			

			for (int l = 0; l < weights.size() - 1; l++)
			{
				if (weight >= weights[l] && weight <= weights[l + 1])   mult = coreUtils.ofMap(weight, weights[l], weights[l + 1], multVals[l], multVals[l + 1]);
			}
			
			tempZ.x = (startNorm.x * (1 - mult)) + (endNorm.x * mult);
			tempZ.y = (startNorm.y * (1 - mult)) + (endNorm.y * mult);
			tempZ.z = (startNorm.z * (1 - mult)) + (endNorm.z * mult);
			tempZ.normalize();

			X = Y ^ tempZ;
			X.normalize();

			Y = tempZ ^ X;
			Y.normalize();

			// add frame
			pFrame = setTransformFromVectors(O, X, Y, tempZ);
			sectionFrames.push_back(pFrame);

			pOnCurve = O;
		}

	}

	ZSPACE_INLINE void zTsSDFSlicer::computePrintBlockSections()
	{
		zFnMesh  fn_sliceMesh(*o_SliceMesh);

		o_sectionGraphs.clear();
		o_sectionGraphs.assign(sectionFrames.size(), zObjGraph());

		zScalarArray scalars;

		int start = 0;
		int end = sectionFrames.size();

		for (int i = start; i < end; i++)
		{
			scalars.clear();
			zPoint O(sectionFrames[i](3, 0), sectionFrames[i](3, 1), sectionFrames[i](3, 2));
			zVector N(sectionFrames[i](2, 0), sectionFrames[i](2, 1), sectionFrames[i](2, 2));

			zVector X(sectionFrames[i](0, 0), sectionFrames[i](0, 1), sectionFrames[i](0, 2));
			X.normalize();

			for (zItMeshVertex v(*o_SliceMesh); !v.end(); v++)
			{
				zPoint P = v.getPosition();
				float minDist_Plane = coreUtils.minDist_Point_Plane(P, O, N);
				scalars.push_back(minDist_Plane);
			}

			zPointArray positions;
			zIntArray edgeConnects;
			zColorArray vColors;
			fn_sliceMesh.getIsoContour(scalars, 0.0, positions, edgeConnects, vColors);

			// create graphs
			zFnGraph tempFn(o_sectionGraphs[i]);
			tempFn.create(positions, edgeConnects);;
		
			tempFn.setEdgeColor(magenta);
		
		}
	}

	ZSPACE_INLINE void zTsSDFSlicer::computeSDF(float printWidth, float neopreneOffset, float raftWidth)
	{
		o_contourGraphs.clear();
		o_contourGraphs.assign(o_sectionGraphs.size(), zObjGraph());	

		o_raftGraphs.clear();
		o_raftGraphs.assign(1, zObjGraph());

		int kStart =  1;
		int kEnd =  o_sectionGraphs.size() -1;


		for (int k = kStart; k < kEnd; k++)
		{
			computeBlockSDF_Internal(k, printWidth, neopreneOffset, false, 0, raftWidth);
		}

	}

	ZSPACE_INLINE void zTsSDFSlicer::computeBlockSDF_Internal(int graphId, float printWidth, float neopreneOffset, bool addRaft, int raftId, float raftWidth)
	{
		if (graphId >= o_sectionGraphs.size())return;

		printf("\n fREP graphID %i ", graphId);

		zFnGraph fnGraph(o_sectionGraphs[graphId]);
		float pWidth = (addRaft) ? printWidth : printWidth;
		float inc3dWidth = (addRaft) ? 0.025 : 0.025;


		zPoint* positions = fnGraph.getRawVertexPositions();

		zTransform t = sectionFrames[graphId];
		fnGraph.setTransform(t, true, false);

		zPoint o(t(3, 0), t(3, 1), t(3, 2));
		zVector n(t(2, 0), t(2, 1), t(2, 2));

		// Transform

		zTransform tLocal;
		tLocal.setIdentity();
		fnGraph.setTransform(tLocal, true, true);

		// field
		zFnMeshScalarField fnField(o_field);

		float offset_outer = 0.5 * pWidth;
		float offset_inner = 1.5 * pWidth;

		// Profile polygon field
		zScalarArray polyField;
		fnField.getScalars_Polygon(polyField, o_sectionGraphs[graphId], false);

		// RESULT FIELDS
		fnField.setFieldValues(polyField);

		zFnGraph fnIsoGraph(o_contourGraphs[graphId]);
		fnField.getIsocontour(o_contourGraphs[graphId], 0.0);

		fnIsoGraph.setEdgeWeight(2);

		// transform back 
		fnGraph.setTransform(t, true, true);
		fnIsoGraph.setTransform(t, true, true);
	}
	
	
	//---- UTILITY METHODS


	ZSPACE_INLINE zTransform zTsSDFSlicer::setTransformFromVectors(zPoint& O, zVector& X, zVector& Y, zVector& Z)
	{
		zTransform out;

		out(0, 0) = X.x; out(0, 1) = X.y; out(0, 2) = X.z; out(0, 3) = 1;
		out(1, 0) = Y.x; out(1, 1) = Y.y; out(1, 2) = Y.z; out(1, 3) = 1;
		out(2, 0) = Z.x; out(2, 1) = Z.y; out(2, 2) = Z.z; out(2, 3) = 1;
		out(3, 0) = O.x; out(3, 1) = O.y; out(3, 2) = O.z; out(3, 3) = 1;


		return out;
	}


	//WIP
	ZSPACE_INLINE void zTsSDFSlicer::CreatePaths(JSON Path, int FieldResolution, float FieldBounds, float* outLeftPoints, float* outRightPoints, int* leftEdgeCOnnects, int* rightEdgeConnects)
	{
		//function body from Vishu
	}



	ZSPACE_INLINE void extCreatePaths(JSON Path, int FieldResolution, float FieldBounds[6], float* outLeftPoints, float* outRightPoints, int* leftEdgeCOnnects, int* rightEdgeConnects)
	{
		//body
		zTsSDSlicer slicer = zTsSDSlicer();
		slicer.CreatePaths(Path, FieldResolution, FieldBounds, outLeftPoints, outRightPoints, leftEdgeCOnnects, rightEdgeConnects);
	}


}