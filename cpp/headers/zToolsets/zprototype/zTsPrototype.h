// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms,
// data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : zSpace Prototype Toolset
//

#ifndef ZSPACE_TS_PROTOTYPE_H
#define ZSPACE_TS_PROTOTYPE_H

#pragma once
#include <headers/base/zSpace_Toolsets.h>

#include <zspace/interface.h>
#include <zspace/zCore/utilities/zUtilsCore.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>

#include <string>
#include <unordered_map>
#include <vector>

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

	/*! \struct PrototypeSettings
	*	\brief Parametric inputs for the zTsPrototype arch-bridge generator.
	*/
	struct PrototypeSettings
	{
		float ColumnWidth = 0.0f;
		float ColumnHeight = 0.0f;
		float ArchHeightCoefficient = 0.0f;
		float BlockSize = 1.0f;

		float BeamTopOffset = 0.0f;
		float ArchTopOffset = 0.0f;
		float ArchBottomOffset = 0.0f;
		float ColumnBottomOffset = 0.0f;

		int SubdivisionCount = 0;
	};

	/*! \struct ColumnPier
	*	\brief A pier, or a zero-width phantom pier standing in for the x = 0 / x = Width outer edges.
	*/
	struct ColumnPier
	{
		float center = 0.0f;
		float halfWidth = 0.0f;
		bool isPhantom = false;
	};

	enum GraphEdgeType
	{
		BeamEdge = 0,
		ColumnEdge = 1,
		Other = 2
	};

	/*! \struct Span
	*	\brief The arch spanning between two adjacent piers - a full ellipse
	*	between two real piers, a half ellipse against an outer edge.
	*/
	struct Span
	{
		float xLeft = 0.0f;
		float xRight = 0.0f;
		float radius = 0.0f;
		float archHeight = 0.0f;
		zPoint center;
		float angleStart = 180.0f;
		float angleEnd = 0.0f;

		int panelCount = 1;

		zIntArray curveVertexIds;
		zIntArray beamTopVertexIds;
	};

	/*! \struct Plane
	*	\brief A plane frame - origin plus its own x/y axes and normal.
	*/
	struct Plane
	{
		zPoint origin;
		zPoint xAxis;
		zPoint yAxis;
		zPoint normal;
	};

	/*! \class zTsPrototype
	*	\brief Builds the panelised mesh topology of a multi-span elliptical
	*	arch bridge from an input "comb" graph (i_Graph) - Width, Height,
	*	the bridge's own plane, and every pier's position are all
	*	auto-detected from it - a pier shaft height, and a constant-thickness
	*	beam.
	*/
	class ZSPACE_TOOLSETS zTsPrototype
	{
	public:
		zTsPrototype() { _Status = "Ready"; }
		explicit zTsPrototype(const PrototypeSettings& inSettings) : zTsPrototype() { i_Settings = inSettings; }
		~zTsPrototype() {}

		//--------------------------
		//---- DATA
		//--------------------------

		PrototypeSettings i_Settings;

		/*!	\brief the input "comb" graph the bridge is generated from.  */
		zObjGraph i_Graph;

		std::string _Status;

		/*!	\brief overall span/height, auto-detected from i_Graph.  */
		float _Width = 0.0f;
		float _Height = 0.0f;

		/*!	\brief the plane the input graph lies in.  */
		Plane _Plane;

		/*!	\brief i_Graph's vertex positions in the plane's own local 2D coordinates.  */
		zPointArray _LocalGraphPositions;

		/*!	\brief classification of each i_Graph edge (Beam/Column/Other).  */
		zIntArray _GraphEdgeType;

		/*!	\brief piers, including the phantom x = 0 / x = Width end piers at index 0 and back()  */
		std::vector<ColumnPier> _Piers;

		/*!	\brief spans, _Spans[i] sits between _Piers[i] and _Piers[i + 1]  */
		std::vector<Span> _Spans;

		/*!	\brief working buffers used to assemble the output mesh.  */
		zPointArray _Positions;
		zIntArray _PCounts;
		zIntArray _PConnects;

		zIntArray _PanelGroupId;
		zIntArray _PierGroupId;

		zIntArray _BlockType;

		zFloatArray _Thickness;

		zPointArray _PositionsPositive;

		zObjMesh o_Mesh;
		zObjMesh o_MeshPositive;

		/*!	\brief the flat mesh split into one small, uniformly-coloured mesh per _PanelGroupId id.  */
		zObjMeshArray o_FaceBlocks;

		zObjMeshArray o_ClosedMeshes;

		//--------------------------
		//---- SET / GET METHODS
		//--------------------------

		void setSettings(const PrototypeSettings& inSettings) { i_Settings = inSettings; }
		const PrototypeSettings& getSettings() const { return i_Settings; }
		PrototypeSettings& getSettings() { return i_Settings; }

		void setStatus(const std::string& inStatus) { _Status = inStatus; }
		const std::string& getStatus() const { return _Status; }

		/*! \brief Runs the full generation pipeline.
		*/
		bool compute();

		zObjMesh& getMesh() { return o_Mesh; }
		const zObjMesh& getMesh() const { return o_Mesh; }

		zObjMesh& getMeshPositive() { return o_MeshPositive; }
		const zObjMesh& getMeshPositive() const { return o_MeshPositive; }

		zObjMeshArray& getFaceBlocks() { return o_FaceBlocks; }
		const zObjMeshArray& getFaceBlocks() const { return o_FaceBlocks; }
		int getFaceBlockCount() const { return int(o_FaceBlocks.size()); }
		zObjMesh& getFaceBlock(int index) { return o_FaceBlocks[index]; }
		const zObjMesh& getFaceBlock(int index) const { return o_FaceBlocks[index]; }

		zObjMeshArray& getClosedMeshes() { return o_ClosedMeshes; }
		const zObjMeshArray& getClosedMeshes() const { return o_ClosedMeshes; }
		int getClosedMeshCount() const { return int(o_ClosedMeshes.size()); }
		zObjMesh& getClosedMesh(int index) { return o_ClosedMeshes[index]; }
		const zObjMesh& getClosedMesh(int index) const { return o_ClosedMeshes[index]; }

		float getWidth() const { return _Width; }
		float getHeight() const { return _Height; }

		const Plane& getPlane() const { return _Plane; }
		const zPoint& getPlaneOrigin() const { return _Plane.origin; }
		const zPoint& getPlaneXAxis() const { return _Plane.xAxis; }
		const zPoint& getPlaneYAxis() const { return _Plane.yAxis; }
		const zPoint& getPlaneNormal() const { return _Plane.normal; }

		int getPierCount() const { return int(_Piers.size()); }
		int getSpanCount() const { return int(_Spans.size()); }
		const ColumnPier& getPier(int index) const { return _Piers[index]; }
		const Span& getSpan(int index) const { return _Spans[index]; }
	};

	/** @}*/

	/** @}*/
}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
#else
#include<source/zToolsets/zprototype/zTsPrototype.cpp>
#endif

#endif
