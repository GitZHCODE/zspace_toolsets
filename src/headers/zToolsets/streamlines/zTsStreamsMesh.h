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

#ifndef ZSPACE_TS_STREAMLINES_MESH_H
#define ZSPACE_TS_STREAMLINES_MESH_H

#pragma once
#include "headers/base/zSpace_Toolsets.h"
#include <headers/zToolsets/streamlines/zTsStreams2D.h>


#include <igl/point_mesh_squared_distance.h>

namespace zSpace
{

	
	/** \addtogroup zToolsets
	*	\brief Collection of tool sets for applications. 
	*  @{
	*/

	/** \addtogroup zTsStreamlines
	*	\brief tool sets for field stream lines.
	*  @{
	*/

	/*! \class zTsStreamsMesh
	*	\brief A streamlines tool set for creating streams on a mesh.
	*	\details Based on evenly spaced streamlines (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.29.9498&rep=rep1&integrationType=pdf)
	*	\since version 0.0.2
	*/

	/** @}*/

	/** @}*/

	class ZSPACE_TOOLSETS zTsStreamsMesh
	{
	protected:
		//--------------------------
		//---- PROTECTED ATTRIBUTES
		//--------------------------

		/*!	\brief core utilities Object  */
		zUtilsCore coreUtils;
		
		/*!	\brief pointer to a mesh object  */
		zObjMesh* meshObj;

		/*!	\brief container of vectors per vertex  */
		vector<zVector> fieldVectors;

		vector<zIntArray> ringNeighbours;

		/*!<\brief seperation distance between stream lines.*/
		double dSep;

		/*!<\brief test seperation distance between stream lines.*/
		double dTest;

		/*!<\brief minimum length of stream.*/
		double minLength;

		/*!<\brief maximum length of stream.*/
		double maxLength;

		/*!<\brief streamType - zForwardbackward / zForward/ zBackward.*/
		zFieldStreamType streamType;

		/*!<\brief timestep for integration.*/
		double dT;

		/*!<\brief integration integrationType - zEuler / zRK4.*/
		zIntergrationType integrationType;

		/*!<\brief angle of stream rotation.*/
		double angle;

		/*!<\brief boolean is true if the backward direction is flipped.*/
		bool flipBackward;

		/*!<\brief vertex matrix for closest point compute using IGL.*/
		MatrixXd mesh_V;

		/*!<\brief face matrix for closest point compute using IGL.*/
		MatrixXi mesh_F;

	public:

		//--------------------------
		//---- PUBLIC ATTRIBUTES
		//--------------------------
				

		/*!	\brief 2 dimensional container of stream positions per field index.  */
		vector<zPointArray> fieldIndex_streamPositions;		

		//--------------------------
		//---- CONSTRUCTOR
		//--------------------------

		/*! \brief Default constructor.
		*
		*	\since version 0.0.1
		*/
		zTsStreamsMesh();

		/*! \brief Overloaded constructor.
		*
		*	\param		[in]	_field			- input vector field 2D.
		*	\since version 0.0.1
		*/
		zTsStreamsMesh(zObjMesh &_meshObj, vector<zVector> &_fieldVectors);

		//--------------------------
		//---- DESTRUCTOR
		//--------------------------

		/*! \brief Default destructor.
		*
		*	\since version 0.0.1
		*/
		~zTsStreamsMesh();
		
		//--------------------------
		//----  SET METHODS
		//--------------------------

		/*! \brief This method sets the separation distance.
		*
		*	\param	[in]	_dSep		- input separation distance.
		*	\since version 0.0.1
		*/
		void setSeperationDistance(double _dSep);

		/*! \brief This method sets the test separation distance.
		*
		*	\param	[in]	_dSep		- input test separation distance.
		*	\since version 0.0.1
		*/
		void setTestSeperationDistance(double _dTest);

		/*! \brief This method sets the minimum length of the streams.
		*
		*	\param	[in]	_minLength		- input maximum length.
		*	\since version 0.0.1
		*/
		void setMinLength(double _minLength);

		/*! \brief This method sets the maximum length of the streams.
		*
		*	\param	[in]	_maxLength		- input maximum length.
		*	\since version 0.0.1
		*/
		void setMaxLength(double _maxLength);

		/*! \brief This method sets the stream type.
		*
		*	\param	[in]	_streamType		- input stream type
		*	\since version 0.0.1
		*/
		void setIntegrationType(zFieldStreamType _streamType);

		/*! \brief This method sets the integration time step.
		*
		*	\param	[in]	_dT		- input time step
		*	\since version 0.0.1
		*/
		void setTimeStep(double _dT);

		/*! \brief This method sets the integration type.
		*
		*	\param	[in]	_integrationType		- input integration type.
		*	\since version 0.0.1
		*/
		void setIntegrationType(zIntergrationType _integrationType);

		/*! \brief This method sets the rotation angle.
		*
		*	\param	[in]	_angle		- input angle
		*	\since version 0.0.1
		*/
		void setAngle(double _angle);

		/*! \brief This method sets the flip backwards boolean.
		*
		*	\param	[in]	_flipBackward		- input flip backwards boolean.
		*	\since version 0.0.1
		*/
		void setFlipBackwards(bool _flipBackward);
		
		//--------------------------
		//----  2D STREAM LINES METHODS
		//--------------------------

		/*! \brief This method creates the stream lines and stores them as a graph.
		*
		*	\param	[out]	streams							- container of streams.
		*	\param	[in]	start_seedPoints				- container of start seed positions. If empty a random position in the field is considered.
		*	\param	[in]	seedStreamsOnly					- generates streams from the seed points only if true.
		*	\since version 0.0.1
		*/
		void createStreams(vector<zStreamLine>& streams, vector<zVector> &start_seedPoints, bool seedStreamsOnly = false, int maxStreams = 10);

		//--------------------------
		//----  2D STREAM LINES METHODS WITH INFLUENCE SCALAR FIELD
		//--------------------------
		
		/*! \brief This method creates the stream lines and stores them as a graph.
		*
		*	\param	[out]	streams							- container of streams.
		*	\param	[in]	start_seedPoints				- container of start seed positions. If empty a random position in the field is considered.
		*	\param	[in]	influenceField					- input scalar field.
		*	\param	[in]	min_Power						- input minimum power value.
		*	\param	[in]	max_Power						- input maximum power value.
		*	\param	[in]	seedStreamsOnly					- generates streams from the seed points only if true.
		*	\since version 0.0.1
		*/
		void createStreams_Influence(vector<zStreamLine>& streams, vector<zVector> &start_seedPoints, zFnMeshField<zScalar>& fnInfluenceField, double min_Power, double max_Power, bool seedStreamsOnly = false);

		//--------------------------
		//---- COMPUTE METHODS
		//--------------------------

		/*! \brief This method computes the closest point on the mesh given an input point.
		*
		*	\param	[in]	inPoint			- input point.
		*	\param	[out]	faceID			- output faceID.
		*	\param	[out]	closestPoint	- output closest point.
		*	\since version 0.0.4
		*/
		void computeClosestPointToMesh(zPoint& inPoint, int& faceID, zPoint& closestPoint);

		bool getFieldValue(zPoint &inPoint, zVector& fieldForce, zVector &axis);

		bool projectToMesh(zObjGraph& inGraph, zIntArray &faceIDs);
		//--------------------------
		//---- PROTECTED METHODS
		//--------------------------
	protected:

		/*! \brief This method creates a single stream line as a graph.
		*
		*	\param	[in]	streamGraph						- stream graph created from the field.
		*	\param	[in]	seedPoint						- input seed point.
		*	\return			bool							- true if the graph is created.
		*	\since version 0.0.1
		*/
		bool createStreamGraph(zObjGraph &streamGraphObj, zVector &seedPoint);

		/*! \brief This method computes the seed points.
		*
		*	\param	[in]	currentStream					- input current stream.
		*	\param	[in]	vertexId						- vertex index in the vertices contatiner of the stream graph.
		*	\param	[in]	seedPoints						- container of seed points.
		*	\since version 0.0.1
		*/
		void getSeedPoints(zStreamLine& currentStream, int vertexId, vector<zVector> &seedPoints);
		
		/*! \brief This method creates a single stream line as a graph based on a influence scalar field.
		*
		*	\param	[in]	streamGraph						- stream graph created from the field.
		*	\param	[in]	seedPoint						- input seed point.
		*	\param	[in]	influenceField					- input scalar field.
		*	\param	[in]	min_Power						- input minimum power value.
		*	\param	[in]	max_Power						- input maximum power value.
		*	\return			bool							- true if the graph is created.
		*	\since version 0.0.1
		*/
		bool createStreamGraph_Influence(zObjGraph &streamGraphObj, zVector &seedPoint, zFnMeshField<zScalar>& fnInfluenceField, double min_Power, double max_Power);

		/*! \brief This method computes the seed points.
		*
		*	\param	[in]	influenceField					- input scalar field.
		*	\param	[in]	currentStream				- input current stream line.
		*	\param	[in]	vertexId						- vertex index in the vertices contatiner of the stream graph.
		*	\param	[in]	min_Power						- input minimum power value.
		*	\param	[in]	max_Power						- input maximum power value.
		*	\param	[in]	seedPoints						- container of seed points.
		*	\since version 0.0.1
		*/
		void getSeedPoints_Influence(zFnMeshField<zScalar>& fnInfluenceField, zStreamLine& currentStream, int vertexId, double min_Power, double max_Power, vector<zVector> &seedPoints);

		//--------------------------
		//----  2D FIELD UTILITIES
		//--------------------------	

		/*! \brief This method checks if the input position is in the bounds of the field.
		*
		*	\param	[in]	inPoint		- input point.
		*	\return			bool		- true if the input position is in bounds.
		*	\since version 0.0.1
		*/
		bool checkFieldBounds(zPoint&inPoint, zPoint &cP , int &faceID);


		/*! \brief This method checks if the input position is a valid stream position.
		*
		*	\param	[in]	inPoint							- input point.
		*	\param	[in]	dTest							- dtest is a percentage of dsep. It is the minimal distance under which the integration of the streamline will be stopped in the current direction.
		*	\return			bool							- true if the input position is valid stream point.
		*	\since version 0.0.1
		*/
		bool checkValidStreamPosition(zVector &inPoint, double &dTest);

		/*! \brief This method checks if the input position is a valid seed position.
		*
		*	\param	[in]	inPoint							- input point.
		*	\param	[in]	dTest							- dtest is a percentage of dsep. It is the minimal distance under which the integration of the streamline will be stopped in the current direction.
		*	\return			bool							- true if the input position is valid stream point.
		*	\since version 0.0.1
		*/
		bool checkValidSeedPosition(zVector &inPoint, double &dSep);

		/*! \brief This method adds the input position to the field stream position container.
		*
		*	\param	[in]	inPoint							- input point.
		*	\param	[in]	dTest							- dtest is a percentage of dsep. It is the minimal distance under which the integration of the streamline will be stopped in the current direction.
		*	\return			bool							- true if the input position is valid stream point.
		*	\since version 0.0.1
		*/
		void addToFieldStreamPositions(zVector &inPoint, int faceId);

		void computeRingNeighbours();
	};

}

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/streamlines/zTsStreamsMesh.cpp>
#endif

#endif