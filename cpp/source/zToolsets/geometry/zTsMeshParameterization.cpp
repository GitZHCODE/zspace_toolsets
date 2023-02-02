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


#include<headers/zToolsets/geometry/zTsMeshParameterization.h>

namespace zSpace
{
	//---- CONSTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsMeshParam::zTsMeshParam()
	{
		

	}


	//---- DESTRUCTOR

	ZSPACE_TOOLSETS_INLINE zTsMeshParam::~zTsMeshParam() {}

	//---- CREATE METHODS


	//--- SET METHODS 

	ZSPACE_TOOLSETS_INLINE void zTsMeshParam::setFromFile(string path, zFileTpye fType)
	{

		bool fileChk = coreUtils.fileExists(path);

		if (!fileChk) return;


		zFnMesh fnParamMesh(o_paramTriMesh);
		fnParamMesh.from(path, fType);

		zFnMesh fnMesh(o_TriMesh);
		fnMesh.from(path, fType);
		fnMesh.computeMeshNormals();

		fnMesh.getMatrices_trimesh(triMesh_V, triMesh_FTris);

		//zPoint* vPositions = fnMesh.getRawVertexPositions();
		//MatrixXd V(fnMesh.numVertices(), 3);

		//// fill vertex matrix
		//for (int i = 0; i < fnMesh.numVertices(); i++)
		//{
		//	V(i, 0) = vPositions[i].x;
		//	V(i, 1) = vPositions[i].y;
		//	V(i, 2) = vPositions[i].z;
		//}

		//triMesh_V = V;

		//// fill triangle matrix
		//MatrixXi FTris(fnMesh.numPolygons(), 3);

		//int nTris = 0;
		//for (zItMeshFace f(o_TriMesh); !f.end(); f++)
		//{
		//	int i = f.getId();

		//	zIntArray fVerts;
		//	f.getVertices(fVerts);

		//	FTris(i, 0) = fVerts[0];
		//	FTris(i, 1) = fVerts[1];
		//	FTris(i, 2) = fVerts[2];
		//}		

		//triMesh_FTris = FTris;



	}

	
	//---- GET METHODS

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsMeshParam::getRawInMesh()
	{
		return &o_TriMesh;

	}

	ZSPACE_TOOLSETS_INLINE zObjMesh* zTsMeshParam::getRawParamMesh()
	{
		return &o_paramTriMesh;
	}

	//--- COMPUTE METHODS 

	ZSPACE_TOOLSETS_INLINE void zTsMeshParam::computeParam_Harmonic()
	{
		Eigen::MatrixXd triMesh_V_uv;

		// Find the open boundary
		Eigen::VectorXi bnd;
		igl::boundary_loop(triMesh_FTris, bnd);

		// Map the boundary to a circle, preserving edge proportions
		Eigen::MatrixXd bnd_uv;
		igl::map_vertices_to_circle(triMesh_V, bnd, bnd_uv);

		// Harmonic parametrization for the internal vertices
		igl::harmonic(triMesh_V, triMesh_FTris, bnd, bnd_uv, 1, triMesh_V_uv);

		// update vertex positions
		zFnMesh fnParamMesh(o_paramTriMesh);
		zPoint* vPositions = fnParamMesh.getRawVertexPositions();
		
		for (int i = 0; i < fnParamMesh.numVertices(); i++)
		{
			vPositions[i].x = triMesh_V_uv(i, 0);
			vPositions[i].y = triMesh_V_uv(i, 1);
			vPositions[i].z = 0;
		}

	}

	ZSPACE_TOOLSETS_INLINE void zTsMeshParam::computeParam_ARAP()
	{
		Eigen::MatrixXd triMesh_V_uv;
		Eigen::MatrixXd initial_guess;

		// Compute the initial solution for ARAP (harmonic parametrization)
		Eigen::VectorXi bnd;
		igl::boundary_loop(triMesh_FTris, bnd);
		Eigen::MatrixXd bnd_uv;
		igl::map_vertices_to_circle(triMesh_V, bnd, bnd_uv);

		igl::harmonic(triMesh_V, triMesh_FTris, bnd, bnd_uv, 1, initial_guess);

		// Add dynamic regularization to avoid to specify boundary conditions
		igl::ARAPData arap_data;
		arap_data.with_dynamics = true;
		Eigen::VectorXi b = Eigen::VectorXi::Zero(0);
		Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0, 0);

		// Initialize ARAP
		arap_data.max_iter = 100;
		// 2 means that we're going to *solve* in 2d
		arap_precomputation(triMesh_V, triMesh_FTris, 2, b, arap_data);


		// Solve arap using the harmonic map as initial guess
		triMesh_V_uv = initial_guess;

		arap_solve(bc, arap_data, triMesh_V_uv);

		// update vertex positions
		zFnMesh fnParamMesh(o_paramTriMesh);
		zPoint* vPositions = fnParamMesh.getRawVertexPositions();

		for (int i = 0; i < fnParamMesh.numVertices(); i++)
		{
			vPositions[i].x = triMesh_V_uv(i, 0);
			vPositions[i].y = triMesh_V_uv(i, 1);
			vPositions[i].z = 0;
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsMeshParam::computeParam_LSCM()
	{
		Eigen::MatrixXd triMesh_V_uv;

		// Fix two points on the boundary
		VectorXi bnd, b(2, 1);
		igl::boundary_loop(triMesh_FTris, bnd);
		b(0) = bnd(0);
		b(1) = bnd(bnd.size() / 2);
		MatrixXd bc(2, 2);
		bc << 0, 0, 1, 0;

		// LSCM parametrization
		igl::lscm(triMesh_V, triMesh_FTris, b, bc, triMesh_V_uv);

		// update vertex positions
		zFnMesh fnParamMesh(o_paramTriMesh);
		zPoint* vPositions = fnParamMesh.getRawVertexPositions();

		for (int i = 0; i < fnParamMesh.numVertices(); i++)
		{
			vPositions[i].x = triMesh_V_uv(i, 0);
			vPositions[i].y = triMesh_V_uv(i, 1);
			vPositions[i].z = 0;
		}
	}

	ZSPACE_TOOLSETS_INLINE void zTsMeshParam::computeGeodesics_Heat(int _vertexID)
	{
		// Precomputation
		igl::HeatGeodesicsData<double> data;
		double t = std::pow(igl::avg_edge_length(triMesh_V, triMesh_FTris), 2);
		
		if (!igl::heat_geodesics_precompute(triMesh_V, triMesh_FTris, t, data))
		{
			cout << "Error: heat_geodesics_precompute failed." << endl;
			return;
		};
		

		Eigen::VectorXd D;
				
		igl::heat_geodesics_solve(data, (Eigen::VectorXi(1, 1) << _vertexID).finished(), D);

		// compute vertex color from Geodesic Distance
		zFnMesh fnTriMesh(o_TriMesh);
		zColor* vColors = fnTriMesh.getRawVertexColors();

		zDomainFloat distanceDomain(D.minCoeff(), D.maxCoeff());
		zDomainColor colDomain(zColor(1,0,0,1), zColor(0, 1, 0, 1));

		for (int i = 0; i < fnTriMesh.numVertices(); i++)
		{
			vColors[i] = coreUtils.blendColor(D(i), distanceDomain, colDomain, zRGB);
		}

		fnTriMesh.computeFaceColorfromVertexColor();

		
	}

}

