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


#include<headers/zToolsets/geometry/zTsPaneling.h>

namespace zSpace
{
	//--zPanel Class

	//---- CONSTRUCTOR

	ZSPACE_INLINE zPanel::zPanel() {
		red = zColor(1, 0, 0, 1);
		yellow = zColor(1, 1, 0, 1);
		green = zColor(0, 1, 0, 1);
		cyan = zColor(0, 1, 1, 1);
		blue = zColor(0, 0, 1, 1);
		magenta = zColor(1, 0, 1, 1);
		grey = zColor(0.5, 0.5, 0.5, 1);
		orange = zColor(1, 0.5, 0, 1);
		tolerance = zDomainDouble(-0.0050, 0.0050);
	}

	//---- DESTRUCTOR

	ZSPACE_INLINE zPanel::~zPanel() {}

	//--- GET METHODS 

	ZSPACE_INLINE zObjMesh* zPanel::getRawPanelMesh()
	{
		return &o_panelMesh;
	}

	ZSPACE_INLINE zObjMesh* zPanel::getRawPanelMesh_Tri()
	{
		return &o_panelMesh_tri;
	}

	ZSPACE_INLINE zCurvatureArray zPanel::getCurvature()
	{
		return curvature;
	}

	ZSPACE_INLINE zPanelType zPanel::getType()
	{
		return panelType;
	}

	//---- COMPUTE METHODS

	ZSPACE_INLINE void zPanel::computeCurvature()
	{
		zFnMesh fnMesh(o_panelMesh_tri);

		zPointArray vPositions;
		fnMesh.getVertexPositions(vPositions);

		MatrixXd V(vPositions.size(), 3);

		curvature.assign(vPositions.size(), zCurvature());

		// fill vertex matrix
		for (int i = 0; i < fnMesh.numVertices(); i++)
		{
			V(i, 0) = vPositions[i].x;
			V(i, 1) = vPositions[i].y;
			V(i, 2) = vPositions[i].z;
		}

		zIntArray pCounts, pConnects;
		fnMesh.getPolygonData(pConnects, pCounts);

		MatrixXi FTris(pCounts.size(), 3);

		for (int i = 0; i < pCounts.size(); i++)
		{
			FTris(i, 0) = pConnects[i * 3 + 0];
			FTris(i, 1) = pConnects[i * 3 + 1];
			FTris(i, 2) = pConnects[i * 3 + 2];
		}

		MatrixXd PD1, PD2;
		VectorXd PV1, PV2;


		igl::principal_curvature(V, FTris, PD1, PD2, PV1, PV2);

		for (int i = 0; i < fnMesh.numVertices(); i++)
		{
			curvature[i].k1 = PV1(i);
			curvature[i].k2 = PV2(i);

			curvature[i].v1[0] = PD1(i, 0);
			curvature[i].v1[1] = PD1(i, 1);
			curvature[i].v1[2] = PD1(i, 2);

			curvature[i].v2[0] = PD2(i, 0);
			curvature[i].v2[1] = PD2(i, 1);
			curvature[i].v2[2] = PD2(i, 2);
		}
	}

	ZSPACE_INLINE void zPanel::computeType()
	{
		//curvature compute types
		//based on k1 & k2 value of center vertex
		
		if (isPanelPlanar())
		{
			panelType = zPlanar;
			printf("\n %s", "zPlanar");
		}
		else if (isPanelConcaveEllipsoid())
		{
			panelType = zConcaveEllipsoid;
			printf("\n %s ", "zConcaveEllipsoid");
		}
		else if (isPanelConcaveCylinder())
		{
			panelType = zConcaveCylinder;
			printf("\n %s ", "zConcaveCylinder");
		}
		else if (isPanelHyperboloid())
		{
			panelType = zHyperboloid;
			printf("\n %s ", "zHyperboloid");
		}
		else if (isPanelConvexCylinder())
		{
			panelType = zConvexCylinder;
			printf("\n %s ", "zConvexCylinder");
		}
		else if (isPanelConvexEllipsoid())
		{
			panelType = zConvexEllipsoid;
			printf("\n %s ", "zConvexEllipsoid");
		}
		else
		{
			panelType = zCustom;
			printf("\n %s ", "zCustom");
		}
	}

	//ZSPACE_INLINE bool zPanel::isPanelPlanar_Vol(float tol)
//{
//	bool out = true;
//	volume calculation
//	zPointArray fCenters;
//	zInt2DArray fTris;
//	zDoubleArray fVols;
//	double panel_vol = 0;
//	zFnMesh fnMesh(o_panelMesh);
//	fnMesh.getMeshFaceVolumes(fTris, fCenters, fVols, true);
//	bool out = true;
//	for (auto& v : fVols)
//	{
//		if (v > tol)
//		{
//			out = false;
//			break;
//		}
//	}
//}

	ZSPACE_INLINE bool zPanel::isPanelPlanar(float tol)
	{
		bool out = true;

		zItMeshVertex v(o_panelMesh, 0);

		zPoint pCenter = v.getPosition();
		zPointArray cornersPts;
		zItMeshHalfEdgeArray cHEdges;
		v.getConnectedHalfEdges(cHEdges);

		for (auto& he : cHEdges)
		{
			cornersPts.push_back(he.getNext().getVertex().getPosition());
		}

		double uA, uB;
		zPoint pA, pB;

		coreUtils.line_lineClosestPoints(cornersPts[0], cornersPts[2], cornersPts[1], cornersPts[3], uA, uB, pA, pB);

		float dist = pA.distanceTo(pCenter);
		if (dist > tol) out = false;

		//k1,k2 addition
		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 > tolerance.min && curvature[vID].k2 > tolerance.min && curvature[vID].k1 < tolerance.max && curvature[vID].k2 < tolerance.max) out = true;

				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);

			}
		}


		return out;
	}

	ZSPACE_INLINE bool zPanel::isPanelConcaveEllipsoid()
	{
		bool out = false;
		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 < tolerance.min && curvature[vID].k2 < tolerance.min) out = true;

				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);

			}
		}
		return out;
	}

	ZSPACE_INLINE bool zPanel::isPanelConcaveCylinder(/*double roundingFactor*/)
	{
		bool out = false;
		

		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 < tolerance.min && curvature[vID].k2< tolerance.max && curvature[vID].k2> tolerance.min) out = true;
				else if (curvature[vID].k1< tolerance.max && curvature[vID].k1> tolerance.min && curvature[vID].k2 < tolerance.min) out = true;
				/*if (curvature[vID].k1 < 0 && curvature[vID].k2 == 0) out = true;
				else if (curvature[vID].k1 == 0 && curvature[vID].k2 < 0) out = true;*/

				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);
			}
		}
		return out;
	}

	ZSPACE_INLINE bool zPanel::isPanelHyperboloid()
	{
		bool out = false;

		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 > tolerance.max && curvature[vID].k2 < tolerance.min) out = true;
				else if (curvature[vID].k1 < tolerance.min && curvature[vID].k2 > tolerance.max) out = true;
				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);
			}
		}
		return out;
	}

	ZSPACE_INLINE bool zPanel::isPanelConvexCylinder()
	{
		bool out = false;

		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 > tolerance.max && (round(curvature[vID].k2) == 0 || round(curvature[vID].k2) == -0)) out = true;
				else if (curvature[vID].k1< tolerance.max && curvature[vID].k1> tolerance.min && curvature[vID].k2 > tolerance.max) out = true;

				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);

			}
		}
		return out;
	}

	ZSPACE_INLINE bool zPanel::isPanelConvexEllipsoid()
	{
		bool out = false;

		for (zItMeshVertex v(o_panelMesh); !v.end(); v++)
		{
			if (!v.onBoundary())
			{
				int vID = v.getId();

				if (curvature[vID].k1 > tolerance.max && curvature[vID].k2 > 0) out = true;

				printf("\nk1,k2| %1.4f %1.4f ", curvature[vID].k1, curvature[vID].k2);

			}
		}
		return out;
	}

	ZSPACE_INLINE void zPanel::color_Type()
	{
		zFnMesh fnMesh(o_panelMesh);

		if (panelType == zPlanar)
		{
			fnMesh.setFaceColor(green);
			//fnMesh.setVertexColor(green);
		}
		else if (panelType == zConcaveEllipsoid)
		{
			fnMesh.setFaceColor(yellow);
			//fnMesh.setVertexColor(yellow);
		}
		else if (panelType == zConcaveCylinder)
		{
			fnMesh.setFaceColor(orange);
			//fnMesh.setVertexColor(orange);
		}
		else if (panelType == zHyperboloid)
		{
			fnMesh.setFaceColor(blue);
			//fnMesh.setVertexColor(blue);
		}
		else if (panelType == zConvexCylinder)
		{
			fnMesh.setFaceColor(magenta);
			//fnMesh.setVertexColor(magenta);
		}
		else if (panelType == zConvexEllipsoid)
		{
			fnMesh.setFaceColor(cyan);
			//fnMesh.setVertexColor(cyan);
		}
		if (panelType == zCustom)
		{
			fnMesh.setFaceColor(grey);
			//fnMesh.setVertexColor(grey);
		}
	}
}

namespace zSpace
{
	//--zFourColorsMesh Class
	//---- CONSTRUCTOR

	ZSPACE_INLINE zFourColorsMesh::zFourColorsMesh() {
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

	ZSPACE_INLINE zFourColorsMesh::~zFourColorsMesh() {}

	//--- SET METHODS 
	
	ZSPACE_INLINE void zFourColorsMesh::setMesh(zObjMesh& o_mesh)
	{
		o_colorMesh = &o_mesh;
	}

	//--- COMPUTE METHODS 
	ZSPACE_INLINE void zFourColorsMesh::computeColors()
	{
		zIntArray inToDual;
		zIntArray dualToIn;

		zFnMesh fnMesh(*o_colorMesh);
		fnMesh.getDualGraph(dualGraph, inToDual, dualToIn, true, false, false);
		const int numF = fnMesh.numPolygons();

		//initialise data
		cols.assign(numF, -1);

		int vNow = 0;
		bool tested = false;
		do
		{
			tested = checkExisted(vNow);
			if (!tested)
			{
				int id;
				int cvCounter = 0;
				zItGraphVertex v(dualGraph, vNow);
				cols[v.getId()] = findColor(v.getId());
				//cout << cols[v.getId()] << endl;

				zIntArray connectedV;
				v.getConnectedVertices(connectedV);
				do
				{
					zItGraphVertex cv(dualGraph, connectedV[cvCounter]);
					if (cols[cv.getId()] != -1)
						cols[cv.getId()] = findColor(cv.getId());

					cvCounter++;
				} while (cvCounter < connectedV.size());
			}
			else vNow++;
		} while (vNow != numF);
	}

	ZSPACE_INLINE void zFourColorsMesh::setMeshColor()
	{
		zFnMesh fnMesh(*o_colorMesh);
		for (zItMeshFace f(*o_colorMesh); !f.end(); f++)
		{
			int id = f.getId();
			if (cols[id] == 0) f.setColor(red);
			else if (cols[id] == 1) f.setColor(blue);
			else if (cols[id] == 2) f.setColor(green);
			else if (cols[id] == 3) f.setColor(magenta);
		}
	}

	//---- PROTECTED UTILITY METHODS
	ZSPACE_INLINE bool zFourColorsMesh::checkExisted(int _id)
	{
		bool found = false;
		for (int i = 0; i < cols.size(); i++)
			if (cols[_id] != -1)
			{
				found = true;
				goto end;
			}
	end:
		return found;
	}

	ZSPACE_INLINE int zFourColorsMesh::findColor(int _id)
	{
		zItGraphVertex v(dualGraph, _id);
		zIntArray connectedV;
		v.getConnectedVertices(connectedV);

		int max = -1;
		int i = 0;
		do
		{
			if (cols[connectedV[i]] - max == 1)
			{
				max++;
				i = 0;
			}
			else
				i++;

		} while (i != connectedV.size());

		max += 1;
		return max;
	}
}
	
namespace zSpace
{
	//---- CONSTRUCTOR
	ZSPACE_INLINE zTsPaneling::zTsPaneling() {
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
	ZSPACE_INLINE zTsPaneling::~zTsPaneling() {}

	//--- SET METHODS
	ZSPACE_INLINE void zTsPaneling::setGuideMesh(zObjMesh& o_mesh)
	{
		//o_guideMesh = o_mesh;
		o_guideMesh = &o_mesh;
	}

	//--- GET METHODS 
	ZSPACE_INLINE zObjMesh* zTsPaneling::getRawGuideMesh()
	{
		//return o_guideMesh;
		return o_guideMesh;
	}

	ZSPACE_INLINE int zTsPaneling::getNumPanels()
	{
		return panels.size();
	}

	ZSPACE_INLINE vector<zPanelType> zTsPaneling::getPanelTypes()
	{
		vector<zPanelType> outTypes;

		for (int i = 0; i < panels.size(); i++)
		{
			outTypes.push_back(panels[i].getType());
		}

		return outTypes;
	}

	ZSPACE_INLINE zPanel* zTsPaneling::getRawPanels(int &numPanels)
	{
		numPanels = panels.size();
		if (panels.size() == 0) return nullptr;
		else return &panels[0];
	}

	//--- CREATE METHODS 	
	ZSPACE_INLINE void zTsPaneling::colorGuideMesh()
	{
		zFnMesh fnMesh(*o_guideMesh);		
		//make all grey
		fnMesh.setFaceColor(grey, true);	
		bool col0 = true;

		for (zItMeshFace f(*o_guideMesh); !f.end(); f++)
		{
			zItMeshFaceArray cFaces;
			f.getConnectedFaces(cFaces);

			bool colors = true;

			bool color01 = true;
			bool color02 = true;
			bool color03 = true;

			for (int i = 0; i < cFaces.size(); i++)
			{
				if (cFaces[i].getColor() == magenta)
				{
					colors = false;
					continue;
				}
			}
			if (colors)f.setColor(magenta);
			else f.setColor(blue);	
		}
	}

	ZSPACE_INLINE void zTsPaneling::colorGuideMeshFour()
	{
		c_mesh.setMesh(*o_guideMesh);
		c_mesh.computeColors();
		c_mesh.setMeshColor();

	/*	for (zItMeshFace f(*o_guideMesh); !f.end(); f++)
		{
			zColor c = f.getColor();
			cout << c.r << endl;
		}*/
	}

	ZSPACE_INLINE void zTsPaneling::createPanels()
	{
		zFnMesh fnGuideMesh(*o_guideMesh);
		//zFnMesh fnGuideMesh(o_guideMesh);
		
		int numPanels = 0;
		zItMeshVertexArray panelVertices;
		//for (zItMeshVertex v(o_guideMesh); !v.end(); v++)
		for (zItMeshVertex v(*o_guideMesh); !v.end(); v++)
		{
			if (v.onBoundary()) continue;

			zItMeshFaceArray cFaces;
			v.getConnectedFaces(cFaces);

			zColor col0 = cFaces[0].getColor();
			bool panelVertex = true;

			for (int i = 1; i < cFaces.size(); i++)
			{
				zColor col = cFaces[i].getColor();
				if (col0 == col) {}
				else panelVertex = false;
			}

			if (panelVertex)
			{
				panelVertices.push_back(v);
			}
		}

		panels.clear();
		panels.assign(panelVertices.size(), zPanel());

		int panelID = 0;
		for (auto& v : panelVertices)
		{
			createPanelMesh(v, panelID);
			panelID++;
		}

		int num=panels.size();
		if(num>0) printf("\n%i %s ", num, "panels created");
		printf("\n%i %s ", fnGuideMesh.numPolygons(), "faces");
	}

	//--- COMPUTE METHODS 
	ZSPACE_INLINE void zTsPaneling::computePanelType()
	{
		printf("\n %s ", "computing panel type");
		for (int i = 0; i < panels.size(); i++)
		{
			panels[i].computeCurvature();
			printf("\n %s ", to_string(i));
			panels[i].computeType();
			panels[i].color_Type();
		}
	}

	//--- EXPORT METHODS 
	ZSPACE_INLINE void zTsPaneling::exportPanelTypes(string folderDir)
	{
		string outFileName = folderDir + "PanelTypes.csv";
		ofstream myfile;
		myfile.open(outFileName.c_str());
		string type;
		if (myfile.fail())
		{
			cout << " error in opening file  " << outFileName.c_str() << endl;
			return;
		}

		myfile << "panelID" << "," << "panelType" << endl;
		for (int i = 0; i < panels.size(); i++)
		{
			if (panels[i].panelType == zPlanar) type = "Planar";
			if (panels[i].panelType == zConcaveEllipsoid) type = "Concave Ellipsoid";
			if (panels[i].panelType == zConcaveCylinder) type = "Concave Cylinder";
			if (panels[i].panelType == zHyperboloid) type = "Hyperboloid";
			if (panels[i].panelType == zConvexCylinder) type = "Convex Cylinder";
			if (panels[i].panelType == zConvexEllipsoid) type = "Convex Ellipsoid";
			if (panels[i].panelType == zCustom)type = "Custom";

			myfile << to_string(i) << "," << type << endl;
		}
		myfile.close();
		cout << " \n file exported : " << outFileName.c_str() << endl;
	}
	
	//---- PROTECTED UTILITY METHODS
	ZSPACE_INLINE void zTsPaneling::createPanelMesh(zItMeshVertex& v, int panelID)
	{
		zPointArray positions;
		zIntArray pCounts, pConnects;

		positions.push_back(v.getPosition());

		zItMeshHalfEdgeArray cHEdges;
		v.getConnectedHalfEdges(cHEdges);

		for (auto& he : cHEdges)
		{
			zItMeshHalfEdge startHE = he.getPrev();
			zItMeshHalfEdge walkHE = he;

			pConnects.push_back(0);

			do
			{
				zPoint p = walkHE.getVertex().getPosition();
				int vID = -1;
				bool chkRepeat = coreUtils.checkRepeatVector(p, positions, vID);

				if (!chkRepeat)
				{
					vID = positions.size();
					positions.push_back(p);
				}

				pConnects.push_back(vID);

				walkHE = walkHE.getNext();

			} while (walkHE != startHE);

			pCounts.push_back(4);
		}

		// quad mesh
		zObjMesh* quad_m = panels[panelID].getRawPanelMesh();
		zFnMesh fnPanelMesh(*quad_m);
		fnPanelMesh.create(positions, pCounts, pConnects);
		//fnPanelMesh.smoothMesh(1);

		//tri mesh
		zObjMesh* tri_m = panels[panelID].getRawPanelMesh_Tri();
		zFnMesh fnPanelTriMesh(*tri_m);
		fnPanelTriMesh.create(positions, pCounts, pConnects);
		//fnPanelTriMesh.smoothMesh(1);
		fnPanelTriMesh.triangulate();
	}
}
