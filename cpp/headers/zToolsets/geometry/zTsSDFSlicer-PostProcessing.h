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

#ifndef ZSPACE_TS_GEOMETRY_SDFSLICER_PP_H
#define ZSPACE_TS_GEOMETRY_SDFSLICER_PP_H



#pragma once

#include <headers/base/zSpace_Toolsets.h>

#include <headers/zCore/base/zExtern.h>

#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>
#include <headers/zInterface/functionsets/zFnParticle.h>

#include <headers/zInterface/functionsets/zFnMeshField.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#include<atlsafe.h>
//#include <comdef.h>
//#include<oaidl.h>

using namespace std;


namespace zSpace
{
	struct GCodeGenerator
	{
		vector<zObjGraph> printPath_Left;
		vector<zObjGraph> printPath_Right;

		vector<zObjGraph> trimPath_Left;
		vector<zObjGraph> trimPath_Right;

	};
	
	


}





namespace zSpace
{


}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/geometry/zTsSDFSlicer.cpp>
#endif

#endif