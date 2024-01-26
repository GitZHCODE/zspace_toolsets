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

#ifndef ZSPACE_EXT_TS_PLANE_H
#define ZSPACE_EXT_TS_PLANE_H



#pragma once

//#include "headers/zToolsets/geometry/zTsSDFSlicer.h"

#include "headers/base/zSpace_Toolsets.h"

#include <headers/zCore/base/zExtern.h>
//#include <headers/zInterface/functionsets/zFnMesh.h>
#include <headers/zInterface/functionsets/zFnGraph.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include<execution>

using namespace std;

namespace zSpace
{
	struct zExtPlane
	{
		float x_X;
		float x_Y;
		float x_Z;
		float x_R;

		float y_X;
		float y_Y;
		float y_Z;	
		float y_R;	

		float n_X;
		float n_Y;
		float n_Z;
		float n_R;
		
		float o_X;
		float o_Y;
		float o_Z;
		float o_R;

		zExtPlane(zTransform plane);

		void updateAttributes(zTransform plane);
	};

}




#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include<source/zToolsets/externalMethods/zExtPlane.cpp>
#endif

#endif