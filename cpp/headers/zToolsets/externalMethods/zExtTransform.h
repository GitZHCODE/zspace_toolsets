//// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
//// data analysis & visualization framework.
////
//// Copyright (C) 2019 ZSPACE 
//// 
//// This Source Code Form is subject to the terms of the MIT License 
//// If a copy of the MIT License was not distributed with this file, You can 
//// obtain one at https://opensource.org/licenses/MIT.
////
//// Author : Heba Eiz <heba.eiz@zaha-hadid.com>
////
//
//#ifndef ZSPACE_EXT_TS_TRANSFORM_H
//#define ZSPACE_EXT_TS_TRANSFORM_H
//
//
//
//#pragma once
//
////#include <headers/zToolsets/geometry/zTsSDFSlicer.h>
//
//#include <headers/base/zSpace_Toolsets.h>
//
//#include <headers/zCore/base/zExtern.h>
////#include <headers/zInterface/functionsets/zFnMesh.h>
//#include <headers/zInterface/functionsets/zFnGraph.h>
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <iostream>
//#include <sstream>
//
//#include<execution>
//
//using namespace std;
//
//namespace zSpace
//{
//	struct zExtTransform
//	{
//		zTransform* _transform;
//		float n00;
//		float n01;
//		float n02;
//		float n03;
//
//		float n04;
//		float n05;
//		float n06;	
//		float n07;	
//
//		float n08;
//		float n09;
//		float n10;
//		float n11;
//		
//		float n12;
//		float n13;
//		float n14;
//		float n15;
//
//		zExtTransform(zTransform* t, bool transpose = false);
//
//		void updateAttributes(zTransform* transform, bool transpose = false);
//		void updateAttributes(bool transpose = false);
//	};
//
//}
//
//
//
//
//#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY)  || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
//// All defined OK so do nothing
//#else
//#include<source/zToolsets/externalMethods/zExtTransform.cpp>
//#endif
//
//#endif