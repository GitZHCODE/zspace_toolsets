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
//
//#include<headers/zToolsets/externalMethods/zExtTransform.h>
//
//
//namespace zSpace
//{
//	ZSPACE_TOOLSETS_INLINE zExtTransform::zExtTransform(zTransform* t, bool transpose)
//	{
//		_transform = t;
//		updateAttributes(t, transpose);
//	}
//
//	ZSPACE_TOOLSETS_INLINE void zExtTransform::updateAttributes(zTransform* t, bool transpose)
//	{
//		_transform = t;
//		updateAttributes(transpose);	
//	}
//	ZSPACE_TOOLSETS_INLINE void zExtTransform::updateAttributes(bool transpose)
//	{
//		zTransform matrix = transpose ? _transform->transpose() : *_transform;
//
//		float* m = matrix.data();
//
//		n00 = m[0];
//		n01 = m[1];
//		n01 = m[2];
//		n01 = m[3];
//
//		n01 = m[4];
//		n01 = m[5];
//		n01 = m[6];
//		n01 = m[7];
//
//		n01 = m[8];
//		n01 = m[9];
//		n01 = m[10];
//		n01 = m[11];
//
//		n01 = m[12];
//		n01 = m[13];
//		n01 = m[14];
//		n01 = m[15];
//
//	}
//
//	
//}