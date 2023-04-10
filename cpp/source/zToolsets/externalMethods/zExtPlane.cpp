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
//#include<headers/zToolsets/externalMethods/zExtPlane.h>
//
//
//namespace zSpace
//{
//	ZSPACE_TOOLSETS_INLINE zExtPlane::zExtPlane(zTransform plane)
//	{
//		zTransform frameTranspose = plane.transpose();
//		float* m = frameTranspose.data();
//
//		x_X = m[0];
//		x_Y = m[1];
//		x_Z = m[2];
//		x_R = m[3];
//
//		y_X = m[4];
//		y_Y = m[5];
//		y_Z = m[6];
//		y_R = m[7];
//
//		n_X = m[8];
//		n_Y = m[9];
//		n_Z = m[10];
//		n_R = m[11];
//
//		o_X = m[12];
//		o_Y = m[13];
//		o_Z = m[14];
//		o_R = m[15];
//	}
//
//	ZSPACE_TOOLSETS_INLINE void zExtPlane::updateAttributes(zTransform plane)
//	{
//			zTransform frameTranspose = plane.transpose();
//			float* m = frameTranspose.data();
//
//			x_X = m[0];
//			x_Y = m[1];
//			x_Z = m[2];
//			x_R = m[3];
//
//			y_X = m[4];
//			y_Y = m[5];
//			y_Z = m[6];
//			y_R = m[7];
//
//			n_X = m[8];
//			n_Y = m[9];
//			n_Z = m[10];
//			n_R = m[11];
//
//			o_X = m[12];
//			o_Y = m[13];
//			o_Z = m[14];
//			o_R = m[15];
//			
//	}
//	
//	
//}