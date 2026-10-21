// This file is part of zspace, a simple C++ collection of geometry data-structures & algorithms, 
// data analysis & visualization framework.
//
// Copyright (C) 2022 ZSPACE 
// 
// This Source Code Form is subject to the terms of the MIT License 
// If a copy of the MIT License was not distributed with this file, You can 
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>; Heba Eiz <heba.eiz@zaha-hadid.com>
//


#ifdef ZSPACE_TOOLSETS
#undef ZSPACE_TOOLSETS
#endif

#if defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY) || defined(ZSPACE_SLICER_BUILD_SHARED)
#define ZSPACE_TOOLSETS  __declspec(dllexport)
#elif defined(ZSPACE_SLICER_USE_SHARED)
#define ZSPACE_TOOLSETS __declspec(dllimport)
#else
#define ZSPACE_TOOLSETS
#endif

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY) || defined(ZSPACE_SLICER_BUILD_SHARED)
#  define ZSPACE_TOOLSETS_INLINE
#else
#  define ZSPACE_TOOLSETS_INLINE inline
#endif

#define ZSPACE_TOOLSETS_EXT extern "C"
