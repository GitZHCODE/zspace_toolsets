#ifndef ZSPACE_TOOLSETS_EXPORT_H
#define ZSPACE_TOOLSETS_EXPORT_H

#pragma once

#if defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
#  if defined(_WIN32)
#    if defined(ZSPACE_TOOLSETS_EXPORTS)
#      define ZSPACE_TOOLSETS __declspec(dllexport)
#    else
#      define ZSPACE_TOOLSETS __declspec(dllimport)
#    endif
#  else
#    define ZSPACE_TOOLSETS
#  endif
#else
#  define ZSPACE_TOOLSETS
#endif

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) || defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
#  define ZSPACE_TOOLSETS_INLINE
#else
#  define ZSPACE_TOOLSETS_INLINE inline
#endif

#endif
