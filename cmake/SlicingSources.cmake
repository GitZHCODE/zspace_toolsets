include_guard(GLOBAL)
# Shared slicer implementation list; platform bridges are added by callers.
function(zspace_slicing_sources output root with_harmonic)
  file(GLOB_RECURSE slicing_sources CONFIGURE_DEPENDS
    "${root}/src/zToolsets/fabrication/slicing/*.cpp")
  if(NOT with_harmonic)
    list(FILTER slicing_sources EXCLUDE REGEX "/slicing/zSlicer/harmonic/")
  endif()
  set(${output} ${slicing_sources} PARENT_SCOPE)
endfunction()
