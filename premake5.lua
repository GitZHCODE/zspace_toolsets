--#############__GENERAL__CONFIGURATION__SETTINGS__#############
local function CommonConfigurationSettings()
    GlobalCommonDefines()

    defines {"IGL_STATIC_LIBRARY"}

--    filter "configurations:Debug"
--        kind "StaticLib"
--        targetname ("%{prj.name}")
--        defines {"ZSPACE_STATIC_LIBRARY"}
--        optimize "Off"
--        warnings "Off"
--        --symbols "On"
--        flags {"MultiProcessorCompile"}
--        buildoptions {"/bigobj"}

    filter "configurations:Debug_DLL*"
        kind "SharedLib"
        targetname ("%{prj.name}")
        defines {"ZSPACE_TOOLSETS_DYNAMIC_LIBRARY",
                 "ZSPACE_DYNAMIC_LIBRARY",
                 "_WINDLL"}
        optimize "Off"
        warnings "Off"
        --symbols "On"
        flags {"MultiProcessorCompile"}
        buildoptions {"/bigobj"}

--    filter "configurations:Release"
--        kind "StaticLib"
--        targetname ("%{prj.name}")
--        defines {"ZSPACE_STATIC_LIBRARY",
--                 "NDEBUG"}
--        optimize "Full"
--        warnings "Off"
--        flags {"LinkTimeOptimization",
--                "MultiProcessorCompile"}
--        buildoptions {"/bigobj"}

    filter "configurations:Release_DLL*"
        kind "SharedLib"
        targetname ("%{prj.name}")
        defines {"ZSPACE_TOOLSETS_DYNAMIC_LIBRARY",
                 "ZSPACE_DYNAMIC_LIBRARY",
                 "NDEBUG",
                 "_WINDLL"}
        optimize "Full"
        warnings "Off"
        flags {"LinkTimeOptimization",
                "MultiProcessorCompile"}

    filter {}
end

ToolsetsIncludeDir = {}
ToolsetsIncludeDir["CORE"] = "%{core_path}src/headers"

--#############__ZSPACE_TOOLSETS__#############
project "zSpace_Toolsets"
    location "projects/zSpace_Toolsets"
    language "C++"
    cppdialect "C++17"

    dependson("zSpace_Interface")

    CommonConfigurationSettings()

    defines{
        "_HAS_STD_BYTE=0",
        "NOMINMAX",
    }

    characterset("MBCS")

    pchheader "zToolsets/ztoolsetspch.h"
    pchsource "src/source/zToolsets/ztoolsetspch.cpp"
    rawforceincludes "zToolsets/ztoolsetspch.h"

    files
    {
        "src/headers/**.h",
        "src/source/**.cpp",
    }

    removefiles {"**/externalMethods/**"}
    removefiles {"**Copy*.cpp"}
    removefiles {"**backup*.cpp"}

    --Exclude From Build
    filter {"files:**SDFSlicer*.* or **KMeans*.* or **OSM*.* or **RobotFab*.* or **Polyhedra*.* or **Mesh2Pix*.* or **Paneling*.* or **Remesh*.* or **SDFBridge*.* or **Spectral*.* or **VariableExtrude*.* or **/pathNetworks/** or **Polytopal*.* or **Spatial*.* or **TopOpt*.*"}
        flags {"ExcludeFromBuild"}
    filter {}

    --###__BASE__###
    includedirs
    {
        "%{ToolsetsIncludeDir.CORE}",
        "src/headers",
    }

    includedirs {prependPath(deps_path, get_include_dirs())}

    libdirs { "%{core_path}bin/%{cfg.buildcfg}"}

    libdirs {prependPath(deps_path, get_lib_dirs())}

    links
    {
        "igl.lib",
        "zSpace_Core.lib",
        "zSpace_Interface.lib",
    }

    --Omniverse is default
    links {get_omniverse_links()}
