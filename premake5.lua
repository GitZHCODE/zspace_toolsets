include ("../zspace_core/premake/delay_load.lua")
include ("../zspace_core/includes.lua")

workspace "zSpace_toolsets"
    filename "zSpace_toolsets"
    architecture "x64"
    configurations {
                    --"Release",
                    "Release_DLL",
                    --"Debug",
                    "Debug_DLL",
                }

project_path = "projects"

core_path = "../ZSPACE_CORE/"

--CUSTOM FUNCTION
function prependPath(path, table)
    for key, value in pairs(table) do 
        table[key] = path .. value
    end
    return table
end

IncludeDir = {}
IncludeDir["CORE"] = "%{core_path}src/headers"
IncludeDir["SRC"] = "src/headers"

--#############__GENERAL__CONFIGURATION__SETTINGS__#############
function CommonConfigurationSettings()

--    filter "configurations:Debug"
--        kind "StaticLib"
--        objdir ("bin-int/%{cfg.buildcfg}")
--        targetdir ("bin/lib/debug/")
--        targetname ("%{prj.name}")
--        defines {"ZSPACE_STATIC_LIBRARY"}
--        optimize "Off"
--        warnings "Off"
--        --symbols "On"
--        flags {"MultiProcessorCompile"}
--        buildoptions {"/bigobj"}

    filter "configurations:Debug_DLL"
        kind "SharedLib"
        objdir ("bin-int/%{cfg.buildcfg}")
        targetdir ("bin/dll/debug")
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
--        objdir ("bin-int/%{cfg.buildcfg}")
--        targetdir ("bin/lib/")
--        targetname ("%{prj.name}")
--        defines {"ZSPACE_STATIC_LIBRARY",
--                 "NDEBUG"}
--        optimize "Full"
--        warnings "Off"
--        flags {"LinkTimeOptimization",
--                "MultiProcessorCompile"}
--        buildoptions {"/bigobj"}

    filter "configurations:Release_DLL"
        kind "SharedLib"
        objdir ("bin-int/%{cfg.buildcfg}")
        targetdir ("bin/dll/")
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

--#########################################
project "zSpace_Toolsets"
    location "%{project_path}/zSpace_Toolsets"
    language "C++"
    cppdialect "C++17"

    CommonConfigurationSettings()

    characterset("MBCS")
    
    defines {"IGL_STATIC_LIBRARY"}

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

    includedirs
    {
        "%{IncludeDir.CORE}",
        "%{IncludeDir.SRC}",
    }

    includedirs {prependPath(core_path, get_include_dirs())}

    libdirs {prependPath(core_path, get_lib_dirs())}

    links
    {
        "igl.lib",
        "zSpace_Core.lib",
        "zSpace_Interface.lib",
    }