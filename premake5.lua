require('vstudio')

include ("../zspace_core/includes.lua")

workspace "zSpace_toolsets"
    filename "zSpace_toolsets"
    architecture "x64"
    configurations {"Debug", "Debug_DLL", "Release", "Release_DLL", "Release_DLL_OV", "Release_Make", "Release_Unreal"}
    startproject "zSpace_Core"

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
IncludeDir["IGL"] = "Dependencies/igl/headers"
IncludeDir["OMNI"] = "Dependencies/omniverse"
IncludeDir["CORE"] = "%{core_path}src/headers"
IncludeDir["SRC"] = "src/headers"

--core_includes = prependPath(core_path, get_include_dirs());

LibDir = {}
LibDir["IGL"] = "Dependencies/igl/build/lib"
LibDir["OMNI"] = "Dependencies/omniverse"
LibDir["CORE"] = "../zspace_core/bin/dll"

--#############__GENERAL__CONFIGURATION__SETTINGS__#############
function CommonConfigurationSettings()
    filter "configurations:Debug"
        kind "StaticLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/lib/debug")
        targetname ("%{prj.name}")
        defines {"ZSPACE_STATIC_LIBRARY",
                "USING_ARMA"}
        symbols "On"

    filter "configurations:Debug_DLL"
        kind "SharedLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/dll/")
        targetname ("%{prj.name}")
        symbols "On"

    filter "configurations:Release"
        kind "StaticLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/lib/")
        targetname ("%{prj.name}")
        defines {"ZSPACE_STATIC_LIBRARY",
                "USING_ARMA"}
        optimize "Full"
        warnings "Off"
        flags {"LinkTimeOptimization"}

    filter "configurations:Release_DLL"
        kind "SharedLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/dll/")
        targetname ("%{prj.name}")
        defines {"ZSPACE_DYNAMIC_LIBRARY", "ZSPACE_TOOLSETS_DYNAMIC_LIBRARY"}
        optimize "Speed"
        warnings "Off"
        flags {"LinkTimeOptimization"}

    filter "configurations:Release_DLL_OV"
        kind "SharedLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/dll/")
        targetname ("%{prj.name}")
        defines {"ZSPACE_DYNAMIC_LIBRARY", "ZSPACE_TOOLSETS_DYNAMIC_LIBRARY"}
        optimize "Speed"
        warnings "Off"
        flags {"LinkTimeOptimization"}
    
    filter "configurations:Release_Unreal"
        kind "StaticLib"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/lib/")
        targetname ("%{prj.name}")
        defines {"ZSPACE_STATIC_LIBRARY"}
        optimize "Speed"
        warnings "Off"
        flags {"LinkTimeOptimization"}
    
    filter {}
end

--#########################################
project "zSpace_Toolsets"
    location "%{project_path}/zSpace_Toolsets"
    language "C++"
    cppdialect "C++17"

    CommonConfigurationSettings()

    characterset("MBCS")
    
    defines {"IGL_STATIC_LIBRARY", "_WINDLL"}

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
        "%{IncludeDir.IGL}",
        "%{IncludeDir.OMNI}",
        "%{IncludeDir.CORE}",
        "%{IncludeDir.SRC}",
    }

    includedirs {prependPath(core_path, get_include_dirs())}

    libdirs
    {
        "%{LibDir.IGL}",
        "%{LibDir.CORE}",
    }

    libdirs {prependPath(core_path, get_lib_dirs())}

    links
    {
        "igl.lib",
        "zSpace_Core.lib",
        "zSpace_Interface.lib",
    }

    filter "configurations:Release_Unreal"
        defines {"ZSPACE_UNREAL_INTEROP"}

    filter "configurations:Release_Make"
        kind "Makefile"
        objdir ("bin-int/%{cfg.architecture}/%{cfg.buildcfg}")
        targetdir ("bin/make")
        targetname ("%{prj.name}")

--#########################################