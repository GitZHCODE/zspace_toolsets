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

    filter "configurations:*Debug*"
        targetname ("%{prj.name}")
        optimize "Off"
        warnings "Off"
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

    filter "configurations:*Release*"
        targetname ("%{prj.name}")
        defines {"NDEBUG"}
        optimize "Full"
        warnings "Off"
        flags {"LinkTimeOptimization",
                "MultiProcessorCompile"}

    filter "configurations:*DLL*"
        kind "SharedLib"
        defines {"ZSPACE_TOOLSETS_DYNAMIC_LIBRARY",
                 "ZSPACE_DYNAMIC_LIBRARY",
                 "_WINDLL"}

    filter {}
end

path_from_toolsets_to_zspace_deps = path.getrelative(zspace_toolsets_path, zspace_deps_path)
path_from_toolsets_to_zspace_core = path.getrelative(zspace_toolsets_path, zspace_core_path)

--#############__ZSPACE_TOOLSETS__#############
project "zSpace_Toolsets"
    location "projects/zSpace_Toolsets/%{wks.name}"
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
        --Manualy include alglib files
        "%{prependPath(path_from_toolsets_to_zspace_deps, get_zspace_include_dirs()).ALGLIB}/alglib/**.cpp"
    }

    removefiles {"**/externalMethods/**"}
    removefiles {"**Copy*.cpp"}
    removefiles {"**backup*.cpp"}

    --Exclude From Build
    filter {"files:**SDFSlicer*.* or **KMeans*.* or **OSM*.* or **RobotFab*.* or **Polyhedra*.* or **Mesh2Pix*.* or **Paneling*.* or **Remesh*.* or **SDFBridge*.* or **Spectral*.* or **VariableExtrude*.* or **/pathNetworks/** or **Spatial*.* or **TopOpt*.*"}
        flags {"ExcludeFromBuild"}
    filter {}

    --###__BASE__###
    includedirs
    {
        "%{path_from_toolsets_to_zspace_core}/src/headers",
        "src/headers",
    }

    -- Add Core include directories
    includedirs {prependPath(path_from_toolsets_to_zspace_deps , get_zspace_include_dirs())}

    -- Add Core lib directories
    libdirs {prependPath(path_from_toolsets_to_zspace_deps , get_zspace_lib_dirs())}

    -- Core build libdirs
    libdirs { "%{path_from_toolsets_to_zspace_core}/bin/%{cfg.buildcfg}/%{cfg.platform}"}

    links
    {
        "igl.lib",
        "zSpace_Core.lib",
        "zSpace_Interface.lib",
    }
    --##############

    --###__Injected__###
        -- Inject included dirs
        includedirs {get_include_dirs_injection()}

        -- Inject lib dirs
        libdirs {get_lib_dirs_injection()}
    --##################

    --###__Omniverse__###
    filter {"platforms:Omniverse or Both"}
        -- Add omniverse includes
        includedirs {prependPath(path_from_toolsets_to_zspace_deps , get_omniverse_includes())}

        -- Add omniverse libdirs
        libdirs {prependPath(path_from_toolsets_to_zspace_deps , get_omniverse_libdirs())}

        -- Add omniverse links
        links {get_omniverse_links()}
    --###################
