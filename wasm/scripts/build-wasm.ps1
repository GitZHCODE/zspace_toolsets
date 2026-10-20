param([string]$CoreDir = "$PSScriptRoot/../../../zspace_core", [string]$EmsdkDir = "$env:USERPROFILE/emsdk")
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path "$PSScriptRoot/../..").Path
$core = (Resolve-Path $CoreDir).Path
$em = Join-Path $EmsdkDir 'upstream/emscripten/em++.exe'
if (!(Test-Path $em)) { throw "Emscripten not found at $em" }
$env:EM_CONFIG = Join-Path $EmsdkDir '.emscripten'
$env:EM_CACHE = "$root/wasm/cache"
$env:TEMP = "$root/wasm/tmp"
$env:TMP = $env:TEMP
$env:EMCC_TEMP_DIR = $env:TEMP
$env:EMSDK_PYTHON = Join-Path $EmsdkDir 'python/3.13.3_64bit/python.exe'
New-Item -ItemType Directory -Force "$root/wasm/out", "$root/wasm/build", $env:TEMP, $env:EM_CACHE | Out-Null
$exclude = @('zFnComputeMesh.cpp','zObjectComputeMesh.cpp','zFnMeshField.cpp','zFnPointField.cpp','zObjectMeshField.cpp','zObjectPointField.cpp','zObjectComputeField2D.cpp')
# Tagged slicing is native-only until its bridge and field/TetGen dependencies
# are explicitly integrated. Keep the existing browser runtime source set stable.
$sources = @(Get-ChildItem "$core/src/zCore", "$core/src/zInterface", "$core/src/zIO", "$root/src/zToolsets", "$root/wasm/bridge" -Recurse -Filter '*.cpp' | Where-Object { $_.Name -notin $exclude -and $_.FullName -notmatch '[\\/]fabrication[\\/]slicing[\\/]' } | ForEach-Object { $_.FullName })
$argsCommon = @('-fexceptions','-include',"$core/wasm/bridge/zspace_wasm_compat.h",'-DZSPACE_STATIC_LIBRARY','-DZSPACE_TOOLSETS_STATIC_LIBRARY',"-I$core/wasm/bridge","-I$core","-I$core/include","-I$core/third_party","-I$core/third_party/depends","-I$root/include",'-std=c++17','-O2')
# Compile objects separately so a toolset edit does not rebuild the entire kernel.
# A header or build-script change invalidates the cache conservatively.
$headerStamp = (Get-ChildItem "$core/include", "$core/src", "$core/wasm/bridge", "$root/include", "$root/src", "$root/wasm/bridge" -Recurse -File | Where-Object { $_.Extension -in @('.h','.hpp') } | Measure-Object LastWriteTimeUtc -Maximum).Maximum
$configFile = "$root/wasm/build/config.txt"
$config = ($argsCommon -join "`n") + "`n$core`n$em`n$((Get-Item $PSCommandPath).LastWriteTimeUtc.Ticks)"
$configChanged = !(Test-Path $configFile) -or (Get-Content $configFile -Raw) -ne $config
$objects = @()
$sha = [System.Security.Cryptography.SHA256]::Create()
foreach ($source in $sources) {
  $hash = [BitConverter]::ToString($sha.ComputeHash([Text.Encoding]::UTF8.GetBytes($source))).Replace('-','').Substring(0,16)
  $object = "$root/wasm/build/$hash.o"
  if ($configChanged -or !(Test-Path $object) -or (Get-Item $object).LastWriteTimeUtc -lt (Get-Item $source).LastWriteTimeUtc -or (Get-Item $object).LastWriteTimeUtc -lt $headerStamp) {
    Write-Host "Compiling $source"
    & $em @argsCommon -c $source -o $object
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
  }
  $objects += $object
}
[IO.File]::WriteAllText($configFile, $config)
# Response file avoids the Windows command-line limit at link time.
# Eigen sparse ordering uses several temporary buffers; 64 KB is insufficient for
# the full Rhino model. Keep overflow checks for an explicit failure on larger inputs.
$arguments = $objects + $argsCommon + @('-sMODULARIZE=1','-sEXPORT_ES6=1','-sENVIRONMENT=web,node','-sSTACK_SIZE=2097152','-sSTACK_OVERFLOW_CHECK=1','-sALLOW_MEMORY_GROWTH=1','-sFORCE_FILESYSTEM=1',"-sEXPORTED_RUNTIME_METHODS=['HEAPF32','HEAPU32','UTF8ToString','FS','ccall']",'-o',"$root/wasm/out/zspace_toolsets.js")
$response = "$root/wasm/build/compile.rsp"
($arguments | ForEach-Object { '"' + $_.Replace('\','/').Replace('"','\"') + '"' }) | Set-Content -Encoding utf8 $response
& $em "@$response"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
