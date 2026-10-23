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
$sources = @(Get-ChildItem "$core/src/zCore", "$core/src/zInterface", "$core/src/zIO", "$root/src/zToolsets", "$root/wasm/bridge" -Recurse -Filter '*.cpp' | Where-Object { $_.Name -notin $exclude } | ForEach-Object { $_.FullName })
$argsCommon = @('-fexceptions','-include',"$core/wasm/bridge/zspace_wasm_compat.h",'-DZSPACE_STATIC_LIBRARY','-DZSPACE_TOOLSETS_STATIC_LIBRARY',"-I$core/wasm/bridge","-I$core","-I$core/include","-I$core/third_party","-I$core/third_party/depends","-I$root/include",'-std=c++17','-O2')
# Compile objects separately so a toolset edit does not rebuild the entire kernel.
# Use compiler-generated .d dependency files so header edits only rebuild affected translation units.
$script:headerStamp = $null
function Get-FallbackHeaderStamp {
  if ($null -eq $script:headerStamp) {
    $script:headerStamp = (Get-ChildItem "$core/include", "$core/src", "$core/wasm/bridge", "$root/include", "$root/src", "$root/wasm/bridge" -Recurse -File | Where-Object { $_.Extension -in @('.h','.hpp') } | Measure-Object LastWriteTimeUtc -Maximum).Maximum
  }
  return $script:headerStamp
}

$configFile = "$root/wasm/build/config.txt"
$config = ($argsCommon -join "`n") + "`n$core`n$em`n$([IO.File]::GetLastWriteTimeUtc($PSCommandPath).Ticks)`n" + ($sources -join "`n")
$configChanged = !(Test-Path $configFile) -or ([IO.File]::ReadAllText($configFile) -ne $config)
$headerTimeCache = @{}

function Get-CachedWriteTimeUtc([string]$path) {
  $t = $headerTimeCache[$path]
  if ($null -eq $t) {
    $t = [IO.File]::GetLastWriteTimeUtc($path)
    $headerTimeCache[$path] = $t
  }
  return $t
}

function Test-DepFileDirty([string]$depFile, [datetime]$objTime) {
  if (![IO.File]::Exists($depFile)) { return ($objTime -lt (Get-FallbackHeaderStamp)) }
  $raw = [IO.File]::ReadAllText($depFile) -replace "\\`r?`n", " "
  $colonIdx = $raw.IndexOf(': ')
  if ($colonIdx -lt 0) { return ($objTime -lt (Get-FallbackHeaderStamp)) }
  $depsText = $raw.Substring($colonIdx + 2)
  $matches = [regex]::Matches($depsText, '(?:[^\s\\]|\\.)+')
  if ($matches.Count -eq 0) { return ($objTime -lt (Get-FallbackHeaderStamp)) }
  foreach ($m in $matches) {
    $depPath = $m.Value -replace '\\ ', ' '
    $depTime = Get-CachedWriteTimeUtc $depPath
    if ($depTime.Year -le 1601 -or $depTime -gt $objTime) { return $true }
  }
  return $false
}

$objects = @()
$dirtyItems = @()
$sha = [System.Security.Cryptography.SHA256]::Create()
foreach ($source in $sources) {
  $hash = [BitConverter]::ToString($sha.ComputeHash([Text.Encoding]::UTF8.GetBytes($source))).Replace('-','').Substring(0,16)
  $object = "$root/wasm/build/$hash.o"
  $depFile = "$root/wasm/build/$hash.d"
  $needsCompile = $configChanged -or ![IO.File]::Exists($object)
  if (!$needsCompile) {
    $objTime = [IO.File]::GetLastWriteTimeUtc($object)
    if ($objTime -lt [IO.File]::GetLastWriteTimeUtc($source)) {
      $needsCompile = $true
    } elseif (Test-DepFileDirty $depFile $objTime) {
      $needsCompile = $true
    }
  }
  if ($needsCompile) {
    $dirtyItems += [pscustomobject]@{ Source = $source; Object = $object; DepFile = $depFile }
  }
  $objects += $object
}

if ($dirtyItems.Count -gt 0) {
  $maxJobs = 8
  if ($env:NUMBER_OF_PROCESSORS) { $maxJobs = [Math]::Min(16, [Math]::Max(1, [int]$env:NUMBER_OF_PROCESSORS)) }
  $pool = [runspacefactory]::CreateRunspacePool(1, $maxJobs)
  $pool.Open()
  $jobs = @()
  try {
    foreach ($item in $dirtyItems) {
      Write-Host "Compiling $($item.Source)"
      $ps = [powershell]::Create().AddScript({
        param($emExe, $commonArgs, $src, $obj, $dep)
        & $emExe @commonArgs -MMD -MF $dep -c $src -o $obj 2>&1
        return [pscustomobject]@{ Source = $src; ExitCode = $LASTEXITCODE }
      }).AddArgument($em).AddArgument($argsCommon).AddArgument($item.Source).AddArgument($item.Object).AddArgument($item.DepFile)
      $ps.RunspacePool = $pool
      $jobs += [pscustomobject]@{ Pipe = $ps; Handle = $ps.BeginInvoke() }
    }
    foreach ($job in $jobs) {
      $output = $job.Pipe.EndInvoke($job.Handle)
      $job.Pipe.Dispose()
      $result = $output | Select-Object -Last 1
      if ($output.Count -gt 1) {
        $output[0..($output.Count - 2)] | ForEach-Object { Write-Host $_ }
      }
      if ($result.ExitCode -ne 0) {
        exit $result.ExitCode
      }
    }
  } finally {
    $pool.Close()
    $pool.Dispose()
  }
}

[IO.File]::WriteAllText($configFile, $config)
$outJs = "$root/wasm/out/zspace_toolsets.js"
$outWasm = "$root/wasm/out/zspace_toolsets.wasm"
# Object compilation may have succeeded before an interrupted/failed link.
$needsLink = $dirtyItems.Count -gt 0 -or $configChanged -or !(Test-Path $outJs) -or !(Test-Path $outWasm)
if (!$needsLink) {
  $outputTime = (@([IO.File]::GetLastWriteTimeUtc($outJs), [IO.File]::GetLastWriteTimeUtc($outWasm)) | Sort-Object | Select-Object -First 1)
  foreach ($object in $objects) {
    if ([IO.File]::GetLastWriteTimeUtc($object) -gt $outputTime) { $needsLink = $true; break }
  }
}
if ($needsLink) {
  # Response file avoids the Windows command-line limit at link time.
  # Eigen sparse ordering uses several temporary buffers; 64 KB is insufficient for
  # the full Rhino model. Keep overflow checks for an explicit failure on larger inputs.
  $arguments = $objects + $argsCommon + @('-sMODULARIZE=1','-sEXPORT_ES6=1','-sENVIRONMENT=web,node','-sSTACK_SIZE=2097152','-sSTACK_OVERFLOW_CHECK=1','-sALLOW_MEMORY_GROWTH=1','-sFORCE_FILESYSTEM=1',"-sEXPORTED_RUNTIME_METHODS=['HEAPF32','HEAPU32','UTF8ToString','FS','ccall']",'-o',$outJs)
  $response = "$root/wasm/build/compile.rsp"
  ($arguments | ForEach-Object { '"' + $_.Replace('\','/').Replace('"','\"') + '"' }) | Set-Content -Encoding utf8 $response
  & $em "@$response"
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
