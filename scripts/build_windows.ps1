param(
    [string]$BuildDir = "build\windows-cmake"
)

if (-not $env:MAX_SDK_ROOT) {
    Write-Error "MAX_SDK_ROOT environment variable must point to the Max SDK"
    exit 1
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = (Resolve-Path (Join-Path $scriptDir '..')).Path
$buildPath = Join-Path $root $BuildDir

if (-not (Test-Path $buildPath)) {
    New-Item -ItemType Directory -Path $buildPath | Out-Null
}

cmake -S $root -B $buildPath -DCMAKE_BUILD_TYPE=Release
cmake --build $buildPath --config Release

Write-Host "Built kbeyond~ into $buildPath"
