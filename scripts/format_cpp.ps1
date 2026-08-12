[CmdletBinding()]
param(
    [switch]$Check
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot

function Find-ClangFormat {
    $command = Get-Command clang-format -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    $userBase = & python -m site --user-base
    if ($LASTEXITCODE -eq 0 -and $userBase) {
        $candidate = Get-ChildItem $userBase -Filter clang-format.exe -File -Recurse -ErrorAction SilentlyContinue |
            Select-Object -First 1 -ExpandProperty FullName
        if ($candidate) {
            return $candidate
        }
    }

    throw 'clang-format 22.x was not found. Install it with: python -m pip install --user clang-format==22.1.8'
}

$formatter = Find-ClangFormat
$version = & $formatter --version
if ($LASTEXITCODE -ne 0 -or $version -notmatch 'clang-format version 22\.') {
    throw "clang-format 22.x is required; found: $version"
}

$files = & git -C $root ls-files -- '*.cpp' '*.h'
if ($LASTEXITCODE -ne 0) {
    throw 'Unable to list tracked C++ files.'
}

foreach ($relativePath in $files) {
    $path = Join-Path $root ($relativePath -replace '/', '\')
    if ($Check) {
        & $formatter --dry-run --Werror --style=file $path
    } else {
        & $formatter -i --style=file $path
    }
    if ($LASTEXITCODE -ne 0) {
        throw "clang-format failed for $relativePath"
    }
}

$action = if ($Check) { 'checked' } else { 'formatted' }
Write-Host "clang-format $action $($files.Count) tracked C++ files."
