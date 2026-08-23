[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^\d+\.\d+\.\d+$')]
    [string]$Version
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$versionHeader = Get-Content (Join-Path $root 'src\gui\version.h') -Raw
$declaredVersion = [regex]::Match($versionHeader, 'MANDEL_VERSION_STRING\s+"([^"]+)"').Groups[1].Value
if ($declaredVersion -ne $Version) {
    throw "Requested version $Version does not match version.h ($declaredVersion)."
}

& $env:ComSpec /c "`"$PSScriptRoot\build_gui_release.bat`""
if ($LASTEXITCODE -ne 0) {
    throw 'Release build failed.'
}

$releaseBuild = Join-Path $root 'build\release'
$exe = Join-Path $releaseBuild 'mandel_gui.exe'
$versionInfo = (Get-Item $exe).VersionInfo
if ($versionInfo.ProductVersion -ne $Version -or $versionInfo.FileVersion -ne "$Version.0") {
    throw "Executable version metadata is invalid: product=$($versionInfo.ProductVersion), file=$($versionInfo.FileVersion)."
}

$vcpkgRoots = @(
    $env:VCPKG_ROOT,
    (Join-Path $env:USERPROFILE 'repo\vcpkg'),
    (Join-Path $env:USERPROFILE 'vcpkg'),
    'C:\vcpkg'
) | Where-Object { $_ }
$vcpkgRoot = $vcpkgRoots | Where-Object { Test-Path (Join-Path $_ 'installed\x64-windows\share\mpfr\copyright') } | Select-Object -First 1
if (-not $vcpkgRoot) {
    throw 'Unable to locate dynamic vcpkg GMP/MPFR license files.'
}

$packageRoot = Join-Path $root 'build\package'
$packageName = "Mandelbrot-Explorer-v$Version-win-x64"
$stage = Join-Path $packageRoot $packageName
if (Test-Path $stage) {
    Remove-Item -LiteralPath $stage -Recurse -Force
}
New-Item -ItemType Directory -Path (Join-Path $stage 'licenses') -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $stage 'third-party-source') -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $stage 'third-party-source\vcpkg-ports') -Force | Out-Null

Copy-Item $exe $stage
Copy-Item (Join-Path $releaseBuild '*.dll') $stage
Copy-Item (Join-Path $root 'README.md') $stage
Copy-Item (Join-Path $root 'LICENSE') $stage
Copy-Item (Join-Path $root 'THIRD_PARTY_NOTICES.md') $stage
Copy-Item (Join-Path $root "docs\releases\v$Version.md") (Join-Path $stage 'RELEASE_NOTES.md')
Copy-Item (Join-Path $vcpkgRoot 'installed\x64-windows\share\gmp\copyright') (Join-Path $stage 'licenses\GMP.txt')
Copy-Item (Join-Path $vcpkgRoot 'installed\x64-windows\share\mpfr\copyright') (Join-Path $stage 'licenses\MPFR.txt')
Copy-Item (Join-Path $vcpkgRoot 'installed\x64-windows\share\asmjit\copyright') (Join-Path $stage 'licenses\AsmJit.txt')

function Copy-CorrespondingSource {
    param(
        [string]$Name,
        [string]$Url,
        [string]$Sha512
    )
    $archive = Join-Path $vcpkgRoot "downloads\$Name"
    if (-not (Test-Path $archive)) {
        Invoke-WebRequest -Uri $Url -OutFile $archive
    }
    $actual = (Get-FileHash $archive -Algorithm SHA512).Hash.ToLowerInvariant()
    if ($actual -ne $Sha512.ToLowerInvariant()) {
        throw "Source archive hash mismatch for $Name."
    }
    Copy-Item $archive (Join-Path $stage 'third-party-source')
}

Copy-CorrespondingSource 'gmp-6.3.0.tar.xz' 'https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz' 'e85a0dab5195889948a3462189f0e0598d331d3457612e2d3350799dba2e244316d256f8161df5219538eb003e4b5343f989aaa00f96321559063ed8c8f29fd2'
Copy-CorrespondingSource 'mpfr-4.2.2.tar.xz' 'https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.2.tar.xz' 'eb9e7f51b5385fb349cc4fba3a45ffdf0dd53be6dfc74932dc01258158a10514667960c530c47dd9dfc5aa18be2bd94859d80499844c5713710581e6ac6259a9'
Copy-Item (Join-Path $vcpkgRoot 'ports\gmp') (Join-Path $stage 'third-party-source\vcpkg-ports') -Recurse
Copy-Item (Join-Path $vcpkgRoot 'ports\mpfr') (Join-Path $stage 'third-party-source\vcpkg-ports') -Recurse
$vcpkgCommit = git -C $vcpkgRoot rev-parse HEAD
if ($LASTEXITCODE -ne 0 -or -not $vcpkgCommit) {
    throw 'Unable to record the vcpkg source revision.'
}
$vcpkgCommit | Set-Content (Join-Path $stage 'third-party-source\vcpkg-commit.txt') -Encoding ascii

$manifestLines = Get-ChildItem -Path $stage -File -Recurse |
    Sort-Object FullName |
    ForEach-Object {
        $relative = $_.FullName.Substring($stage.Length).TrimStart('\').Replace('\', '/')
        "$(Get-FileHash $_.FullName -Algorithm SHA256 | Select-Object -ExpandProperty Hash)  $relative"
    }
$manifestLines | Set-Content (Join-Path $stage 'RELEASE-MANIFEST.sha256') -Encoding ascii

$archiveName = "mandelbrot-explorer-v$Version-win-x64.zip"
$archive = Join-Path $root "build\$archiveName"
$checksum = "$archive.sha256"
Remove-Item -LiteralPath $archive -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath $checksum -Force -ErrorAction SilentlyContinue
Compress-Archive -Path $stage -DestinationPath $archive -CompressionLevel Optimal
$archiveHash = Get-FileHash $archive -Algorithm SHA256
"$($archiveHash.Hash)  $([IO.Path]::GetFileName($archive))" | Set-Content $checksum -Encoding ascii

Write-Host "Package: $archive"
Write-Host "SHA-256: $($archiveHash.Hash)"
