# Third-party notices

The Mandelbrot Explorer source code is licensed under the MIT License. Release
binaries also use the following third-party components:

| Component | Version | Linking | License |
| --- | --- | --- | --- |
| GNU MP (GMP) | 6.3.0 | Dynamic | LGPL-3.0-or-later or GPL-2.0-or-later |
| GNU MPFR | 4.2.2 | Dynamic | LGPL-3.0-or-later |
| AsmJit | 2025-10-13 | Dynamic | zlib License |
| Microsoft Visual C++ and OpenMP runtimes | MSVC redistributable version used by the build | Dynamic | Microsoft Visual Studio license |

The Windows release keeps GMP and MPFR as replaceable DLLs beside the
executable. Their complete license notices are included in the package's
`licenses` directory. The exact upstream source archives and vcpkg build
ports/patches used for these DLLs are included under `third-party-source`.
Upstream project pages are:

- GMP: <https://gmplib.org/>
- MPFR: <https://www.mpfr.org/>
- AsmJit: <https://github.com/asmjit/asmjit>

Microsoft runtime files are copied from the Visual C++ Redistributable installed
with the compiler and are redistributed under Microsoft's applicable license.
