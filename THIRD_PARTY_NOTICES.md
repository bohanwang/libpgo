# Third-party dependency notices

libpgo uses the dependencies below. The authoritative license text and
copyright notices are included in the upstream source distribution identified
by the corresponding `CMakeModules/third-party/` declaration. This table
records the exact source selected by the 0.0.4 build so that source, binary,
and wheel audits can reproduce the same license set.

| Dependency | 0.0.4 source pin | License | Build scope |
| --- | --- | --- | --- |
| Alembic | 1.8.9, SHA-256 `8835cc0cd2324510252e9e5b7412dca70305b05a43212b1416a7dfeb85219565` | BSD-3-Clause | `PGO_ENABLE_ALEMBIC` |
| ARPACK-NG | commit `8a4ede774e68d6950e7a680a3d94dca8603be80c` | BSD-3-Clause | Optional utility; not in default wheels |
| argparse | 3.1, SHA-256 `3e5a59ab7688dcd1f918bc92051a10564113d4f36c3bbed3ef596c25e519a062` | MIT | Core |
| autodiff | 1.1.2, SHA-256 `86f68aabdae1eed214bfbf0ddaa182c78ea1bb99e4df404efb7b94d30e06b744` | MIT | Core |
| backward-cpp | commit `0bfd0a07a61551413ccd2ab9a9099af3bad40681` | MIT | Core |
| Boost | 1.85.0, SHA-256 `0a9cc56ceae46986f5f4d43fe0311d90cf6d2fa9028258a95cab49ffdacf92ad` | BSL-1.0 | Python/full |
| CCCL | 3.1.3, SHA-256 `ca8cbb65bb9f0d8d9734d5f70a193f243e088cb5f4594c7d40b09fae34ae3c18` | Apache-2.0 WITH LLVM-exception | CUDA utility; not in default wheels |
| Ceres Solver | 2.2.0, SHA-256 `12efacfadbfdc1bbfa203c236e96f4d3c210bed96994288b3ff0c8e7c6f350d4` | BSD-3-Clause | Python/full |
| CGAL | 6.0.1, SHA-256 `7b0bb231d57261491722b7a0950f8026e17a08ae4a93315495cfb9d91faa31e3` | Component-specific GPL-3.0-or-later/LGPL-3.0-or-later/BSL-1.0 or commercial terms | Python/full |
| cuCollections | commit `1a1e640179e85139765a27d2d376e02628b2ccbc` | Apache-2.0 | CUDA utility; not in default wheels |
| Eigen | 3.4.0, SHA-256 `8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72` | MPL-2.0 with component-specific notices | Core/fallback |
| fmt | 11.1.4, SHA-256 `7e85cbf6125a76daa0f83cd9240eff863d988aca68cd5f66c01ff7b59fa886b6` | MIT | Core |
| Geogram | 1.9.0, SHA-256 `f2b51adf05fc8599893032c79866b4f2ff29326f810dcde649adf205896b76ab` | BSD-3-Clause | Full build |
| GLFW | 3.4, SHA-256 `b5ec004b2712fd08e8861dc271428f048775200a2df719ccf575143ba749a3e9` | Zlib | `PGO_UI` |
| GLM | 1.0.1, SHA-256 `9f3174561fd26904b23f0db5e560971cbf9b3cbda0b280f04d5c379d03bf234c` | MIT | `PGO_UI` |
| GMP/GMPXX | manylinux or Homebrew system package; tracked approved DLLs on Windows | LGPL-3.0-or-later or GPL-2.0-or-later, depending on use | CGAL runtime; bundled into repaired wheels |
| GoogleTest | 1.15.2, SHA-256 `f179ec217f9b3b3f3c6e8b02d3e7eda997b49e4ce26d6b235c9053bec9c0bf9f` | BSD-3-Clause | Tests only |
| Imath | 3.2.2, SHA-256 `b4275d83fb95521510e389b8d13af10298ed5bed1c8e13efd961d91b1105e462` | BSD-3-Clause | Statically linked into Alembic-enabled builds |
| libigl | commit `477e15a3d566a21f415aa5ee62992b12a836b01b` | MPL-2.0 | Python/full; copyleft modules are not enabled |
| libcmaes | commit `4598001f8174a0f2b90fc8983f819c8d047af26f` | Apache-2.0 or LGPL-3.0-or-later | Optional utility; not in default wheels |
| nlohmann/json | 3.11.3, SHA-256 `0d8ef5af7f9794e3263480193c491549b2ba6cc74bb018906202ada498a79406` | MIT | Core |
| NLopt | 2.10.0, SHA-256 `7bae1edca94c104a72b8d1126770c95963cfaf4530ad880295aae2ecdaa928e1` | LGPL-2.1-or-later | Optional utility; not in default wheels |
| MPFR | manylinux or Homebrew system package; tracked approved DLLs on Windows | LGPL-3.0-or-later | CGAL runtime; bundled into repaired wheels |
| oneTBB | Intel PyPI 2022.3.1 on Linux/Windows; Homebrew package on macOS arm64 | Apache-2.0 | External build input; bundled into repaired wheels |
| oneMKL | Intel PyPI `mkl-devel` 2025.3.1 on Linux/Windows | Intel Simplified Software License | Dynamic LP64/TBB runtime and CPU dispatch libraries bundled into repaired wheels |
| pybind11 | 2.12.0, SHA-256 `bf8f242abd1abcd375d516a7067490fb71abd79519a282d22b6e4d19282185a7` | BSD-3-Clause | Python |
| spdlog | 1.15.2, SHA-256 `d91ab0e16964cedb826e65ba1bed5ed4851d15c7b9453609a52056a94068c020` | MIT | Core |
| SuiteSparse | 7.10.3, SHA-256 `d4600765554133fb3c0a830ace87ff89225a01ffe02f4012ecfdcf7c778dbcc0` | Component-specific BSD/GPL/LGPL terms | Python/full |

Additional bundled sources under `third-party/` retain their upstream license
files in the source tree. The standalone release-wheel repair step copies its
native runtime closure into the wheel: GMP, GMPXX, MPFR, TBB, and, on
Linux/Windows, oneMKL with the TBB threading and CPU dispatch libraries.
Imath is built from the pinned source above and linked statically. macOS links
the system Accelerate framework instead of MKL.

The fetched upstream license files control if this summary differs from an
upstream distribution. Release packaging must include this notice file and
must re-run the wheel content and linkage audit whenever a source pin changes.

## Release licensing gate

The former unnecessary `libiglInterface` dependency on
`igl_copyleft::cgal` has been removed from the default build. The Python
extension still links `pgo_c_static`, which exposes CGAL-backed remeshing, and
SuiteSparse and CGAL use component-specific terms. The binary distribution
therefore cannot be treated as an MIT-only artifact merely because libpgo's
own source license is MIT.

Before publishing a 0.0.4 wheel, the release owner must confirm the applicable
GPL/LGPL/MPL obligations, the distribution license and corresponding-source
offer, and the complete license texts/notices that accompany the artifact.
This is a release gate, not legal advice. P4.6 must retain the final linked
component report with the wheel evidence.
