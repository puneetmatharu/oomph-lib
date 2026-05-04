# AGENTS.md

Guidance for automated coding agents working in this repository. This file
applies to the whole tree.

## Project Overview

`oomph-lib` is a C++17 finite-element library for multi-physics problems. The
root project builds the libraries in `src/`, installs headers/libraries into a
local `install/` directory by default, and treats the many `demo_drivers/`
programs as both examples and validation tests. Documentation is generated from
Doxygen-style text and source comments under `doc/`.

Important entry points:

- `README.md`: canonical build, install, test, demo-driver, documentation, and
  developer guidance.
- `CONTRIBUTING.md`: GitHub workflow, formatting, and pre-commit guidance.
- `doc/coding_conventions/coding_conventions.txt`: detailed C++ style and
  design conventions.
- `.clang-format` and `.cmake-format.py`: formatting specifications.
- `CMakeLists.txt`, `CMakePresets.json`, and `cmake/`: build system machinery.
- `oomph_build.py`: wrapper for third-party libraries, root build, install, and
  optional docs.

## Repository Layout

- `src/`: library code, grouped by functionality. Each subdirectory usually has
  headers, sources, and a local `CMakeLists.txt` that calls
  `oomph_library_config(...)`.
- `demo_drivers/`: example and validation drivers. Most leaf directories contain
  one or more `.cc` files, `CMakeLists.txt`, `validate.sh`, and `validata/`.
- `doc/`: Doxygen documentation. Its structure broadly mirrors
  `demo_drivers/`.
- `cmake/`: project-specific CMake modules. Functions with the `oomph_` prefix
  are local project APIs, not built-in CMake commands.
- `external_distributions/`: CMake machinery for fetching/building third-party
  libraries.
- `external_src/`: bundled third-party or compatibility source. Treat this as
  vendor code unless a task explicitly targets it.
- `scripts/`: project helper scripts, including validation and test-filtering
  utilities.
- `user_drivers/`: area for user driver examples outside `demo_drivers/`.

## Build And Test

Prefer Ninja for all CMake builds.

Prerequisites include CMake 3.24+, Ninja, Python 3, C/C++/Fortran compilers,
and OpenBLAS. MPI, Doxygen, LaTeX, and additional third-party libraries are
needed for some configurations.

Useful commands:

```bash
# Full helper-script build and local install. This can take a long time.
python3 oomph_build.py

# Reuse existing helper-script directories during iterative work.
python3 oomph_build.py --reuse-tpl --reuse-oomph

# Manual root configure/build/install, after third-party libraries exist.
cmake -G Ninja -B build $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt)
cmake --build build
cmake --install build

# Developer-oriented debug build.
cmake -G Ninja -B build-debug \
  $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt) \
  -DCMAKE_BUILD_TYPE=Debug \
  -DOOMPH_ENABLE_PARANOID=ON \
  -DOOMPH_ENABLE_RANGE_CHECKING=ON
cmake --build build-debug
cmake --install build-debug
```

`PARANOID` and `RANGE_CHECKING` are development checks. Do not define these
macros directly in driver code; configure the library with
`OOMPH_ENABLE_PARANOID` and `OOMPH_ENABLE_RANGE_CHECKING` instead.

Run validation through demo-driver CTest builds after the library is installed:

```bash
# All demo-driver tests. Expensive.
cd demo_drivers
cmake -G Ninja -B build
cd build
ctest -j4

# Re-run failed tests verbosely.
ctest --rerun-failed --output-on-failure

# A targeted demo-driver test.
cd demo_drivers/poisson/one_d_poisson
cmake -G Ninja -B build
cd build
ctest --output-on-failure

# Filter tests by TEST_NAME from a top-level demo_drivers build.
ctest -R 'poisson.one_d_poisson' --output-on-failure

# Generate a ctest command for drivers containing a keyword.
cd demo_drivers
../scripts/filter_cmake_tests.py --root . --keywords 'QAdvectionDiffusionElement'
```

If `oomph-lib` is installed outside the default `install/` directory or
`/usr/local`, pass `-Doomphlib_ROOT=/absolute/install/prefix` when configuring
demo drivers, or set the `oomphlib_ROOT` environment variable.

Build documentation only when relevant; it is slow and requires Doxygen and
LaTeX:

```bash
cd doc
cmake -G Ninja -B build
cmake --build build
```

## Formatting

- C++ formatting is governed by `.clang-format`.
- CMake formatting is governed by `.cmake-format.py`.
- The project uses 2-space indentation and 80-column formatting.
- Braces follow Allman style.
- Do not sort includes as a drive-by change; `.clang-format` has
  `SortIncludes: false`.
- Avoid blanket-formatting `demo_drivers/`; it can disturb documentation
  examples.
- Avoid blanket-formatting `src/generic/`; some C-like compatibility headers are
  deliberately excluded.
- Treat `external_src/` as vendor code and preserve its local style unless the
  task is specifically about that code.

Manual formatting examples:

```bash
clang-format -i --style=file path/to/file.cc path/to/file.h
cmake-format --config-files .cmake-format.py --in-place path/to/CMakeLists.txt
```

Optional pre-commit hooks are configured in `.pre-commit-config.yaml`. The CI
formatting workflow also runs clang-format over `src/` C++ files, excluding
known C/compatibility files, and cmake-format over CMake files.

## C++ Style And Design

Follow `doc/coding_conventions/coding_conventions.txt` for detailed rules.
Key conventions:

- Source files use `.h` and `.cc`. Use `.template.cc` for non-inline member
  functions of templated classes that must be included by users.
- Put includes at the start of files.
- All library additions belong in namespace `oomph`.
- Do not write `using namespace std;` globally in headers. It is acceptable in
  driver codes, `.cc` files, or narrow function scopes.
- Class names are `UpperCamelCase` with no underscores.
- General variables are lowercase with underscores.
- Private data starts with a capital letter, e.g. `Boundary_node_pt`.
- Public accessors use the lowercase equivalent, e.g. `boundary_node_pt()`.
- Pointer variables and pointer accessors end in `_pt`.
- Accessors for container entries are singular, e.g. `element_pt(e)`.
- Count accessors use an `n` prefix without an underscore, e.g. `nelement()`.
- Template parameters are uppercase, e.g. `DIM`. Expose them, when needed, as
  `static const` members named `TEMPLATE_PARAMETER_DIM`.
- Prefer descriptive names even when they are long.
- Use `unsigned` for counters and sizes unless negative values are meaningful.
- Pass function arguments by reference; use `const` references for input
  arguments.
- Return output arguments at the end of argument lists. Time-level arguments come
  first.
- Use `const` wherever applicable and provide const/non-const overloads for
  subscript-like access.
- Avoid direct public data. Keep data private or protected and provide accessors
  or explicit state-changing methods.
- Boolean state should use `enable_*`/`disable_*` mutators and `is_*` or `has_*`
  queries rather than reference-returning boolean accessors.
- Every class should either delete copy construction/assignment or implement
  them correctly.
- Avoid macros except established project macros such as
  `OOMPH_EXCEPTION_LOCATION` and `OOMPH_CURRENT_FUNCTION`.
- Inline simple get/set functions in headers, but avoid inline functions that
  introduce circular dependencies.
- Precede functions/classes with Doxygen `///` comments bounded by `// =====`
  separator lines, matching surrounding code.
- Put member data declarations on separate documented lines so Doxygen handles
  them correctly.

When adding diagnostics in development checks, put sanity checks under
`#ifdef PARANOID` where appropriate and throw `OomphLibError` with
`OOMPH_CURRENT_FUNCTION` and `OOMPH_EXCEPTION_LOCATION`. Use `OomphLibWarning`
for warnings that should follow project formatting.

## CMake Conventions

- Use project helper functions rather than open-coded CMake when an existing
  helper fits.
- In `src/<library>/CMakeLists.txt`, define `HEADERS`, `SOURCES`, `LIBNAME`,
  and `LINKLIBS`, then call `oomph_library_config(...)`.
- If adding files to an existing library, update that library's
  `CMakeLists.txt`.
- If adding a new `src/` library directory, add it to the parent
  `src/CMakeLists.txt` `SUBDIRS` list.
- In demo-driver directories, use `oomph_add_executable(...)` and
  `oomph_add_test(...)`.
- Demo-driver `TEST_NAME` values should mirror the path below `demo_drivers/`
  with dots, for example `poisson.one_d_poisson`.
- If adding a new demo-driver directory, update the nearest parent
  `CMakeLists.txt` `SUBDIRS` list so top-level builds enter it.
- Keep verbose entry/exit messages in local `CMakeLists.txt` files consistent
  with the directory name.
- Do not pass compiler macros with `-DCMAKE_CXX_FLAGS`; use
  `OOMPH_EXTRA_COMPILE_DEFINES` or the relevant project option.

## Demo Drivers And Validation Data

Demo drivers are part of the test suite. When changing driver behavior, update
`validate.sh`, `validata/`, and the local `oomph_add_test(...)` declaration
together.

Validation guidance:

- Tests run in a `Validation` directory and compare selected output against
  compressed files in `validata/` using `scripts/fpdiff.py`.
- Avoid shell wildcards in `validate.sh` when concatenating output; wildcard
  expansion order can vary by platform.
- Sort output derived from unordered STL containers before comparing it.
- For drivers using `triangle`, prefer robust integral data or norms over
  element-by-element output because generated meshes can vary by machine or
  optimisation level.
- Keep validation output representative but small enough for stable self-tests.

## Documentation

Documentation is generated with Doxygen. Source and tutorial comments often use
Doxygen commands, `///` comments, and `.txt` pages under `doc/`.

- Keep docs and demo drivers in sync; tutorial docs often describe concrete
  driver code.
- Add new tutorials under `doc/` using the same structural pattern as existing
  documentation and demo-driver directories.
- If documentation changes affect generated examples or indexes, run a local doc
  build when practical.
- Do not hand-edit generated HTML or build artifacts.

## Workflow Expectations For Agents

- Read the relevant local docs and neighboring code before editing.
- Keep changes narrowly scoped and preserve existing numerical behavior unless
  the task asks for a behavioral change.
- Do not rewrite formatting across broad areas as part of an unrelated change.
- Do not modify vendored third-party code unless explicitly required.
- Prefer targeted builds/tests based on touched areas; full self-tests and docs
  are expensive.
- If a command depends on missing third-party installs or a long full build,
  state that clearly and provide the exact command that should be run.
- Before finalizing, check `git diff --check` and inspect the final diff.
