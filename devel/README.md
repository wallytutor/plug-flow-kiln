# Development

Full automation of builds is fragile and not recommended. Here we provide guidelines to build locally the main project dependencies for development of the rotary kiln model. Main libraries are compiled locally while tools are recovered from MSYS2 repositories with `pacman`.

It is recommended to work with the latest versions of all packages and some of the commands below might break when an update is performed. Also use a fresh [MSYS2 environment](https://www.msys2.org/) to get reproducible builds. Once installation is finished, open a new MSYS2 window and run `pacman -Syyu`, accepting the following propmpts. After window has automatically closed, open a `mingw64` terminal and start following the recommendations given here.

**NOTE:** it is important that the next steps be followed from a `mingw64` terminal rather than in a `msys2` one because otherwise the installed packages and generated executable will not be portable.

The following snippet provides the minimum environment for building the package.

```bash
# General tools used during differnt build stages.
pacman -S     \
    autoconf  \
    automake  \
    binutils  \
    diffutils \
    git       \
    grep      \
    libtool   \
    make      \
    patch     \
    pkg-config 

# Install compiler suit for the platform.
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran

# It is easier to proceed with platform LAPACK and Metis.
pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-metis

# For `mpi.h` to be available you need to install this.
pacman -S mingw-w64-x86_64-msmpi

# The version from MSYS2 will not work, install the following.
pacman -S mingw-w64-x86_64-cmake

# Required for compilation of Cantera.
pacman -S mingw-w64-x86_64-boost \
          mingw-w64-x86_64-python \
          mingw-w64-x86_64-cython \
          mingw-w64-x86_64-python-pip \
          mingw-w64-x86_64-python-setuptools \
          mingw-w64-x86_64-python-numpy \
          mingw-w64-x86_64-python-ruamel-yaml \
          mingw-w64-x86_64-python-ruamel.yaml.clib \
          mingw-w64-x86_64-python-pytest \
          mingw-w64-x86_64-python-pytest-cov \
          mingw-w64-x86_64-scons \
          mingw-w64-x86_64-eigen3 \
          mingw-w64-x86_64-fmt \
          mingw-w64-x86_64-gtest
```

Other `mingw64` compatible packages can be searched for [here](https://packages.msys2.org/package/?repo=mingw64).

In the following steps it is assumed packages are cloned in the user home directory inside MSYS install folder. This avoids cluttering the present project and gives an easily accessible path to builds. All the next steps assume the following variable was exported:

```bash
export FAKEROOT=${HOME}/fakeroot/
```

## Build MUMPS

This solver is required by both Ipopt and CasADi so it has to be built prior to those software. In [Ipopt documentation](https://coin-or.github.io/Ipopt/INSTALL.html) it is recommended to use a COIN-OR tested version of the software. The following should enable the build of the package.

```bash
cd ${HOME}

git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps

./get.Mumps

CC=/mingw64/bin/gcc CXX=/mingw64/bin/g++ ./configure --prefix=${FAKEROOT}

make
make install
```

## Build Ipopt

Once MUMPS has finished we can start with Ipopt. The reference build instructions is provided [here](https://coin-or.github.io/Ipopt/INSTALL.html). Notice that the patched version of MUMPS compiled above has some extra includes that are not installed because they are used only for the compilation of Ipopt, so the configuration option `--with-mumps-cflags` has to be extended as given below.

```bash
cd ${HOME}

git clone https://github.com/coin-or/Ipopt.git
cd Ipopt && git checkout releases/3.14.12

mkdir build && cd build

CC=/mingw64/bin/gcc CXX=/mingw64/bin/g++ ../configure   \
    --prefix=${FAKEROOT}                              \
    --with-mumps-cflags=-I${FAKEROOT}/include/coin-or/mumps -I${HOME}/ThirdParty-Mumps \
    --with-mumps-lflags=-L${FAKEROOT}/bin -lcoinmumps-3

make
make test
make install
```

## Build CasADi

The reference build instructions is provided [here](https://github.com/casadi/casadi/wiki/InstallationLinux).

```bash
cd ${HOME}

git clone https://github.com/casadi/casadi.git
cd casadi && git checkout 3.6.1

mkdir build && cd build

# Because it seems CMake is not searching on MUMPS_INCLUDE_DIR.
ln -s ${FAKEROOT}/include/coin-or/mumps/ /mingw64/include/mumps

CC=/mingw64/bin/gcc CXX=/mingw64/bin/g++ cmake ..                   \
    -DCMAKE_INSTALL_PREFIX:PATH=${FAKEROOT}                       \
    -DWITH_THREAD_MINGW:BOOL=ON                                     \
    -DFORTRAN_REQUIRED:BOOL=true                                    \
    -DWITH_IPOPT:BOOL=ON                                            \
    -DWITH_LAPACK:BOOL=ON                                           \
    -DWITH_MUMPS:BOOL=ON                                            \
    -DIPOPT_INCLUDE_DIRS=${FAKEROOT}/include/coin-or              \
    -DIPOPT_LIBRARY_DIRS=${FAKEROOT}/lib                          \
    -DIPOPT_LIBRARIES=libipopt                                    \
    -DMUMPS_INCLUDE_DIR:PATH=${FAKEROOT}/include/coin-or          \
    -DMUMPS_LIB_dmumps_seq:PATH=${FAKEROOT}/bin/libcoinmumps-3.dll

cmake --build . --config Release
cmake --install .
```

## Build Cantera

```bash
git clone --recursive https://github.com/Cantera/cantera.git
cd cantera && git checkout tags/v2.5.0 && git submodule update

# Because scons build does not provide a way to configure the name of
# Python library, this workaround has to be run before build of Cantera.
ln -s /mingw64/lib/libpython3.10.dll.a /mingw64/lib/libpython310.dll.a

rm -rf cantera.conf && scons clean && scons build \
    CXX='/mingw64/bin/g++' \
    CC='/mingw64/bin/gcc' \
    FORTRAN='/mingw64/bin/gfortran' \
    FORTRANFLAGS='-O3' \
    VERBOSE='no' \
    extra_inc_dirs='/mingw64/include/eigen3' \
    extra_lib_dirs='/mingw64/lib' \
    prefix=$FAKEROOT/ \
    target_arch='amd64' \
    toolchain='mingw' \
    f90_interface='default' \
    python_package='full' \
    system_eigen='y' \
    system_fmt='n' \
    system_yamlcpp='n' \
    system_sundials='n' \
    sundials_include='/mingw64/include' \
    sundials_libdir='/mingw64/lib' \
    blas_lapack_libs='blas,cblas,lapack' \
    blas_lapack_dir='/mingw64/lib' \
    lapack_names='lower' \
    lapack_ftn_trailing_underscore='yes' \
    lapack_ftn_string_len_at_end='yes' \
    googletest='submodule' \
    use_pch='no' \
    cxx_flags='-g -Wextra -O3 --std=c++17' \
    cc_flags='' \
    thread_flags='-pthread' \
    optimize='yes' \
    optimize_flags='-O3 -Wno-inline' \
    debug='no' \
    warning_flags='-Wall' \
    renamed_shared_libraries='yes' \
    versioned_shared_library='no' \
    use_rpath_linkage='no' \
    layout='standard' \
    fast_fail_tests='no'
```

The `scons build` failed several times and we could not it working out-of-the box. Currently Cantera seems to be incompatyble with `mingw64` version of Sundials, what produced problems with finding Boost or even testing the Fortran compiler, so it was decided to go with their downloaded version of Sundials and build locally. Another major issue was the expansion of paths from POSIX style to Windows style. When running configuration some paths are converted and others not. There also seems to be some issues with order of commands. The following snippet of `cantera.conf` file generated by `scons`. The main difficulty was to get Eigen being found and the final linkage of the program, what required to escape the path in `extra_inc_dirs` and `extra_lib_dirs`, respectively.

```ini
prefix = '<replace-by-the-expanded-FAKEROOT-value>'
CXX = 'C:/msys64/mingw64/bin/g++'
CC = 'C:/msys64/mingw64/bin/gcc'
FORTRAN = 'C:/msys64/mingw64/bin/gfortran'
python_package = 'full'
system_eigen = 'y'
system_fmt = 'n'
system_yamlcpp = 'n'
system_sundials = 'n'
sundials_include = 'C:/msys64/mingw64/include'
sundials_libdir = 'C:/msys64/mingw64/lib'
blas_lapack_libs = 'blas,cblas,lapack'
blas_lapack_dir = 'C:/msys64/mingw64/lib'
googletest = 'submodule'
use_pch = False
cxx_flags = '-g -Wextra -O3 --std=c++17'
thread_flags = '-pthread'
debug = False
extra_inc_dirs = 'C:\\msys64\\mingw64\\include\\eigen3'
extra_lib_dirs = 'C:\\msys64\\mingw64\\lib'
use_rpath_linkage = False
layout = 'standard'
```

Once everything is fine, finish with the following:

```bash
scons test
scons install
```

**NOTE:** Cantera dynamic linking it not working with the current compilation. Only static linking works and produce huge executables.

## System organization

Not all packages place their outputs in the desired folders. Next we make a quick organisation for the ease of generating configuration files for this project. There is room for improvement in the following snippet because it was done on the fly.

```bash
# CMake files
mv "${FAKEROOT}/casadi/cmake/" "${FAKEROOT}/share/"
mv "${FAKEROOT}/lib/cmake/tinyxml2/" "${FAKEROOT}/share/cmake/"

# DLL's
mv "${FAKEROOT}/casadi/*.dll" "${FAKEROOT}/bin/"
mv "${FAKEROOT}/casadi/*.exe" "${FAKEROOT}/bin/"
mv "${FAKEROOT}/lib/*.dll" "${FAKEROOT}/bin/"

# Static libraries
mv "${FAKEROOT}/casadi/*.dll.a" "${FAKEROOT}/lib/"

# Includes
mv "${FAKEROOT}/casadi/include/casadi" "${FAKEROOT}/include/"

# Package configuration
mv "${FAKEROOT}/casadi/pkgconfig/casadi.pc" "${FAKEROOT}/lib/pkgconfig/"

rm -rf "${FAKEROOT}/lib/cmake/"
rm -rf "${FAKEROOT}/casadi/"
```

## Testing

To test the functioning of the libraries, you can now place the [source.bat](source.bat) script provided here inside `$FAKEROOT` and copy one of the Ipopt-dependent executables generated inside `$HOME/casadi/build` (*e.g.* the `rosenbrock.exe`) to that same location. Running the source script from a Windows terminal will set the path so the file can be executed and all dependencies are confirmed to be in the path.

## To-do's

- [ ] Add builds of all third-parties from [here](https://github.com/orgs/coin-or-tools).
- [ ] Check warning related to CasADi's `CASADI_EXPORT` definition.
- [ ] Confirm in a computer without MSYS2 that the distribution still works.
- [ ] Compile Sundials instead of accepting default Cantera version.
- [ ] Add licenses of all packages to a single folder.
