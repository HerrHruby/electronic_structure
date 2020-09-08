cmake_minimum_required(VERSION 3.17)
project(HF)

set(CMAKE_CXX_STANDARD 17)

set(BUILD_1 main)
set(SOURCE_FILES_1 main.cpp)
add_executable(HF main.cpp include/integrals.h src/integrals.cpp src/basis_function.cpp include/basis_function.h src/contracted_integrals.cpp include/contracted_integrals.h src/xyzfile_reader.cpp include/xyzfile_reader.h src/scf.cpp include/STO3G_data.h include/scf.h)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Arb)

target_link_libraries (HF Eigen3::Eigen Arb)



