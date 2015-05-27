cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")
set(CMAKE_C_COMPILER "/usr/bin/mpicc")

project(mpi_matrices)

find_package(MPI REQUIRED)
if ( MPI_FOUND )
        include_directories( ${MPI_INCLUDE_PATH} )
endif( MPI_FOUND )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/dense_matrix.cpp
    src/dense_matrix.h
    src/generator.cpp
    src/generator.h
    src/direction.h
    src/sparse_matrix_elem.h
    src/main.cpp
    src/mpimatrix.cpp
    src/mpimatrix_op.cpp
    src/mpimatrix_lu.cpp
    src/mpimatrix.h
    src/sparse_vector.cpp
    src/sparse_vector_op.cpp
    src/sparse_vector.h
    src/sparse_matrix.cpp
    src/sparse_matrix_op.cpp
    src/sparse_matrix.h
    src/mpimatrix_cg.cpp)

add_executable(mpi_matrices ${SOURCE_FILES} src/mpimatrix_prec.cpp)

target_link_libraries(mpi_matrices ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(mpi_matrices PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(mpi_matrices PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()