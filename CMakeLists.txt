cmake_minimum_required(VERSION 3.1)
project(Parametrize)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(TBB REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG}")

include_directories(src)
include_directories(3rd/pcg32)
include_directories(3rd/pss)
include_directories(3rd/glm)
include_directories(3rd/eigen)
include_directories(${GLM_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIRS})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${TBB_INCLUDE_DIRS})

set(
    parametrize_SRC
    src/adjacent-matrix.cpp
    src/adjacent-matrix.hpp
    src/compare-key.hpp
    src/config.hpp
    src/dedge.cpp
    src/dedge.hpp
    src/disajoint-tree.hpp
    src/dset.hpp
    src/field-math.hpp
    src/flow.hpp
    src/gldraw.cpp
    src/gldraw.hpp
    src/hierarchy.cpp
    src/hierarchy.hpp
    src/loader.cpp
    src/loader.hpp
    src/main.cpp
    src/merge-vertex.cpp
    src/merge-vertex.hpp
    src/optimizer.cpp
    src/Optimizer.cu
    src/optimizer.hpp
    src/parametrizer.cpp
    src/parametrizer.hpp
    src/post-solver.cpp
    src/post-solver.hpp
    src/serialize.hpp
    src/subdivide.cpp
    src/subdivide.hpp
)

add_executable(
    parametrize
    ${parametrize_SRC}
)
target_link_libraries(parametrize ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TBB_LIBRARIES})
