cmake_minimum_required(VERSION 3.1)
project(controlsim)

find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(glm REQUIRED)
find_package(Eigen3 REQUIRED)

# find_package(yaml-cpp REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(YAML_CPP REQUIRED yaml-cpp)


set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fsanitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

include_directories(src)
include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${GLM_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

set(
    controlsim_SRC
    src/config.hpp
    src/common.hpp
    src/geometry.hpp
    src/main.cpp
    src/material.cpp
    src/material.hpp
    src/displaywidget.cpp
    src/displaywidget.hpp
    src/mainwidget.hpp
    src/mainwidget.cpp
    src/simulator.cpp
    src/simulator.hpp
    src/mesh.cpp
    src/mesh.hpp
    src/cloth.cpp
    src/cloth.hpp
)
qt5_add_resources(
    controlsim_RES
    res/resource.qrc
)

add_executable(
    controlsim
    ${controlsim_SRC}
    ${controlsim_RES}
)
target_link_libraries(controlsim Qt5::Widgets ${YAML_CPP_LIBRARIES})
