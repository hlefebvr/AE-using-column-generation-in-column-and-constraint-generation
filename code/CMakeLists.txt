cmake_minimum_required(VERSION 3.22)
project(ccg_with_nested_cg)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

set(USE_GUROBI ON)
set(USE_EIGEN OFF)
set(USE_MOSEK OFF)
set(USE_ROOT OFF)

FetchContent_Declare(
        idol
        GIT_REPOSITORY https://github.com/hlefebvr/idol.git
        GIT_TAG        v0.8.1-alpha # v0.8.0-alpha # origin/main
)

FetchContent_MakeAvailable(idol)

add_subdirectory(algorithms)
add_subdirectory(FLP)
add_subdirectory(JSP)

