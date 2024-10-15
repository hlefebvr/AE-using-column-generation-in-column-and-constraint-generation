cmake_minimum_required(VERSION 3.22)
project(ccg_with_nested_cg)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

set(USE_GUROBI ON)
set(USE_EIGEN OFF)
set(USE_MOSEK OFF)
set(USE_ROOT OFF)
set(USE_MIBS ON)
set(USE_OSI ON)
set(USE_CLP ON)
set(USE_CPLEX ON)

FetchContent_Declare(
        idol
        GIT_REPOSITORY https://github.com/hlefebvr/idol.git
        GIT_TAG        origin/main # v0.8.1-alpha
)

FetchContent_MakeAvailable(idol)

add_subdirectory(algorithms)
add_subdirectory(FLP)
add_subdirectory(JSP)
add_subdirectory(MJSP)
