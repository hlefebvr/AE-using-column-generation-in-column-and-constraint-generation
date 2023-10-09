
add_executable(FLP_solve
        main.cpp
        FacilityLocationProblem.cpp
        FacilityLocationProblem.h
        Instance.cpp
        Instance.h)

target_link_libraries(FLP_solve PUBLIC idol_ccg)

if (WITH_DEBUG_INFO)
    message(STATUS "Activated debug info")
    target_compile_definitions(FLP_solve PRIVATE WITH_DEBUG_INFO)
    target_compile_options(FLP_solve PRIVATE -g)
    target_compile_options(FLP_solve PRIVATE -rdynamic)
endif()
