add_executable(MJSP_solve main.cpp Instance.cpp Instance.h
        JobSchedulingProblem.cpp
        JobSchedulingProblem.h)

target_link_libraries(MJSP_solve PUBLIC idol_ccg)

add_executable(MJSP_generate
        generate.cpp
        Instance.cpp
        Instance.h
)
