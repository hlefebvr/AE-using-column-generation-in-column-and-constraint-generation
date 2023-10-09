add_executable(JSP_solve main.cpp Instance.cpp Instance.h
        JobSchedulingProblem.cpp
        JobSchedulingProblem.h)

target_link_libraries(JSP_solve PUBLIC idol_ccg)

add_executable(JSP_generate
        generate.cpp
        Instance.cpp
        Instance.h
)
