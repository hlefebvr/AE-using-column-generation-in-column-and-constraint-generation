add_executable(JSP2_solve main.cpp
        JobSchedulingProblem.cpp
        JobSchedulingProblem.h
        ../JSP/Instance.cpp
        ../JSP/Instance.h)

target_link_libraries(JSP2_solve PUBLIC idol_ccg)

add_executable(JSP2_generate
        generate.cpp
)
