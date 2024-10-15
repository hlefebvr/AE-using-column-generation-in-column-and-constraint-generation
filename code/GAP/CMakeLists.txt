
add_executable(GAP_solve
        main.cpp
        GeneralAssignmentProblem.cpp
        GeneralAssignmentProblem.h
        Instance.cpp
        Instance.h)

target_link_libraries(GAP_solve PUBLIC idol_ccg)
