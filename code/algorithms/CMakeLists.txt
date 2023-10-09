add_library(idol_ccg STATIC
        ColumnAndConstraintGeneration.cpp
        ColumnAndConstraintGeneration.h
        ColumnAndConstraintGenerator.h
    )

target_link_libraries(idol_ccg PUBLIC idol)