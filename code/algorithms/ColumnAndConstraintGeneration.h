//
// Created by henri on 07/04/23.
//

#ifndef CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATION_H
#define CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATION_H

#include "ColumnAndConstraintGenerator.h"

class ColumnAndConstraintGeneration {
    ColumnAndConstraintGenerator& m_generator;
    const double m_time_limit;
public:
    explicit ColumnAndConstraintGeneration(ColumnAndConstraintGenerator& t_generator, double t_time_limit = 3600);

    Solution::Primal solve(double t_std_phase_time_limit = std::numeric_limits<double>::max(), const std::string& t_tag = "");
};

#endif //CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATION_H
