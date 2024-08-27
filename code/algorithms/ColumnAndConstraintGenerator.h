//
// Created by henri on 07/04/23.
//

#ifndef CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATOR_H
#define CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATOR_H

#include "idol/modeling.h"

using namespace idol;

class ColumnAndConstraintGenerator {
    friend class ColumnAndConstraintGeneration;
public:
    virtual ~ColumnAndConstraintGenerator() = default;
protected:
    [[nodiscard]] virtual Model create_master_problem() = 0;

    virtual void set_large_scale_optimizer(Model& t_master) = 0;

    virtual void set_default_optimizer(Model& t_master) = 0;

    virtual Solution::Primal compute_initial_scenario() = 0;

    virtual Solution::Primal compute_worst_case_scenario(const Model& t_master, const Solution::Primal& t_first_stage_solution, double t_time_limit) = 0;

    virtual void add_scenario_to_master_problem(Model& t_master,
                                                const Solution::Primal& t_worst_case_scenario,
                                                unsigned int t_iteration) = 0;
};

#endif //CCG_WITH_NESTED_CG_COLUMNANDCONSTRAINTGENERATOR_H
