//
// Created by henri on 15.08.23.
//

#ifndef CCG_WITH_NESTED_CG_JOBSCHEDULINGPROBLEM_H
#define CCG_WITH_NESTED_CG_JOBSCHEDULINGPROBLEM_H

#include "../algorithms/ColumnAndConstraintGenerator.h"
#include "Instance.h"
#include "idol/optimizers/mixed-integer-optimization/wrappers/Gurobi/Gurobi.h"

class JobSchedulingProblem : public ColumnAndConstraintGenerator {
    idol::Env m_env;
    const Instance& m_instance;
    const std::vector<std::vector<JobOccurrence>> m_job_occurrences;
    const std::vector<std::vector<double>> m_big_M;
    const double m_Gamma;
    const double m_percentage_increase = .3;
    Annotation<Ctr> m_annotation;

    bool m_use_heuristic = true;
    unsigned int m_parallel_pricing = 1;

    idol::Var m_theta;
    idol::Vector<idol::Var, 1> m_x;
    idol::Vector<idol::Var, 1> m_xi;

public:
    JobSchedulingProblem(const Instance& t_instance, double t_Gamma);

    [[nodiscard]] double Gamma() const { return m_Gamma; }

    void use_heuristic(bool t_value) { m_use_heuristic = t_value; }

    void use_parallel_pricing(unsigned int t_concurrency) { m_parallel_pricing = t_concurrency; }

protected:
    Model create_master_problem() override;

    void set_large_scale_optimizer(Model &t_master) override;

    void set_default_optimizer(Model &t_master) override;

    Solution::Primal compute_initial_scenario() override;

    void add_scenario_to_master_problem(Model &t_master,
                                        const Solution::Primal &t_worst_case_scenario,
                                        unsigned int t_iteration) override;

    Solution::Primal
    compute_worst_case_scenario(const Model &t_master, const Solution::Primal &t_first_stage_solution, double t_time_limit) override;

    [[nodiscard]] Gurobi create_gurobi() const;
};


#endif //CCG_WITH_NESTED_CG_JOBSCHEDULINGPROBLEM_H
