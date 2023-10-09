//
// Created by henri on 09/03/23.
//

#ifndef CCG_WITH_NESTED_CG_FACILITYLOCATIONPROBLEM_H
#define CCG_WITH_NESTED_CG_FACILITYLOCATIONPROBLEM_H

#include "Instance.h"
#include "../algorithms/ColumnAndConstraintGenerator.h"
#include "optimizers/solvers/gurobi/Gurobi.h"

using namespace idol;

class FacilityLocationProblem : public ColumnAndConstraintGenerator {
    Env m_env;
    const Instance& m_instance;
    const double m_Gamma;
    Annotation<Ctr> m_annotation;

    bool m_use_heuristic = true;
    unsigned int m_parallel_pricing = 5;

    Var m_theta;
    Vector<Var, 1> m_x;
    Vector<Var, 1> m_xi;

    [[nodiscard]] double cost(unsigned int t_facility, unsigned int t_customer) const;
public:
    explicit FacilityLocationProblem(const Instance& t_instance, double t_Gamma);

    void use_heuristic(bool t_value) { m_use_heuristic = t_value; }

    void use_parallel_pricing(unsigned int t_concurrency) { m_parallel_pricing = t_concurrency; }

    [[nodiscard]] double Gamma() const { return m_Gamma; }
protected:
    Model create_master_problem() override;

    Solution::Primal compute_worst_case_scenario(const Model& t_master, const Solution::Primal &t_first_stage_solution, double t_time_limit) override;

    void add_scenario_to_master_problem(Model &t_master,
                                        const Solution::Primal &t_worst_case_scenario,
                                        unsigned int t_iteration) override;

    void set_large_scale_optimizer(Model &t_master) override;

    void set_default_optimizer(Model &t_master) override;

    Solution::Primal compute_initial_scenario() override;

    [[nodiscard]] Gurobi create_gurobi() const;
};


#endif //CCG_WITH_NESTED_CG_FACILITYLOCATIONPROBLEM_H
