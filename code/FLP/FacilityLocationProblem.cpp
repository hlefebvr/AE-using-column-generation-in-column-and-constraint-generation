//
// Created by henri on 09/03/23.
//

#include "FacilityLocationProblem.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/Optimizers_DantzigWolfeDecomposition.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/BranchAndBound.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/node-selection-rules/factories/BestBound.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/DantzigWolfeDecomposition.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/branching-rules/factories/MostInfeasible.h"
#include "idol/optimizers/mixed-integer-optimization/callbacks/heuristics/IntegerMaster.h"
#include "idol/optimizers/mixed-integer-optimization/callbacks/cutting-planes/LazyCutCallback.h"
#include "idol/optimizers/mixed-integer-optimization/callbacks/cutting-planes/KnapsackCover.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/infeasibility-strategies/FarkasPricing.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/stabilization/Neame.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/logs/Info.h"

FacilityLocationProblem::FacilityLocationProblem(const Instance &t_instance, double t_Gamma)
    : m_instance(t_instance),
      m_Gamma(t_Gamma),
      m_theta(m_env, -1e14, 0, Continuous, "theta"),
      m_annotation(m_env, "decomposition", MasterId) {

    const unsigned int n_facilities = m_instance.n_facilities();

    m_x = Var::make_vector(m_env, Dim<1>(n_facilities), 0., 1., Binary, "x");
    m_xi = Var::make_vector(m_env, Dim<1>(n_facilities), 0., 1., Binary, "xi");

}

Model FacilityLocationProblem::create_master_problem() {

    const unsigned int n_facilities = m_instance.n_facilities();

    Model result(m_env);

    result.add(m_theta);
    result.add_vector<Var, 1>(m_x);

    result.set_obj_expr(m_theta + idol_Sum(i, Range(n_facilities), m_instance.fixed_cost(i) * m_x[i]));

    return result;
}

void FacilityLocationProblem::add_scenario_to_master_problem(Model &t_master,
                                                             const Solution::Primal &t_worst_case_scenario,
                                                             unsigned int t_iteration) {

    const unsigned int n_facilities = m_instance.n_facilities();
    const unsigned int n_customers = m_instance.n_customers();

    const auto is_disrupted = [&](unsigned int t_facility) -> bool {
        return t_worst_case_scenario.get(m_xi[t_facility]) > .5;
    };

    auto y = Var::make_vector(m_env, Dim<2>(n_facilities, n_customers), 0, 1, Binary, "y_" + std::to_string(t_iteration));

    if (t_master.optimizer().is<Optimizers::BranchAndBound<DefaultNodeInfo>>()) {

        const auto& branch_and_bound = t_master.optimizer().as<Optimizers::BranchAndBound<DefaultNodeInfo>>();
        const auto& dantzig_wolfe = branch_and_bound.relaxation().optimizer().as<Optimizers::DantzigWolfeDecomposition>();

        for (auto i  : Range(n_facilities)) {
            for (auto j : Range(n_customers)) {
                y[i][j].set(dantzig_wolfe.formulation().decomposition_by_variable(), t_iteration);
            }
        }

    }

    t_master.add_vector<Var, 2>(y);

    // Add objective constraint
    t_master.add_ctr(m_theta >=
                         idol_Sum(i, Range(n_facilities),
                                  idol_Sum(j, Range(n_customers), cost(i,j) * y[i][j])));

    // Add assignment constraints
    for (auto j : Range(n_customers)) {
        Ctr assignment_j(m_env, idol_Sum(i, Range(n_facilities), y[i][j]) <= 1);
        assignment_j.set(m_annotation, t_iteration);
        t_master.add(assignment_j);
    }

    // Add capacity constraints
    for (auto i : Range(n_facilities)) {

        Ctr capacity_i(m_env, idol_Sum(j, Range(n_customers), m_instance.demand(j) * y[i][j]) <= m_instance.capacity(i));
        capacity_i.set(m_annotation, t_iteration);
        t_master.add(capacity_i);

        t_master.add_ctr(idol_Sum(j, Range(n_customers), m_instance.demand(j) * y[i][j]) <= m_instance.capacity(i) * m_x[i]);

    }

    // Adding interdiction
    /*
    for (auto i : Range(n_facilities)) {
        for (auto j : Range(n_customers)) {
            t_master.add_ctr(y[i][j] <= m_x[i]);
        }
    }
     */

    // Add interdiction constraints
    for (auto i : Range(n_facilities)) {

        for (auto j: Range(n_customers)) {

            if (is_disrupted(i)) {
                t_master.set_var_ub(y[i][j], 0);
                // Alternatively: t_master.remove(y[i][j]); but not implemented yet in DantzigWolfeDecomposition
            }

        }

    }

}

double FacilityLocationProblem::cost(unsigned int t_facility, unsigned int t_customer) const {
    return std::floor(m_instance.demand(t_customer) * m_instance.per_unit_transportation_cost(t_facility, t_customer) - m_instance.profit(t_customer));
}

void FacilityLocationProblem::set_large_scale_optimizer(Model &t_master) {

    const double time_limit = t_master.optimizer().get_remaining_time();

    t_master.use(
            BranchAndBound()
                    .with_node_optimizer(
                            DantzigWolfeDecomposition(m_annotation)
                                    .with_master_optimizer(
                                            create_gurobi().with_continuous_relaxation_only(true)
                                    )
                                    .with_default_sub_problem_spec(
                                            DantzigWolfe::SubProblem()
                                                .add_optimizer(create_gurobi())
                                                .with_max_column_per_pricing(20)
                                                .with_column_pool_clean_up(300, .66)
                                    )
                                    .with_infeasibility_strategy(DantzigWolfe::FarkasPricing())
                                    .with_hard_branching(true)
                                    .with_dual_price_smoothing_stabilization(DantzigWolfe::Neame(.3))
                                    .with_max_parallel_sub_problems(m_parallel_pricing)
                                    .with_logger(Logs::DantzigWolfe::Info().with_frequency_in_seconds(5))
                                    .with_logs(true)
                    )
                    .with_logs(true)
                    .with_branching_rule(MostInfeasible(m_x.begin(), m_x.end()))
                    .with_node_selection_rule(BestBound())
                    .conditional(m_use_heuristic, [this](auto& x){
                        x.add_callback(
                                Heuristics::IntegerMaster()
                                        .with_optimizer(create_gurobi())
                                        .with_integer_columns(false)
                        );
                    })
                    .with_time_limit(time_limit)
    );

}

void FacilityLocationProblem::set_default_optimizer(Model &t_master) {
    t_master.use(create_gurobi());
}

Solution::Primal FacilityLocationProblem::compute_worst_case_scenario(const Model& t_master,
                                                                      const Solution::Primal &t_first_stage_solution,
                                                                      double t_time_limit) {

    const unsigned int n_facilities = m_instance.n_facilities();
    const unsigned int n_customers = m_instance.n_customers();

    Model model(m_env);

    // Create master problem for separation
    model.set_obj_sense(Maximize);
    model.add_vector<Var, 1>(m_xi);
    auto pi = model.add_var(-Inf, 0, Continuous, Column(1), "pi");
    auto budget = model.add_ctr(idol_Sum(i, Range(n_facilities), m_xi[i]) <= m_Gamma);

    // Create separation problem
    Model separation(m_env);

    auto y = separation.add_vars(Dim<2>(n_facilities, n_customers), 0, 1, Binary);

    // Capacities
    for (auto i : Range(n_facilities)) {
        separation.add_ctr(idol_Sum(j, Range(n_customers), m_instance.demand(j) * y[i][j]) <= m_instance.capacity(i) * t_first_stage_solution.get(m_x[i]));
    }

    // Assignments
    for (auto j : Range(n_customers)) {
        separation.add_ctr(idol_Sum(i, Range(n_facilities), y[i][j]) <= 1);
    }

    // Create interdiction cut
    auto cut = pi <= idol_Sum(i, Range(n_facilities),
                                          idol_Sum(j, Range(n_customers),
                                                   cost(i,j) * !y[i][j]
                                                   - cost(i,j) * !y[i][j] * m_xi[i]
                                          )
    );

    model.use(
            create_gurobi().with_lazy_cut(true)
                .add_callback(
                        LazyCutCallback(separation, std::move(cut))
                            .with_separation_optimizer(create_gurobi())
                    )
                .with_time_limit(t_time_limit)
            );

    model.optimize();

    if (model.get_status() != Optimal) {
        Solution::Primal result;
        result.set_status(model.get_status());
        result.set_reason(model.get_reason());
        return result;
    }

    auto solution = save_primal(model);

    solution.set_objective_value(solution.objective_value() - t_first_stage_solution.get(m_theta));

    return solution;
}

Gurobi FacilityLocationProblem::create_gurobi() const {
    return std::move(Gurobi().with_thread_limit(1));
}

Solution::Primal FacilityLocationProblem::compute_initial_scenario() {
    return {};
}
