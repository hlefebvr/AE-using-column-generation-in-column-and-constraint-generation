//
// Created by henri on 15.08.23.
//

#include "JobSchedulingProblem.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/DantzigWolfeDecomposition.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/BranchAndBound.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/node-selection-rules/factories/DepthFirst.h"
#include "idol/optimizers/mixed-integer-optimization/callbacks/heuristics/IntegerMaster.h"
#include "idol/optimizers/mixed-integer-optimization/branch-and-bound/branching-rules/factories/MostInfeasible.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/Optimizers_DantzigWolfeDecomposition.h"
#include "idol/optimizers/mixed-integer-optimization/callbacks/cutting-planes/LazyCutCallback.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/infeasibility-strategies/FarkasPricing.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/stabilization/Neame.h"
#include "idol/optimizers/mixed-integer-optimization/dantzig-wolfe/logs/Info.h"
#include "idol/modeling/bilevel-optimization/LowerLevelDescription.h"
#include "idol/optimizers/bilevel-optimization/wrappers/MibS/MibS.h"

JobSchedulingProblem::JobSchedulingProblem(const Instance &t_instance, double t_Gamma)
    : m_instance(t_instance),
      m_Gamma(t_Gamma),
      m_job_occurrences(m_instance.compute_job_occurrences()),
      m_big_M(Instance::compute_big_M(m_job_occurrences)),
      m_annotation(m_env, "decomposition", MasterId),
      m_theta(m_env, -Inf, Inf, Continuous, "theta"),
      m_x(idol::Var::make_vector<1>(m_env, Dim<1>(m_instance.n_machines()), 0, 1, Binary, "x")),
      m_xi(idol::Var::make_vector<1>(m_env, Dim<1>(m_instance.n_jobs()), 0, 1, Binary, "xi"))
    {

}

Model JobSchedulingProblem::create_master_problem() {

    const auto n_machines = m_instance.n_machines();

    idol::Model result(m_env);

    result.add(m_theta);
    result.add_vector<idol::Var, 1>(m_x);

    Expr objective = m_theta;
    for (unsigned int i = 0; i < n_machines; ++i) {
        objective += m_instance.machine(i).fixed_cost * m_x[i];
    }

    result.set_obj_expr(objective);

    return result;
}

Solution::Primal JobSchedulingProblem::compute_initial_scenario() {
    return {};
}

void JobSchedulingProblem::set_default_optimizer(Model &t_master) {

    t_master.use(create_gurobi());

}

void JobSchedulingProblem::set_large_scale_optimizer(Model &t_master) {

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
                    .with_subtree_depth(1)
                    .with_branching_rule(MostInfeasible(m_x.begin(), m_x.end()))
                    .with_node_selection_rule(DepthFirst())
                    .conditional(m_use_heuristic, [this](auto& x){
                        x.add_callback(
                                Heuristics::IntegerMaster()
                                        .with_optimizer(create_gurobi())
                                        .with_integer_columns(true)
                        );
                    })
                    .with_time_limit(time_limit)
    );

}

Gurobi JobSchedulingProblem::create_gurobi() const {
    return std::move(Gurobi().with_thread_limit(1));
}

void
JobSchedulingProblem::add_scenario_to_master_problem(Model &t_master,
                                                     const Solution::Primal &t_worst_case_scenario,
                                                     unsigned int t_iteration) {

    const unsigned int n_machines = m_instance.n_machines();
    const auto n_jobs = m_instance.n_jobs();

    const auto is_attacked = [&](unsigned int t_job_index) -> bool {
        return t_worst_case_scenario.get(m_xi[t_job_index]) > .5;
    };

    auto U = idol::Var::make_vector(m_env, Dim<1>(n_jobs), 0, 1, Binary, "U");

    std::vector<std::vector<Var>> y;
    y.reserve(n_machines);

    std::vector<std::vector<Var>> t;
    t.reserve(n_machines);

    for (unsigned int i = 0; i < n_machines; ++i) {
        const auto &n_job_occurrences = m_job_occurrences[i].size();
        y.push_back(idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, 1, Binary, "y"));
        t.push_back(idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, Inf, Continuous, "t"));
    }


    if (t_master.optimizer().is<idol::Optimizers::BranchAndBound<DefaultNodeInfo>>()) {

        auto &branch_and_bound = t_master.optimizer().as<idol::Optimizers::BranchAndBound<DefaultNodeInfo>>();

        if (branch_and_bound.relaxation().optimizer().is<idol::Optimizers::DantzigWolfeDecomposition>()) {

            const auto &dantzig_wolfe = branch_and_bound.relaxation().optimizer().as<idol::Optimizers::DantzigWolfeDecomposition>();
            const auto &var_annotation = dantzig_wolfe.formulation().decomposition_by_variable();

            for (unsigned int j = 0; j < n_jobs; ++j) {
                U[j].set(var_annotation, t_iteration);
            }

            for (unsigned int i = 0; i < n_machines; ++i) {
                for (auto &y_ik: y[i]) {
                    y_ik.set(var_annotation, t_iteration);
                }
                for (auto &t_ik: t[i]) {
                    t_ik.set(var_annotation, t_iteration);
                }
            }

        }

    }

    t_master.add_vector<idol::Var, 1>(U);
    t_master.add_vector<idol::Var, 2>(y);
    t_master.add_vector<idol::Var, 2>(t);

    // Compute sum of y_k over G_j for each job j and each machine i
    std::vector<std::vector<Expr<Var, Var>>> sum_y_k(n_machines, std::vector<Expr<Var, Var>>(n_jobs, 0));
    for (unsigned int i = 0; i < n_machines; ++i) {
        const auto &n_job_occurrences = m_job_occurrences[i].size();
        for (unsigned int k = 0; k < n_job_occurrences; ++k) {
            sum_y_k[i][m_job_occurrences[i][k].parent->index] += y[i][k];
        }
    }

    // Objective Function
    Expr objective;
    for (unsigned int j = 0; j < n_jobs; ++j) {
        objective += m_instance.job(j).weight * U[j];
    }
    for (unsigned int i = 0; i < n_machines; ++i) {
        for (unsigned int j = 0; j < n_jobs; ++j) {
            objective -= m_instance.job(j).profit[i] * sum_y_k[i][j];
        }
    }
    t_master.add_ctr(m_theta >= objective);

    // Linking constraints
    for (unsigned int i = 0; i < n_machines; ++i) {
        t_master.add_ctr(idol_Sum(j, Range(n_jobs), sum_y_k[i][j]) <= m_x[i]);
    }

    // Assignment constraints
    for (unsigned int j = 0; j < n_jobs; ++j) {
        t_master.add_ctr(U[j] + idol_Sum(i, Range(n_machines), sum_y_k[i][j]) == 1);
    }

    // Deadlines
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 0; k < n_job_occurrences; ++k) {
            Ctr c(m_env, t[i][k] <= m_job_occurrences[i][k].deadline);
            c.set(m_annotation, t_iteration);
            t_master.add(c);
        }

    }

    // Non-overlapping constraints
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 1; k < n_job_occurrences; ++k) {

            const auto &job_occurrence = m_job_occurrences[i][k];
            const double processing_time = job_occurrence.parent->processing_time[i];

            Ctr c(m_env, t[i][k] - t[i][k - 1] - processing_time * y[i][k] >= 0);
            c.set(m_annotation, t_iteration);
            t_master.add(c);
        }
    }

    // Disjunctive + release dates
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 0; k < n_job_occurrences; ++k) {

            const auto &job_occurrence = m_job_occurrences[i][k];
            const double processing_time = job_occurrence.parent->processing_time[i];

            Ctr c(m_env, t[i][k] - processing_time * y[i][k] - m_big_M[i][k] * y[i][k] >= job_occurrence.release_date - m_big_M[i][k]);
            c.set(m_annotation, t_iteration);
            t_master.add(c);
        }

    }

    // Scenario constraints
    for (unsigned int j = 0; j < n_jobs; ++j) {
        if (is_attacked(j)) {
            t_master.add_ctr(U[j] == 0);
        }
    }

}

Solution::Primal JobSchedulingProblem::compute_worst_case_scenario(const Model &t_master,
                                                                   const Solution::Primal &t_first_stage_solution,
                                                                   double t_time_limit) {

    for (const auto& [var, val] : t_first_stage_solution) {
        if (var.name().find("x") == 0) {
            std::cout << var.name() << " = " << val << std::endl;
        }
    }

    const auto n_jobs = m_instance.n_jobs();
    const auto n_machines = m_instance.n_machines();

    idol::Model model(m_env);
    Bilevel::LowerLevelDescription description(m_env);

    auto U = idol::Var::make_vector(m_env, Dim<1>(n_jobs), 0, 1, Binary, "U");

    std::vector<std::vector<Var>> y;
    y.reserve(n_machines);

    std::vector<std::vector<Var>> t;
    t.reserve(n_machines);

    for (unsigned int i = 0; i < n_machines; ++i) {
        const auto &n_job_occurrences = m_job_occurrences[i].size();
        y.push_back(idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, 1, Binary, "y_" + std::to_string(i)));
        t.push_back(idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, Inf, Continuous, "t_" + std::to_string(i)));
    }

    model.add_vector<Var, 1>(U);
    model.add_vector<Var, 2>(y);
    model.add_vector<Var, 2>(t);

    // Compute sum of y_k over G_j for each job j and each machine i
    std::vector<std::vector<Expr<Var, Var>>> sum_y_k(n_machines, std::vector<Expr<Var, Var>>(n_jobs, 0));
    for (unsigned int i = 0; i < n_machines; ++i) {
        const auto &n_job_occurrences = m_job_occurrences[i].size();
        for (unsigned int k = 0; k < n_job_occurrences; ++k) {
            sum_y_k[i][m_job_occurrences[i][k].parent->index] += y[i][k];
        }
    }

    // Objective Function
    Expr objective;
    for (unsigned int j = 0; j < n_jobs; ++j) {
        objective += m_instance.job(j).weight * U[j];
    }
    for (unsigned int i = 0; i < n_machines; ++i) {
        for (unsigned int j = 0; j < n_jobs; ++j) {
            objective -= m_instance.job(j).profit[i] * sum_y_k[i][j];
        }
    }
    model.set_obj_expr(-1. * objective);
    description.set_follower_obj_expr(objective);

    // Linking constraints
    for (unsigned int i = 0; i < n_machines; ++i) {
        model.add_ctr(idol_Sum(j, Range(n_jobs), sum_y_k[i][j]) <= std::round(t_first_stage_solution.get(m_x[i])));
    }

    // Assignment constraints
    for (unsigned int j = 0; j < n_jobs; ++j) {
        model.add_ctr(U[j] + idol_Sum(i, Range(n_machines), sum_y_k[i][j]) == 1);
    }

    // Deadlines
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 0; k < n_job_occurrences; ++k) {
            model.add_ctr(t[i][k] <= m_job_occurrences[i][k].deadline);
        }

    }

    // Non-overlapping constraints
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 1; k < n_job_occurrences; ++k) {

            const auto &job_occurrence = m_job_occurrences[i][k];
            const double processing_time = job_occurrence.parent->processing_time[i];

            model.add_ctr(t[i][k] - t[i][k - 1] - processing_time * y[i][k] >= 0);
        }
    }

    // Disjunctive + release dates
    for (unsigned int i = 0; i < n_machines; ++i) {

        const auto &n_job_occurrences = m_job_occurrences[i].size();

        for (unsigned int k = 0; k < n_job_occurrences; ++k) {

            const auto &job_occurrence = m_job_occurrences[i][k];
            const double processing_time = job_occurrence.parent->processing_time[i];

            model.add_ctr(t[i][k] - processing_time * y[i][k] - m_big_M[i][k] * y[i][k] >= job_occurrence.release_date - m_big_M[i][k]);
        }

    }

    for (const auto& ctr : model.ctrs()) {
        description.make_follower_ctr(ctr);
    }
    for (const auto& var : model.vars()) {
        description.make_follower_var(var);
    }

    model.add_vector<Var, 1>(m_xi);

    model.add_ctr(idol_Sum(j, Range(n_jobs), m_xi[j]) <= m_Gamma);

    for (unsigned int j = 0; j < n_jobs; ++j) {
        auto c = model.add_ctr(U[j] <= 1 - m_xi[j]);
        description.make_follower_ctr(c);
    }

    model.use(Bilevel::MibS(description).with_logs(false));

    model.optimize();

    auto result = save_primal(model);

    result.set_objective_value(-result.objective_value() - t_first_stage_solution.get(m_theta));

    for (const auto& [var, val] : result) {
        if (var.name().find("xi") == 0) {
            std::cout << var.name() << " = " << val << std::endl;
        }
    }

    return result;

}
