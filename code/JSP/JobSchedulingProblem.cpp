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
#include "idol/optimizers/bilevel-optimization/wrappers/MibS/MibS.h"

JobSchedulingProblem::JobSchedulingProblem(const Instance &t_instance, double t_Gamma)
    : m_instance(t_instance),
      m_Gamma(t_Gamma),
      m_job_occurrences(m_instance.compute_job_occurrences()),
      m_big_M(Instance::compute_big_M(m_job_occurrences)),
      m_annotation(m_env, "decomposition", MasterId),
      m_theta(m_env, -Inf, Inf, Continuous, "theta"),
      m_x(idol::Var::make_vector<1>(m_env, Dim<1>(m_instance.n_jobs()), 0, 1, Binary, "x")),
      m_xi(idol::Var::make_vector<1>(m_env, Dim<1>(m_instance.n_jobs()), 0, 1, Binary, "xi"))
    {

}

Model JobSchedulingProblem::create_master_problem() {

    const auto n_jobs = m_instance.n_jobs();

    idol::Model result(m_env);

    result.add(m_theta);
    result.add_vector<idol::Var, 1>(m_x);

    result.set_obj_expr(m_theta + idol_Sum(j, Range(n_jobs), m_instance.job(j).weight * m_x[j]));

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

    const auto n_jobs = m_instance.n_jobs();
    const auto n_job_occurrences = m_job_occurrences.size();

    const auto is_attacked = [&](unsigned int t_job_index) -> bool {
        return t_worst_case_scenario.get(m_xi[t_job_index]) > .5;
    };

    auto y = idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, 1, Binary, "y");
    auto t = idol::Var::make_vector(m_env, Dim<1>(n_job_occurrences), 0, Inf, Continuous, "t");

    if (t_master.optimizer().is<idol::Optimizers::BranchAndBound<DefaultNodeInfo>>()) {

        auto& branch_and_bound = t_master.optimizer().as<idol::Optimizers::BranchAndBound<DefaultNodeInfo>>();

        if (branch_and_bound.relaxation().optimizer().is<idol::Optimizers::DantzigWolfeDecomposition>()) {

            const auto &dantzig_wolfe = branch_and_bound.relaxation().optimizer().as<idol::Optimizers::DantzigWolfeDecomposition>();
            const auto &var_annotation = dantzig_wolfe.formulation().decomposition_by_variable();

            for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
                y[k].set(var_annotation, t_iteration);
                t[k].set(var_annotation, t_iteration);
            }

        }

    }

    t_master.add_vector<Var, 1>(y);
    t_master.add_vector<Var, 1>(t);

    // Compute sum of y_k over G_j for each job j
    std::vector<Expr<Var, Var>> sum_y_k(n_jobs);
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        sum_y_k[m_job_occurrences[k].parent->index] += y[k];
    }

    // Objective function
    t_master.add_ctr(m_theta >= idol_Sum(j, Range(n_jobs), - (m_instance.job(j).weight + m_instance.job(j).profit) * sum_y_k[j]));

    // Linking constraints
    for (unsigned int j = 0 ; j < n_jobs ; ++j) {
        t_master.add_ctr(sum_y_k[j] <= m_x[j]);
    }

    // GUB constraint
    for (unsigned int j = 0 ; j < n_jobs ; ++j) {
        Ctr c(m_env, sum_y_k[j] <= 1);
        c.set(m_annotation, t_iteration);
        t_master.add(c);
    }

    // Deadlines
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        Ctr c(m_env, t[k] <= m_job_occurrences[k].deadline);
        c.set(m_annotation, t_iteration);
        t_master.add(c);
    }

    // Non-overlapping
    for (unsigned int k = 1 ; k < n_job_occurrences ; ++k) {

        const auto& job_occurrence = m_job_occurrences[k];
        const double factor = 1 + is_attacked(job_occurrence.parent->index) * m_percentage_increase;
        const double processing_time = factor * job_occurrence.parent->processing_time;

        Ctr c(m_env, t[k] - t[k-1] - processing_time * y[k] >= 0);
        c.set(m_annotation, t_iteration);
        t_master.add(c);
    }

    // Disjunctive + release dates
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {

        const auto& job_occurrence = m_job_occurrences[k];
        const double factor = 1 + is_attacked(job_occurrence.parent->index) * m_percentage_increase;
        const double processing_time = factor * job_occurrence.parent->processing_time;

        Ctr c(m_env, t[k] - processing_time * y[k] - m_big_M[k] * y[k] >= job_occurrence.release_date - m_big_M[k]);
        c.set(m_annotation, t_iteration);
        t_master.add(c);
    }

}

Solution::Primal JobSchedulingProblem::compute_worst_case_scenario(const Model &t_master,
                                                                   const Solution::Primal &t_first_stage_solution,
                                                                   double t_time_limit) {
    const auto n_jobs = m_instance.n_jobs();
    const auto n_job_occurrences = m_job_occurrences.size();

    idol::Model model(m_env);

    // Create master problem for separation
    model.set_obj_sense(Maximize);
    model.add_vector<Var, 1>(m_xi);
    auto pi = model.add_var(-Inf, n_jobs, Continuous, Column(1), "pi");
    auto budget = model.add_ctr(idol_Sum(i, Range(n_jobs), m_xi[i]) <= m_Gamma);

    // Create separation problem
    Model separation(m_env);

    auto y = separation.add_vars(Dim<1>(n_job_occurrences), 0, 1, Binary, "y");
    auto z = separation.add_vars(Dim<1>(n_job_occurrences), 0, 1, Binary, "z");
    auto t = separation.add_vars(Dim<1>(n_job_occurrences), 0, Inf, Continuous, "t");

    // Compute sum of y_k over G_j for each job j
    std::vector<Expr<Var, Var>> sum_y_k(n_jobs);
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        sum_y_k[m_job_occurrences[k].parent->index] += y[k];
    }

    // Linking constraints
    for (unsigned int j = 0 ; j < n_jobs ; ++j) {
        separation.add_ctr(sum_y_k[j] <= std::round(t_first_stage_solution.get(m_x[j])));
    }

    // Deadlines
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        separation.add_ctr(t[k] <= m_job_occurrences[k].deadline);
    }

    // Packing
    for (unsigned int k = 1 ; k < n_job_occurrences ; ++k) {
        const auto& job_occurrence = m_job_occurrences[k];
        const double p_k = job_occurrence.parent->processing_time;
        const double tau_k = m_percentage_increase * p_k;
        separation.add_ctr(t[k] - t[k-1] - p_k * y[k] - tau_k * z[k] >= 0);
    }

    // Disjunctive + release
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        const auto& job_occurrence = m_job_occurrences[k];
        const double r_k = job_occurrence.parent->release_date;
        const double p_k = job_occurrence.parent->processing_time;
        const double tau_k = m_percentage_increase * p_k;
        separation.add_ctr(t[k] - p_k * y[k] - tau_k * z[k] - m_big_M[k] * y[k] >= r_k - m_big_M[k]);
    }

    // z <= y
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        separation.add_ctr(z[k] <= y[k]);
    }

    // Compute sum of z_k over G_j for each job j
    std::vector<Constant> sum_z_k_constant(n_jobs);
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        sum_z_k_constant[m_job_occurrences[k].parent->index] += !z[k];
    }

    // Compute sum of y_k over G_j for each job j
    std::vector<Constant> sum_y_k_constant(n_jobs);
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        sum_y_k_constant[m_job_occurrences[k].parent->index] += !y[k];
    }

    Expr<Var, Var> rhs;
    for (unsigned int j = 0 ; j < n_jobs ; ++j) {

        const auto& job = m_instance.job(j);
        const auto cost = - (job.weight + job.profit);

        rhs += cost * sum_y_k_constant[j];
        rhs -= cost * (m_xi[j] * sum_y_k_constant[j]);
        rhs += cost * (m_xi[j] * sum_z_k_constant[j]);

    }

    auto cut = pi <= rhs;

    model.use(
            create_gurobi()
                    .with_lazy_cut(true)
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

    auto result = save_primal(model);

    result.set_objective_value(result.objective_value() - t_first_stage_solution.get(m_theta));

    return result;

}
