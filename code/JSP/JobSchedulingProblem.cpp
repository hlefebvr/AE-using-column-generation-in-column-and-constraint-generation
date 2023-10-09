//
// Created by henri on 15.08.23.
//

#include "JobSchedulingProblem.h"
#include "optimizers/dantzig-wolfe/DantzigWolfeDecomposition.h"
#include "optimizers/branch-and-bound/BranchAndBound.h"
#include "optimizers/branch-and-bound/node-selection-rules/factories/DepthFirst.h"
#include "optimizers/column-generation/IntegerMasterHeuristic.h"
#include "optimizers/branch-and-bound/branching-rules/factories/MostInfeasible.h"
#include "optimizers/dantzig-wolfe/Optimizers_DantzigWolfeDecomposition.h"
#include "optimizers/callbacks/LazyCutCallback.h"

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

    t_master.use(create_gurobi()
                         .with_log_level(Info, Black)
    );

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
                                    .with_pricing_optimizer(
                                            create_gurobi()
                                                    .with_max_n_solution_in_pool(20)
                                                    .with_infeasible_or_unbounded_info(true)
                                                    //.with_log_level(Info, Black)
                                    )
                                    .with_farkas_pricing(true)
                                    .with_branching_on_master(false)
                                    .with_parallel_pricing_limit(m_parallel_pricing)
                                    //.with_column_pool_clean_up(300, .33)
                                    .with_max_columns_per_pricing(20)
                                    .conditional(m_use_non_optimal_pricing, [this](auto& x){
                                        x.with_non_optimal_pricing_phase(10, .1);
                                    }, [this](auto& x) {
                                        x.with_dual_price_smoothing_stabilization(.4);
                                    })
                                    .with_log_frequency(1)
                                    .with_log_level(Info, Yellow)
                    )
                    .with_log_frequency(1)
                    .with_subtree_depth(1)
                    .with_log_level(Info, Blue)
                    .with_branching_rule(MostInfeasible(m_x.begin(), m_x.end()))
                    .with_node_selection_rule(DepthFirst())
                    .conditional(m_use_heuristic, [this](auto& x){
                        x.with_callback(
                                IntegerMasterHeuristic()
                                        .with_optimizer(create_gurobi())
                                        .with_integer_columns(true)
                        );
                    })
                    .with_time_limit(time_limit)
    );

}

Gurobi JobSchedulingProblem::create_gurobi() const {
    return std::move(Gurobi()
            .with_presolve(false)
            .with_thread_limit(1)
            .with_external_param(GRB_DoubleParam_Heuristics, 0.)
            .with_external_param(GRB_IntParam_Cuts, 0)
            .with_external_param(GRB_IntParam_RINS, 0)
            .with_external_param(GRB_IntParam_Aggregate, 0)
            .with_external_param(GRB_IntParam_Symmetry, 0)
            .with_external_param(GRB_IntParam_Disconnected, 0)
            .with_external_param(GRB_IntParam_VarBranch, 2)
    );
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

    if (t_master.optimizer().is<idol::Optimizers::BranchAndBound<NodeInfo>>()) {

        auto& branch_and_bound = t_master.optimizer().as<idol::Optimizers::BranchAndBound<NodeInfo>>();

        if (branch_and_bound.relaxation().optimizer().is<idol::Optimizers::DantzigWolfeDecomposition>()) {

            const auto &dantzig_wolfe = branch_and_bound.relaxation().optimizer().as<idol::Optimizers::DantzigWolfeDecomposition>();
            const auto &var_annotation = dantzig_wolfe.var_annotation();

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
        separation.add_ctr(sum_y_k[j] <= t_first_stage_solution.get(m_x[j]));
    }

    // Deadlines
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        separation.set_var_ub(t[k], m_job_occurrences[k].deadline);
    }

    // Packing
    for (unsigned int k = 1 ; k < n_job_occurrences ; ++k) {
        const auto& job_occurrence = m_job_occurrences[k];
        const int p_k = job_occurrence.parent->processing_time;
        const int tau_k = m_percentage_increase * p_k;
        separation.add_ctr(t[k] - t[k-1] - p_k * y[k] - tau_k * z[k] >= 0);
    }

    // Disjunctive + release
    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        const auto& job_occurrence = m_job_occurrences[k];
        const int r_k = job_occurrence.parent->release_date;
        const int p_k = job_occurrence.parent->processing_time;
        const int tau_k = m_percentage_increase * p_k;
        separation.add_ctr(t[k] - p_k * y[k] - tau_k * z[k] >= r_k - m_big_M[k]);
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
                    .with_callback(
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
