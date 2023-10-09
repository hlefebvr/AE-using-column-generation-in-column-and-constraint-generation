//
// Created by henri on 07/04/23.
//
#include "ColumnAndConstraintGeneration.h"

ColumnAndConstraintGeneration::ColumnAndConstraintGeneration(ColumnAndConstraintGenerator &t_generator, double t_time_limit)
    : m_generator(t_generator), m_time_limit(t_time_limit) {

}

idol::Solution::Primal ColumnAndConstraintGeneration::solve(double t_std_phase_time_limit, const std::string& t_tag) {

    idol::Solution::Primal first_stage_solution;

    enum Event { MasterSolved, AdversarialSolved };

    auto master_problem = m_generator.create_master_problem();
    const auto sense = master_problem.get_obj_sense();

    std::optional<SolutionStatus> adversarial_unexpected_status;
    Timer total_timer, master_timer, adversarial_timer;
    double best_obj = sense == Minimize ? Inf : -Inf;
    double best_bound = sense == Minimize ? -Inf : Inf;
    unsigned int iteration = 0;
    bool large_scale_phase = false;

    const auto log = [&](Event t_event) {
        // if (t_tag == "inner") { return; }
        std::cout
                << "<Tag=" << t_tag << "> "
                << "<Iteration=" << iteration << "> "
                << "<Type=" << t_event << "> "
                << "<TimeTotal=" << std::setw(10) << total_timer.count() << "> "
                << "<TimeMastr=" << std::setw(10) << master_timer.count() << "> "
                << "<TimeAdver=" << std::setw(10) << (t_event == AdversarialSolved ? adversarial_timer.count() : 0) << "> "
                << "<Status=" << std::setw(10) << master_problem.get_status() << "> "
                << "<Status=" << std::setw(10) << master_problem.get_reason() << "> "
                << "<BestBound=" << std::setw(10) << best_bound << "> "
                << "<BestObj="   << std::setw(10) << best_obj << "> "
                << "<RelGap="    << std::setw(10) << relative_gap(best_bound, best_obj) * 100 << "%> "
                << "<ABsGap="    << std::setw(10) << absolute_gap(best_bound, best_obj) << "> "
                << std::endl;
    };

    const auto stopping_condition = [&]() {
        return relative_gap(best_bound, best_obj) < 1e-5
                || total_timer.count() >= m_time_limit
                || master_problem.get_status() != Optimal;
    };

    const auto set_time_limit = [&]() {

        const double remaining_time = std::max(1e-3, m_time_limit - total_timer.count());

        if (!large_scale_phase) {
            master_problem.optimizer().set_param_time_limit(std::min(remaining_time, t_std_phase_time_limit));
        } else {
            master_problem.optimizer().set_param_time_limit(remaining_time);
        }

    };

    m_generator.set_default_optimizer(master_problem);

    Solution::Primal initial_scenario = m_generator.compute_initial_scenario();
    m_generator.add_scenario_to_master_problem(master_problem, initial_scenario, iteration);

    ++iteration;

    total_timer.start();
    do {

        //master_problem.optimizer().set_best_bound_stop(best_obj);

        set_time_limit();

        master_timer.start();

        master_problem.optimize();

        if (master_problem.get_status() != Optimal && master_problem.get_reason() == TimeLimit) {

            std::cout << "******** LARGE SCALE PHASE ********" << std::endl;

            large_scale_phase = true;
            m_generator.set_large_scale_optimizer(master_problem);

            set_time_limit();
            //master_problem.optimizer().set_best_bound_stop(master_problem.get_best_obj());
            master_problem.optimize();

        }

        master_timer.stop();

        if (stopping_condition()) { break; }

        best_bound = master_problem.get_best_obj();
        first_stage_solution = save_primal(master_problem);

        if (master_problem.get_status() == Feasible && master_problem.get_reason() == ObjLimit) {
            best_bound = best_obj;
        }

        log(MasterSolved);

        adversarial_timer.start();
        const auto worst_case_scenario = m_generator.compute_worst_case_scenario(master_problem, first_stage_solution, std::max(.0, m_time_limit - total_timer.count()));
        adversarial_timer.stop();

        if (worst_case_scenario.status() != Optimal) {
            adversarial_unexpected_status = worst_case_scenario.status();
            break;
        }

        const double iter_obj = first_stage_solution.objective_value() + worst_case_scenario.objective_value();

        if (sense == Minimize) {
            best_obj = std::min(best_obj, iter_obj);
        } else {
            best_obj = std::max(best_obj, iter_obj);
        }

        if (stopping_condition()) { break; }

        m_generator.add_scenario_to_master_problem(master_problem, worst_case_scenario, iteration);

        log(AdversarialSolved);

        ++iteration;

    } while (true);
    total_timer.stop();

    std::cout << "result,"
              << t_tag << ','
              << (t_std_phase_time_limit < m_time_limit ? std::to_string(t_std_phase_time_limit) : "inf") << ','
              << (t_std_phase_time_limit < m_time_limit ? "CG" : "STD") << ','
              << master_problem.get_status() << ','
              << master_problem.get_reason() << ','
              << large_scale_phase << ','
              << iteration << ','
              << total_timer.count() << ','
              << master_timer.cumulative_count() << ','
              << adversarial_timer.cumulative_count() << ','
              << best_bound << ','
              << best_obj << ','
              << relative_gap(best_bound, best_obj) * 100 << ','
              << absolute_gap(best_bound, best_obj) << ','
              ;

    if (adversarial_unexpected_status.has_value()) {
        std::cout << adversarial_unexpected_status.value();
    }
    std::cout << ',';

    first_stage_solution.set_objective_value(best_obj);

    return first_stage_solution;
}
