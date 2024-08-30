//
// Created by henri on 15.08.23.
//
#include <iostream>
#include "Instance.h"
#include "idol/errors/Exception.h"
#include "JobSchedulingProblem.h"
#include "../algorithms/ColumnAndConstraintGeneration.h"

int main(int t_argc, const char** t_argv) {

    using namespace idol;

    if (t_argc != 5) {
        throw Exception("Parameters: <instance_file> <std_phase_time_limit> <gamma> <with_heuristic>");
    }

    const std::string filename = t_argv[1];
    const double std_phase_time_limit = std::stof(t_argv[2]);
    const double Gamma = std::stof(t_argv[3]);
    const bool with_heuristic = std::string_view (t_argv[4]) == "true";

    auto instance = Instance::from_file(filename);

    std::cout << "Parsed " << filename << "." << std::endl;

    JobSchedulingProblem problem(instance, Gamma);

    ColumnAndConstraintGeneration ccg(problem, 10800);

    problem.use_parallel_pricing(1);

    problem.use_heuristic(with_heuristic);

    ccg.solve(std_phase_time_limit, filename.substr(filename.find_last_of('/') + 1));

    std::cout << with_heuristic << ','
              << instance.n_jobs() << ','
              << problem.Gamma() << ','
              << std::endl;



    return 0;
}
