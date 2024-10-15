//
// Created by henri on 09/03/23.
//
#include <iostream>
#include "GeneralAssignmentProblem.h"
#include "../algorithms/ColumnAndConstraintGeneration.h"

int main(int t_argc, char** t_argv) {

    if (t_argc != 5) {
        throw Exception("Parameters: <instance_file> <std_phase_time_limit> <gamma> <with_heuristic>");
    }

    const std::string filename = t_argv[1];
    const double std_phase_time_limit = std::stof(t_argv[2]);
    const double Gamma = std::stof(t_argv[3]);
    const bool with_heuristic = std::string_view (t_argv[4]) == "true";

    const Instance instance(Problems::FLP::read_instance_1991_Cornuejols_et_al(filename));

    GeneralAssignmentProblem problem(instance, Gamma);

    ColumnAndConstraintGeneration ccg(problem, 10800);

    problem.use_parallel_pricing(1);

    problem.use_heuristic(with_heuristic);

	try {

    ccg.solve(std_phase_time_limit, filename);

	} catch (GRBException& err) {
		std::cout << err.getMessage() << std::endl;
		__throw_exception_again;
	}

    std::cout << with_heuristic << ','
              << instance.n_facilities() << ','
              << instance.n_customers() << ','
              << problem.Gamma() << ','
              << std::endl;

    return 0;
}
