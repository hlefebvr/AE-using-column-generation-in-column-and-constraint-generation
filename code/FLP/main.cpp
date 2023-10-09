//
// Created by henri on 09/03/23.
//
#include <iostream>
#include "FacilityLocationProblem.h"
#include "../algorithms/ColumnAndConstraintGeneration.h"

#ifdef WITH_DEBUG_INFO
#include <cstdio>
#include <execinfo.h>
#include <csignal>
#include <cstdlib>
#include <unistd.h>

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 20);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}
#endif

int main(int t_argc, char** t_argv) {

#ifdef WITH_DEBUG_INFO
    signal(SIGSEGV, handler);
    signal(SIGABRT, handler);
#endif

    if (t_argc != 5) {
        throw Exception("Parameters: <instance_file> <std_phase_time_limit> <gamma> <with_heuristic>");
    }

    const std::string filename = t_argv[1];
    const double std_phase_time_limit = std::stof(t_argv[2]);
    const double Gamma = std::stof(t_argv[3]);
    const bool with_heuristic = std::string_view (t_argv[4]) == "true";

    const Instance instance(Problems::FLP::read_instance_1991_Cornuejols_et_al(filename));

    FacilityLocationProblem problem(instance, Gamma);

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

    // Root::run();

    return 0;
}
