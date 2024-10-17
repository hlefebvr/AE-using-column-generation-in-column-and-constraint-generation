//
// Created by henri on 16.08.23.
//
#include <iostream>
#include <fstream>
#include "Instance.h"

int main(int t_argc, const char** t_argv) {

    const std::string destination_folder = "/home/henri/Research/AE-using-column-generation-in-column-and-constraint-generation/code/JSP/data";
    const auto n_instances_per_config = 10;
    const auto n_values = { 15, 20 };
    const auto R_values = { 1, 5, 10, 20 };
    const auto D_values = { 1, 5, 10, 20 };

    for (auto N : n_values) {
        for (auto R : R_values) {
            for (auto D : D_values) {

                for (unsigned int i = 0 ; i < n_instances_per_config ; ++i) {

                    auto instance = Instance::generate(N, R, D);

                    auto filename = "instance_N" + std::to_string(N) + "_R" + std::to_string(R) + "_D" + std::to_string(D) + "__" + std::to_string(i) + ".txt";

                    std::ofstream file(destination_folder + "/" + filename);

                    if (file.fail()) {
                        throw std::runtime_error("unable to open destination file.");
                    }

                    file << instance << "\n\n\n# Generated with parameters N = " << N << ", R = " << R << " and D = "
                         << D;

                    file.close();

                }

            }
        }
    }

    return 0;
}
