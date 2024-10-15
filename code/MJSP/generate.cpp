//
// Created by henri on 16.08.23.
//
#include <iostream>
#include <fstream>
#include "Instance.h"

int main(int t_argc, const char** t_argv) {

    const std::string destination_folder = "/home/henri/Research/AE-using-column-generation-in-column-and-constraint-generation/code/MJSP/data";
    const auto n_instances_per_config = 5;
    const auto n_values = { 20, 30, 40 };
    const auto m_values = { 3, 4, 5 };
    const auto R_values = { 1, 5, 10, 20 };
    const auto D_values = { 1, 5, 10, 20 };

    for (auto M : m_values) {
        for (auto N: n_values) {
            for (auto R: R_values) {
                for (auto D: D_values) {

                    for (unsigned int q = 0; q < n_instances_per_config; ++q) {

                        auto instance = Instance::generate(N, M, R, D);

                        auto filename =
                                "instance_M" + std::to_string(M) +  "_N" + std::to_string(N) + "_R" + std::to_string(R) + "_D" + std::to_string(D) +
                                "__" + std::to_string(q) + ".txt";

                        std::ofstream file(destination_folder + "/" + filename);

                        if (file.fail()) {
                            throw std::runtime_error("unable to open destination file.");
                        }

                        file << instance << "\n\n\n# Generated with parameters M = " << M << ", N = " << N << ", R = " << R
                             << " and D = "
                             << D;

                        file.close();

                    }

                }
            }
        }
    }

    return 0;
}
