//
// Created by henri on 09/03/23.
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include "Instance.h"

Instance::Instance(idol::Problems::FLP::Instance &&t_instance_base)
        : idol::Problems::FLP::Instance(std::move(t_instance_base)) {

    compute_profits();

}

void Instance::compute_profits() {

    const unsigned int n_facilities = this->n_facilities();
    const unsigned int n_customers = this->n_customers();

    m_profits.clear();
    m_profits.resize(n_customers);

    for (unsigned int j = 0 ; j < n_customers ; ++j) {

        std::vector<double> costs;
        costs.reserve(n_facilities);
        for (unsigned int i = 0 ; i < n_facilities ; ++i) {
            costs.emplace_back(per_unit_transportation_cost(i,j));
        }

        std::sort(costs.begin(), costs.end());

        const double median_cost = costs[ std::floor(costs.size() / 2.) ];

        set_profit(j, idol::Problems::FLP::Instance::demand(j) * median_cost * 4);

    }

}

double Instance::fixed_cost(unsigned int t_i) const {
    return std::floor(m_scaling_factor * idol::Problems::FLP::Instance::fixed_cost(t_i));
}

double Instance::capacity(unsigned int t_i) const {
    return std::floor(m_scaling_factor * idol::Problems::FLP::Instance::capacity(t_i));
}

double Instance::demand(unsigned int t_j) const {
    return std::floor(m_scaling_factor * idol::Problems::FLP::Instance::demand(t_j));
}
