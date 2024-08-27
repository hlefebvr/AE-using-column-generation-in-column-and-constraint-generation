//
// Created by henri on 09/03/23.
//

#ifndef CCG_WITH_NESTED_CG_INSTANCE_H
#define CCG_WITH_NESTED_CG_INSTANCE_H

#include "idol/problems/facility-location-problem/FLP_Instance.h"
#include <cmath>

class Instance : public idol::Problems::FLP::Instance {
    double m_scaling_factor = 10;
    std::vector<double> m_profits;

    void compute_profits();
public:
    explicit Instance(idol::Problems::FLP::Instance&& t_instance_base);

    [[nodiscard]] double fixed_cost(unsigned int t_i) const override;

    [[nodiscard]] double capacity(unsigned int t_i) const override;

    [[nodiscard]] double demand(unsigned int t_j) const override;

    [[nodiscard]] double profit(unsigned int t_j) const { return std::floor(m_scaling_factor * m_profits[t_j]); }
    void set_profit(unsigned int t_j, double t_value) { m_profits[t_j] = t_value; }
};


#endif //CCG_WITH_NESTED_CG_INSTANCE_H
