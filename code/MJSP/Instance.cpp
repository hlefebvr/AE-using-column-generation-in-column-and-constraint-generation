//
// Created by henri on 15.08.23.
//

#include <list>
#include <algorithm>
#include <limits>
#include <fstream>
#include <random>
#include <iomanip>
#include <cassert>
#include <iostream>
#include "Instance.h"

Instance::Instance(unsigned int t_n_jobs, unsigned int t_n_machines) {

    m_jobs.reserve(t_n_jobs);
    for (unsigned int i = 0 ; i < t_n_jobs ; ++i) {
        m_jobs.emplace_back(i, t_n_machines);
    }

    m_machines.reserve(t_n_machines);
    for (unsigned int i = 0 ; i < t_n_machines ; ++i) {
        m_machines.emplace_back(i);
    }

}

std::vector<std::vector<JobOccurrence>>
Instance::compute_job_occurrences() const {

    const unsigned int n_machines = this->n_machines();

    std::vector<std::vector<JobOccurrence>> result;
    result.reserve(n_machines);

    for (unsigned int i = 0 ; i < n_machines ; ++i) {
        result.emplace_back(compute_job_occurrences(i));
    }

    return result;
}

std::vector<JobOccurrence> Instance::compute_job_occurrences(unsigned int t_machine_id) const {

    // Create job occurrences
    std::list<JobOccurrence> tmp_result;
    for (unsigned long int i = 0, n_jobs = m_jobs.size() ; i < n_jobs ; ++i) {

        const Job& job_i = m_jobs.at(i);

        for (unsigned long int j = i + 1 ; j < n_jobs ; ++j) {

            const Job& job_j = m_jobs.at(j);

            if (job_i.deadline > job_j.deadline && job_i.release_date + job_i.processing_time[t_machine_id] + job_j.processing_time[t_machine_id] < job_j.deadline) {

                tmp_result.emplace_back(job_i, job_i.release_date, job_j.deadline);

            }

        }

        tmp_result.emplace_back(job_i, job_i.release_date, job_i.deadline);

    }

    // Copy job occurrences into vector
    std::vector<JobOccurrence> result;
    result.reserve(tmp_result.size());
    std::copy(tmp_result.begin(), tmp_result.end(), std::back_inserter(result));

    // Sort job occurrences by deadline
    std::sort(result.begin(), result.end(), [](const JobOccurrence& t_a, const JobOccurrence& t_b) { return t_a.deadline < t_b.deadline; });

    return result;
}

std::vector<double> Instance::compute_big_M(const std::vector<JobOccurrence> &t_job_occurrence) {

    const unsigned int n_job_occurrences = t_job_occurrence.size();

    std::vector<double> result;
    result.reserve(n_job_occurrences);

    for (unsigned int k = 0 ; k < n_job_occurrences ; ++k) {
        double min = std::numeric_limits<double>::max();
        for (unsigned int l = k + 1 ; l < n_job_occurrences ; ++l) {
            const double r_l = t_job_occurrence.at(l).release_date;
            if (min > r_l) {
                min = r_l;
            }
        }
        if (min == std::numeric_limits<double>::max()) {
            min = 0;
        }
        result.emplace_back(std::max<double>(min + 1, 1) - 1);
    }

    return result;
}

Instance Instance::from_file(const std::string &t_filename) {

    std::ifstream file(t_filename);

    if (file.fail()) {
        throw std::runtime_error("Unable to open instance file.");
    }

    unsigned int n_jobs;
    unsigned int n_machines;
    double placeholder;

    file >> n_machines;
    file >> n_jobs;

    Instance result(n_jobs, n_machines);

    for (unsigned int i = 0 ; i < n_machines ; ++i) {

        auto& machine = result.machine(i);

        file >> placeholder;
        machine.fixed_cost = placeholder;

    }

    for (unsigned int j = 0 ; j < n_jobs ; ++j) {

        auto& job = result.job(j);

        file >> placeholder;
        job.weight = placeholder;

        file >> placeholder;
        job.release_date = placeholder;

        file >> placeholder;
        job.deadline = placeholder;

        for (unsigned int i = 0 ; i < n_machines ; ++i) {
            file >> placeholder;
            job.processing_time[i] = placeholder;
        }

        for (unsigned int i = 0 ; i < n_machines ; ++i) {
            file >> placeholder;
            job.profit[i] = placeholder;
        }

    }

    return result;
}

Instance Instance::generate(unsigned int t_n_jobs, unsigned int t_n_machines, unsigned int t_R, unsigned int t_D) {

    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<double> processing_time_distribution(1,100);
    std::uniform_real_distribution<double> weight_distribution(1,100);
    std::uniform_real_distribution<double> profit_distribution(1,100);
    std::uniform_real_distribution<double> fixed_cost_distribution(1, 100);
    std::uniform_real_distribution<double> release_date_distribution(0,t_n_jobs * t_R);

    Instance result(t_n_jobs, t_n_machines);

    for (unsigned int i = 0 ; i < t_n_machines ; ++i) {
        auto& machine = result.machine(i);
        machine.fixed_cost = fixed_cost_distribution(generator);
    }

    for (unsigned int j = 0 ; j < t_n_jobs ; ++j) {

        auto& job = result.job(j);

        double max_processing_time = 0;
        for (unsigned int i = 0 ; i < t_n_machines ; ++i) {
            job.profit[i] = profit_distribution(generator);
            job.processing_time[i] = processing_time_distribution(generator);
            max_processing_time = std::max(max_processing_time, job.processing_time[i]);
        }
        job.weight = weight_distribution(generator);
        job.release_date = release_date_distribution(generator);

        std::uniform_real_distribution<double> deadline_distribution(0,t_n_jobs * t_D);

        job.deadline = job.release_date + max_processing_time + deadline_distribution(generator);

    }

    return result;
}

std::vector<std::vector<double>>
Instance::compute_big_M(const std::vector<std::vector<JobOccurrence>> &t_job_occurrence) {

    std::vector<std::vector<double>> result;
    result.reserve(t_job_occurrence.size());

    for (const auto& machine_job_occurrences : t_job_occurrence) {
        result.emplace_back(compute_big_M(machine_job_occurrences));
    }

    return result;
}

JobOccurrence::JobOccurrence(const Job &t_parent,
                             int t_release_date,
                             int t_deadline)
                             : parent(&t_parent),
                               release_date(t_release_date),
                               deadline(t_deadline)
                             {

}

Job::Job(unsigned int t_index, unsigned int t_n_machines)
    : index(t_index),
      processing_time(t_n_machines, 0),
      profit(t_n_machines, 0) {

}

std::ostream &operator<<(std::ostream &t_os, const std::vector<double>& t_vector) {

    for (const auto& element : t_vector) {
        t_os << element << '\t';
    }

    return t_os;
}

std::ostream &operator<<(std::ostream &t_os, const Instance& t_instance) {

    const unsigned int n_jobs = t_instance.n_jobs();
    const unsigned int n_machines = t_instance.n_machines();

    t_os << n_machines << '\n';
    t_os << n_jobs << '\n';

    t_os << std::setprecision(10);

    for (unsigned int i = 0 ; i < n_machines ; ++i) {
        const auto& machine = t_instance.machine(i);
        t_os << machine.fixed_cost << '\n';
    }

    for (unsigned int i = 0 ; i < n_jobs ; ++i) {
        const auto& job = t_instance.job(i);
        t_os << job.weight << '\t';
        t_os << job.release_date << '\t';
        t_os << job.deadline << '\t';
        t_os << job.processing_time << '\t';
        t_os << job.profit << '\t';
        t_os << '\n';
    }

    return t_os;
}