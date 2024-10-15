//
// Created by henri on 15.08.23.
//

#include <list>
#include <algorithm>
#include <limits>
#include <fstream>
#include <random>
#include <iomanip>
#include <iostream>
#include "Instance.h"

Instance::Instance(unsigned int t_n_jobs) {

    m_jobs.reserve(t_n_jobs);
    for (unsigned int i = 0 ; i < t_n_jobs ; ++i) {
        m_jobs.emplace_back(i);
    }

}

std::vector<JobOccurrence> Instance::compute_job_occurrences() const {

    // Create job occurrences
    std::list<JobOccurrence> tmp_result;
    for (unsigned long int i = 0, n_jobs = m_jobs.size() ; i < n_jobs ; ++i) {

        const Job& job_i = m_jobs.at(i);

        for (unsigned long int j = i + 1 ; j < n_jobs ; ++j) {

            const Job& job_j = m_jobs.at(j);

            if (job_i.deadline > job_j.deadline && job_i.release_date + job_i.processing_time + job_j.processing_time < job_j.deadline) {

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
    double placeholder;

    file >> n_jobs;

    Instance result(n_jobs);

    for (unsigned int i = 0 ; i < n_jobs ; ++i) {

        auto& job = result.job(i);

        file >> placeholder;
        job.weight = placeholder;

        file >> placeholder;
        job.release_date = placeholder;

        file >> placeholder;
        job.deadline = placeholder;

        file >> placeholder;
        job.processing_time = placeholder;

        file >> placeholder;
        job.profit = placeholder;

        file >> placeholder;
        job.outsourcing_cost = placeholder;
    }

    return result;
}

Instance Instance::generate(unsigned int t_n_jobs, unsigned int t_R, unsigned int t_D) {

    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<double> processing_time_distribution(1,100);
    std::uniform_real_distribution<double> weight_distribution(1,100);
    std::uniform_real_distribution<double> profit_distribution(1,100);
    std::uniform_real_distribution<double> outsourcing_cost_distribution(1,100);
    std::uniform_real_distribution<double> release_date_distribution(0,t_n_jobs * t_R);

    Instance result(t_n_jobs);

    for (unsigned int j = 0 ; j < t_n_jobs ; ++j) {

        auto& job = result.job(j);

        job.processing_time = processing_time_distribution(generator);
        job.weight = weight_distribution(generator);
        job.release_date = release_date_distribution(generator);
        job.profit = profit_distribution(generator);
        job.outsourcing_cost = outsourcing_cost_distribution(generator);

        std::uniform_real_distribution<double> deadline_distribution(0,t_n_jobs * t_D);

        job.deadline = job.release_date + job.processing_time + deadline_distribution(generator);

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

Job::Job(unsigned int t_index) : index(t_index) {

}

std::ostream &operator<<(std::ostream &t_os, const Instance& t_instance) {

    const unsigned int n_jobs = t_instance.n_jobs();

    t_os << n_jobs << '\n';

    t_os << std::setprecision(10);

    for (unsigned int i = 0 ; i < n_jobs ; ++i) {
        const auto& job = t_instance.job(i);
        t_os << job.weight << '\t';
        t_os << job.release_date << '\t';
        t_os << job.deadline << '\t';
        t_os << job.processing_time << '\t';
        t_os << job.profit << '\t';
        t_os << job.outsourcing_cost << '\n';
    }

    return t_os;
}