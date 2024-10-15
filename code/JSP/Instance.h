//
// Created by henri on 15.08.23.
//

#ifndef CCG_WITH_NESTED_CG_INSTANCE_H
#define CCG_WITH_NESTED_CG_INSTANCE_H


#include <vector>
#include <string>

struct Job {
    const unsigned int index;
    double release_date = 0;
    double deadline = 0;
    double processing_time = 0;
    double weight = 0;
    double profit = 0;
    double outsourcing_cost = 0;
public:
    explicit Job(unsigned int t_index);

    Job(const Job&) = default;
    Job(Job&&) = default;

    Job& operator=(const Job&) = delete;
    Job& operator=(Job&&) = delete;
};

struct JobOccurrence {
    const Job* parent;
    double release_date;
    double deadline;

    JobOccurrence(const Job& t_parent, int t_release_date, int t_deadline);

    JobOccurrence(const JobOccurrence&) = default;
    JobOccurrence(JobOccurrence&&) = default;

    JobOccurrence& operator=(const JobOccurrence&) = default;
    JobOccurrence& operator=(JobOccurrence&&) = default;
};

class Instance {
    std::vector<Job> m_jobs;
public:
    explicit Instance(unsigned int t_n_jobs);

    [[nodiscard]] unsigned int n_jobs() const { return m_jobs.size(); }

    Job& job(unsigned int t_index) { return m_jobs.at(t_index); }

    [[nodiscard]] const Job& job(unsigned int t_index) const { return m_jobs.at(t_index); }

    [[nodiscard]] std::vector<JobOccurrence> compute_job_occurrences() const;

    static std::vector<double> compute_big_M(const std::vector<JobOccurrence>& t_job_occurrence);

    static Instance from_file(const std::string& t_filename);

    static Instance generate(unsigned int t_n_jobs, unsigned int t_R, unsigned int t_D);
};

std::ostream &operator<<(std::ostream &t_os, const Instance& t_instance);

#endif //CCG_WITH_NESTED_CG_INSTANCE_H
