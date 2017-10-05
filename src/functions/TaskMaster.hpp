#ifndef TASKMASTER_HPP
#define TASKMASTER_HPP

#include <iostream>
#include <thread>
#include <mutex>
#include <functional>
#include <exception>

#include "TaskWorker.hpp"

template <typename TaskResult, typename TaskParam>
class TaskMaster {
public:
    TaskMaster();
    void runJobs(std::function<TaskResult(const TaskParam &)> const& F, std::vector<TaskParam> parameters);
    std::vector<TaskResult> getResults();

private:
    int threadPoolSize;
    bool running;
    int batchSize;
    std::vector<TaskWorker<TaskResult, TaskParam>*> workers;

    int currentId;
    std::mutex mtx;
    std::vector<TaskParam> params;
    std::vector<Worker_Result_Wrap<TaskResult>*> resultPtrs;

};

template <typename TaskResult, typename TaskParam>
TaskMaster<TaskResult, TaskParam>::TaskMaster()
    : threadPoolSize(8), running(false) {

}

template <typename TaskResult, typename TaskParam>
void TaskMaster<TaskResult, TaskParam>::runJobs(std::function<TaskResult(const TaskParam &)> const& F, std::vector<TaskParam> parameters) {

    if (running) {
        throw std::runtime_error("Task master already running");
    }

    running = true;
    batchSize = parameters.size();
    currentId = 0;
    params = parameters;

    workers.resize(0);
    resultPtrs.resize(batchSize);

    for (int i = 0; i < threadPoolSize; ++i) {
        workers.push_back(new TaskWorker<TaskResult, TaskParam>(F, params, resultPtrs, mtx, currentId));
    }
}

template <typename TaskResult, typename TaskParam>
std::vector<TaskResult> TaskMaster<TaskResult, TaskParam>::getResults() {

    if (!running) {
        throw std::runtime_error("Task master is not running");
    }

    running = false;

    for (int i = 0; i < threadPoolSize; ++i) {
        workers[i]->join();
    }

    std::vector<TaskResult> results;
    results.reserve(batchSize);
    for (int i = 0; i < batchSize; ++i) {
        results.push_back(resultPtrs[i]->result);
        delete resultPtrs[i];
    }

    return results;
}

#endif // TASKMASTER_HPP