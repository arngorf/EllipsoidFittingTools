#ifndef TASKWORKER_HPP
#define TASKWORKER_HPP

#include <iostream>
#include <thread>
#include <mutex>
#include <functional>
#include <exception>

template <typename TaskResult>
class Worker_Result_Wrap {
public:
    Worker_Result_Wrap(TaskResult result);
    TaskResult result;
};

template <typename TaskResult>
Worker_Result_Wrap<TaskResult>::Worker_Result_Wrap(TaskResult result) : result(result) {

}

template <typename TaskResult, typename TaskParam>
class TaskWorker {
public:
    TaskWorker(std::function<TaskResult(const TaskParam &)> const &F, std::vector<TaskParam> &parameters, std::vector<Worker_Result_Wrap<TaskResult>*> &results, std::mutex &mtx, int &currentId);
    void run();
    void join();

private:
    std::function<TaskResult(const TaskParam &)> const &F;
    std::vector<TaskParam> &parameters;
    std::vector<Worker_Result_Wrap<TaskResult>*> &results;
    std::mutex &mtx;
    int &currentId;
    const int batchSize;

    std::thread workerThread;

};

template <typename TaskResult, typename TaskParam>
TaskWorker<TaskResult, TaskParam>::TaskWorker(std::function<TaskResult(const TaskParam &)> const &F, std::vector<TaskParam> &parameters, std::vector<Worker_Result_Wrap<TaskResult>*> &results, std::mutex &mtx, int &currentId)
    : F(F), parameters(parameters), results(results), mtx(mtx), currentId(currentId), batchSize(parameters.size()) {
    workerThread = std::thread(&TaskWorker::run, this);
}

template <typename TaskResult, typename TaskParam>
void TaskWorker<TaskResult, TaskParam>::run() {

    int workerJobId;

    while (true) {

        mtx.lock();

        if (currentId >= batchSize) {
            mtx.unlock();
            break;
        }

        // Get new job to do here.
        workerJobId = currentId;
        ++currentId;
        mtx.unlock();

        // Do the job
        Worker_Result_Wrap<TaskResult> *resultWrap;
        TaskResult result = F(parameters[workerJobId]);
        resultWrap = new Worker_Result_Wrap<TaskResult>(result);

        //mtx.lock();
        results[workerJobId] = resultWrap;
        //mtx.unlock();
    }
}

template <typename TaskResult, typename TaskParam>
void TaskWorker<TaskResult, TaskParam>::join() {
    workerThread.join();
}

#endif // TASKWORKER_HPP