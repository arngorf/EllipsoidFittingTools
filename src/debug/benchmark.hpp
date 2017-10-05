/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include "../unistd.h" /* POSIX flags */
#include "../time.h"   /* clock_gettime(), time() */
#include "../sys/time.h"   /* gethrtime(), gettimeofday() */
#include <vector>
#include <iostream>

class KeeperOfTimes {
public:
    KeeperOfTimes();
    int addVar(std::string name);
    void addTime(int id, double time);
    void printTimes();

private:
    std::vector<std::string> names;
    std::vector<double> times;
};

/**
 * Returns the real time, in seconds, or -1.0 if an error occurred.
 *
 * Time is measured since an arbitrary and OS-dependent start time.
 * The returned real time is only useful for computing an elapsed time
 * between two calls to this function.
 */
double getRealTime();
double getCPUTime();
#endif