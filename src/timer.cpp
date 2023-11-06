#include "timer.hpp"
#include<iomanip>

Timer::Timer()
{
    start();
}

void Timer::start()
{
    st = std::chrono::high_resolution_clock::now();
}

void Timer::stop_and_log(char const *label)
{
    ed = std::chrono::high_resolution_clock::now();
    elapsed = ed - st;
    t = elapsed.count();

    std::cout << "\n" << label << ":\n";
    std::cout << "    time elapsed : " << std::fixed << std::setprecision(3) << t << "\n";
}