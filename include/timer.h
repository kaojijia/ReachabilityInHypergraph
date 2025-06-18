

#ifndef DP2H_TIMER_H
#include <string>
#include <map>
#include <chrono>
#include <ratio>
#include <iostream>

namespace cx
{

    class Timer
    {
    public:
        std::map<std::string, std::chrono::high_resolution_clock::time_point> eventMap;
        std::map<std::string, unsigned long long> durationMap;

        void StartTimer(std::string eventName);
        unsigned long long EndTimer(std::string eventName);
        void EndTimerAndPrint(std::string eventName);
        void StopTimerAddDuration(std::string eventName);
        void PrintDuration(std::string eventName);
    };

}

#define DP2H_TIMER_H

#endif
