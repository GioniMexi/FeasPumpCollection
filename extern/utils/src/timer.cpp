/**
 * @file timer.cpp
 * @brief CPU Stopwatch for benchmarks
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2019
 */

#include <ctime>
#include <cstdlib>

#include "utils/timer.h"

namespace dominiqs {

StopWatch& gStopWatch()
{
	static StopWatch theStopWatch;
	return theStopWatch;
}


std::string currentDateTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[100];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d %X (%z)", &tstruct);
	return buf;
}

} // namespace dominiqs
