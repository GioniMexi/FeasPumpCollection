/**
 * @file timer.h
 * @brief CPU Stopwatch for benchmarks
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2019
 */

#ifndef TIMER_H
#define TIMER_H

#include "asserter.h"
#include <string>
#include <chrono>


namespace dominiqs {

/**
 * Benchmarking class
 * Simulates a simple stopwatch for wallclock time
 */

class StopWatch
{
public:
	/** default constructor */
	StopWatch(bool autoStart = false);
	/** start stopwatch */
	void start();
	/** stop stopwatch */
	void stop();
	/** reset stopwatch */
	void reset();

	/** @return partial elapsed time in seconds */
	double getPartial() const;
	/** @return total elapsed time in seconds */
	double getTotal() const;
	/** @return elapsed time in seconds */
	double getElapsed() const;
private:
	std::chrono::high_resolution_clock::time_point w_begin;
	std::chrono::high_resolution_clock::time_point w_end;
	double w_total;
};

inline StopWatch::StopWatch(bool autoStart) : w_total(0.0)
{
	if (autoStart) start();
}

inline void StopWatch::start()
{
	w_begin = std::chrono::high_resolution_clock::now();
}

inline void StopWatch::stop()
{
	w_end = std::chrono::high_resolution_clock::now();
	w_total += getPartial();
}

inline void StopWatch::reset()
{
	w_total = 0.0;
}

inline double StopWatch::getPartial() const
{
	return std::chrono::duration<double>(w_end - w_begin).count();
}

inline double StopWatch::getTotal() const
{
	return w_total;
}

inline double StopWatch::getElapsed() const
{
	auto w_now = std::chrono::high_resolution_clock::now();
	return std::chrono::duration<double>(w_now - w_begin).count();
}

/**
 * Global Stopwatch
 */

StopWatch& gStopWatch();

/**
 * Automatic stopwatch stopper (RAII principle)
 */

class AutoStopWatch
{
public:
	AutoStopWatch(StopWatch* c) : theStopWatch(c), pending(true)
	{
		DOMINIQS_ASSERT( theStopWatch );
		theStopWatch->start();
	}
	~AutoStopWatch() { stop(); }
	void stop()
	{
		if (pending)
		{
			theStopWatch->stop();
			pending = false;
		}
	}
private:
	StopWatch* theStopWatch;
	bool pending;
};

/**
 * Automatic global stopwatch stopper (RAII principle)
 */

class GlobalAutoStopWatch
{
public:
	GlobalAutoStopWatch() : pending(true) { gStopWatch().start(); }
	~GlobalAutoStopWatch() { stop(); }
	void stop()
	{
		if (pending)
		{
			gStopWatch().stop();
			pending = false;
		}
	}
private:
	bool pending;
};

/**
 * Get current date/time as std::string, format is YYYY-MM-DD.HH:mm:ss
 */
std::string currentDateTime();

} // namespace dominiqs

#endif /* TIMER_H */
