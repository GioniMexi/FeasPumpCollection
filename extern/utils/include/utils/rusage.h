/**
 * @file rusage.h
 * @brief Resource Usage class
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2011 Domenico Salvagnin
 */

#ifndef RUSAGE_H
#define RUSAGE_H

#include <cstdlib>
#include <sys/resource.h>

namespace dominiqs {

/**
 * Get Average Load of the system
 */

inline double getLoadAvg()
{
	double res = 0.0;
	getloadavg(&res, 1);
	return res;
}

/**
 * Resource Usage Class
 */

class ResourceUsage
{
public:
	/** default constructor */
	ResourceUsage(bool autoStart = false);
	/** start */
	void start();
	/** stop */
	void stop();
	/** reset */
	void reset();
	/** get total user time in seconds */
	double getUserTime() const;
	/** get total system time in seconds */
	double getSysTime() const;
	/** get maximum occupied memory in MB */
	double getMaxMemory() const;
	/** get number of pagefaults without I/O */
 	int getMinorPageFaults() const;
	/** get number of pagefaults with I/O */
 	int getMajorPageFaults() const;
	/** get number of times process was completely swapped out */
	int getSwapCount() const;
private:
	// rusage structs
	rusage begin, end;
	// resources
	double userTime; //< user time (in seconds)
	double sysTime; //< system time (in seconds)
	double maxMemory; //< max memory (in MB)
	int pageFaults; //< # of page faults without I/O
	int pageFaultsIO; //< # of page faults with I/O
	int numSwap; //< # of times process was completely swapped out
	// constants
	static constexpr double KB_PER_MB = 1024.0;
	static constexpr double MICROSEC_TO_SEC = 1.0e-6;
};

ResourceUsage::ResourceUsage(bool autoStart)
{
	reset();
	if (autoStart) start();
}

inline void ResourceUsage::start()
{
	getrusage(RUSAGE_SELF, &begin);
}

inline void ResourceUsage::stop()
{
	getrusage(RUSAGE_SELF, &end);
	// time
	double t1 = ((double)begin.ru_utime.tv_sec + MICROSEC_TO_SEC*((double)begin.ru_utime.tv_usec));
	double t2 = ((double)end.ru_utime.tv_sec + MICROSEC_TO_SEC*((double)end.ru_utime.tv_usec));
	userTime += (t2 - t1);
	t1 = ((double)begin.ru_stime.tv_sec + MICROSEC_TO_SEC*((double)begin.ru_stime.tv_usec));
	t2 = ((double)end.ru_stime.tv_sec + MICROSEC_TO_SEC*((double)end.ru_stime.tv_usec));
	sysTime += (t2 - t1);
	// memory
	maxMemory = end.ru_maxrss / KB_PER_MB;
#if defined(__MACH__) && defined(__APPLE__)
	// on Mac, ru_maxrss is in bytes and not in kilobytes
	maxMemory /= KB_PER_MB;
#endif
	pageFaults += (end.ru_minflt - begin.ru_minflt);
	pageFaultsIO += (end.ru_majflt - begin.ru_majflt);
	numSwap += (end.ru_nswap - begin.ru_nswap);
}

inline void ResourceUsage::reset()
{
	userTime = 0.0;
	sysTime = 0.0;
	maxMemory = 0.0;
	pageFaults = 0;
	pageFaultsIO = 0;
	numSwap = 0;
}

double ResourceUsage::getUserTime() const
{
	return userTime;
}

double ResourceUsage::getSysTime() const
{
	return sysTime;
}

double ResourceUsage::getMaxMemory() const
{
	return maxMemory;
}

int ResourceUsage::getMinorPageFaults() const
{
	return pageFaults;
}

int ResourceUsage::getMajorPageFaults() const
{
	return pageFaultsIO;
}

int ResourceUsage::getSwapCount() const
{
	return numSwap;
}

} // namespace dominiqs

#endif /* RUSAGE_H */
