/**
 * @file machine_utils.h
 * @brief Machine utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2009-2012 Domenico Salvagnin
 */

#ifndef MACHINE_UTILS_H
#define MACHINE_UTILS_H

#include <stddef.h> //< for size_t and ptrdiff_t
#include <unistd.h> //< for sysconf
#include <cstdlib>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include <emmintrin.h>
#endif // __APPLE__

#ifndef NDEBUG
#include <exception> //< for terminate
#include <iostream> //< for clog
#endif // NDEBUG

#include "asserter.h"

namespace dominiqs {

/**
 * Default memory alignment
 */

static const size_t MEM_ALIGNMENT = 16;

/**
 * Allocate @param size bytes of memory
 * The returned block is guaranteed to be 16-aligned (to play nicely with SSE2)
 * @return a pointer to the allocated block
 */

inline void* mallocSSE2(std::size_t size)
{
#ifdef __APPLE__
	void* result = std::malloc(size);
#else
	void* result = _mm_malloc(size, MEM_ALIGNMENT);
#endif
	DOMINIQS_ASSERT( (reinterpret_cast<ptrdiff_t>(result) % MEM_ALIGNMENT) == 0 );
	return result;
}

/**
 * Free a 16-aligned memmory block allocated by mallocSSE2
 */

inline void freeSSE2(void* ptr)
{
#ifndef NDEBUG
	if ((reinterpret_cast<ptrdiff_t>(ptr) % MEM_ALIGNMENT) != 0)
	{
		std::clog << "freeSSE2 called with bad pointer address: " << ptr << std::endl;
		std::terminate();
	}
#endif // NDEBUG
	#ifdef __APPLE__
		std::free(ptr);
	#else
		_mm_free(ptr);
	#endif
}


namespace implementation
{
	/**
	 * Computes the numbers of bytes needed to store an array of
	 * count elements of type T, such that the total is aligned
	 * according to align
	 */

	template<typename T>
	size_t bytesNeeded(size_t count, size_t align = MEM_ALIGNMENT)
	{
		size_t ret = sizeof(T) * count;
		return ((ret + MEM_ALIGNMENT - 1) & ~(MEM_ALIGNMENT - 1));
	}

	/**
	 * Computes the contribution (in bytes) of an array of
	 * count elements of type T to the size of block
	 * This is the "base case" for the variadic template allocByteAdder
	 */

	template<typename T>
	size_t allocByteAdder(T** pointer, size_t count) {
		return bytesNeeded<T>(count);
	}

	/**
	 * Computes the size of a memory block needed to hold a list
	 * of arrays of types T_1,..T_k and sizes count_1,...,count_k.
	 */

	template<typename T, typename... Args>
	size_t allocByteAdder(T** pointer, size_t count, Args... args) {
		return allocByteAdder(pointer, count) + allocByteAdder(args...);
	}

	/**
	 * Sets *pointer to point to the proper offset within the memory block
	 * This is the "base case" for the variadic template allocOffsetSet
	 */

	template<typename T>
	void allocOffsetSet(void** mem, T** pointer, size_t count) {
		*pointer = static_cast<T*>(*mem);
		*mem = static_cast<unsigned char*>(*mem) + bytesNeeded<T>(count);
	}

	/**
	 * Sets a list of pointers to point to the proper offset
	 * within the memory block.
	 */

	template<typename T, typename... Args>
	void allocOffsetSet(void** mem, T** pointer, size_t count, Args... args) {
		allocOffsetSet(mem, pointer, count);
		allocOffsetSet(mem, args...);
	}
}

/**
 * Allocate a single contiguous block of memory which is split into several
 * arrays, each of type T_i and with count_i elements.
 *
 * Usage: calling allocateBlock(T_1** p1, size_t count1,...,T_k** pk, size_t countk)
 * allocate a single memory block big enough (and properly aligned) to hold the k arrays.
 * Each pointer p will then be set to the proper offset within the memory block
 *
 * @return the address of the allocated memory block
 */
template<typename... Args>
void* allocateBlock(Args... args)
{
	size_t tot = implementation::allocByteAdder(args...);
	void* mem = mallocSSE2(tot);
	void* itr = mem;
	implementation::allocOffsetSet(&itr, args...);
	return mem;
}

/**
 * Free the memory block pointed by ptr
 */

inline void freeBlock(void* ptr)
{
	freeSSE2(ptr);
}


/**
 * Return the amount of physical memory of the system in bytes
 */

inline double getPhysicalMemory()
{
#ifdef __APPLE__
	uint64_t memsize;
	size_t len = sizeof(memsize);
	int mib[2] = { CTL_HW, HW_MEMSIZE };
	if ((sysctl(mib, 2, &memsize, &len, 0, 0) == 0) && (len == sizeof(memsize))) return (double)memsize;
	return 0;
#else
	double res = sysconf(_SC_PAGESIZE) * (double)sysconf(_SC_PHYS_PAGES);
	if (res >= 0) return res;
	return 0;
#endif // __APPLE__
}


/**
 * Return the amount of available memory of the system in bytes
 */

inline double getAvailableMemory()
{
#ifdef __APPLE__
	int memsize;
	size_t len = sizeof(memsize);
	int mib[2] = { CTL_HW, HW_USERMEM };
	if ((sysctl(mib, 2, &memsize, &len, 0, 0) == 0) && (len == sizeof(memsize))) return (double)memsize;
	return 0;
#else
	double res = sysconf(_SC_PAGESIZE) * (double)sysconf(_SC_AVPHYS_PAGES);
	if (res >= 0) return res;
	return 0;
#endif // __APPLE__
}

/**
 * Return the number of physical cores in the system
 * This number may be the same as the number of logical cores,
 * depending on hyper-threading.
 * The returned number is accured on Mac platforms (less so on
 * Linux, where we resort to sysconf)
 */

inline int getNumPhysicalCores()
{
#ifdef __APPLE__
	int physicalCores;
	size_t len = sizeof(physicalCores);
	sysctlbyname("hw.physicalcpu", &physicalCores, &len, 0, 0);
	return physicalCores;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif // __APPLE__
}

/**
 * Return the number of logical cores in the system
 */

inline int getNumLogicalCores()
{
	return sysconf(_SC_NPROCESSORS_ONLN);
}

} // namespace dominiqs

#endif /* MACHINE_UTILS_H */
