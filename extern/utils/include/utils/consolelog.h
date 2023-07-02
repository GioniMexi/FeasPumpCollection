/**
 * @file consolelog.h
 * @brief Console logging functions
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2023
 */

#ifndef CONSOLELOG_H
#define CONSOLELOG_H

#include <fmt/format.h>
#include <fmt/color.h>
#include <unistd.h>


template<typename ...Args>
void consoleLog(Args&&... args)
{
	fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleInfo(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::green), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleWarn(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::yellow), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleError(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::red), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


#ifdef DEBUG_LOG

enum class DebugLevel {
	None = 0,
	Minimal = 1,
	Normal = 2,
	Verbose = 3,
	VeryVerbose = 4,
	All = 5
};

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL DebugLevel::Normal
#endif //< DEBUG_LEVEL

template<typename ...Args>
void consoleDebug(DebugLevel level, Args&&... args)
{
	if (level > DEBUG_LEVEL)  return;
	consoleLog(std::forward<Args>(args)...);
}

#else //< DEBUG_LOG

#define consoleDebug(level, ...)  do {} while (false)

#endif //<DEBUG_LOG

#endif /* CONSOLELOG_H */
