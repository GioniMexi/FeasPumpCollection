/**
 * @file app.cpp
 * @brief Basic Application
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2011 Domenico Salvagnin
 */

#include "utils/app.h"
#include "utils/fileconfig.h"
#include "utils/consolelog.h"
#include <signal.h>

namespace dominiqs
{

int UserBreak = 0;

static void userSignalBreak(int signum)
{
	UserBreak = 1;
}

App::App() : parseDone(false) {}

void App::addExtension(const std::string& ext)
{
	extensions.push_back(ext);
}

bool App::parseArgsAndConfig(int argc, char const *argv[])
{
	// parse command line arguments and check usage
	args.parse(argc, argv);
	if (!checkUsage()) return false;
	// read config
	mergeConfig(args, gConfig());
	parseDone = true;
	return true;
}

int App::run()
{
	if (!parseDone) return -1;
	int retValue = 0;

	// read config
	readConfig();
	// Ctrl-C handling
	UserBreak = 0;
	void (*previousHandler)(int) = ::signal(SIGINT, userSignalBreak);
	try
	{
		startup();
		exec();
		shutdown();
	}
	catch(std::exception& e)
	{
		retValue = -1;
		consoleError(e.what());
	}
	// restore signal handler
	::signal(SIGINT, previousHandler);
	return retValue;
}

} // namespace dominiqs
