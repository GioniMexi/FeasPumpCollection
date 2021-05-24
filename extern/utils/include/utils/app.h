/**
 * @file app.h
 * @brief Basic Application
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2011 Domenico Salvagnin
 */

#ifndef APP_H
#define APP_H

#include <stdint.h>
#include "args_parser.h"

namespace dominiqs
{

/**
* Ctrl-C handling
*/

extern int UserBreak;

/**
 * @brief Provides a basic commmand-line application logic
 *
 * The purpose of this class is to centralize common stuff such as:
 * - Ctrl-C handling
 * - parsing command line options / config files
 * - exception handling (main try/catch block)
 *
 * Users must derived from this class to implement custom behaviour.
 * In particular, derived classes MUST implement:
 * - checkUsage()
 * - exec()
 *
 * Derived classes MAY also override:
 * - readConfig()
 * - startup() / shutdown()
 */

class App
{
public:
	App();
	/** destructor */
	virtual ~App() {}
	/**
	 * @brief add the extension of a type of input files
	 */
	void addExtension(const std::string& ext);
	/**
	* @brief parse command line arguments and config files
	*
	* @param argc number of command line arguments (including program name)
	* @param argv vector of command line arguments
	* @return true if successful, false otherwise
	*/
	bool parseArgsAndConfig(int argc, char const *argv[]);
	/**
	* @brief read config params: to be implemented in derived classes
	*/
	virtual void readConfig() {}
	/**
	* @brief application overall logic
	*
	* DO NOT to override in derived classes: derive exec() instead!!!
	* @return 0 if successful, -1 otherwise
	*/
	int run();
protected:
	/**
	* @brief check if required command line parameters are there
	* MUST be implemented in derived classes
	* @return true if successful, false otherwise
	*/
	virtual bool checkUsage() = 0;
	/**
	* @brief additional startup work
	* to be implemented in derived classes (optional)
	*/
	virtual void startup() {}
	/**
	* @brief all custom application logic is here
	* MUST be implemented in derived classes
	*/
	virtual void exec() = 0;
	/**
	* @brief additional shutdown work
	* to be implemented in derived classes (optional)
	*/
	virtual void shutdown() {}
	/* data */
	ArgsParser args;
	std::vector<std::string> extensions;
private:
	bool parseDone;
};

} // namespace dominiqs

#endif /* APP_H */
