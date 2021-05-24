/**
 * @file args_parser.cpp
 * @brief Command line parsing utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2007-2014
 */

#include "utils/args_parser.h"
#include "utils/fileconfig.h"
#include "utils/asserter.h"
#include "utils/str_utils.h"


namespace dominiqs {

ArgsParser::ArgsParser() {}

void ArgsParser::reset()
{
	binary = "";
	input.clear();
	config.clear();
	overrides.clear();
}

bool ArgsParser::parse(int argc, char const* argv[])
{
	binary = std::string(argv[0]);
	std::vector<std::string> tokens;
	for (int i = 1; i < argc; i++) tokens.emplace_back(argv[i]);
	return parse(tokens);
}


bool ArgsParser::parse(const std::vector<std::string>& tokens)
{
	std::ostringstream out;
	enum class ParserStatus {DEFAULT, CONFIG, OVERRIDE};
	ParserStatus status = ParserStatus::DEFAULT;
	for (auto tok: tokens)
	{
		if (status == ParserStatus::CONFIG)
		{
			config.push_back(tok);
			status = ParserStatus::DEFAULT;
		}
		else if (status == ParserStatus::OVERRIDE)
		{
			overrides.push_back(tok);
			status = ParserStatus::DEFAULT;
		}
		else
		{
			// status == DEFAULT
			if (tok == "-c" || tok == "--config") status = ParserStatus::CONFIG;
			else if (tok == "-C") status = ParserStatus::OVERRIDE;
			else if (tok.find_first_of('=') != std::string::npos) overrides.push_back(tok);
			else input.push_back(tok);
		}
	}
	return true;
}

void mergeConfig(const ArgsParser& args, FileConfig& config)
{
	/* load config files */
	for (auto conf: args.config) config.load(conf);
	/* overrides */
	for (auto over: args.overrides)
	{
		size_t eqIdx = over.find_first_of('=');
		DOMINIQS_ASSERT( eqIdx != std::string::npos );
		std::string entry(over.begin(), over.begin() + eqIdx);
		std::string value(over.begin() + eqIdx + 1, over.end());
		config.set(entry, value);
	}
}

} // namespace dominiqs
