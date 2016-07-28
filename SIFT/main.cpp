#include "boost\program_options.hpp"
#include "boost\filesystem.hpp"
extern "C" {
#include "vl/generic.h"
}

namespace p_opts = boost::program_options;
namespace fsys = boost::filesystem;

const enum RETURN_CODES
{
	// Error codes
	REQUIRED_OPT_NOT_SET = -1,
	INVALID_INPUT = -2,
	INVALID_OUTPUT = -3,

	// Success codes
	SUCCESSFUL_RUN = 0
};

// Holds the description of the program options
p_opts::options_description p_opts_desc("Arguments description");

// Holds positional description of program options
p_opts::positional_options_description p_opts_pos;

// Array used to map variable names with their corresponding values
p_opts::variables_map p_opts_vm;

/// <summary>Sets up the program options and parses the argv vector accordingly
/// </summary>
/// <param name="argc">the arguments count</param>
/// <param name="argc">the standard c++ arguments vector</param>
void setup_p_opts(int argc, char ** argv)
{
	p_opts_desc.add_options()
		("help", "Print this help message")
		("action", p_opts::value<std::string>()->required(), "[extract|match] Which action to perform")
		("input,i", p_opts::value<std::vector<std::string>>()
			->default_value(std::vector<std::string>(), "current working directory"),
			"Absolute or relative path to the input image file or files, or a directory containing image files")
		("output,o", p_opts::value<std::string>()
			->default_value(std::string(), "input path directory"),
			"Absolute or relative path to an existing directory where filtered images are saved");
	p_opts_pos.add("action", 1).add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}


int main(int argc, char **  argv)
{
	VL_PRINT("Hello world!\n");
	return SUCCESSFUL_RUN;
}