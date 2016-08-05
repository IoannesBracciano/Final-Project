#include <iostream>
#include "../Filter/generic.hpp"

extern "C" {
#include "vl/sift.h"
}


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
		//("action", p_opts::value<std::string>()->required(),
		//	"[extract|match] Which action to perform")
		("input,i", p_opts::value<std::vector<std::string>>()
			->default_value(std::vector<std::string>(), "current working directory"),
			"Absolute or relative path to the input image file or files, "
			"or a directory containing image files")
		("output,o", p_opts::value<std::string>()
			->default_value(std::string(), "input path directory"),
			"Absolute or relative path to an existing directory "
			"where extracted descriptor files are saved");
	p_opts_pos//.add("action", 1)
		.add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}


/** The next portion of code was ported from the matlab interface of the vlfeat library
as there was no suitable code for matching between two image descriptors in C API
and was copied and pasted almost "as is".
See http://stackoverflow.com/questions/26925038/how-to-use-vlfeat-sift-matching-function-in-c-code
*/

/// <summary>Holds matched descriptor pairs</summary>
typedef struct {
	/// <summary>Descriptor of first image</summary>
	int k1;
	/// <summary>Descriptor of second image</summary>
	int k2;
	/// <summary>Score of matching</summary>
	double score;
} Pair;

/// <summary>Compares the descriptors of two different images and matches the pairs</summary>
/// <param name="pairs">pointer to the vector to be filled with the matched pairs.
///                     Needs to be big enough to hold all possible pairs</param>
/// <param name="descr1">pointer to the descriptor matrix of the first image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="descr2">pointer to the descriptor matrix of the second image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="K1">Number of descriptors of first image</param>
/// <param name="K2">Number of descriptors of second image</param>
/// <param name="ND">Descriptor size (=128 for SIFT)</param>
/// <param name="thresh">Ratio test threshold value (=1.5)</param>
Pair * compare(
	Pair *pairs,
	const float *descr1,
	const float *descr2,
	int K1,
	int K2,
	int ND = 128,
	float thresh = 1.5
	)
{
	int k1, k2;

	/* Loop over 1st image descr. */
	for (k1 = 0; k1 < K1; ++k1, descr1 += ND) {
		float best = FLT_MAX;
		float second_best = FLT_MAX;
		int bestk = -1;

		/* Loop over 2nd image descr. and find the 1st and 2nd closest descr. */
		for (k2 = 0; k2 < K2; ++k2, descr2 += ND) {
			int bin;
			float acc = 0;

			/* Compute the square L2 distance between descriptors */
			for (bin = 0; bin < ND; ++bin) {
				float delta = descr1[bin] - descr2[bin];
				acc += delta*delta;
				if (acc >= second_best)
					break;
			}

			if (acc < best) {
				second_best = best;
				best = acc;
				bestk = k2;
			}
			else if (acc < second_best) {
				second_best = acc;
			}
		}

		/* Rewind */
		descr2 -= ND*K2;

		/* Record the correspondence if the best descr. passes the ratio test */
		if (thresh * best < second_best && bestk != -1) {
			pairs->k1 = k1;
			pairs->k2 = bestk;
			pairs->score = best;
			pairs++;
		}
	}

	return pairs;
}

///
/// ----- Main entry point -----
///
int main(int argc, char **  argv)
{
	// Try to set up program options
	try
	{
		setup_p_opts(argc, argv);

		// Deal with help option before notifying the options variable map
		// to avoid conflict with required variables
		if (p_opts_vm.count("help"))
		{
			std::cout << p_opts_desc << std::endl;
			return SUCCESSFUL_RUN;
		}

		// Notify the options variable map,
		// throws an exception if required arguments are missing
		p_opts::notify(p_opts_vm);
	}
	// Show error message if required options are indeed missing
	catch (p_opts::required_option& e)
	{
		std::cout << e.what() << std::endl;
		return REQUIRED_OPT_NOT_SET;
	}

	/* Program was called correctly, OK to continue */
}