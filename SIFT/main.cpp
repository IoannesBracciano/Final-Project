#include <iostream>
#include "../Filter/generic.hpp"
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"

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
			"where extracted descriptor files are saved")
		("peak-thresh,t", p_opts::value<double>()->default_value(0.03),
			"Minimun amount of contrast to accept a keypoint."
			"Increasing this value will result in fewer keypoints.")
		("edge-thresh,T", p_opts::value<double>()->default_value(10),
			"Edge rejection threshold. Increasing this will result is fewer keypoints.")
		("octaves", p_opts::value<int>()->default_value(-1),
			"Number of octaves")
		("levels", p_opts::value<int>()->default_value(3),
			"Number of levels per octave");
	p_opts_pos//.add("action", 1)
		.add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}

void extract(fsys::path & input, std::vector<vl_sift_pix*> & descriptors,
	int octaves, int levels, double peak_thresh, double edge_thresh)
{
	// Open image file
	cv::Mat I = cv::imread(input.string(), cv::IMREAD_GRAYSCALE);
	                                       // SIFT works only on gray scale images
	// Normalize pixel values in the range [0,1]
	I.convertTo(I, CV_32F, 1.0 / 255.5);

	// Create a new SIFT detector
	VlSiftFilt * detector = vl_sift_new(I.cols, I.rows, octaves, levels, 0);

	// Setting peak and edge thresholds
	vl_sift_set_peak_thresh(detector, peak_thresh);
	vl_sift_set_edge_thresh(detector, edge_thresh);

	// Return if there is no octaves to calculate
	if (vl_sift_process_first_octave(detector, (float*)I.data) == VL_ERR_EOF)
	{
		vl_sift_delete(detector);
		return;
	}

	// Extract the keypoints of each octave
	do
	{
		const VlSiftKeypoint * keypoints;
		vl_sift_detect(detector);
		keypoints = vl_sift_get_keypoints(detector);

		// Compute the descriptors of all keypoints in the current octave
		for (int i = 0; i < vl_sift_get_nkeypoints(detector); ++i)
		{
			// First calculate the orientation of the keypoint
			// (library supports up to 4 orientations for each keypoint)
			double angles[4];
			int nangles = vl_sift_calc_keypoint_orientations(detector, angles, keypoints + i);

			// For each orientation construct a different keypoint descriptor
			for (int j = 0; j < nangles; ++j)
			{
				vl_sift_pix descriptor[128];
				vl_sift_calc_keypoint_descriptor(detector, descriptor, keypoints + i, angles[j]);
				descriptors.push_back(descriptor);
			}
		}
	} while (vl_sift_process_next_octave(detector) != VL_ERR_EOF);
	                                   //           A
	                                   //  One God  | damn day
	                                   // to figure | out I was
	                                   //     doing | that wrong
	                                   //  ;@@@@@@@@@@@@@@@@@@

	// Delete SIFT detector used for this image
	vl_sift_delete(detector);
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

	// Getting input paths to resolve them
	std::vector<std::string> input;
	// If user has not specified an input
	// we consider it to be the current working directory
	if (p_opts_vm["input"].defaulted())
	{
		input.push_back(fsys::current_path().string());
	}
	else
	{
		input = p_opts_vm["input"].as<std::vector<std::string>>();
	}
	// Create another vector that will hold all input files after resolving the paths
	std::vector<fsys::path> files = std::vector<fsys::path>();

	// Try to resolve the input paths
	try
	{
		for (int i = 0; i < input.size(); ++i)
		{
			resolveInput(input[i], files);
		}
	}
	// Report any filesystem errors and terminate
	catch (fsys::filesystem_error & e)
	{
		std::cout << e.what() << std::endl;
		return INVALID_INPUT;
	}
	// Report runtime errors thrown by resolveInput()
	catch (std::runtime_error & e)
	{
		std::cout << e.what() << std::endl;
		return INVALID_INPUT;
	}

	// Resolving output
	fsys::path output;
	// If user providing an output, extract it
	if (!p_opts_vm["output"].defaulted())
	{
		output = p_opts_vm["output"].as<std::string>();
		if (!fsys::exists(output) || !fsys::is_directory(output))
		{
			std::cout << "Output must be an existing directory" << std::endl;
			return INVALID_OUTPUT;
		}
	}

	// Extracting action
	// std::string action = p_opts_vm["action"].as<std::string>();

	// Extracting peak thresh
	double peak_thresh = p_opts_vm["peak-thresh"].as<double>();
	// Extracting edge threshold
	double edge_thresh = p_opts_vm["edge-thresh"].as<double>();
	// Extracting octaves
	int octaves = p_opts_vm["octaves"].as<int>();
	// Extracting levels of each octave
	int levels = p_opts_vm["levels"].as<int>();

	// Used to save the user choice when he is asked whether to 
	// override an already existing file or not.
	char replace_file = 0;

	// Looping through each input image file to extract its feature keypoints
	for (int i = 0; i < files.size(); ++i)
	{
		// Vector that holds the image descriptors
		std::vector<vl_sift_pix*> descriptors = std::vector<vl_sift_pix*>();

		// When user has not provided a specific output path
		// save filtered images in the same directory as the original ones
		if (p_opts_vm["output"].defaulted())
		{
			output = files[i].parent_path();
		}

		// Output file is named after the input image
		// but with a different extension to reflect that it is
		// a Scale Invariant Feature Descriptors file
		fsys::path new_filename = files[i].filename().replace_extension(".sifd");

		// If filename exists we are asking the user if he wants to replace the file
		if (fsys::exists(output / new_filename))
		{
			char perform_on_all;

			if (replace_file == 0)
			{
				std::cout << (output / new_filename) << " exists. Do you want to replace it? (y/n, default no):";
				do { std::cin >> replace_file; } while (replace_file == '\n');

				while (replace_file != 'n' && replace_file != 'N' &&
					replace_file != 'y' && replace_file != 'Y')
				{
					std::cout << "Didn't get that, type 'y' for yes or 'n' for no:";
					do { std::cin >> replace_file; } while (replace_file == '\n');
				}

				std::cout << "Perform this action for all existing files? (y/n, default no):";
				do { std::cin >> perform_on_all; } while (replace_file == '\n');

				while (perform_on_all != 'n' && perform_on_all != 'N' &&
					perform_on_all != 'y' && perform_on_all != 'Y')
				{
					std::cout << "Didn't get that, type 'y' for yes or 'n' for no:";
					do { std::cin >> perform_on_all; } while (replace_file == '\n');
				}
			}

			if (replace_file == 'n' || replace_file == 'N')
			{
				std::cout << "Skipping " << files[i].filename() << std::endl;
				if (perform_on_all == 'n' || perform_on_all == 'N')
				{
					replace_file = 0;
				}
				continue;
			}
		}

		std::cout << "Extracting keypoints from " << files[i].filename() << "...";
		extract(files[i], descriptors, octaves, levels, peak_thresh, edge_thresh);
		std::cout << " Extracted " << descriptors.size() << " keypoints" << std::endl;

		// Create a binary file and write the descriptors in raw format
		fsys::ofstream ofs((output / new_filename), std::ios::out | std::ios::binary);
		for (std::vector<vl_sift_pix*>::iterator it = descriptors.begin(); it != descriptors.end(); ++it)
		{
			vl_sift_pix * descriptor = *it;
			ofs.write(reinterpret_cast<const char *>(descriptor), sizeof(vl_sift_pix) * 128);
		}
		ofs.close();
	}

	return SUCCESSFUL_RUN;
}