/**
 * Filter.exe
 * usage: Filter.exe filter [input] [-o output] [-S number] [-s number] [-gd]
 *
 * Filters images using either a fast approximation of bilateral filter
 * or spatial-range joint mean shift. 
 * Part of the final project for the dept. of Computer Science,
 * Faculty of Applied Sciences, T.E.I. of Athens
 *
 * Author: Ioannes Bracciano <john.bracciano@gmail.com>
 * All rights reserved
*/

#include "generic.hpp"
#include "FBLF.h"

/*
 * Holds the description of the program options
 */
p_opts::options_description p_opts_desc("Arguments description");

/*
 * Holds positional description of program options
 */
p_opts::positional_options_description p_opts_pos;

/*
 * Array used to map variable names with their corresponding values
 */
p_opts::variables_map p_opts_vm;


/// <summary>Sets up the program options and parses the argv vector accordingly
/// </summary>
/// <param name="argc">the arguments count</param>
/// <param name="argc">the standard c++ arguments vector</param>
void setup_p_opts(int argc, const char ** argv)
{
	p_opts_desc.add_options()
		("help",			"Print this help message")
		("filter,f", p_opts::value<std::string>()->required(), 
							"[bilateral|meanshift] Name of filter to be used")
		("input,i", p_opts::value<std::vector<std::string>>()
								->default_value(std::vector<std::string>(), "current working directory"),
							"Absolute or relative path to the input image file or files, or a directory containing image files")
		("output,o", p_opts::value<std::string>()
								->default_value(std::string(), "input path directory"),
							"Absolute or relative path to an existing directory where filtered images are saved")
		("sp,S", p_opts::value<double>()->default_value(10.0),
							"Value of spatial domain gaussian sigma parameter (bilateral), or spatial window radius (meanshift)")
		("sr,s", p_opts::value<double>()->default_value(10.0),
							"Value of range domain gaussian sigma parameter (bilateral), or color window radius (meanshift)")
		("gray-scale,g",	"Perform the filtering on gray scale")
		("debug,d",			"Enable the debug mode");
	p_opts_pos.add("filter", 1).add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}


/// <summary>
/// Opens the input image file to the memory, applies the bilateral filter on it,
/// and saves the filtered image in a new file in the directory speciied by output
/// </summary>
/// <param name="input">the input file path</param>
/// <param name="output">the output file path</param>
/// <param name="sigma_space">the gaussian parameter sigma in the spatial domain</param>
/// <param name="sigma_range">the gaussian parameter sigma in the range domain</param>
/// <param name="gray_scale">if set, filtering is applied on a gray scaled version of the input image</param>
/// <returns>the time of execution of the filtering in seconds</returns>
double apply_bilateral(const fsys::path & input, const fsys::path & output,
	double sigma_space, double sigma_range, bool gray_scale)
{
	std::cout << "Applying bilateral filter on " << input.filename() << "... ";
	double exec_time_s = 0.0; // Execution time in seconds

	if (gray_scale)
	{
		cv::Mat src = cv::imread(input.string(), cv::IMREAD_GRAYSCALE);
		cv::Mat dst = Mat(src.rows, src.cols, CV_8U);
		exec_time_s = fastBilateralFilter(src, dst, sigma_space, sigma_range);
		std::cout << exec_time_s << " seconds" << std::endl;
		imwrite(output.string(), dst);
	}
	else
	{
		cv::Mat src = cv::imread(input.string(), cv::IMREAD_UNCHANGED);
		cv::Mat dst = Mat(src.rows, src.cols, CV_8U);
		cv::Mat bgr[3]; cv::split(src, bgr);
		std::vector<cv::Mat> res;

		for (uint8_t i = 0; i < 3; ++i) {
			res.push_back(Mat(src.rows, src.cols, CV_8U));
			exec_time_s += fastBilateralFilter(bgr[i], res.at(i), sigma_space, sigma_range);
		}
		std::cout << exec_time_s << " seconds" << std::endl;

		merge(res, dst);
		imwrite(output.string(), dst);
	}
	
	return exec_time_s;
}


/// <summary>
/// Opens the input image file to the memory, applies the joint mean shift technique on it,
/// and saves the filtered image in a new file in the directory speciied by output
/// The mean shift algorithm is found in the opencv library and it is actually a pyramid
/// implementation of the algorith. Setting the number of pyramid levels to 1, the algorithm
/// performs the mean shift technique to the original sized image only
/// </summary>
/// <param name="input">the input file path</param>
/// <param name="output">the output file path</param>
/// <param name="sp">the spatial window radius</param>
/// <param name="sr">the color window radius</param>
/// <param name="gray_scale">if set, filtering is applied on a gray scaled version of the input image</param>
/// <returns>the time of execution of the filtering in seconds</returns>
double apply_meanshift(const fsys::path & input, const fsys::path & output,
	double sp, double sr, bool gray_scale)
{
	std::cout << "Applying mean shift on " << input.filename() << "... ";

	cv::Mat src = cv::imread(input.string(), cv::IMREAD_COLOR), dst;
	double time_start = (double)cv::getTickCount();
	cv::pyrMeanShiftFiltering(src, dst, sp, sr);	// Pyramid levels is defaulted to 1
	double exec_time_s = ((double)getTickCount() - time_start) / getTickFrequency();
	std::cout << exec_time_s << " seconds" << std::endl;

	if (gray_scale)
	{
		cv::cvtColor(dst, dst, cv::COLOR_BGR2GRAY);
	}

	imwrite(output.string(), dst);

	return exec_time_s;
}

///
/// --- Main entry point ---
///
int main(int argc, const char** argv) {

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
	// Show error message if command line syntaxt is wrong
	catch (p_opts::invalid_command_line_syntax & e)
	{
		std::cout << e.what() << std::endl;
		return INVALID_SYNTAX;
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

	// Extracting filter option
	std::string filter = p_opts_vm["filter"].as<std::string>();

	// Extracting sigma options
	double sigma_space = p_opts_vm["sp"].as<double>();
	double sigma_range = p_opts_vm["sr"].as<double>();

	// Exctracting gray scale option
	bool gray_scale = p_opts_vm.count("gray-scale");

	// Keeps the total time it took to process the images (in seconds)
	double t_total = 0.0;

	// Used to save the user choice when he is asked whether to 
	// override an already existing file or not.
	char replace_file = 0;

	// Number of files that already existed and skipped
	unsigned int skipped_files_count = 0;

	for (int i = 0; i < files.size(); ++i)
	{
		// When user has not provided a specific output path
		// save filtered images in the same directory as the original ones
		if (p_opts_vm["output"].defaulted())
		{
			output = files[i].parent_path();
		}

		// Constructing the output filename
		fsys::path new_filename = fsys::path(
			files[i].filename().stem().string() +				// Removing file extension
			"_filtered")										// Appending filter to the old file name
			.replace_extension(files[i].extension());			// Replacing extension

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
				std::cout << "Skipping " << new_filename << std::endl;
				++skipped_files_count;
				if (perform_on_all == 'n' || perform_on_all == 'N')
				{
					replace_file = 0;
				}
				continue;
			}
		}

		if (filter == "bilateral")
		{
			t_total += apply_bilateral(files[i], (output / new_filename),
				sigma_space, sigma_range, gray_scale);
		}
		else if (filter == "meanshift")
		{
			t_total += apply_meanshift(files[i], (output / new_filename),
				sigma_space, sigma_range, gray_scale);
		}
	}

	std::cout << "Total execution time (image processing only): " << t_total << " seconds" << std::endl
		<< "Mean execution time per image: " << t_total / ( files.size() - skipped_files_count ) << " seconds" << std::endl;

	return SUCCESSFUL_RUN;

}