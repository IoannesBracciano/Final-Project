#include <iostream>
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"
#include "boost\program_options.hpp"
#include "boost\filesystem.hpp"

extern "C" {
#include "vl/sift.h"
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
		("action", p_opts::value<std::string>()->required(),
			"[extract|match] Which action to perform")
		("input,i", p_opts::value<std::vector<std::string>>()
			->default_value(std::vector<std::string>(), "current working directory"),
			"Absolute or relative path to the input image file or files, or a directory containing image files")
		("output,o", p_opts::value<std::string>()
			->default_value(std::string(), "input path directory"),
			"Absolute or relative path to an existing directory where filtered images are saved")
		("peak-thresh,t", p_opts::value<double>()->default_value(0.03),
			"Minimun amount of contrast to accept a keypoint."
			"Increasing this value will result in fewer keypoints.")
		("edge-thresh,T", p_opts::value<double>()->default_value(3.5),
			"Edge rejection threshold. Increasing this will result is fewer keypoints.")
		("octaves", p_opts::value<int>()->default_value(-1),
			"Number of octaves")
		("levels", p_opts::value<int>()->default_value(3),
			"Number of levels per octave");
	p_opts_pos.add("action", 1).add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}

/// <summary>
/// Checks if given path is a supported image file.
/// The list of supported OpenCV 3 image files, was taken from
/// http://docs.opencv.org/3.0-beta/modules/imgcodecs/doc/reading_and_writing_images.html#imread
/// </summary>
/// <param name="path">the path of the file to be examined</param>
bool is_supported_image_file(fsys::path path)
{
	return
		// Windows bitmaps
		path.extension() == ".bmp" || path.extension() == ".dib" ||
		// JPEG files
		path.extension() == ".jpg" || path.extension() == ".jpeg" || path.extension() == ".jpe" ||
		// JPEG 2000 files
		path.extension() == ".jp2" ||
		// Portable Network Graphics
		path.extension() == ".png" ||
		// WebP
		path.extension() == ".webp" ||
		// Portable image format
		path.extension() == ".pbm" || path.extension() == ".pgm" || path.extension() == ".ppm" ||
		// Sun rasters
		path.extension() == ".sr" || path.extension() == ".ras" ||
		// TIFF files
		path.extension() == ".tiff" || path.extension() == ".tif";
}


/// <summary>
/// Resolves the input program option,
/// populating a vector with the paths of supported files that are resolved
/// </summary>
/// <param name="path">the path to be resolved.
/// if it points to a directory, all supported files are extracted
/// </param>
/// <param name="files">a reference to the vector where extracted paths are kept</param>
/// <exception cref="std::runtime_error">
/// Thrown when input path does not exist or is neither a file nor a directory
/// </exception>
void resolveInput(const fsys::path & path, std::vector<fsys::path> & files)
{
	//std::cout << "Input: " << path << std::endl;
	if (fsys::exists(path))									// If path exists
	{
		if (fsys::is_regular_file(path))					// and is a regular file
		{
			if (is_supported_image_file(path))				// and it is a supported image file
			{
				files.push_back(path);						// add it to the files vector
			}
			else // Unsupported file or not an image file at all
			{
				throw std::runtime_error("Unsupported file. For a list of supported image files see: "
					"http://docs.opencv.org/3.0-beta/modules/imgcodecs/doc/reading_and_writing_images.html#imread");
			}
		}
		else if (fsys::is_directory(path))					// else if path is a directory
		{
			fsys::directory_iterator end_itr;
			for (fsys::directory_iterator itr(path); itr != end_itr; ++itr)
			{												// Iterate through all files one level down
				if (fsys::is_regular_file(itr->path()))
				{
					if (is_supported_image_file(itr->path()))
					{										// and if any is a supported image file
						files.push_back(itr->path());		// add it to the files vector
					}
				}
			}
		}
		else // Input path is neither a file nor a directory
		{
			throw std::runtime_error("Input path is neither a file nor a directory");
		}
	}
	else // Input path does not exist
	{
		throw std::runtime_error("Input path does not exist");
	}
}

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

	// Extracting action
	std::string action = p_opts_vm["action"].as<std::string>();
	// Extracting peak thresh
	double peak_thresh = p_opts_vm["peak-thresh"].as<double>();
	// Extracting edge threshold
	double edge_thresh = p_opts_vm["edge-thresh"].as<double>();
	// Extracting octaves
	int octaves = p_opts_vm["octaves"].as<int>();
	// Extracting levels of each octave
	int levels = p_opts_vm["levels"].as<int>();

	// Vector that holds the image descriptors
	std::vector<vl_sift_pix*> descriptors = std::vector<vl_sift_pix*>();

	// Looping through each input image file to extract its feature keypoints
	for (int i = 0; i < files.size(); ++i)
	{
		// Open image file
		cv::Mat I = cv::imread(files[i].string(), cv::IMREAD_GRAYSCALE);
		                                          // SIFT works only on gray scale images
		// Normalize pixel values in the range [0,1]
		I.convertTo(I, CV_32F, 1.0 / 255.5);

		std::cout << "Extracting keypoints from " << files[i].filename() << "...";

		// Create a new SIFT detector
		VlSiftFilt * detector = vl_sift_new(I.cols, I.rows, octaves, levels, 0);

		// If first octave cannot be calculated then there is probably
		// some sort of error with the image file, so continue to the next one
		if (vl_sift_process_first_octave(detector, (float*)I.data) == VL_ERR_EOF)
		{
			vl_sift_delete(detector);
			continue;
		}

		// Extract the keypoints of each octave
		do
		{
			const VlSiftKeypoint * keypoints;
			vl_sift_detect(detector);
			keypoints = vl_sift_get_keypoints(detector);

			// Compute the descriptors of all keypoints in the current octave
			for (int j = 0; j < vl_sift_get_nkeypoints(detector); ++j)
			{
				// First calculate the orientation of the keypoint
				// (library supports up to 4 orientations for each keypoint)
				double angles[4];
				int nangles = vl_sift_calc_keypoint_orientations(detector, angles, keypoints + j);

				// For each orientation construct a different keypoint descriptor
				for (int k = 0; k < nangles; ++k)
				{
					vl_sift_pix descriptor[128];
					vl_sift_calc_keypoint_descriptor(detector, descriptor, keypoints + j, angles[k]);
					// Save the descriptor
					descriptors.push_back(descriptor);
				}
			}
		} while (vl_sift_process_next_octave(detector) == VL_ERR_EOF);

		// Delete SIFT detector used for this image
		vl_sift_delete(detector);

		std::cout << " Extracted " << descriptors.size() << " keypoints" << std::endl;
		getchar();
	}

	return SUCCESSFUL_RUN;
}