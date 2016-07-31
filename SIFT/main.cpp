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
		//("action", p_opts::value<std::string>()->required(),
		//	"[extract|match] Which action to perform")
		("input,i", p_opts::value<std::vector<std::string>>()
			->default_value(std::vector<std::string>(), "current working directory"),
			"Absolute or relative path to the input image file or files, or a directory containing image files")
		("output,o", p_opts::value<std::string>()
			->default_value(std::string(), "input path directory"),
			"Absolute or relative path to an existing directory where filtered images are saved")
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

/// <summary>
/// Checks if given path is a supported image file.
/// The list of supported OpenCV 3 image files, was taken from
/// http://docs.opencv.org/3.0-beta/modules/imgcodecs/doc/reading_and_writing_images.html#imread
/// </summary>
/// <param name="path">the path of the file to be examined</param>
bool is_supported_image_file(const fsys::path & path)
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
	                                   // ;@@@@@@@@@@@@@@@@@@@

	// Delete SIFT detector used for this image
	vl_sift_delete(detector);
}

/** The next portion of code was ported from the matlab interface of the vlfeat library
    as there was no suitable code for matching between two image descriptors and was copied
	and pasted almost "as is".
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
	int ND=128,
	float thresh=1.5
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