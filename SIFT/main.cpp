#include <iostream>
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"
#include "boost/regex.hpp"
#include "../Filter/generic.hpp"
#include "log.hpp"

extern "C" {
#include "vl/sift.h"
}

// Holds the description of the program options
p_opts::options_description p_opts_desc("Arguments description");

// Holds positional description of program options
p_opts::positional_options_description p_opts_pos;

// Array used to map variable names with their corresponding values
p_opts::variables_map p_opts_vm;

// Regular expressions for query and database image file names
boost::regex query_images{ "q_img\([0-9]+\).*" };
boost::regex base_images{ "db_img\([0-9]+\).*" };

/// <summary>Sets up the program options and parses the argv vector accordingly
/// </summary>
/// <param name="argc">the arguments count</param>
/// <param name="argc">the standard c++ arguments vector</param>
void setup_p_opts(int argc, char ** argv)
{
	p_opts_desc.add_options()
		("help", "Print this help message")
		("input,i", p_opts::value<std::vector<std::string>>()
			->default_value(std::vector<std::string>(), "current working directory"),
			"Absolute or relative path to the input image or descriptors file(s), "
			"or a directory containing image or descriptor files")
		("output,o", p_opts::value<std::string>()
			->default_value(std::string(), "input path directory"),
			"Absolute or relative path to an existing directory "
			"where extracted descriptor files are saved")
		("match,m", p_opts::value<std::vector<std::string>>(),
			"A path to a descriptor file or a directory containing descriptor files "
			"that the input image descriptors will be matched against")
		("log,l", p_opts::value<std::string>(),
			"A path to a directory where the log will be output")
		("peak-thresh,t", p_opts::value<double>()->default_value(0.01),
			"Minimun amount of contrast to accept a keypoint."
			"Increasing this value will result in fewer keypoints")
		("edge-thresh,T", p_opts::value<double>()->default_value(10),
			"Edge rejection threshold. Increasing this will result is fewer keypoints")
		("octaves", p_opts::value<int>()->default_value(-1),
			"Number of octaves")
		("levels", p_opts::value<int>()->default_value(3),
			"Number of levels per octave");
	p_opts_pos
		.add("input", -1);

	// Try to store the command line arguments into the options variable map,
	// but do not notify yet
	p_opts::store(p_opts::command_line_parser(argc, argv)
		.options(p_opts_desc)
		.positional(p_opts_pos).run(), p_opts_vm);
}

/// <summary>
/// Checks if given path is a supported descriptors file.
/// Supported descriptors files are created using SIFT.exe
/// and they are simply 128 floating point numbers (default SIFT descriptor size)
/// written sequentially in binary form
/// </summary>
/// <param name="path">the path of the file to be examined</param>
/// <returns>true if file passes the filter, false if not</returns>
bool descriptors_file_filter(const fsys::path & path)
{
	return path.extension() == ".sifd";
}

/// <summary>
/// Combined image and descriptors file filter
/// </summary>
/// <param name="path">the path of the file to be examined</param>
/// <returns>true if file passes the filter, false if not</returns>
bool extra_supported_file_filter(const fsys::path & path)
{
	return image_file_filter(path) || descriptors_file_filter(path);
}

/// <summary>
/// Loads a descriptors file into memory
/// </summary>
/// <param name="path">a valid path to the descriptors file</param>
/// <param name="vector">a pointer to the first block of memory where the descriptors are loaded</param>
/// <returns>number of descriptors in the vector</returns>
int load_descriptors_file(const fsys::path & path, vl_sift_pix * & vector)
{
	// Getting the file size in bytes
	uintmax_t file_size = fsys::file_size(path);
	// Creating a vector of the same size to store the descriptors
	vector = (vl_sift_pix*)malloc(file_size);
	// Opening the file in binary mode for reading
	fsys::ifstream file(path, std::ios::in | std::ios::binary);
	// Reading the file into the previously created vector
	file.read(reinterpret_cast<char*>(vector), file_size);
	// Closing file stream
	file.close();

	return file_size / (128*sizeof(vl_sift_pix));
}

/// <summary>
/// Opens the image found under path, extracts its SIFT descriptors
/// and saves the into the given descriptors vector
/// </summary>
/// <param name="octaves">number of octaves (SIFT specific parameter)</param>
/// <param name="levels">number of levels per octave (SIFT specific parameter)</param>
/// <param name="peak_thresh">minimun amount of contrast to accept a keypoint (SIFT specific parameter)</param>
/// <param name="edge_thresh">edge rejection threshold (SIFT specific parameter)</param>
/// <returns>execution time in seconds</returns>
double extract(fsys::path & path, std::vector<vl_sift_pix*> & descriptors,
	int octaves, int levels, double peak_thresh, double edge_thresh)
{
	// Open image file
	cv::Mat I = cv::imread(path.string(), cv::IMREAD_GRAYSCALE);
	                                       // SIFT works only on gray scale images
	// Normalize pixel values in the range [0,1]
	I.convertTo(I, CV_32F, 1.0 / 255.5);

	// Start clock
	int64 t = cv::getTickCount();

	// Create a new SIFT detector
	VlSiftFilt * detector = vl_sift_new(I.cols, I.rows, octaves, levels, 0);

	// Setting peak and edge thresholds
	vl_sift_set_peak_thresh(detector, peak_thresh);
	vl_sift_set_edge_thresh(detector, edge_thresh);

	// Return if there is no octaves to calculate
	if (vl_sift_process_first_octave(detector, (float*)I.data) == VL_ERR_EOF)
	{
		vl_sift_delete(detector);
		// Return -1 to indicate failure
		return -1;
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
				vl_sift_pix * descriptor = (vl_sift_pix*) malloc(sizeof(vl_sift_pix)*128);
				vl_sift_calc_keypoint_descriptor(detector, descriptor, keypoints + i, angles[j]);
				descriptors.push_back(descriptor);
			}
		}
	} while (vl_sift_process_next_octave(detector) != VL_ERR_EOF);

	// Calculate execution time in seconds
	double dt = double(cv::getTickCount() - t) / double(cv::getTickFrequency());

	// Delete SIFT detector used for this image
	vl_sift_delete(detector);

	return dt;
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
/// <param name="pairs">pointer to a vector to be filled with the matched pairs.
///                     Needs to be big enough to hold all possible pairs</param>
/// <param name="descr1">pointer to the descriptor matrix of the first image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="descr2">pointer to the descriptor matrix of the second image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="K1">Number of descriptors of first image</param>
/// <param name="K2">Number of descriptors of second image</param>
/// <param name="ND">Descriptor size (=128 for SIFT)</param>
/// <param name="thresh">Ratio test threshold value (=1.5)</param>
/// <returns>a pointer pointing to the end of the vector conaitning the matched pairs</returns>
Pair * compare(Pair * pairs, const float *descr1, const float *descr2,
	int K1, int K2, int ND = 128, float thresh = 1.5)
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
/** End of ported code section */


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
	catch (p_opts::required_option & e)
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
			if (p_opts_vm.count("match"))
			{
				// Input is also allowed to be a descriptors file if a match path is given
				resolvePath(input[i], files, extra_supported_file_filter);
			}
			else
			{
				// Input must be a supported image file only to extract its descriptors
				resolvePath(input[i], files);
			}
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

	// Resolving match files (if option was commanded)
	std::vector<fsys::path> match_files = std::vector<fsys::path>();
	if (p_opts_vm.count("match"))
	{
		// Extracting match program option
		std::vector<std::string> match = p_opts_vm["match"].as<std::vector<std::string>>();
		for (int i = 0; i < match.size(); ++i)
		{
			resolvePath(match[i], match_files, descriptors_file_filter);
			                                   // match files can only be descriptors files
		}
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

	// Resolving log path
	fsys::path log;
	if (p_opts_vm.count("log"))
	{
		log = p_opts_vm["log"].as<std::string>();
		if (!fsys::exists(log) || !fsys::is_directory(log))
		{
			std::cout << "Specify an existing directory for the log file" << std::endl;
			return INVALID_LOG;
		}
	}

	// Extracting peak thresh
	double peak_thresh = p_opts_vm["peak-thresh"].as<double>();
	// Extracting edge threshold
	double edge_thresh = p_opts_vm["edge-thresh"].as<double>();
	// Extracting octaves
	int octaves = p_opts_vm["octaves"].as<int>();
	// Extracting levels of each octave
	int levels = p_opts_vm["levels"].as<int>();

	// Time and statistical metrics
	unsigned int descr_mean_count = 0;
	double exec_mean_time = 0;

	// Looping through each input file
	for (int i = 0; i < files.size(); ++i)
	{
		// Program was called to match descriptor files
		if (p_opts_vm.count("match"))
		{
			// Open the log file if specified in the program options
			if (p_opts_vm.count("log"))
			{
				try
				{
					log_init(log / "log.txt");
				}
				catch (std::exception e)
				{
					std::cout << "Could not open log file: " << std::endl
						<< e.what() << std::endl;

					return INVALID_LOG;
				}
			}

			// Declaring pointers to the descriptors to be matched
			vl_sift_pix * descr1 = nullptr, * descr2 = nullptr;
			int n_descr1 = 0, n_descr2 = 0;

			// Opening input file
			if (image_file_filter(files[i]))		    // If input file is image
			{
				std::vector<vl_sift_pix*> descriptors = std::vector<vl_sift_pix*>();
				extract(files[i], descriptors, octaves, levels, peak_thresh, edge_thresh);
				n_descr1 = descriptors.size();
				descr1 = (vl_sift_pix*)malloc(n_descr1 * 128 * sizeof(vl_sift_pix));
				for (int j = 0; j < n_descr1; ++j)
				{
					for (int k = 0; k < 128; ++k)
					{
						descr1[j * 128 + k] = descriptors[j][k];
					}
				}
			}
			else if (descriptors_file_filter(files[i])) // else if input file is a descriptors file
			{
				// load descriptors from file
				n_descr1 = load_descriptors_file(files[i], descr1);
			}

			if (p_opts_vm.count("log"))
			{
				if (boost::regex_match(files[i].filename().string(), query_images))
				{
					log_write("Compairing query image " + files[i].filename().string() + " against...");
				}
				else if (boost::regex_match(files[i].filename().string(), base_images))
				{
					log_write("Compairing database image " + files[i].filename().string() + " against...");
				}
			}

			// Match metrics
			unsigned int max_pairs = 0;
			fsys::path best_match;
			double total_time = 0;

			

			// Looping through files to match the input descriptors against
			for (int j = 0; j < match_files.size(); ++j)
			{
				// load descriptors from matching file
				n_descr2 = load_descriptors_file(match_files[j], descr2);

				std::cout << "Compairing " << files[i].filename() <<
					" against " << match_files[j].filename() << "..." << std:: endl;

				if (p_opts_vm.count("log"))
				{
					log_write("\t- " + match_files[j].filename().string() + ": ", false);
				}

				// compare the two descriptors vectors
				Pair * pairs = (Pair*)malloc(sizeof(Pair) * (n_descr1 > n_descr2 ? n_descr1 : n_descr2));
				Pair * pairs_first = pairs;
				double tic = cv::getTickCount();
				Pair * pairs_end = compare(pairs, descr1, descr2, n_descr1, n_descr2);
				double tac = double(cv::getTickCount() - tic) / double(cv::getTickFrequency());
				int n_pairs = pairs_end - pairs_first;

				std::cout << "\t" << n_pairs << " pairs found in " << tac << " seconds" << std::endl;

				if (p_opts_vm.count("log"))
				{
					log_write(boost::lexical_cast<std::string>(n_pairs) + " pairs found in " +
						boost::lexical_cast<std::string>(tac) + " seconds");
				}

				if (n_pairs > max_pairs)
				{
					max_pairs = n_pairs;
					best_match = match_files[j];
				}

				total_time += tac;

				delete[] pairs_first;
				delete[] descr2;
			}

			if (match_files.size() > 1)
			{
				std::cout << "Done in " << total_time << " seconds" << std::endl;
				std::cout << "\t*** Best match: " << best_match.filename()  << " ***" << std::endl;

				if (p_opts_vm.count("log"))
				{
					log_write("Done in " + boost::lexical_cast<std::string>(total_time) + " seconds");
					log_write("Best match: " + best_match.filename().string());
					log_end();
				}
			}

			delete[] descr1;
		}
		
		else // Program was called to extract descriptors from images
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
				if (!replace_file_dialog(output / new_filename))
					continue;

			std::cout << "Extracting keypoints from " << files[i].filename() << "...";
			double t = extract(files[i], descriptors, octaves, levels, peak_thresh, edge_thresh);
			std::cout << " Extracted " << descriptors.size() << " keypoints " <<
				"in "<< t << " seconds" << std::endl;

			// Updating metrics
			descr_mean_count += descriptors.size() / files.size();
			exec_mean_time += t / files.size();

			// Create a binary file and write the descriptors in raw binary format
			fsys::ofstream ofs((output / new_filename), std::ios::out | std::ios::binary);
			for (std::vector<vl_sift_pix*>::iterator it = descriptors.begin(); it != descriptors.end(); ++it)
			{
				vl_sift_pix * descriptor = *it;
				ofs.write(reinterpret_cast<const char *>(descriptor), sizeof(vl_sift_pix) * 128);
			}
			ofs.close();
		}
	}

	// Log
	if (!p_opts_vm.count("match"))
	{
		std::cout << "Mean descriptors count per image: " << descr_mean_count << std::endl;
		std::cout << "Mean extraction time per image: " << exec_mean_time << std::endl;
	}

	return SUCCESSFUL_RUN;
}