#include <iostream>
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"
#include "boost/regex.hpp"
#include "../Filter/generic.hpp"

#include "sift.hpp"

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
				// Extract descriptors
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
				if (get_asset_type(files[i]) == QUERY_ASSET)
				{
					log_write("Compairing query asset " + get_asset_id(files[i]) + " against...");
				}
				else if (get_asset_type(files[i]) == DB_ASSET)
				{
					log_write("Compairing database asset " + get_asset_id(files[i]) + " against...");
				}
				else
				{
					log_write("Compairing asset " + files[i].filename().string() + " against...");
				}
			}

			/** TEMPORARY FILE FOR EASE OF COPYING TO EXCEL - TO BE REMOVED **/
			//fsys::ofstream tmplog(log / "excel.txt", std::ios::out);

			// Match metrics
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
				Pair * pairs_begin = (Pair*)malloc(sizeof(Pair) * (n_descr1 > n_descr2 ? n_descr1 : n_descr2));
				double tic = cv::getTickCount();
				Pair * pairs_end = compare(pairs_begin, descr1, descr2, n_descr1, n_descr2);
				double tac = double(cv::getTickCount() - tic) / double(cv::getTickFrequency());

				register_asset_descr_pairs(get_asset_id(files[i]), get_asset_id(match_files[j]),
					pairs_begin, pairs_end);

				unsigned int n_pairs = get_asset_descr_pairs(get_asset_id(files[i]), get_asset_id(match_files[j])).size();

				std::cout << "\t" << n_pairs << " pairs found in " << tac << " seconds" << std::endl;
				//tmplog << boost::lexical_cast<std::string>(n_pairs) << std::endl;

				if (p_opts_vm.count("log"))
				{
					log_write(boost::lexical_cast<std::string>(n_pairs) + " pairs found in " +
						boost::lexical_cast<std::string>(tac) + " seconds");
				}

				total_time += tac;

				delete[] pairs_begin;
				delete[] descr2;
			}

			if (match_files.size() > 1)
			{
				std::pair<std::string, std::string> best_match = get_best_match(get_asset_id(files[i]));
				std::cout << "Done in " << total_time << " seconds" << std::endl;
				std::cout << "\t*** Best asset match pair: " << best_match.first << " => " << best_match.second << " ***" << std::endl;
				register_asset_match(best_match.first, best_match.second);

				if (p_opts_vm.count("log"))
				{
					log_write("Done in " + boost::lexical_cast<std::string>(total_time) + " seconds");
					log_write("Best asset match pair: " + best_match.first + " => " + best_match.second);
					log_end();
				}
			}

			//tmplog.close();

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
	if (p_opts_vm.count("match"))
	{
		std::cout << "Correct matches percentage: " << get_correct_match_percentage() << "%" << std::endl;
	}
	else
	{
		std::cout << "Mean descriptors count per image: " << descr_mean_count << std::endl;
		std::cout << "Mean extraction time per image: " << exec_mean_time << std::endl;
	}

	return SUCCESSFUL_RUN;
}