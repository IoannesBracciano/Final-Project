#pragma once

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

namespace p_opts = boost::program_options;
namespace fsys = boost::filesystem;

/// <summary>Enumeration of program return codes</summary>
const enum RETURN_CODES
{
	// User did not provide a required program option
	REQUIRED_OPT_NOT_SET = -1,
	// Input path is invalid
	INVALID_INPUT = -2,
	// Output path is invalid
	INVALID_OUTPUT = -3,

	// Program terminated successfully
	SUCCESSFUL_RUN = 0
};

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
void resolveInput(const fsys::path & path, std::vector<fsys::path> & files,
	bool(*file_ext_filter)(const fsys::path&) = is_supported_image_file)
{
	//std::cout << "Input: " << path << std::endl;
	if (fsys::exists(path))									// If path exists
	{
		if (fsys::is_regular_file(path))					// and is a regular file
		{
			if (file_ext_filter(path))						// and it is a supported image file
			{
				files.push_back(path);						// add it to the files vector
			}
			else // Unsupported file or not an image file at all
			{
				throw std::runtime_error("Unsupported file extension");
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

///<summary>
/// Presents a dialog asking the user if he want to replace
/// the already existing file or not. User can choose to perform
/// his selected action on all files, and the dialog won't repeat
/// itself whenever the function is called again
/// </summary>
/// <param name="path">the path to the already existing file</param>
bool replace_file_dialog(const fsys::path & path)
{
	static char replace_file = 0;	// making this static to remember user's action
                                    // in case he selects to perform it on all files
	static char perform_on_all = 0;

	if (replace_file == 0)
	{
		std::cout << (path) << " exists. Do you want to replace it? (y/n, default no):";
		do { std::cin >> replace_file; } while (replace_file == '\n');

		while (replace_file != 'n' && replace_file != 'N' &&
			replace_file != 'y' && replace_file != 'Y')
		{
			std::cout << "Didn't get that, type 'y' for yes or 'n' for no:";
			do { std::cin >> replace_file; } while (replace_file == '\n');
		}

		std::cout << "Perform this action for all existing files? (y/n, default no):";
		do { std::cin >> perform_on_all; } while (perform_on_all == '\n');

		while (perform_on_all != 'n' && perform_on_all != 'N' &&
			perform_on_all != 'y' && perform_on_all != 'Y')
		{
			std::cout << "Didn't get that, type 'y' for yes or 'n' for no:";
			do { std::cin >> perform_on_all; } while (perform_on_all == '\n');
		}
	}

	if (replace_file == 'n' || replace_file == 'N')
	{
		if (perform_on_all == 'n' || perform_on_all == 'N')
		{
			replace_file = 0;
		}
		return false;
	}

	return true;
}