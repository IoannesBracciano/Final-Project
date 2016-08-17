/* Dataset related utilities */
#pragma once
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

/// Represents the types of assets in the dataset.
/// The database is a collection of image and SIFT descriptor files.
/// SIFT descriptors are extracted from query images (QUERY_ASSETs)
/// and are matched against a database of other descriptors
/// already extracted from a different set of images (DB_ASSETs)
typedef enum { QUERY_ASSET, DB_ASSET, UNKNOWN } AssetType;

/// Regular expression to match query asset filenames
const boost::regex query_asset_filename_pattern{ "q_img\\([0-9]+\\).*" };
/// Regular expression to match database asset filenames
const boost::regex db_asset_filename_pattern{ "db_img\\([0-9]+\\).*" };
// Regular expression to extract asset ids (the numeric part of their filename)
const boost::regex asset_id_extractor{ "([0-9]+)" };


/// <summary>Finds the asset type based on its filename</summary>
/// <param name="path">a valid path to the file</param>
/// <returns>the type of the asset</returns>
AssetType get_asset_type(const boost::filesystem::path & path)
{
	if (path.has_filename())
	{
		if (boost::regex_match(path.filename().string(), query_asset_filename_pattern))
		{
			return QUERY_ASSET;
		}
		else if (boost::regex_match(path.filename().string(), db_asset_filename_pattern))
		{
			return DB_ASSET;
		}
	}

	// If no filename was found or no pattern was matched
	return UNKNOWN;
}


/// <summary>Extracts the asset id from its filename</summary>
/// <param name="path">a valid path to the file</param>
/// <returns>the id of the asset</returns>
std::string get_asset_id(const boost::filesystem::path & path)
{
	if (path.has_filename())
	{
		std::string filename = path.filename().string();
		// Search for the id in the filename
		boost::smatch matches;
		boost::regex_search(filename, matches, asset_id_extractor);

		if (matches.size() > 1)
		{
			return matches[1].str();
		}
	}

	return "";
}