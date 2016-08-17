/* Statistical analysis tools */
#pragma once
#include "boost/functional/hash.hpp"

/// <summary>
/// Custom key hasher for string pairs. Standard library
/// does not provide one by default
/// </summary>
/*
struct StringPairKeyHasher
{
	std::size_t operator()(const std::pair<std::string, std::string> & key) const
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, boost::hash_value(key.first));
		boost::hash_combine(seed, boost::hash_value(key.second));
		return seed;
	}
};
*/

typedef std::unordered_map <std::pair<std::string, std::string>, std::vector<Pair>,
	boost::hash<std::pair<std::string, std::string>>> AssetDescrPairsMap;

/// Holds the pairs found between two assets
AssetDescrPairsMap asset_descr_pairs;

/// Short name for iterator that iterates over asset descriptor pairs
typedef AssetDescrPairsMap::iterator pairs_iterator;

/// Holds the registered matched assets
std::unordered_map <std::string, std::string> asset_matches;

/// Short name for iterator that iterates over asset matches
typedef std::unordered_map<std::string, std::string>::iterator matches_iterator;

/// <summary>
/// Registers the SIFT descriptor pairs found between two assets
/// </summary>
/// <param name="id1">the id of the first asset</param>
/// <param name="id1">the id of the second asset</param>
/// <param name="pairs_begin">pointer to the beginning of an array of asset descriptor pairs</param>
/// <param name="pairs_end">pointer after the last element of the asset descriptor pairs array</param>
void register_asset_descr_pairs(const std::string id1, const std::string id2,
		const Pair * pairs_begin, const Pair * pairs_end)
{
	asset_descr_pairs[std::pair<std::string, std::string>(id1, id2)] =
		std::vector<Pair>(pairs_begin, pairs_end);
}

const std::vector<Pair> & get_asset_descr_pairs(std::string id1, std::string id2)
{
	return asset_descr_pairs[std::pair<std::string, std::string>(id1, id2)];
}

/// <summary>Finds the pair of assets with the highest number of descriptor pairs</summary>
/// <returns>a pair of asset ids with the highest number of descriptor pairs</returns>
const std::pair<std::string, std::string> get_best_match(const std::string & id1)
{
	std::pair<std::string, std::string> best_match;
	unsigned int max_pairs_count = 0;

	// For every registered asset descriptor pairs
	for (pairs_iterator it = asset_descr_pairs.begin(); it != asset_descr_pairs.end(); ++it)
	{
		// Check if the assets have the highest number of pairs,
		// and declare them as best match if yes. This could be changed
		// to match different criteria apart from highest number of pairs
		// (maybe pass a pointer to a criteria function as argument?)
		std::pair<std::pair<std::string, std::string>, std::vector<Pair>> pairs = *it;
		if (pairs.first.first == id1 &&
				pairs.second.size() > max_pairs_count)
		{
			best_match = pairs.first;
			max_pairs_count = pairs.second.size();
		}
	}

	return best_match;
}

/// <summary>
/// Registers a match between (query) asset with id1
/// and (database) asset with id2. The ids must correspond
/// to the numerical part of the asset filenames and there
/// is no check whether they are correct or not
/// </summary>
/// <param name="id1">the id of the first asset</param>
/// <param name="id1">the id of the second asset</param>
void register_asset_match(const std::string id1, const std::string id2)
{
	asset_matches[id1] = id2;
}

/// <summary>Counts correct registered matches, according to mappings</summary>
/// <returns>the correct match count</returns>
unsigned int get_correct_match_count()
{
	unsigned int counter = 0;

	// For every registered match
	for (matches_iterator it = asset_matches.begin(); it != asset_matches.end(); ++it)
	{
		std::pair<std::string, std::string> match = *it;
		// Check if there is a corresponding mapping and increment the counter
		if (maps(match.first, match.second)) ++counter;
	}

	return counter;
}

/// <summary>
/// Calculates the percentage of correct matches as the ratio of the number
/// of correct matches found to the number of all matches registered
/// </summary>
/// <returns>the correct match percentage</returns>
double get_correct_match_percentage()
{
	return (double(get_correct_match_count()) / double(asset_matches.size())) * 100;
}