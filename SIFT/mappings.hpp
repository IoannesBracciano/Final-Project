/* Asset mappings */
#pragma once
#include <unordered_map>

/// Holds all the mappings between the query assets and the database assets
/// using the numerical part of their filenames as identifiers
const std::unordered_map <std::string, std::list<std::string>> mappings =
		std::unordered_map <std::string, std::list<std::string>>({
			{ "1",	{ "12" } },
			{ "2",	{ "13","14" } },
			{ "3",	{ "50" } },
			{ "4",	{ "59" } },
			{ "5",	{ "55","58" } },
			{ "6",	{ "57" } },
			{ "7",	{ "62" } },
			{ "8",	{ "63","64" } },
			{ "9",	{ "66" } },
			{ "10", { "67" } },
			{ "11", { "69","70" } },
			{ "12", { "71","72" } },
			{ "13", { "73" } },
			{ "14", { "76" } },
			{ "15", { "77" } },
			{ "16", { "82" } },
			{ "17", { "86" } },
			{ "18", { "86" } },
			{ "19", { "85" } },
			{ "20", { "89" } },
			{ "21", { "92","94" } },
			{ "22", { "94" } },
			{ "23", { "97" } },
			{ "24", { "96" } },
			{ "25", { "105","106" } },
			{ "26",	{ "107" } },
			{ "27",	{ "109","111" } },
			{ "28",	{ "108" } },
			{ "29",	{ "110","112" } },
			{ "30",	{ "116" } },
			{ "31",	{ "120" } },
			{ "32",	{ "124","125" } },
			{ "33",	{ "128" } },
			{ "34",	{ "129" } },
			{ "35",	{ "135" } },
			{ "36",	{ "138" } },
			{ "37",	{ "140" } }
		});

/// <summary>Checks whether or not the first asset id maps the second</summary>
/// <param name="id1">the id of the first asset</param>
/// <param name="id2">the id of the second asset</param>
/// <returns>true if there is a match in mappings, false othrewise</returns>
bool maps(const std::string & id1, const std::string & id2)
{
	std::list<std::string> ids = mappings.at(id1);
	if (!ids.empty() &&
		std::find(ids.begin(), ids.end(), id2) != ids.end())
	{
		return true;
	}
	else return false;
}
