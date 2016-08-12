#pragma once
#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/date_time.hpp"

// Holds the log stream
static boost::filesystem::ofstream logofstr;

/// <summary>Initiates the log file</summary>
/// <param name="path">the path to the log file</param>
/// <param name="stream">an uninitialized output file stream</param>
void log_init(const boost::filesystem::path & path)
{
	logofstr.open(path);
	logofstr << boost::posix_time::to_iso_string(
		boost::posix_time::second_clock::local_time())
		<< "\tsift.exe Init log file" << std::endl;
}

/// <summary>Writes the content in the log file</summary>
/// <param name="content">the content to be written to the file</param>
/// <param name="endl">if true, a line ending is appended after the content</param>
void log_write(const std::string content, bool endl = true)
{
	if (!logofstr.is_open()) return;
	logofstr << content;
	if (endl) logofstr << std::endl;
}

void log_end()
{
	logofstr << std::endl  << "End of log\t" << boost::posix_time::to_iso_string(
		boost::posix_time::second_clock::local_time()) << std::endl;
	logofstr.close();
}