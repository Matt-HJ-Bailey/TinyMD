#pragma once
#include <string>
#include <vector>
#include <unordered_map>

// Thanks to user Evan Teran on StackExchange for
// this string split method
template<typename Out>
void split(const std::string &s, char delim, Out result);

std::vector<std::string> split(const std::string &s, char delim);

std::unordered_map<std::string, std::string> read_in_from_file(std::string fname);