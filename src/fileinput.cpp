#include "fileinput.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Thanks to user Evan Teran on StackExchange for
// this string split method
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

std::unordered_map<std::string, std::string> read_in_from_file(std::string fname) {
    std::ifstream inFile;
    std::unordered_map<std::string, std::string> parameterMap;
    inFile.open(fname);
    
    if (!inFile)  {
        std::cerr << "Unable to open " << fname << "\n";
        exit(1);
    }

    std::string tempLine;
    while(!inFile.eof()) {
        inFile >> tempLine;
        std::vector<std::string> tempLineSplit = split(tempLine, ':');
        parameterMap.insert({tempLineSplit[0], tempLineSplit[1]});
    }

    inFile.close();
    return parameterMap;
}