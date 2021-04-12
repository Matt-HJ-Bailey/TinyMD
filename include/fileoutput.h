#pragma once
#include <string>
#include <Eigen/Dense>

void dumpToFile(Eigen::ArrayXXd &array, const std::string fname, const int step);

void clearOutputFiles(std::string fname);