#include <string>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include "fileoutput.h"


void dumpToFile(Eigen::ArrayXXd &array, const std::string fname, const int step) {
    std::ofstream outfile;
    outfile.open(fname, std::ios_base::app);
    outfile << array.rows() <<"\n";
    outfile << "=====STEP " << step << "=====\n";
    outfile << std::scientific;
    outfile << std::setprecision(6);
    for (int i=0; i < array.rows(); ++i) {
        outfile << "Ar\t" << array(i, 0) << "\t" << array(i, 1) << "\t" << array(i, 2) << "\n";
    }
}

void clearOutputFiles(std::string fname) {
    std::ofstream tmpfile;
    tmpfile.open(fname, std::ios::out | std::ios::trunc);
    tmpfile.close();
}