#pragma once
#include <Eigen/Dense>

class Box
{
    public:
        Box(double, double, double);
        Eigen::Vector3d sides;
        int dimensions;
        double volume();
};