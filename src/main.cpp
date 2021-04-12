#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <iomanip>
#include <Eigen/Dense>
#include <cmath>
#include <random>
#include "constants.h"
#include "simulationbox.h"
#include "fileinput.h"
#include "fileoutput.h"
#include "thermostats.h"

using namespace Eigen;

class LennardJones
{
    private:
         Vector3d cutoffadjust;
        double sigma;
        double epsilon;
        double cutoff;         
    public:
        LennardJones(double, double, double);
        Vector3d getForce(Vector3d);
};

LennardJones::LennardJones(double insigma, double inepsilon, double incutoff=10000) {
    // Constructs the Lennard Jones class, converting to atomic units as we go
    // We expect the sigma to be in angstrom,
    // The epsilon to be in Hartrees
    // and the cutoff to also be in angstrom
 
    sigma = insigma * constants::bohrAng;
    epsilon = inepsilon * constants::boltzHar;
    cutoff = incutoff * constants::bohrAng;
    double cutoff_pair_virial = 24 * epsilon / pow(cutoff, 2);
    double cutoff_sigma_6 = pow(sigma, 6) / pow(cutoff, 6);
    double cutoff_sigma_12 = pow(cutoff_sigma_6, 2);
    double cutoffadjust_value = cutoff_pair_virial * (2 * cutoff_sigma_12 - cutoff_sigma_6);
    cutoffadjust(0) = cutoffadjust_value;
    cutoffadjust(1) = cutoffadjust_value;
    cutoffadjust(2) = cutoffadjust_value;
}

Vector3d LennardJones::getForce(Vector3d distance)
{
    /* Class function for the Lennard Jones model that takes in
     * a radial distance and outputs the force vector,
     * and the energy in relation to that */
    double radiusSquared=distance.squaredNorm();
    if (radiusSquared > pow(cutoff,2)){
        // We're outside the range we wanted to look at
        // so just return no force and carry on
        Vector3d force(0.0, 0.0, 0.0);
        return force;
    }
    double sigma6 = pow(sigma, 6) / pow(radiusSquared, 3);
    double sigma12 = pow(sigma6, 2);
    //Allen and Tildesey eq 2.59-2.63, then eq 5.3
    double pair_virial = 24.0 * epsilon / radiusSquared;
    Vector3d force = (distance * pair_virial * (2 * sigma12 - sigma6)) + cutoffadjust;
    return force;
}

void implement_boundary_conditions(ArrayXXd &positions, const Box &inputBox) {
    // Mutates the input array to fix any atoms that are outside the
    // boundary conditions
    // Fails when an atom is more than twice the boundary away,
    // but something else has gone wrong then

    for (int atom=0; atom < positions.rows(); ++atom) {
        for (int dim=0; dim < positions.cols(); ++dim) {
            if (positions(atom, dim) < 0.0) {
                positions(atom, dim) += inputBox.sides(dim);
            }
            else if (positions(atom, dim) > inputBox.sides(dim)) {
                positions(atom, dim) -= inputBox.sides(dim);
            }
        }
    }
}

ArrayXXd populate_velocities(const Box &input_box, const ArrayXd &masses, double temperature){
    // Populates the velocity array with velocities sampled
    // randomly from a Gaussian distribution with mean 0 and
    // variance 1.
    
    // This is then rescaled to match the temperature before
    // it is used


    std::random_device device;
    std::mt19937 mt_rand(device());
    
    // This is the normal distribution generator
    double mean = 0.0;
    double var = 1.0;
    std::normal_distribution<double> gaussian(mean, var);

    ArrayXXd velocities = ArrayXXd::Zero(masses.rows(), input_box.dimensions);
    // Use the number of rows in the masses as a proxy for 
    // number of atoms
    for (int atom=0; atom < masses.rows(); ++atom) {
        for (int dim=0; dim < input_box.dimensions; ++dim) {
            velocities(atom, dim) = gaussian(mt_rand);
        }
    }
    
    auto rescaler = VelocityRescale();
    rescaler.apply(velocities, masses, temperature);
    return velocities;
}

ArrayXXd populate_positions(const Box &inputBox, int num_x, int num_y, int num_z) {
    /* Fills the position array with atoms evenly spaced in the box
     * starting from (0, 0) */
    ArrayXXd positions(num_x * num_y * num_z, inputBox.dimensions);
    double x_increment = inputBox.sides(0) / static_cast<double>(num_x);
    double y_increment = inputBox.sides(1) / static_cast<double>(num_y);
    double z_increment = inputBox.sides(2) / static_cast<double>(num_z);
    int counter = 0;
    for (int row=0; row < num_x; ++row) {
        for (int col=0; col < num_y; ++col) {
            for (int plane=0; plane < num_z; ++plane) {
                positions(counter, 0) = row * x_increment;
                positions(counter, 1) = col * y_increment;
                positions(counter, 2) = plane * z_increment;
                ++counter;
            }
        }
    }
    return positions;
}

Vector3d implement_minimum_image(const Vector3d distance, const Box &box) {
    auto new_distance = distance;
    for (int dim=0; dim < box.dimensions; ++dim) {
        if (distance(dim) <= -0.5 * box.sides(dim)) {
            new_distance(dim) += box.sides(dim);
        } else if (distance(dim) > 0.5 * box.sides(dim)) {
            new_distance(dim) -= box.sides(dim);
        }
    }   
    return new_distance;
}

ArrayXXd calculate_forces(ArrayXXd &atoms_array, const Box &box, LennardJones &potential_model) {
    /* Calculates the force according to a Lennard-Jones potential
     * from one atom to its nearest image, and returns
     * an array of these forces */
    ArrayXXd forces = ArrayXXd::Zero(atoms_array.rows(), atoms_array.cols());
    for (int i=0; i < atoms_array.rows(); ++i) {
        for (int j=0; j<i; ++j) {
            Vector3d distance = atoms_array.row(i) - atoms_array.row(j);
            distance = implement_minimum_image(distance, box);
            Vector3d instant_force(potential_model.getForce(distance));

            for (int dim=0; dim < box.dimensions; ++dim) {
                forces(i, dim) += instant_force(dim);
                forces(j, dim) -= instant_force(dim);
            }
        }
    }
    return forces;
}

void velocityVerlet(ArrayXXd &positions, ArrayXXd &velocities, double timestep, ArrayXd &masses,
                    const Box &box, LennardJones &potentialModel) {
    ArrayXXd accelerations = calculate_forces(positions, box, potentialModel).colwise() / masses;
    ArrayXXd positions_next = positions 
                              + (velocities * timestep)
                              + (0.5 * accelerations * pow(timestep, 2));
    ArrayXXd accelerations_next = calculate_forces(positions_next, box, potentialModel).colwise() / masses;
    ArrayXXd velocities_next = velocities 
                               + 0.5 * ( accelerations + accelerations_next) * timestep;
    velocities = velocities_next;
    implement_boundary_conditions(positions_next, box);
    positions = positions_next;
}

int main() {
    auto parameter_map = read_in_from_file("input.inpt");
    LennardJones potentialModel(std::stod(parameter_map["LJ_sigma"]),
                                std::stod(parameter_map["LJ_epsilon"]),
                                std::stod(parameter_map["LJ_cutoff"]));
    Box simulation_box(std::stod(parameter_map["length_x"]),
                       std::stod(parameter_map["length_y"]),
                       std::stod(parameter_map["length_z"]));
    int num_x = stoi(parameter_map["num_x"]);
    int num_y = stoi(parameter_map["num_y"]);
    int num_z = stoi(parameter_map["num_z"]);
    double temperature = std::stod(parameter_map["temperature"]);
    double timestep = std::stod(parameter_map["timestep"]);
    int rescaleStep=std::stoi(parameter_map["num_rescale"]);
    int outputStep=std::stoi(parameter_map["num_output"]);
    int numSteps=std::stoi(parameter_map["num_steps"]);
    
    ArrayXXd atoms = populate_positions(simulation_box, num_x, num_y, num_z);
    ArrayXd masses = ArrayXd::Zero(atoms.rows(), 1);
    masses += std::stod(parameter_map["mass"]);
    
    ArrayXXd velocities = populate_velocities(simulation_box, masses, temperature);
 
    clearOutputFiles(parameter_map["position_file"]);

    for (int step=0; step < numSteps; ++step) {
        velocityVerlet(atoms, velocities, timestep, masses, simulation_box, potentialModel);
        if (step % rescaleStep == 0) {
            rescale_velocities(velocities, masses, temperature);
        }

        if (step % outputStep == 0) {
            dumpToFile(atoms, parameter_map["position_file"], step);
            std::cout << "Step " << step << std::endl;
        }
    }
}
