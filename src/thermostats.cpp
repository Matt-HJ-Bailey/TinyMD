#include <Eigen/Dense>
#include <cmath>
#include <random>
#include "thermostats.h"
#include "constants.h"

std::random_device g_device;
std::mt19937 g_mt_rand(g_device());
std::uniform_real_distribution<double> distribution(0.0, 1.0);
  
Thermostat::Thermostat(double aim_temp) {
    m_aim_temp = aim_temp;
}

VelocityRescale::VelocityRescale(int steps, double aim_temp) {
    m_steps = steps;
    m_aim_temp = aim_temp;
}

Berendsen::Berendsen(double timescale, double aim_temp) {
    m_timescale = timescale;
    m_aim_temp = aim_temp;
}

Andersen::Andersen(double freq, double timestep, double aim_temp) {
    // The user provides us with a timestep in ps
    // and a collision frequency in ps^-1
   
    m_cutoff = freq * timestep;
    m_aim_temp = aim_temp;
}

void Berendsen::apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses) {
    auto velocities_squared = (velocities * velocities).rowwise().sum();
    auto effective_temp = (0.5 * masses * velocities_squared).sum() / constants::boltzHar;
    auto delta_T = effective_temp - m_aim_temp / m_timescale;
}
      
void VelocityRescale::apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses) {
    auto velocities_squared = (velocities * velocities).rowwise().sum();
    auto translational_energy = (0.5 * masses * velocities_squared).sum();
    double thermal_energy = constants::boltzHar * m_aim_temp;
    double velocity_factor = sqrt(thermal_energy / translational_energy);
    velocities *= velocity_factor;
}

void Andersen::apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses) {
    //Implements an Andersen thermostat. This keeps the system at
    //constant temperature by the following method:
    //    1) For each atom, calculate a "collision chance" in the last
    //       timestep. If this chance > collision frequency * timestep,
    //       then progress to 2
    //    2) Simulate a "collision" with an external heat bath. This
    //       randomly selects the velocity from a Maxwell-Boltzmann distribution


    for (int atom = 0; atom < velocities.rows(); ++atom) {
        double random_num = distribution(g_mt_rand);
        
        if (random_num > m_cutoff) {
            // The lucky atom is hit!
            double mean = 0.0;
            double variance = std::sqrt(constants::boltz * m_aim_temp / masses(atom));
            std::normal_distribution<double> gaussian(mean, variance);
            for (int dim=0; dim < velocities.cols(); ++dim) {
                velocities(atom, dim) = gaussian(g_mt_rand);
            }
        }
    }
}