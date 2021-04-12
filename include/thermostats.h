#pragma once

#include <Eigen/Dense>

class Thermostat{
    private:
        double m_aim_temp;
    public:
        Thermostat(double aim_temp);
        virtual void apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses, double temperature) = 0;
};

class VelocityRescale: public Thermostat {
    private:
        int m_steps;
        double m_aim_temp;
    public:
        VelocityRescale(int steps, double aim_temp);
        void apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses, double temperature);
};

class Andersen: public Thermostat {
    private:
        double m_cutoff;
        double m_aim_temp;
    public:
        Andersen(double freq, double timestep, double aim_temp);
        void apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses, double temperature);
};

class Berendsen: public Thermostat {
    private:
        double m_timescale;
        double m_aim_temp;
    public:
        Berendsen(double timescale, double aim_temp);
        void apply(Eigen::ArrayXXd &velocities, const Eigen::ArrayXd &masses, double temperature);
};        