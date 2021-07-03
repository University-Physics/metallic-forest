
/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "Body.h"
#include "animate.h"

std::string filename(int n);

int main(int argc, char *argv[])
{
    Body Molecule[N];
    CRandom rand(42);
    double mu = 0, sigma = 0.1;
    double move_x, move_y;
    Vector3D move;
    double tdibujo = 0;

    double x, y, z, vx, vy, vz, x0 = 2 * Lx, y0 = 2 * Ly;

    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);

    // Initial configuration
    for (int k = 0; k < N; k++) // Run through every molecule
    {
        // Initial positions in cubic lattice
        x = dx + (k % Nx) * dx + x0;
        y = dy + ((k / Nx) % Ny) * dy + y0;
        z = 0;

        // Initial null velocities
        vx = 0;
        vy = 0;
        vz = 0;

        Molecule[k].init(x, y, z, vx, vy, vz, m0);
    }

    start_animation(argc);

    for (int t = 0; t < 10000; t++)
    {
        if (t % 20 == 0)
        {
            begin_frame(argc);
            for (int k = 0; k < N; k++)
                Molecule[k].print();
            end_frame(argc);
        }
        for (int k = 0; k < N; k++)
        {

            move_x = rand.gauss(mu, sigma);
            move_y = rand.gauss(mu, sigma);
            move.load(move_x, move_y, 0);

            Molecule[k].setR(Molecule[k].getR() + move);
        }
    }

    return 0;
}

std::string filename(int n)
{
    std::string name;
    std::stringstream n_s;
    n_s << n;
    name = "./" + n_s.str() + ".csv";
    return name;
}
