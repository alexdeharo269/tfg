#ifndef ISING_SIMULATOR_H
#define ISING_SIMULATOR_H

#include <vector>
#include <string>
#include <random>

// Encapsulate simulation parameters.
struct SimulationParams {
    double temperature;      // e.g., 0.01, 0.2, 0.5, …
    double avg_mag;          // e.g., 0.5, 0.66, 0.75, …
    unsigned int n;          // System size
};

class IsingSimulator {
public:
    // Constructor: pass simulation parameters plus connectivity info and a seed.
    IsingSimulator(const SimulationParams &params,
                   unsigned int R, unsigned int Rg,
                   double pp,
                   const int t_max, 
                   unsigned int seed);

    // Runs the full simulation: initialization, simulation loop, and writing output.
    void runSimulation(long long unsigned sim_id);

    // Output simulation parameters --> aquí es donde tengo que meter los demás que salen del analisis de quimeras.
    double init_mag;
    double meanDomainLength;
    int frozentime;

private:
    // Simulation parameters.
    SimulationParams params;
    unsigned int R;    // connectivity parameter
    unsigned int Rg;   // second connectivity parameter
    double pp;         // reaction-diffusion parameter (example)
    const int t_max;   // Maximum simulation steps
    unsigned int seed;

    

    // Simulation state.
    std::vector<int> init_ring;
    std::vector<int> ring;
    std::vector<std::vector<int>> matk; // connectivity matrix for Kawasaki/Glauber
    std::vector<std::vector<int>> matg; // connectivity matrix for reaction-diffusion

    // Helper functions.
    void initializeRing();
    void buildConnectivity();
    void simulationLoop(const std::string &dataFilename);
    void simulationStep(
        std::mt19937_64 &generator,
        std::uniform_int_distribution<unsigned int> &i_dist,
        std::uniform_real_distribution<double> &r_dist
    );
    double calculateMDL();
    void writeMetadata(const std::string &filename);

    int campokawa(unsigned int i1, unsigned int i2)
    {

        int sum = 0;
        for (unsigned int k = 0; k < params.n; k++)
        {
            sum = sum + ring[k] * matk[i1][k] - ring[k] * matk[i2][k];
        }
        sum = sum + matk[i2][i1] * ring[i1] - matk[i1][i2] * ring[i2];
        sum *= (ring[i1] - ring[i2]);
        return sum;
    }
    int campoglaub(unsigned int i)
    {
        int sum = 0;
        for (unsigned int j = 0; j < params.n; j++)
        {
            sum += matg[i][j] * ring[j];
        }
        return sum * 2 * ring[i];
    }
};

#endif // ISING_SIMULATOR_H
