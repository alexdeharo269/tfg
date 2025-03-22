#include "Ising_Simulator.h"
#include <cstdio>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <atomic>

//Los ficheros los hace bien pero hay un error en la dinamica. comprobar como esta haciendo los parametros de output. 
//Comparar simulaciones clave con el otro. 

// Constructor: allocate memory for ring and connectivity matrices.
IsingSimulator::IsingSimulator(const SimulationParams &params,
                               unsigned int R, unsigned int Rg,
                               double pp, const int t_max,
                               unsigned int seed)
    :params(params), R(R), Rg(Rg), pp(pp), t_max(t_max), seed(seed)
{
    init_ring.resize(params.n, 0);
    ring.resize(params.n, 0);
    matk.assign(params.n, std::vector<int>(params.n, 0));
    matg.assign(params.n, std::vector<int>(params.n, 0));
}
// Initialize the ring spins based on avg_mag.
void IsingSimulator::initializeRing() {
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (unsigned int i = 0; i < params.n; i++) {
        double s = dist(gen);
        init_ring[i] = (s < params.avg_mag) ? 1 : -1;
    }
    ring = init_ring;
}
// Initialize the ring spins ensuring an equal number of up (+1) and down (-1) spins,
// and then shuffling them randomly.
void IsingSimulator::initializeRingBalanced() {
    init_ring.resize(params.n);
    // Fill first half with +1 and second half with -1.
    std::fill(init_ring.begin(), init_ring.begin() + params.n / 2, 1);
    std::fill(init_ring.begin() + params.n / 2, init_ring.end(), -1);
    // Shuffle the spins randomly.
    std::mt19937_64 gen(seed);
    std::shuffle(init_ring.begin(), init_ring.end(), gen);
    ring = init_ring;
}

// Build the connectivity matrices for matk and matg.
void IsingSimulator::buildConnectivity() {
    for (unsigned int i = 0; i < params.n; i++) {
        // Build matk (neighbors within R).
        for (unsigned int j = i + 1; j <= (i + R ); j++) {
            int idx = j % params.n;
            matk[i][idx] = 1;
            matk[idx][i] = 1;
        }
        // Build matg (neighbors within R+Rg excluding those already in matk).
        for (unsigned int j = i + 1; j <= (i + R + Rg) ; j++) {
            int idx = j % params.n;
            matg[i][idx] = 1-matk[i][idx];
            matg[idx][i] = 1-matk[idx][i];
        }
    }
}

// A single simulation step (here simplified; insert your Kawasaki/Glauber/other logic).
void IsingSimulator::simulationStep(std::mt19937_64 &generator,
                                    std::uniform_int_distribution<unsigned int> &i_dist,
                                    std::uniform_real_distribution<double> &r_dist)
{
    for (unsigned int try_ = 0; try_ < params.n; try_++)
    {

        unsigned int i1, i2;
        i1 = i_dist(generator);
        i2 = i_dist(generator);

        double xr = r_dist(generator);

        int delta;
        int mov;
        if (xr < pp)
        {
            delta = campokawa(i1, i2);
            mov = 0;
        }
        else
        { // reac

            delta = campoglaub(i1);
            mov = 1;
        }
        double p, ji;
        p = exp(-delta / params.temperature);
        if (p > 1)
        {
            p = 1;
        }
        //  Cambio de sitio (o no),
        ji = r_dist(generator);
        if (ji < p)
        {
            if ((mov == 0) && (matk[i1][i2] == 1))
            {
                // fprintf(stdout, "error");//swaps++;
                std::swap(ring[i1], ring[i2]);
            }
            if ((mov == 1) && (matg[i1][i2] == 1))
            {
                // flips++;
                // fprintf(stdout, "Reaccion");
                ring[i1] = -ring[i1];
            }
        }
    }
}

// Calculate the mean domain length (MDL) by counting domain boundaries.
double IsingSimulator::calculateMDL() {
    double boundaries = 0;
    for (unsigned int i = 1; i < params.n; i++) {
        boundaries += std::abs(ring[i] - ring[i - 1]) / 2;
    }
    boundaries += std::abs(ring[params.n - 1] - ring[0]) / 2;
    return (boundaries > 0) ? static_cast<double>(params.n) / boundaries : params.n;
}

// Run the simulation loop, writing the spin configurations to a file.
void IsingSimulator::simulationLoop(const std::string &dataFilename)
{
    int t_lim = 1000, t = 0;
    double length_tm1 = 0, length_t = 0;
    frozentime = 0;
    
    char filename[256];
    sprintf(filename, "%s", dataFilename.c_str());
    FILE *flips_data = fopen(filename, "w");

    if (!flips_data) {
        perror("Error opening simulation data file");
        exit(EXIT_FAILURE);
    }
    // Set up a local random generator.
    std::mt19937_64 local_generator(seed);

    std::uniform_int_distribution<unsigned int> i_dist(0, params.n - 1);
    std::uniform_real_distribution<double> r_dist(0.0, 1.0);
    
    // Main simulation loop.
    for (t = 0; t < t_lim; t++) {
        simulationStep(local_generator, i_dist, r_dist);
        
        // Write current configuration (mapping {-1,1} to {0,1}).
        for (unsigned int i = 0; i < params.n; i++) {
            fprintf(flips_data, "%i ", (ring[i] + 1) / 2);
        }
        fprintf(flips_data, "\n");
        
        // Update MDL and track last time it changed.
        length_t = calculateMDL();
        if (length_t != length_tm1) {
            length_tm1 = length_t;
            frozentime = t;
        }
        if ((t == t_lim - 1) && ((t_lim - frozentime) < 100) && (t < t_max - 1000))
        {
            t_lim += 1000;
        }
    }
    fclose(flips_data);
    meanDomainLength = calculateMDL();
}


// Run the complete simulation: initialize, build connectivity, run the loop, and write outputs.
void IsingSimulator::runSimulation(long long unsigned sim_id) {
    // Initialize spins.
    initializeRing();
    double init_sum = std::accumulate(init_ring.begin(), init_ring.end(), 0.0);
    init_mag = init_sum / params.n;
    

    // Build connectivity matrices.
    buildConnectivity();
    /*
    for (unsigned i = 0; i < params.n; i++)
    {
        for (unsigned j = 0; j < params.n; j++)
        {
            fprintf(stdout, "%i", matk[i][j]);
        }
        fprintf(stdout, "\n");
    }*/
    // Prepare output filenames (here embedding R and temperature).
    char dataFilename[256], metaFilename[256];
    sprintf(dataFilename, "./ising_data/simulations/%llu.dat", sim_id);
    
    
    // Run simulation loop.
    simulationLoop(dataFilename);
    
    // Write out the metadata file.
    //writeMetadata(metaFilename, frozentime, meanDomainLength);
    
    fprintf(stdout, "ID: %llu, T=%.2f, Avg_mag=%.2f, N=%u, R=%u, Rg=%u, In_mag=%.2f, MDL=%.2f, Frozen=%i\n",
            sim_id, params.temperature, params.avg_mag,
            params.n, R, Rg,init_mag, meanDomainLength, frozentime);
}


int main()
{
    auto start = std::chrono::system_clock::now();

    // Define parameter sets.
    std::vector<double> temperatures = {0.00001};
    std::vector<double> avg_mags = {0.6,0.65,0.7,0.75,0.8,0.85};

    // Connectivity parameters.
    unsigned int R_init = 9, r_jump = 5, r_size = 10;
    unsigned int Rg_init = 0, rg_jump = 50, rg_size = 1;

    // Fixed simulation parameters.
    unsigned int n = 200;
    double pp = 1.00; // example reaction-diffusion parameter
    const int t_max=5000;
    unsigned int seed1 = 1944243;

    FILE *masterMeta = fopen("./ising_data/simulation_index.csv", "w");

    // Generate a simulation ID or filename
    std::atomic<unsigned long long> KeyCounter(0);
    // Loop over parameter combinations.
    for (const auto &temp : temperatures)
    {
        for (const auto &mag : avg_mags)
        {
            // Pack simulation parameters.
            SimulationParams params;
            params.temperature = temp;
            params.avg_mag = mag;
            params.n = n;

            // Define connectivity R values.
            std::vector<unsigned int> R_vals(r_size);
            for (unsigned i = 0; i < r_size; i++)
            {
                R_vals[i] = R_init + i * r_jump + 1;
            }

            // Parallelize over connectivity values.
            #pragma omp parallel for schedule(dynamic)
            for (unsigned r_iter = 0; r_iter < R_vals.size(); r_iter++)
            {
                unsigned int R = R_vals[r_iter];
                // For simplicity, choose the first Rg value.
                unsigned int Rg = (rg_size > 0) ? (Rg_init + 0 * rg_jump + 1) : 0;
                // Create an instance of the simulator with its own seed.
                // No esta iterando sobre Rg

                unsigned long long sim_id = KeyCounter.fetch_add(1, std::memory_order_relaxed);

                IsingSimulator simulator(params, R, Rg, pp,t_max, seed1 + r_iter);
                simulator.runSimulation(sim_id);

                // Append a row to the master metadata file (make thread-safe).
                #pragma omp critical
                {
                    fprintf(masterMeta, "%llu,%.2f,%.2f,%u,%u,%u,%.2f,%.2f,%i\n",
                            sim_id, params.temperature, params.avg_mag,
                            params.n, R, Rg, simulator.init_mag, simulator.meanDomainLength, simulator.frozentime);
                }
            }
        }
    }

    auto end = std::chrono::system_clock::now();
    printf("\nTotal execution time: %lli s\n",
           std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
    return EXIT_SUCCESS;
}