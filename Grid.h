
#include <vector>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
using namespace std; 

class Grid{
    private:
        using index_t = vector<int>::size_type; // Define index_t as the size type of vector<int>
        vector<vector<int>> grid;
        index_t size;
        mt19937_64 generator;
        uniform_real_distribution<double> r_distribution;

        public:
            // Constructor
            Grid(index_t n, unsigned seed) : size(n), generator(seed), r_distribution(0.0, 1.0)
        {
            grid.resize(n, vector<int>(n));
        }

        
        float magnetization(int domain,float init_mg) const
        // Magnetización del grid
        // Domains: completa -->3, abajo-->1, arriba-->2
        // Es const al final para poder acceder a ella dentro de otras funciones.
        {
            int plus = floor((1 + init_mg / 2))*size;
            int l_lim; 
            int u_lim; 
            if (domain==0){l_lim=0;u_lim=plus;}else if(domain==1){l_lim=plus;u_lim=size;}else if(domain==3){l_lim=0;u_lim=size;}
            double sum = 0;
            for (int i = 0 ; i < static_cast<int>(size) ; i++)
            {
                for (int j =static_cast<int>(l_lim); j <  static_cast<int>(u_lim)  ; j++)
                {
                    sum += grid[static_cast<index_t>(i)][static_cast<index_t>(j)];
                }
            }
            return float(sum / pow(size, 2));
        }
        
        struct domains
        {
            float ur;
            float ul;
            float dr;
            float dl;
        };
        
        struct region{
            float up;
            float down;
        };

        float initialize(const string &init_options,float init_mg)
        // Queremos iniciarla con magnetización inicial nula. Varios enfoques:
        // uqr, opposite sign up, pongo el signo contrario en el j+1 cuando inicialize un j;
        // ubd, until bound, magnetizo hasta 0 con un cierto error.
        {
            if (init_options == "ubd")
            {
                for (index_t i = 0; i < size; i++)
                {
                    grid[i][size - 1] = -1;
                    grid[i][0] = 1;
                    for (index_t j = 1; j < size - 1; j++)
                    {
                        double s;
                        s = r_distribution(generator);
                        if (s < 0.5)
                        {
                            grid[i][j] = 1;
                        }
                        else
                        {
                            grid[i][j] = -1;
                        }
                    }
                }
                return magnetization(3,init_mg);
                
            }
            else if (init_options == "uqr")
            {
                for (index_t i = 0; i < size; i++)
                {
                    grid[i][size - 1] = -1;
                    grid[i][0] = 1;
                    if(i%2==0)
                    {
                        for (index_t j = 1; j < size - 1; j++)
                        {
                            if (j % 2 == 0)
                            {
                                grid[i][j] = -1;
                            }
                            else
                            {
                                grid[i][j] = 1;
                            }
                        }
                    }
                    else{
                        for (index_t j = 1; j < size - 1; j++)
                        {
                            if (j % 2 == 0)
                            {
                                grid[i][j] = 1;
                            }
                            else
                            {
                                grid[i][j] = -1;
                            }
                        }
                    }
                }
                return 0.0f;
            }
            else if(init_options=="mgt"){
                //Vamos a poner dominios uniformes
                int plus_bound=floor((1+init_mg/2)*size);
                for (index_t i=0; i<size;i++){
                    for (index_t j=0; j<plus_bound; j++){
                        grid[i][j]=1;
                    }
                    for (index_t j = plus_bound; j < size; j++){
                        grid[i][j]=-1;
                    }
                }
            }
            return 0;
        }

        double energy(){ //Energy of the configuration
            double e=0;
            index_t i,j;
            int left, right;  //Para poder calcular los vecinos izquierdo y derecho teniendo en cuenta CC continuas en x


            j=0;
            for(i=0;i<size;i++){ 
                e+=grid[i][1]+2;
            }

            for(j=1;j<size;j++){i=0;
                left=grid[size-1][j];
                e+=grid[i][j]*(grid[i+1][j]+left+grid[i][j-1]+grid[i][j+1]);    //Ahoramismo se esta saliendo para j
                for(i=1; i<size-1;i++){
                    left=grid[i-1][j]; right=grid[i+1][j];
                    e+=grid[i][j]*(left+right+grid[i][j-1]+grid[i][j+1]);
                }i=size-1;                                                          //Se podría comprobar que aquí i llega hasta size-1.
                right=grid[0][j];
                e+=grid[i][j]*(grid[i][j-1]+grid[i][j+1]+grid[i-1][j]+right);               
            }j=size-1;

            for(i=0;i<size;i++){
                e+=2-grid[i][size-2];
            }
            e=-0.5*e;

            return e;
        }

        int get(int i, int j){
            return grid[static_cast<index_t>(i)][static_cast<index_t>(j)];
        }
        
        void set(int i, int j, int val){
            grid[static_cast<index_t>(i)][static_cast<index_t>(j)]=val;
        } 

        

         
        region magnet_regions(){
            float init_mgt;
            region hemisferio;
            hemisferio.up = magnetization(0, init_mgt);
            hemisferio.down = magnetization(1, init_mgt);
            return hemisferio;
        }

};
