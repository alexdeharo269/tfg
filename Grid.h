
#include <vector>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
using namespace std; 

class Ring{
    private:
        using index_t = vector<int>::size_type; // Define index_t as the size type of vector<int>
        vector<int> ring;
        index_t size;
        mt19937_64 generator;
        uniform_real_distribution<double> r_distribution;

    public:
            // Constructor
            Ring(index_t n, unsigned seed) : size(n), generator(seed), r_distribution(0.0, 1.0)
        {
            ring.resize(n);
        }

        
        float magnetization() const
        // Magnetización del ring
        // Domains: completa -->3, abajo-->1, arriba-->2
        // Es const al final para poder acceder a ella dentro de otras funciones.
        {
            double sum = 0;
            for (int i = 0 ; i <static_cast<int>(size) ; i++)
            {
                sum += ring[static_cast<index_t>(i)];
            }
            return float(sum/static_cast<int>(size));
        }

        float initialize(const string &init_options)
        // Queremos iniciarla con magnetización inicial nula. Varios enfoques:
        // uqr, opposite sign up, pongo el signo contrario en el j+1 cuando inicialize un j;
        // ubd, until bound, magnetizo hasta 0 con un cierto error.
        // brr, barrera, todos los ups metidos en una barrera cuadadrada centrada. El único que ahoramismo está adaptado a 1D
        {
            if (init_options == "brr")
            {
                //
                index_t plusrange=static_cast<index_t>(floor(static_cast<int>(size/10)));  //NOT CHECKED

                for (index_t i = 0; i < size; i++)
                {
                    if ((i>size/2-plusrange) & (i<size/2+plusrange)){ring[i]=1;}else{ring[i]=-1;}
                }
            }        
            return 0;
        }

        double energy(){ //Energy of the configuration
            double e =0; //E+=-(1/2)*spin(i)*(spin(i+1)+spin(i-1))
            index_t i;
            int left, right;  //Para poder calcular los vecinos izquierdo y derecho teniendo en cuenta CC continuas en x

            for(i=0;i<size;i++){
                left=ring[i-1]; right=ring[i+1];
                if(i==0){left=ring[size];}
                if(i==size){right=ring[0];}     
                e+=ring[i]*(left+right);
            }
            e=-0.5*e;
            return e;
        }

        int get(int i){
            return ring[static_cast<index_t>(i)];
        }
        
        void set(int i, int val){
            ring[static_cast<index_t>(i)]=val;
        }
};
