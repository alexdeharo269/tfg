// PENDIENTE:
// Comprobar la inicialización
// La idea es hacer difusión (Kawasaki) hasta R y reacción (Glauber) en el resto
// Condiciones de contorno periodicas
#include <vector>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
#include "Grid.h"
using namespace std;
using index_t = vector<int>::size_type; // Define index_t as the size type of vector<int>

const int n = 1000;            // Make n of type index_t
const int t_max = 500;        // Tiempo
const int R = n/10;             // Alcance

float temperature_vals[] = {0.5f,0.8f,1.2f,1.7f,1.0f,5.0f}; // TEMPERATURAS ENTERAS DE 1 A 10

double p_val(int i1, int i2, float temperature);
void print(FILE *out,int t);

unsigned seed1 = 1946402705; // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64

//Tenemos varias distribuciones de numeros aleatorios
uniform_int_distribution<int> i_distribution(0, n - 1);     //Distribución random para los i
uniform_real_distribution<double> r_distribution(0., 1.);   // initialize the distribution r_distribution

//Constructor
Ring ring(n, seed1);

int main()
{

    int t_size = sizeof(temperature_vals) / sizeof(temperature_vals[0]); // ¿Cuantas temperaturas hay?
    float temperature;
    char filename[31]; // Restringido a 31 caracteres, entonces solo caben en el nombre temperaturas menores de 10 y con 1 cifra decimal
    
    //Flips-data
    FILE *magnetization_data;
    magnetization_data = fopen("./resultados/magnet_domains.dat", "w");

    for (int t_iter = 0; t_iter < t_size; t_iter++)
    {
        // Inicializo el ring
        ring.initialize("brr");


        //fprintf(stdout, "\n \n \n Magnet 0: %f", ring.magnetization());

        temperature = temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%.1f.dat", temperature); // si no funciona %i probar %d
        fprintf(stdout, "\n Archivo: %s",filename);

        FILE *flips_data = fopen(filename, "w");

        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            t_iter++;
        }
        
        // Defino i1 e i2 para la difusión y i para la reacción 
        //int i;
        int i1,i2;
        double p;
        double ji;
        int t=0;
        for ( t = 0; t < t_max; t++) // va en pasos Monte Carlo de n intentos de intercambio de posición o de spin flip
        {
            for (int try_ = 0; try_ <n; try_++)
            {
                //Una pareja de números, uno al azar y otro un vecino al azar pero dentro del alcance.
                i1 = i_distribution(generator);
                uniform_int_distribution<int> vecino_distribution(i1-R, i1+R);
                i2 = vecino_distribution(generator);
                
                if(i2<0){i2+=n;}else if (i2>n-1){i2-=n;}
                
                p= p_val(i1,i2, temperature);
                // Cambio de sitio (o no)
                ji = r_distribution(generator);
                if (ji > p)     //No está mal, hay que tener en cuenta que al ser grid global realmente ya he cambiado los siios al calucular Eprime
                {
                    int aux;
                    aux = ring.get(i1);
                    ring.set(i1,ring.get(i2));
                    ring.set(i2,aux);
                }
                
            }
            print(flips_data,t);
        }
        fclose(flips_data);
    }
    fclose(magnetization_data);
    

    return 0;
}


double p_val(int i1, int i2, float temperature)
{
    double p;
    double E,Eprime;
    int aux;
    //Energría de la configuración
    E=ring.energy();
    //Cambio por su vecino aleatorio y vuelvo a calcular la energía
    aux = ring.get(i1);
    ring.set(i1, ring.get(i2));
    ring.set(i2, aux);
    Eprime=ring.energy();
    //Calculo la probabiliad de cambio de posición como
    p=exp(-(Eprime-E)/temperature);//Función de partición Z con constante de Boltzman=1
    if (p >1){p=1;};
    return p;
}


void print(FILE *out,int t)   //Solo para imprimir el grid
{
    //fprintf(out,"%i ",t);
    for (int i = 0; i < n; i++)
    {
        int bin;
        if (ring.get(i)==-1){bin=0;}else{bin=1;}
        fprintf(out, "%i ",bin);
    }
    fprintf(out, "\n");

    fflush(out);
};


