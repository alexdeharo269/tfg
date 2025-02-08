#include <vector> 
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
//#include "Grid.h"
using namespace std;
using index = vector<int>::size_type; // Define index_t as the size type of vector<int>

const int n = 20;            // Make n of type index_t
const int t_max = 1000;        // Tiempoñ
const int R = 5;             // Alcance
const double alphadif=1.0;     //Probablidad de difusión
const double alphareac=0.0000;

float temperature_vals[] = {0.000001f,0.2f,0.5f,1.0f}; // TEMPERATURAS ENTERAS DE 1 A 10

double p_val(int i1, int i2, float temperature);
double p_val_2(index i1, index i2, float temperature);
double p_val_reac(int i1, float temperature);
double energy(index i1, index i2);
float magnetization();
void initialize();
//void print(FILE *out);

unsigned seed1 = 19442448; // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64

//Tenemos varias distribuciones de numeros aleatorios
uniform_int_distribution<int> i_distribution(0, n - 1);     //Distribución random para los i

uniform_real_distribution<double> r_distribution(0., 1.);   // initialize the distribution r_distribution

//Constructor
//Ring ring(n, seed1,R);
vector<int> ring(static_cast<index>(n),0);
vector<vector<int>> mat(n,vector<int>(static_cast<index>(n),0));

index s = static_cast<index>(n);
index r = static_cast<index>(R);

int main()
{

    int t_size = sizeof(temperature_vals) / sizeof(temperature_vals[0]); // ¿Cuantas temperaturas hay?
    float temperature;
    char filename[33]; // Restringido a 33 caracteres, entonces solo caben en el nombre temperaturas menores de 10 y con 2 cifra decimal
    
    //Flips-data
    //FILE *magnetization_data;
    //magnetization_data = fopen("./resultados/magnet_domains.dat", "w");

    for (int t_iter = 0; t_iter < t_size; t_iter++)
    {
        // Inicializo el ring
        initialize();
        

        fprintf(stdout, "\n \n \n Magnet 0: %f", magnetization());

        temperature = temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%.3f.dat", temperature); // si no funciona %i probar %d
        fprintf(stdout, "\n Archivo: %s",filename);

        FILE *flips_data = fopen(filename, "w");

        //if (flips_data == NULL)
        //{
        //    perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
        //    t_iter++;
        //}
        for(index i=0;i<n;i++){
            for(index j=0;j<n;j++){
                fprintf(stdout,"%i",mat[i][j]);
            }

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
                //uniform_int_distribution<int> vecino_dif_distribution(i1-R, i1+R);
                //i2 = vecino_dif_distribution(generator);
                i2= i_distribution(generator);


                index j1=static_cast<index>(i1);
                index j2 = static_cast<index>(i2);

                if ((i2>=(i1-R))&&(i2<=(i1+R))){

                    //if(i2<0){i2+=n;}else if (i2>n-1){i2-=n;}

                    p = p_val_2(j1, j2, temperature);
                    //if (p>1){
                    //    p=1.0;}
                    // Cambio de sitio (o no)
                    ji = r_distribution(generator);
                    if (ji < p)
                    {
                        int aux;
                        aux = ring[j1];
                        ring[j1]= ring[j2];
                        ring[j2]= aux;
                    }
                }
                
                else //if ((i2 > (i1 - R)) & (i2 < (i1 + R)))
                {
                    p = p_val_reac(i2, temperature);
                    p=p*alphareac;

                    // Cambio de spin
                    ji = r_distribution(generator);
                    if (ji < p)
                    {

                        ring[i2]=-ring[i2];
                    }
                }
            }

            for (int i = 0; i < n; i++)
            {

                fprintf(flips_data, "%i ", (ring[static_cast<index>(i)] + 1) / 2);
            }
            fprintf(flips_data, "\n");
        }
        fprintf(stdout, "\n \n \n Magnet final: %f", magnetization());
        fclose(flips_data);
    }
    //fclose(magnetization_data);

    return 0;
    
}

/*
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
    p=exp(-alphadif*(Eprime-E)/temperature);//Función de partición Z con constante de Boltzman=1
    if (p >1){p=1;};
    return p;
}*/
//Calculo con el Hamiltoniano de joaquin, tengo que calcular en el hamiltoniano no solo la interacción con los vecinos
//sino la interaccion hasta R, redefiniendo la matriz de conectiviad. Ademas no multiplico por -0.5ino por la diferencia entre 
//espines = 2 o -2.

double p_val_2(index i1, index i2, float temperature)
{
    double p;
    double deltaE;
    //Energría de la configuración
    //deltaE=energy(i1,i2);
    deltaE=1.0;
    //Calculo la probabiliad de cambio de posición como
    p=exp(-alphadif*(deltaE)/temperature);//Función de partición Z con constante de Boltzman=1
    if (p >1){p=1;};
    return p;
}
/*
double p_val_reac(int i1, float temperature){


    double exponential = exp(ring.energy_spin(i1) / temperature);
    double p;
    if (exponential < 1.)
    {
        p = exponential;
    }
    else
    {
        p = 1;
    } // minimo
    return p;
}*/

void initialize(){
    for (index i = 0; i < n; i++)
    {
        double s;
        s = r_distribution(generator);
        if (s < 0.5)
        {
            ring[i] = 1;
        }
        else
        {
            ring[i] = -1;
        }
        for(index j=i+1;j<(i+r);j++){
            
            
            //static_cast<int>(j)=static_cast<int>(j)%n;
            if(j<0){j+=s;}
            mat[i][j]=1;
            mat[j][i]=1;
        }
    }


};



float magnetization()
// Magnetización del ring
// Domains: completa -->3, abajo-->1, arriba-->2
// Es const al final para poder acceder a ella dentro de otras funciones.
{
    double sum = 0;
    for (index i = 0; i < static_cast<index>(n); i++)
    {
        sum += ring[i];
    }
    return float(sum / n);
};

double energy(index i1, index i2){
    double e = 0;
    index j;
    for (j = i1 - r; j < i1 + r; j++)
    {
        if (j != i1)
        {
            e += ring[((j % s) + s) % s];
            //if(j<0){j+=s;}else if(j>s-1){j-=s;}
            //if((j<0)||(j>=s)){fprintf(stdout,"bad");}
            //e+=ring[j];
        }

    }
    for (j = i2 - r; j < i2 + r; j++)
    {
        if (j != i2)
        {
            //if(j<0){j+=s;}else if(j>s-1){j-=s;}
            //if((j<0)||(j>=s)){fprintf(stdout,"bad");}
            //e-=ring[j];
            e -= ring[((j % s)+s)%s];
        }
    }

    return e * (ring[i1] - ring[i2]);
}
