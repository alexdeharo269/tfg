// usar la libreria vector para vectores grandes
// grid cerrado (definir función vecino o añadir un marco al grid copiando los valores de los bordes)
// muchos pasos monte carlo L^2.
// g++ programa.cpp -o Pprograma nohup ./Program
//Comprobar el otro inicialization mode. HACER COPIA DEL GRID, INICIALIZAR ANTES DEL BUCLE TEMPERATURA PARA PODER COMPARAR BIEN 
//LAS SIMULACIONES A DISTINTAS T. HACERLO REINICIALIZANDO EL GRID EN LA COPIA. VER COMO COPIAR EL GRID ENTERO.
//Inicalización en modo QR funciona bien, probar modo continuo para el apartado 0;
//En concreto comprobar también que magnet 0 en los fprintf(stdout) corresponden mas o menos con init_magnet
//Esto hace que para este analisis tengamos qe redefinir los dominios de magnetización
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

const int n = 32;            // Make n of type index_t
const int t_max = 20;             // Tiempo
float temperature_vals[] = {0.2f,0.3f,0.4f,0.5f,0.6f,1.0f,5.0f}; // TEMPERATURAS ENTERAS DE 1 A 10
float max_asymetry=0.01f;
float init_mgt=0.5;

double p_val(int i1, int i2, int j1, int j2, float temperature);
void print(FILE *out);

unsigned seed1 = 1946402705; // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64

//Tenemos varias distribuciones de numeros aleatorios
uniform_int_distribution<int> i_distribution(0, n - 1);     //Distribución random para los i
uniform_int_distribution<int> j_distribution(1, n - 2);     //Distribución random para los j (respetar extremos en y fijos)
uniform_real_distribution<double> r_distribution(0., 1.);   // initialize the distribution r_distribution

//Constructor
Grid grid(n, seed1);



int main()
{

    int t_size = sizeof(temperature_vals) / sizeof(temperature_vals[0]); // ¿Cuantas temperaturas hay?
    float temperature;
    char filename[31]; // Restringido a 31 caracteres, entonces solo caben en el nombre temperaturas menores de 10 y con 1 cifra decimal
    

    
    

    /*Apartado 2
    Vamos a calcular la magnetización por dominios para cada temperatura:
    El archivo de salida contiene en cada columna el promedio sobre los pasos montecarlo de la magnetización para cada dominio.
    Distitas filas corresponden a distintas columnas. Realmente es calcular sobre las mitades superior e inferior, pero haciendo 
    4 cuadrantes realmente tengo más información, aunque no debería ser relevante estadísticamente.
    */
    FILE *magnetization_data;
    magnetization_data = fopen("./apartados/magnet_domains.dat", "w");

    /*Apartado 4
    Calculamos la energía por partícula en función de la temperatura. Esto es equivalente a calcular grid.energy()/n^2 al final de cada 
    paso Montecarlo (aka una medida) y promediarlo sobre t_max.
    */
    FILE *epp;
    epp=fopen("./apartados/energy_per_particle.dat","w");

    /*Apartado 6:
    Calor específico a partir de las fluctuaciones de energía. Teniendo en cuenta que <E> es la suma sobre todas las medidas 
    del hamiltoniano y divididas por el nº de medidas (t_max), se ve que <E^2> es la suma de las energías al cuadrado.
    Estos apartados, al ser simulaciones Monte Carlo, creo que en principio no pasa nada si en vez de empeze en tiempo=0
    empiezo en tiempo=100, las diferencias deberían ser mayores. Luego ya es calcular la temperatura critica
    */

    FILE *CV_dat;
    CV_dat=fopen("./apartados/specific_heat.dat","w");
    
    /*Apartado 7:
    Calculo de la susceptibilidad a partir de las fluctuaciones de magnetización en cada dominio. ¿No falta volver a 
    dividir por t_max en la formula? Mismo planteamiento que en el apartado 7
    */

    FILE *Chi_dat;
    Chi_dat = fopen("./apartados/chi.dat", "w");


    /*Apartado 8: 
    Quiero partir de una magnetización no nula
    */

    for (int t_iter = 0; t_iter < t_size; t_iter++)
    {
        // Inicializo el grid
        
        for (int tries=0; tries<10;tries ++){
            float asymetry=grid.initialize("mgt",init_mgt);
            if (abs(asymetry)<max_asymetry){
                break;                         
            }
            else{fprintf(stdout,"\n Intento %i, asimetría: %f",tries+1,asymetry);}
        }


        fprintf(stdout, "\n \n \n Magnet 0: %f", grid.magnetization(3,init_mgt));

        temperature = temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%.1f.dat", temperature); // si no funciona %i probar %d
        fprintf(stdout, "\n Archivo: %s",filename);

        FILE *flips_data = fopen(filename, "w");

        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            t_iter++;
        }
        
        // Inicializo los resultados de la magnetización por dominios
        Grid::region result;
        result.up=0;result.down=0;
        int i1,i2,j1,j2;
        double p;
        double ji;
        double energy; double m_up; double m_down;
        float sum_m_up=0; float sum_m_down=0; 
        double sum_epp=0;
        double sum_eppsrd=0;
        double sum_msqrd_up=0; double sum_msqrd_down=0;
        print(flips_data);
        for (int t = 0; t < t_max; t++) // va en pasos Monte Carlo de n^2 intentos de spin flip
        {
            for (int try_ = 0; try_ < pow(n, 2); try_++)
            {
                //Dos parejas
                i1 = i_distribution(generator);i2 = i_distribution(generator);j1=j_distribution(generator);j2=j_distribution(generator);
                p= p_val(i1,i2,j1,j2, temperature);
                // Cambio de sitio (o no)
                ji = r_distribution(generator);
                if (ji > p)     //No está mal, hay que tener en cuenta que al ser grid global realmente ya he cambiado los siios al calucular Eprime
                {
                    int aux;
                    aux = grid.get(i1,j1);
                    grid.set(i1,j1,grid.get(i2,j2));
                    grid.set(i2,j2,aux);
                }
                
            }
            fprintf(flips_data, "\n");
            print(flips_data);
            
            //Aqui no lo estoy haciendo con el valor absoluto, aunque luego lo pruebo.
            energy=grid.energy();
            sum_epp+=energy/pow(n,2);
            sum_eppsrd+=pow(energy/n,2);

            result = grid.magnet_regions();
            m_up = result.up;
            m_down = result.down;
            sum_m_up += m_up;
            sum_m_down += m_down;
            sum_msqrd_up+=pow(m_up,2); sum_msqrd_down+=pow(m_down,2);
            //fprintf(stdout, "eep: %lf \n", sum_epp);
            
        }
        
        fprintf(magnetization_data,"%.1f, %f, %f \n",
        temperature, sum_m_up/t_max, sum_m_down/t_max);
        fprintf(epp,"%.1f, %lf \n",temperature,sum_epp/t_max );
        fprintf(CV_dat,"%.1f, %lf \n",temperature,(sum_eppsrd-pow(n*sum_epp,2))/pow(t_max,2));
        fprintf(Chi_dat, "%lf, %lf, %lf  \n", temperature,
                (sum_msqrd_up - pow(abs(sum_m_up), 2)) / pow(t_max, 2), // Aquí me la estoy jugando poniendo el abs a la magnetización, formula del apartado 7
                (sum_msqrd_down - pow(abs(sum_m_down), 2)) / pow(t_max, 2));
        fclose(flips_data);
    }
    fclose(magnetization_data);
    fclose(epp);
    fclose(CV_dat);
    fclose(Chi_dat);

    return 0;
}


double p_val(int i1, int i2, int j1, int j2, float temperature)
{
    double p;
    double E,Eprime;
    int aux;
    //Energría de la configuración
    E=grid.energy();
    //Cambio por su vecino aleatorio y vuelvo a calcular la energía
    aux = grid.get(i1, j1);
    grid.set(i1, j1, grid.get(i2, j2));
    grid.set(i2, j2, aux);
    Eprime=grid.energy();
    //Calculo la probabiliad de cambio de posición como
    p=exp(-(Eprime-E)/temperature);
    if (p >1){p=1;};
    return p;
}


void print(FILE *out)   //Solo para imprimir el grid
{
    for (int i = 0; i < n; i++)
    {
        fprintf(out, "%i", grid.get(i,0));
        for (int j = 1; j < n; j++)
        {
            fprintf(out, ",%i  ", grid.get(i,j));
        }
        fprintf(out, "\n");
    }
    fflush(out);
};


