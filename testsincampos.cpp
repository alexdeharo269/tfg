#include<stdio.h>
#include<vector>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
using namespace std;
//using index = vector<int>::size_type;

//pararelizar
//computo en cloud
//calcular campo, modelasr dinamica glauber.

unsigned int n = 50;
unsigned int R = 10;
unsigned int Rg=0;
const int t_max = 10000; // Tiempo
const double alphadif = 1.0; // Probablidad de difusión
const double alphareac = 0.0000000001;
float temperature_vals[] = {0.01f, 0.2f, 0.5f, 1.0f,2.0f,5.0f,9.0f}; // TEMPERATURAS ENTERAS DE 1 A 10

vector<int> ring(n, 0);                                         
vector<vector<int>> matk(n, vector<int>(n, 0));   //matriz de conectividad
vector<vector<int>> matg(n, vector<int>(n, 0));  // matriz de conectividad


vector<int>init_ring(n,0);
vector<int> frozenring(n,0);
vector<int> hg(n, 0);    // campo dinamica glauber

unsigned seed1 = 1944244;   // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64
// Tenemos varias distribuciones de numeros aleatorios
uniform_int_distribution<unsigned int> i_distribution(0, n - 1); // Distribución random para los i
uniform_real_distribution<double> r_distribution(0., 1.); // initialize the distribution r_distribution

void initialize();
double campokawa(unsigned int i1, unsigned int i2);
void campoglaub();
int count = 0;

int main(){
    
    initialize();
    //Print a la matriz de conectividad

    double sum = 0;
    for (unsigned int i = 0; i < n; i++)
    {
        sum += ring[i];
    }
    fprintf(stdout, " \n Magnetizacion inicial: %.4f\n", sum/n);
    
    /*for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            fprintf(stdout, "%i", matg[i][j]);

        }
        fprintf(stdout,"\n");
    }*/
    
    int t_size = sizeof(temperature_vals) / sizeof(temperature_vals[0]); // ¿Cuantas temperaturas hay?
    float temperature;
    char filename[35]; // Restringido a 35 caracteres, entonces solo caben en el nombre temperaturas menores de 10 y con 2 cifra decimal

    for (int t_iter = 0; t_iter < t_size; t_iter++){
        temperature = temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%.3f.dat", temperature); // si no funciona %i probar %d
        fprintf(stdout, "\n Archivo: %s",filename);
        FILE *flips_data = fopen(filename, "w");
        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            t_iter++;
        }

        //Vuelvo al primer anillo que inicilcé:
        for(unsigned int i=0;i<n;i++){
            ring[i]=init_ring[i];
        }

        unsigned int i1, i2;
        double p;
        double ji;
        int t = 0;
        int d = 0;
        count=0;
        

        
        for (t = 0; t < t_max; t++) // va en pasos Monte Carlo de n intentos de intercambio de posición o de spin flip
        {
            for(unsigned int i=0;i<n;i++){
                frozenring[i]=ring[i];
            }
            campoglaub(); 

            for (unsigned int try_ = 0; try_ < n; try_++)
            {
                i1 = i_distribution(generator);
                i2 = i_distribution(generator);
                /*
                int s;
                if(i1>i2){
                    s=i1;
                    i1=i2;
                    i2=s;
                }*/

                //Para trabajar con índices unsigned teniendo en cuenta la condición circular.
                unsigned low_bound=(i1>=R)? (i1-R):(i1+n-R);
                unsigned up_bound=(i1+R)%n; 
                unsigned int low_bound_reac=(i1>=R+Rg)?(i1-R-Rg):(i1-R-Rg+n);
                unsigned int up_bound_reac=(i1+R+Rg)%n;
                
                //Este es el rango de difusión
                if (((i2 >= low_bound) && (i2 <= up_bound)) ||
                    ((up_bound < low_bound) && ((i2 >= low_bound) || (i2 < up_bound))))
                {

                    
                    p = exp(-alphadif * campokawa(i1,i2) / temperature);
                    if (p >= 1)
                    {
                        //count++;
                        p = 1;
                    }
                    //  Cambio de sitio (o no),
                    ji = r_distribution(generator);
                    if (ji < p)
                    {
                        d = d + 1;
                        swap(frozenring[i1],frozenring[i2]);
                    }
                    //continue;
                }
                else if (!(((i2 >= low_bound) && (i2 <= up_bound)) ||
                           ((up_bound < low_bound) && ((i2 >= low_bound) || (i2 < up_bound)))))
                {
                    if (((i2 >= low_bound_reac) && (i2 <= up_bound_reac)) ||
                        ((up_bound_reac < low_bound_reac) && ((i2 >= low_bound_reac) || (i2 < up_bound_reac))))
                    {
                    // fprintf(stdout,"\n");
                    double pg = exp(-alphareac*(hg[i1]) / temperature);
                    if (pg>1){pg=1;}
                    // Cambio de spin
                    ji = r_distribution(generator);
                    if (ji < pg)
                    {
                        count++;
                        frozenring[i1] = -frozenring[i1];
                    }
                    }
                }
            }
            for(unsigned int i=0;i<n;i++){
                ring[i]=frozenring[i];
            }
            

            for (unsigned int i = 0; i < n; i++)
            {

                fprintf(flips_data, "%i ", (ring[i] + 1) / 2);
            }
            fprintf(flips_data, "\n");
        }
        fclose(flips_data);
        fprintf(stdout,"  %i  %i",d,count);
        
    }
}




void initialize(){
    double s;

    for (unsigned int i = 0; i < n; i++)
    {
        s = r_distribution(generator);
        if (s < 0.5)
        {
            init_ring[i] = 1;
        }
        else
        {
            init_ring[i] = -1;
        }

        for(unsigned int j=i+1;j<=(i+R);j++){
            matk[i][j % n] = 1;
            matk[j%n][i]=1;
        }

        for(unsigned int j=i+1;j<=(i+R+Rg);j++){
            matg[i][j%n]=1-matk[i][j%n];
            matg[j % n][i] = 1 - matk[j%n][i];
        }
    }


}



double campokawa(unsigned int i1,unsigned int i2){
    
    int sum=0;
    for (unsigned int k=0;k<n;k++){
        sum = sum + ring[k] * matk[i1][k] - ring[k] * matk[i2][k];
    }
    sum=sum+matk[i2][i1]*ring[i1]-matk[i1][i2]*ring[i2];
    sum *= (ring[i1] - ring[i2]);
    return sum;
}

void campoglaub(){
    for(unsigned int i=0;i<n;i++){
        int sum = 0;
        for (unsigned int j=0;j<=i;j++){
            sum+=matg[i][j]*ring[j];
        }
        hg[i]=sum*2*ring[i];
    }    
}