// usar la libreria vector para vectores grandes
//grid cerrado (definir función vecino o añadir un marco al grid copiando los valores de los bordes)
// muchos pasos monte carlo L^2.
// g++ programa.cpp -o Pprograma nohup ./Program

#include<vector>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
using namespace std;
using index_t = vector<int>::size_type; // Define index_t as the size type of vector<int>

const index_t n = 20; // Make n of type index_t
const int t_max = 2; // Tiempo 
int temperature_vals[] = {1,2}; //TEMPERATURAS ENTERAS DE 1 A 10
vector<vector<int>> grid(n, vector<int>(n));


void initialize_grid();
double p_val(index_t i, index_t j,int temperature);
void print(FILE *out, FILE *out2,float m,int t);
float m_number(); // magnetización

unsigned seed1 = 1946402705;  //semilla generador de números aleatorios
mt19937_64 generator(seed1);  // generador  mt19937_64

uniform_int_distribution<index_t> i_distribution(0, n-1); // initialize the distribution i_distribution (from which we want to draw random numbers)
uniform_real_distribution<double> r_distribution(0., 1.); // initialize the distribution r_distribution

int main(){
    

    int t_size = sizeof(temperature_vals) / sizeof(temperature_vals[0]); // ¿Cuantas temperaturas hay?
    int temperature;
    char filename[29]; // Restringido a 28 caracteres, entonces solo caben en el nombre temperaturas menores de 10
    
    FILE *magnetization_data;
    magnetization_data = fopen("magnet_data.dat", "w");
    FILE *magnet_T;
    magnet_T=fopen("magnet_vs_t.dat","w");

    for (int t_iter=0; t_iter<t_size;t_iter++){
        temperature=temperature_vals[t_iter];
        sprintf(filename, "./ising_data/ising_dataT%i.dat",temperature); // si no funciona %i probar %d

        FILE *flips_data = fopen(filename, "a");

        if (flips_data == NULL)
        {
            perror("Error opening file: ");  //Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            t_iter++;
        }

        initialize_grid();
        index_t i;
        index_t j;
        double p;
        double ji;
        float m;
        float sum_m=0;
        for (int t = 0; t < t_max; t++) // va en pasos Monte Carlo de n^2 intentos de spin flip
        {
            for (int try_ = 0; try_ < pow(n, 2); try_++)
            {
                //while (i>n-1 || j>n-1)
                
                    i = i_distribution(generator); // escojo dos puntos al azar.
                    j = i_distribution(generator);
                
                
                p = p_val(i, j,temperature);
                
                // Cambio de spin
                ji = r_distribution(generator);
                if (ji < p)
                {
                    grid[i][j] = -grid[i][j];
                }
            }
            fprintf(flips_data, "\n");
            m = m_number();
            sum_m+=abs(m);         ////////////////////////////////
            if (t==0){
                fprintf(magnetization_data,"%.4f",m); //Formatear la salida de magnetizacion
            }
            print(flips_data, magnetization_data, m,t); 
        }
        fprintf(magnetization_data, "\n");
        sum_m=sum_m/t_max;         ////////////////////////////////   
        fprintf(magnet_T,"%i %.5f \n",temperature,sum_m);
        fclose(flips_data);        
    }
    fclose(magnetization_data);
    fclose(magnet_T);

    return 0;
}

void initialize_grid()
{
    for(index_t i=0;i<n;i++){   
        for(index_t j=0;j<n;j++){
            double s;

            s=r_distribution(generator);
            if(s<0.5){
                grid[i][j] = 1;
            }
            else
            {
                grid[i][j] =-1;
            }
        }

    }
}

double p_val(index_t i, index_t j,int temperature){
    double p;
    double deltaE; double exponential;
    int left, down, up, right;
    

    
    if (i == 0){left = grid[n - 1][j];}else{left = grid[i - 1][j];}
    if (i == n - 1){right = grid[0][j];}else{right = grid[i + 1][j];}
    if (j == 0){down = grid[i][n - 1];}else{down = grid[i][j - 1];}
    if (j == n - 1){up = grid[i][0];}else{up = grid[i][j + 1];}
    
    deltaE=2*grid[i][j]*(int(left)+int(right)+int(up)+int(down));
    deltaE=0.5;
    exponential=exp(-deltaE/temperature);

    if (exponential<1.){
        p=exponential;
        }
    else{p=1;}  //minimo
    return p;
}

float m_number(){
    index_t i; index_t j;
    double sum=0; 
    double sumj;
    for(i=0; i<n;i++){
        sumj=0;
        for(j=0;j<n;j++){
            sumj+=grid[i][j];
        }
        sum+=sumj;
    }
    return float(sum/pow(n,2));
}

void print(FILE* out,FILE *out2,float m,int t){
    for(index_t i=0;i<n;i++){
        fprintf(out, "%i", grid[i][0]);
        for (index_t j=1; j < n; j++)
        {
            fprintf(out, ",%i  ",grid[i][j]);
        }
        fprintf(out,"\n");
    }
    fflush(out);
    if(t!=0)
    {
        fprintf(out2, ",%.5f", m);} // Para formatear el archivo de salida.
    }
