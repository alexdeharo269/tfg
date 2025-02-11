// En la dinamica glauber hay problemas de contorno, ver C:\Users\Ale\Desktop\UGR FISICA\TFG\TFG_things\ising_data\glaub_T0\raw\ising_dataR12.dat
//Comprobar por qué están y si están también afectando esta dinámica. De todas formas la dinamica glauber pura
//no parece muy interesante para este tipo de analisis, quizás probando con un dominio central y con magnetización 0?

// Comprobar como se comoportan el analisis cuando se unen los dominios
// Comprobar con 4 dominios
// Analizar y entender el dbscan
// para el fihero de centroids, quiero que compruebe cada 100 pasos que el número de centroids no ha cambiado, si cambia hay un
// merge, y ya no podemos trackear el úlitmo paso (number changes) de la misma manera


/*
-Quimera normal, estática
-4 quimeras que se mueve
-Quimeras que se juntan
-Aparecerán más casos con la dinamica de glauber y reac dif, pero con esos ya le puedo presentar a joaquin y hacer figuras buenas.
*/

/*
Hacer un initialize_shuffle
*/



#include<stdio.h>
#include<vector>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
#include <omp.h>
#include<numeric>
#include<cstdlib>
using namespace std;

double pp = 1;

unsigned int n =256;
unsigned int Rg=0;
const int t_max = 100000;// Tiempo
unsigned int R_init=10;
long long unsigned int r_size=10;
float temperature_vals[] = {0.01f, 0.2f, 0.5f, 1.0f, 2.0f, 5.0f, 9.0f, 10.0f}; // TEMPERATURAS ENTERAS DE 1 A 10
double temperature=2.00;

vector<unsigned int>R_vals(r_size,0);
 
vector<int>init_ring(n,0);



unsigned seed1 = 1944243;   // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64
//random_device{}(); para generar numeros realmente aleatorios que paso como semillas
//Como están los threads es lo más correcto porque cada simulación llega exactamente a la misma secuencia del 
//generador de numeros pseudoaleatorios.


void initialize();
double campokawa(unsigned int i1, unsigned int i2, vector<int> &ring, vector<vector<int>> &matk);
void initialize_dominio_central(unsigned ancho);
double mdl(vector<int> &ring);

int main()
{
    auto start = std::chrono::system_clock::now();
    

    //initialize_dominio_central(64);
    initialize();
    //Print a la matriz de conectividad

    double sum = 0;
    for (unsigned int i = 0; i < n; i++)
    {
        sum += init_ring[i];
    }
    fprintf(stdout, " \n Magnetizacion inicial: %.4f.\n Size: %i", sum/n,n);
    /*
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            fprintf(stdout, "%i", matk[i][j]);

        }
        fprintf(stdout,"\n");
    }*/

    FILE *r_data = fopen("./ising_data/kawa_T0/r_data.dat", "w");
    
    fprintf(stdout, "\nR_size= %llu", r_size);
    for (unsigned i = 0; i < r_size; i++)
    {
        R_vals[i] = R_init+i + 1;
    }

    char filename[256]; // Restringido a 42 caracteres, entonces solo caben en el nombre de R hasta 100.
    unsigned int R;

    #pragma omp parallel for shared(R_vals ,init_ring, seed1) private(filename, R)
    for (unsigned r_iter = 0; r_iter < r_size; r_iter++){
        R = R_vals[r_iter];
        vector<vector<int>> matk(n, vector<int>(n, 0)); // matriz de conectividad
        vector<vector<int>> matg(n, vector<int>(n, 0)); // matriz de conectividad
        for(unsigned int i=0;i<n;i++){
            for (unsigned int j = i + 1; j <= (i + R); j++)
            {
                matk[i][j % n] = 1;
                matk[j % n][i] = 1;
            }
        }

        double mean_domain_length;
        double length_t=0;
        double length_tm1=0;
        int frozentime=0;

        sprintf(filename, "./ising_data/kawa_T0/raw/ising_dataR%iT%1f.dat", R, temperature);
        
        FILE *flips_data = fopen(filename, "w");
        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            r_iter++;
        }


        //Thread-safe rnadonm number generator
      //mt19937_64 local_generator(seed1 + static_cast<unsigned>(omp_get_thread_num())+ r_iter);
        mt19937_64 local_generator(seed1 + r_iter);

        uniform_int_distribution<unsigned int> i_distribution(0, n - 1); // Distribución random para los i
        uniform_real_distribution<double> r_distribution(0., 1.);        // initialize the distribution r_distribution

        vector<int> frozenring(n, 0);
        vector<int> ring = init_ring;

        //Vuelvo al primer anillo que inicilcé:
        ring=init_ring;

        //fprintf(stdout, "\n%f", pp);

        unsigned int i1, i2;
        //double p;
        //double ji;
        int t = 0;
        //int d = 0;
        //int count=0;
        unsigned int low_bound;
        unsigned int up_bound;
        //unsigned int low_bound_reac;
        //unsigned int up_bound_reac;
        bool inR=false;
        //bool inRg=false;
        //int t_lim=t_max;
        int t_lim=1000;
        

        for (t = 0; t < t_lim; t++) // va en pasos Monte Carlo de n intentos de intercambio de posición o de spin flip
        {
            frozenring = ring;
            
            /*
            if(t==t_max/2){
                //fprintf(stdout, "\n");
                for (int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        if(hk[i][j]!=hk[j][i]){
                            //fprintf(stdout, "%i", hk[i][j]);
                        }
                        
                    }
                    //fprintf(stdout,"\n");
                }
            }*/

            for (unsigned int try_ = 0; try_ < n / ((2 * R * pp / n + 2 * Rg * (1 - pp) / n)); try_++)
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

                // Para trabajar con índices unsigned teniendo en cuenta la condición circular.
                low_bound = (i1 >= R) ? (i1 - R) : (i1 + n - R);
                up_bound = (i1 + R) % n;
                //low_bound_reac = (i1 >= R + Rg) ? (i1 - R - Rg) : (i1 - R - Rg + n);
                //up_bound_reac = (i1 + R + Rg) % n;
                // Calculo a partir de i1 e i2 si está en rango R o Rg
                inR = (((i2 >= low_bound) && (i2 <= up_bound)) ||
                       ((up_bound < low_bound) && ((i2 >= low_bound) || (i2 < up_bound))));
               

                double delta=campokawa(i1,i2,ring,matk);
                double p =std::min(exp(-delta / temperature),1.0);
                
                double ji = r_distribution(generator);
                if ((ji<p)&&(inR)){
                    std::swap(frozenring[i1], frozenring[i2]);
                }   

                /*
                int mov;
                if(xr<pp){
                    //dif
                    //fprintf(stdout, "%f", xr);
                    delta=alphadif*campokawa(i1,i2,ring,matk); mov=0;
                }
                else{//reac
                    delta=alphareac*hg[i1]; mov=1;
                }*/
                
                //En el limite de temperatura 0 sin reacción:
                //if(delta<=0){swap(frozenring[i1], frozenring[i2]);}
                /*
                p = exp(-delta / temperature);
                
                if (delta <= 0)
                {
                    // count++;
                    p = 1;
                }
                //  Cambio de sitio (o no),
                ji = r_distribution(generator);
                if (ji < p)
                {
                    if((mov==0)&&(inR)){
                        d = d + 1;
                        swap(frozenring[i1], frozenring[i2]);
                    }
                    if((mov==1)&&(inRg)&&(!inR)){
                        count++;
                        frozenring[i1] = -frozenring[i1];
                    }
                }*/
            }

            ring = frozenring;
            

            for (unsigned int i = 0; i < n; i++)
            {

                fprintf(flips_data, "%i ", (ring[i] + 1) / 2);
            }
            fprintf(flips_data, "\n");

            length_t = mdl(ring);
            if(length_t!=length_tm1){length_tm1=length_t;frozentime=t;}
            
            if((t==t_lim-1)&((t_lim-frozentime)<100)&(t<t_max-1000)){
                t_lim+=1000;
            }

            
        }
        fclose(flips_data);
        //fprintf(stdout,"\nT=%.1f  ex=%i  sf=%i  thread=%i",temperature,d,count,omp_get_thread_num());
        mean_domain_length = mdl(ring);

        //fprintf(stdout, "\nR=%i  MDL=%.2f, thread=%i,  Tr=%i", R, mean_domain_length,
        //        omp_get_thread_num(), frozentime);
        fprintf(stdout, "\nR=%i  MDL=%.2f, Frozen time=%i", R, mean_domain_length,
                frozentime);
        #pragma omp critical
        {
            fprintf(r_data, "%i %.2f %i\n", R, mean_domain_length, frozentime);
        }

        double sumf = 0;
        for (unsigned int i = 0; i < n; i++)
        {
            sumf += ring[i];
        }
        if (sumf != sum)
        {
            fprintf(stdout, "Conservación magnetización violada R=%u, %f!=%f",r_iter, sum / n, sumf / n);
        }
    }
    fclose(r_data);
    auto end = std::chrono::system_clock::now();
    fprintf(stdout, "\nExec time: %lli s", chrono::duration_cast<std::chrono::seconds>(end - start).count());
}

void initialize(){
    double s;
    uniform_real_distribution<double> r_distribution(0., 1.);        // initialize the distribution r_distribution

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
    }
}

double campokawa(unsigned int i1, unsigned int i2, vector<int> &ring, vector<vector<int>> &matk)
{

    int sum = 0;
    for (unsigned int k = 0; k < n; k++)
    {
        sum = sum + ring[k] * matk[i1][k] - ring[k] * matk[i2][k];
    }
    sum = sum + matk[i2][i1] * ring[i1] - matk[i1][i2] * ring[i2];
    sum *= (ring[i1] - ring[i2]);
    return sum;
}



double mdl(vector<int> &ring){
    //Que cuente a partir del primer cambio para que no ponga dominios extra por las condicioens de contorno
    double k=0;
    for(unsigned j=1; j<n; j++){
        k+=abs(ring[j]-ring[j-1])/2;
    }
    k+=abs(ring[n-1]-ring[0])/2;
    
    return n/k;
}

// Vamos a añadir la función de inicializacion con un dominio central gordo y el resto aleatorio
void initialize_dominio_central(unsigned ancho)
{
    double s;
    uniform_real_distribution<double> r_distribution(0., 1.); // initialize the distribution r_distribution

    for (unsigned i = n / 2 - ancho / 2; i < n / 2 + ancho / 2; i++)
    {
        init_ring[i] = 1;
    }

    for (unsigned int i = 0; i < n / 2 - ancho / 2; i++)
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
    }
    for (unsigned int i = n / 2 + ancho / 2; i < n; i++)
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
    }
}


