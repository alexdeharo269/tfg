
// Comprobar como se comoportan el analisis cuando se unen los dominios
// para el fihero de centroids, quiero que compruebe cada 100 pasos que el número de centroids no ha cambiado, si cambia hay un
// merge, y ya no podemos trackear el úlitmo paso (number changes) de la misma manera



//Con temperatura puede no llegar nunca al estado estacionario, aunque puedo correr la simulación como mucho hasta el tiempo que
//la simulación a T=0 se congela o el doble, por ejemplo.

//Flowchart de mi allgoritmo. La idea es hacer ver que para dinamicas de equilibrio mientras se respeten las probabiidades no pasa na, 
// el problema es en dinamicas de no equilibrio. Aun así me queda la duda de si tendría que coger i2 entre 1 y N o entre i-R y y+R.
//Antes, cuando hacía la actualización secuencial del anllo, había un junte de dominios mucho mas rapido, que llegamos a hablarlo. 
//No hay referencia a kawasaki con parallel updating en Landau and Binkler book 

//Hay que ir metiendo y comprobando como son las quimeras con reaccion difusión.



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

double pp=0.9;
unsigned int n =1000;
unsigned int Rg=10;
const int t_max = 1000;      // Tiempo
unsigned int R_init=20;
long long unsigned int r_size=1;
unsigned int r_jump=5;
bool Zero_Temperature=false;
float temperature_vals[] = {0.01f, 0.2f, 0.5f, 1.0f, 2.0f, 5.0f, 9.0f, 10.0f}; // TEMPERATURAS ENTERAS DE 1 A 10

double temperature=0.5;
const double avg_mag=0.50;

vector<unsigned int>R_vals(r_size,0);
 
vector<int>init_ring(n,0);



unsigned seed1 = 1944243;   // semilla generador de números aleatorios
mt19937_64 generator(seed1); // generador  mt19937_64
//random_device{}(); para generar numeros realmente aleatorios que paso como semillas

typedef void (*UpdateFunc)(uniform_int_distribution<unsigned int> &, uniform_real_distribution<double> &,mt19937_64 &, vector<int> &, 
    vector<vector<int>> &, int, bool, double);

void initialize();
int campokawa(unsigned int i1, unsigned int i2, vector<int> &ring, vector<vector<int>> &matk);
void initialize_dominio_central(unsigned ancho);
double mdl(vector<int> &ring);
void initialize_shuffle();
void Kawasaki_Step(uniform_int_distribution<unsigned int> &i_distribution , uniform_real_distribution<double> &r_distribution ,
                   mt19937_64 &gen, vector<int> &ring, vector<int> &frozen_ring, vector<vector<int>> &matk,
                   unsigned int R, bool zero_Temperature, double T);

void Glauber_Step(uniform_int_distribution<unsigned int> &i_distribution,
                  uniform_real_distribution<double> &r_distribution,
                  mt19937_64 &gen, vector<int> &ring, vector<int> &frozen_ring,
                  vector<vector<int>> &matg,
                  bool zero_Temperature, double T);

int campoglaub(vector<int> &ring, vector<vector<int>> &matg, unsigned int i1);

void Kawasaki_Step_sequential(uniform_int_distribution<unsigned int> &i_distribution, uniform_real_distribution<double> &r_distribution,
                              mt19937_64 &gen, vector<int> &ring, vector<vector<int>> &matk,
                            bool zero_Temperature, double T);

void Reac_dif_Step_sequential(uniform_int_distribution<unsigned int> &i_distribution, uniform_real_distribution<double> &r_distribution,
                              mt19937_64 &gen, vector<int> &ring, vector<vector<int>> &matk,vector<vector<int>> &matg, double Pp, unsigned N,
                            bool zero_Temperature, double T);

int main()
{
    auto start = std::chrono::system_clock::now();
    

    //initialize_dominio_central(64);
    initialize();
    double sum = 0;
    for (unsigned int i = 0; i < n; i++)
    {
        sum += init_ring[i];
    }
    fprintf(stdout, " \n Magnetizacion inicial: %.4f.\n Size: %i\n", sum/n,n);


    FILE *r_data = fopen("./ising_data/reac_dif_T0/r_data.dat", "w");
    
    fprintf(stdout, "\nR_size= %llu", r_size);
    for (unsigned i = 0; i < r_size; i++)
    {
        R_vals[i] = R_init+i*r_jump;
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
            for (unsigned int j = i + 1; j <= (i + R + Rg); j++)
            {
                matg[i][j % n] = 1 - matk[i][j % n];
                matg[j % n][i] = 1 - matk[j % n][i];
            }
        }
        
        
        /*
        fprintf(stdout, "\n");

        for (unsigned int i=0;i<n;i++){
            for(unsigned int j=0;j<n;j++){
            fprintf(stdout,"%i",matg[i][j]);}
            fprintf(stdout,"\n");
            
        }*/

        double mean_domain_length;
        double length_t=0;
        double length_tm1=0;
        int frozentime=0;

        sprintf(filename, "./ising_data/reac_dif_T0/raw/ising_dataR%iT%.1f.dat", R, temperature);
        
        FILE *flips_data = fopen(filename, "w");
        if (flips_data == NULL)
        {
            perror("Error opening file: "); // Esto se asegura de que si hay algun error en la temperatura pase a la siguiente.
            r_iter++;
        }
        mt19937_64 local_generator(seed1 + r_iter);

        uniform_int_distribution<unsigned int> i_distribution(0, n - 1); // Distribución random para los i
        uniform_real_distribution<double> r_distribution(0., 1.);        // initialize the distribution r_distribution

        vector<int> frozenring(n, 0);
        vector<int> ring = init_ring;

        int t = 0;
        int t_lim=1000;
        
        ring=init_ring;

        for (t = 0; t < t_lim; t++) // va en pasos Monte Carlo de n intentos de intercambio de posición o de spin flip
        {
            //frozenring = ring;

            //for (unsigned int try_ = 0; try_ < n / ((2 * R * pp / n + 2 * Rg * (1 - pp) / n)); try_++)
            //Kawasaki_Step_sequential(i_distribution, r_distribution,
            //                   local_generator, ring, matk, Zero_Temperature, temperature);
            
            Reac_dif_Step_sequential(i_distribution, r_distribution,
                                          local_generator,ring, matk,  matg, pp,  n,
                                          Zero_Temperature, temperature);
            
            //Glauber_Step(i_distribution,r_distribution,local_generator,ring,frozenring,matg,Zero_Temperature,temperature);

            // ring = frozenring;

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
        mean_domain_length = mdl(ring);
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
        }if (sumf != sum){fprintf(stdout, "Conservación magnetización violada R=%u, %f!=%f",r_iter, sum / n, sumf / n);}
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
        if (s < avg_mag)
        {
            init_ring[i] = 1;
        }
        else
        {
            init_ring[i] = -1;
        }
    }
}

int campokawa(unsigned int i1, unsigned int i2, vector<int> &ring, vector<vector<int>> &matk)
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

void Kawasaki_Step(uniform_int_distribution<unsigned int> &i_distribution, uniform_real_distribution<double> &r_distribution,
                   mt19937_64 &gen, vector<int> &ring, vector<int> &frozen_ring, vector<vector<int>> &matk, 
                    bool zero_Temperature, double T)
{
    for (unsigned int try_ = 0; try_ < n; try_++)

    {
        unsigned int i1, i2;
        i1=i_distribution(gen);
        i2=i_distribution(gen);

        double delta = campokawa(i1, i2, ring, matk);
        if(zero_Temperature){
            if((delta<=0)&&(matk[i1][i2]==1)&&(ring[i1]!=ring[i2]))
            {
                if (i1==i2){ fprintf(stdout, "Autoconexion");}
                std::swap(frozen_ring[i1], frozen_ring[i2]);
            }
        }
        else{
            double p =std::min(exp(-delta / T),1.0);
            double ji = r_distribution(gen);
            if ((ji < p) && (matk[i1][i2] == 1) && (ring[i1] != ring[i2]))
            {
                std::swap(frozen_ring[i1], frozen_ring[i2]);
            }
        }
    }
}

void Glauber_Step(uniform_int_distribution<unsigned int> &i_distribution,
                  uniform_real_distribution<double> &r_distribution,
                  mt19937_64 &gen, vector<int> &ring, vector<int> &frozen_ring,
                  vector<vector<int>> &matg, 
                  bool zero_Temperature, double T)

{
    vector<int> hg(n, 0);
    for (unsigned int try_ = 0; try_ < n; try_++)
    {
        unsigned int i1, i2;
        i1 = i_distribution(gen);

        i2=i_distribution(gen);
        double delta = campoglaub(ring, matg, i1);
        if (zero_Temperature)
        {
            if ((delta <= 0) && (matg[i1][i2] == 1))
            {
                ring[i1] = -ring[i1];
            }
        }
        else
        {
            double p = std::min(exp(-delta / T), 1.0);
            double ji = r_distribution(gen);
            if ((ji < p) && (matg[i1][i2] == 1))
            {
                ring[i1] = -ring[i1];
            }
        }
    }
}

void Kawasaki_Step_sequential(uniform_int_distribution<unsigned int> &i_distribution, uniform_real_distribution<double> &r_distribution,
                   mt19937_64 &gen, vector<int> &ring, vector<vector<int>> &matk,
                    bool zero_Temperature, double T)
{
    for (unsigned int try_ = 0; try_ < n; try_++)

    {
        unsigned int i1, i2;
        i1 = i_distribution(gen);
        i2 = i_distribution(gen);


        double delta = campokawa(i1, i2, ring, matk);
        if (zero_Temperature)
        {
            if((delta<=0)&&(matk[i1][i2]==1)&&(ring[i1]!=ring[i2]))
            {
                if (i1 == i2)
                {
                    fprintf(stdout, "Autoconexion");

                }

                std::swap(ring[i1], ring[i2]);
            }
        }
        else
        {
            double p = std::min(exp(-delta / T), 1.0);
            double ji = r_distribution(gen);
            if ((ji < p) && (matk[i1][i2] == 1) && (ring[i1] != ring[i2]))
            {
                std::swap(ring[i1], ring[i2]);
            }
        }
    }
}

void Reac_dif_Step_sequential(uniform_int_distribution<unsigned int> &i_distribution, uniform_real_distribution<double> &r_distribution,
                              mt19937_64 &gen, vector<int> &ring, vector<vector<int>> &matk, vector<vector<int>> &matg, double Pp, unsigned N,
                             bool zero_Temperature, double T)
{
    for (unsigned int try_ = 0; try_ < N; try_++)
    {

        
        unsigned int i1, i2;
        i1 = i_distribution(gen);
        i2 = i_distribution(gen);

        double xr = r_distribution(gen);

        int delta;
        int mov;
        if (xr < Pp)
        {
            delta = campokawa(i1, i2, ring, matk);
            mov = 0;
        }
        else
        { // reac
            
            delta = campoglaub(ring, matg, i1);
            mov = 1;
        }
        double p, ji;
        p = exp(-delta / T);
        if (p >1)
        {
            p = 1;
        }
        //  Cambio de sitio (o no),
        ji = r_distribution(generator);
        if (ji < p)
        {
            if ((mov == 0) && (matk[i1][i2] == 1))
            {
                //fprintf(stdout, "error");//swaps++;
                swap(ring[i1], ring[i2]);
            }
            if ((mov == 1) && (matg[i1][i2] == 1))
            {
                // flips++;
                //fprintf(stdout, "Reaccion");
                ring[i1] = -ring[i1];
            }
        }
    }
}

int campoglaub(vector<int> &ring, vector<vector<int>> &matg, unsigned int i )
{
    int sum = 0;
    for (unsigned int j = 0; j < n; j++)
    {
        sum += matg[i][j] * ring[j];
    }
    return  sum * 2 * ring[i];
}
