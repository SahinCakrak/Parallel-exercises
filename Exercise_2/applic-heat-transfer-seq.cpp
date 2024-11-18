#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <chrono>


void output_T(int j, std::vector<double> T, int nx, double dt) //outputs temperature and current time of simulation
{
    double t = j*dt;
    printf("the current time of simulation is: %f \n", t );
    std::cout<<"The current temperature is: \n";
    for(int i=0; i<nx;i++)
        std::cout<< T[i]<<" ";

}

void copy_T(int nx, std::vector<double> T, std::vector<double> Tcopy)   //helpful for copying temperature in a temporary temperature 
{
    for(int i=0; i<nx; i++)
        Tcopy[i]=T[i];

}

double get_2nd_der(int i, double dx, std::vector<double> T) //for the second spatial derivative, necessary for temp calculation
{
    double d2Tdx;
    d2Tdx = (T[i+1]-2*T[i]+T[i-1])/(dx*dx);

    return d2Tdx;
}

int main(){
    int s = 1000;
    int nx = 1000;
    int nt = 1000000;  
    double lambda = 2;
    double T0 = 0.7;
    double T1 = 0;
    double T2 = 1;
    double dt;
    double dx;
    double tmax = 100000;
    /*
    std::cout<<"initialize width of wall \n";  //initialize all variables
    std::cin>>s;
    std::cout<<"initialize thermal conductivity \n";
    std::cin>>lambda;
    std::cout<<"initialize simulated time \n";
    std::cin>>tmax;
    std::cout<<"initialize T1\n";
    std::cin>>T1;
    std::cout<<"initialize T2\n";
    std::cin>>T2;
    std::cout<<"initialize the constant initial temperature T0";
    std::cin>>T0;
    std::cout<<"initialize number of points in space nx \n";
    std::cin>>nx;
    std::cout<<"initialize number of time steps nt \n";
    std::cin>>nt;
    */
    std::vector<double> T;
    std::vector<double> Tcopy;
    T.resize(nx,T0);
    Tcopy.resize(T.size());
    T[0] = T1;
    T[nx-1] = T2;
    dx = s/(double) nx; 
    dt = tmax/nt; 
    double cfl = dt/dx;
    double dfl = lambda*dt/(dx*dx);
    std::cout<<"the Courant-Friedrichs-Lewy-number is: "<<cfl<<"\n";
    std::cout<<"the diffusive cfl number is: "<<dfl<<"\n";
    
    output_T(0, T, nx,  dt);
    auto start = std::chrono::high_resolution_clock::now(); //Laufzeit-Messung Beginn
    for(int j=0; j<nt; j++){
        std::copy(T.begin(), T.end(), Tcopy.begin());
            
            for(int i=1; i<nx-1; i++){
                T[i] = Tcopy[i]+ lambda*dt*get_2nd_der( i, dx,  Tcopy);
            }
        
    }
    auto stop = std::chrono::high_resolution_clock::now(); //Ende der Laufzeit-Messung
    output_T( tmax , T, nx,dt);
    std::cout<<" \n";
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout<< "The duration of the program was "<<duration.count()<< " seconds"<<std::endl;
    return 0;
}