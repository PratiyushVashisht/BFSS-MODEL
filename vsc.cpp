/******************************************************************************
Metropolis Algorithm for Gaussian Model <X> , <X^2> and their respective Monte Carlo errors using Metropolis sampling.
*******************************************************************************/

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <random>
//#include <matplotlib.h>
//namespace plt = matplotlibcpp;
using namespace std;

int main()
{
    cout.precision(6);
    cout.setf(ios::fixed|ios::showpoint);
    
    int i,n,N;
    double EPS;     // Machine epsilon
    double u,x,dx,x_new;
    double I ,stder;
    
    
    static ofstream f_data;
    static ofstream n_data;
    f_data.precision(4);
    n_data.precision(2);
    
    double x_val = 0.0 , x_sq_val = 0.0;
    double x_val_e =0.0 , x_sq_val_e = 0.0;
    double avg_x_val =0.0 , avg_x_sq_val =0.0;
    double std_err_x_val =0.0,std_err_x_sq_val = 0.0;
    
    //initialize x value
    x=0.0;
    //simulation parameters
    N=100;       //No. of samples
    EPS= 0.75; //Metropolis step size
    
     f_data.open("data1.txt");
        if (f_data.bad())
        {cout <<"Failed to open data file\n"<<flush;
        }
    n_data.open("data2.txt");
        if (n_data.bad())
        {cout <<"Failed to open data file\n"<<flush;
        }
        
    //loop
    for (i=0;i<N;i++)
    {
        
        dx = (double)rand()/RAND_MAX-0.5;
        cout<< dx<<endl; //random jump in x
        x_new = x +EPS*dx; // Proposed value of x_new
        
                // acceptance Or Rejection Metropolis update with weight exp(−x ^2)
        u= (double)rand()/RAND_MAX;
        
    
        if (u < exp(-(x_new*x_new - x*x)))
        x= x_new;
        
        // for <X> = 
        x_val= x_val+x;
        x_val_e = x_val_e +x*x;
        
        x_sq_val = x_sq_val +x*x;
        x_sq_val_e = x_sq_val_e +(x*x)*(x*x);
        
        f_data<< x_val<<"\t"<<x_val_e<<"\t"<<x_sq_val<<"\t"<<x_sq_val_e<<endl;
        n_data<< x_val/i<<"\t"<<x_sq_val /i<<"\t"<<i<<endl;
    }
    
    avg_x_val = x_val/N ;
    avg_x_sq_val = x_sq_val /N ;
    
    
    
    //standard error 
    std_err_x_val =sqrt((x_val_e/N - pow(avg_x_val,2) )/N);
    std_err_x_sq_val =sqrt((x_sq_val_e/N - pow(avg_x_sq_val,2) )/N);
    
    cout<<"<x>:"<<avg_x_val<<"\t"
    <<std_err_x_val<<endl;
    cout<<"<x^2>:"<<avg_x_sq_val<<"\t"
    <<std_err_x_sq_val<<endl;

  
    
    
    return 0;
    
}