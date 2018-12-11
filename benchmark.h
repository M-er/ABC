#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <math.h>

#include "TypeDef.h"

#define abss(a)     (a<0 ? (-a) : a)
//#define max(a , b)  (a>b ?   a  : b)
//#define min(a , b)  (a>b ?   b  : a)

using namespace std;

double mejor_fitness(int func_number);
double evaluar_fitness(Sol &x,  int func_number, int dim);
double Shifted_Sphere( int dim , double* x );
double Schwefel_Problem( int dim , double* x );
double Shifted_Rosenbrock( int dim , double* x );
double Shifted_Rastrigin( int dim , double* x );
double Shifted_Griewank( int dim , double* x );
double Shifted_Ackley( int dim , double* x );
double Shifted_Weierstrass( int dim , double* x );
//funciones para F7 CEC 2008
double Fast_Fractal (int dim, double* x);
double twist(double y);
double fractal1D(double x);
double doubledip(double x, double c, double s);

void inicializacion (Sol &lista, int funcion, int dim);
void reinicia (Sol &lista, int funcion, int dim, int itera_total, int itera);
int max_min (int funcion);
//int dimension (int funcion);
void controlar_limite (Sol &lista, int funcion, int dim);
double tiempo(time_t inicial);
int Random_Entero(int min, int max);
double Random_Real(double min, double max);

#endif
