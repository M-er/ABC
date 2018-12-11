
#include "benchmark.h"
#include "data.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

static double pi;
static double e;

double numEval;

double mejor_fitness(int func_number){
	return f_bias[func_number];
}
double evaluar_fitness(Sol &x,  int func_number, int dim){
//double func(int func_number , int dim , double* x){
	pi = acos(-1.0);
	e = exp(1.0);

    if(dim>1000) dim = 1000;

	numEval++;
    switch(func_number){
        case 1 : return Shifted_Sphere(dim , x );
            break;
        case 2 : return Schwefel_Problem(dim , x );
            break;
        case 3 : return Shifted_Rosenbrock(dim , x );
            break;
        case 4 : return Shifted_Rastrigin(dim , x );
            break;
        case 5 : return Shifted_Griewank(dim , x );
            break;
        case 6 : return Shifted_Ackley(dim , x );
            break;
		case 7 : //return Fast_Fractal(dim , x );
				for (int i = 0; i < dim; i++)
					x[i]=0.12323;
				cout<<Fast_Fractal(dim , x )<<endl;
				exit(0);
			break;
        default :
            printf("Error : Function number out of range\n");
            exit(0);
            break;
    }
}

double Shifted_Sphere( int dim , double* x ){
    int i;
    double z;
    double F = 0;
    for(i=0;i<dim;i++){
        z = x[i] - sphere[i];
        F += z*z;
    }
		
    return F + f_bias[0];
}

double Schwefel_Problem( int dim , double* x ){
    int i;
    double z;
    double F = abss(x[0] - schwefel[0]);
    for(i=1;i<dim;i++){
    	  z = x[i] - schwefel[i];
        F = (F>abss(z) ?   F  : abss(z)); 
		//F=max(F , abss(z));
    }
    return F + f_bias[1]; 
}

double Shifted_Rosenbrock( int dim , double* x ){
    int i;
    double z[dim];
    double F = 0.0;

    for(i=0;i<dim;i++) z[i] = x[i] - rosenbrock[i] + 1;   

    for(i=0;i<dim-1;i++){    
        F = F + 100*double( pow(double(pow(z[i],2)-z[i+1]) , 2) ) + double(pow((z[i]-1) , 2));
    }
    return F + f_bias[2]; 
}


double Shifted_Rastrigin( int dim , double* x )
{
    int i;
    double z;
    double F = 0;
    for(i=0;i<dim;i++){  
        z = x[i] - rastrigin[i];
        F = F + ( pow(z,2) - 10*cos(2*pi*z) + 10);
    }
    return F + f_bias[3]; 
}

double Shifted_Griewank( int dim , double* x ){
    int i;
    double z;
    double F1 = 0;
    double F2 = 1;
    for(i=0;i<dim;i++){       
        z = x[i] - griewank[i];
        F1 = F1 + ( pow(z,2) / 4000 );
        F2 = F2 * ( cos(z/sqrt(i+1)));

    }
    return (F1 - F2 + 1 + f_bias[4]); 
}

double Shifted_Ackley( int dim , double* x ){
    int i;
    double z;
    double Sum1 = 0;
    double Sum2 = 0;
    double F = 0;
    for(i=0;i<dim;i++){   
        z = x[i]- ackley[i];
        Sum1 = Sum1 + pow(z , 2 );
        Sum2 = Sum2 + cos(2*pi*z);
    }
    F = (-20)*exp((-0.2)*sqrt(Sum1/dim)) -exp(Sum2/dim) + 20 + e + f_bias[5];

    return F; 
}

double Fast_Fractal(int dim, double* x){
	cout<<"Fast_Fractal\n";
	double F=0;
	for(int i=0;i<dim;i++){ 
		int index = int(i%dim + 1);
		double auxFrac = x[i]+twist(x[index]);
		F = F+fractal1D(auxFrac);
	}
	return F;

}

double twist(double y){
cout<<"twist\n";
	double T = 0;
	T = 4 * (y*y*y*y - 2*y*y*y + y*y);
	return T;
}

double fractal1D(double x){
cout<<"fractal1D\n";
	double F1D = 0;
	double limJ = 0;
	int ran2 = 0;
	for(int k=1; k<=3; k++){
		limJ = pow(2,(k-1));
		for(int j=1; j<=limJ; j++){
			ran2 = Random_Entero(0,2);
			for(int i=1; i<=ran2; i++){
				double ran1 = Random_Real(0,1);
				double s = (1 / (limJ * (2 - ran1)));
				F1D = F1D + doubledip(x,ran1,s);
			}
		}
	}
}

double doubledip(double x, double c, double s){
cout<<"doubledip\n";
	double dd = 0;
	double p1,p2,p3,termino;
	if((0.5<x) && (x<0.5)){
		termino = (x-c);
		p1 = (-6144) * pow(termino,6);
		p2 = 3088 * pow(termino,4);
		p3 = -392 * pow(termino,2);
		dd = (p1+p2+p3+1) * s;
		return dd;
	}
	else
	return 0.0;
}

/*----------------------Funcion de Inicializacion---------------------------*/
void inicializacion (Sol &lista, int funcion, int dimf)
{
	switch ( funcion )
	{
		/*----------F1 CEC2008-----------*/
		case 1:  
		{	
			
			
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-100,100);
			}
		}	
		break;
        
		/*----------F2 CEC2008-----------*/
		case 2:
		{	
					
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-100,100);
			}
		}
		break; 
		
		/*----------F3 CEC2008-----------*/
		case 3:
		{	
			
			
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-100,100);
			}
		}
		break; 
		
		/*----------F4 CEC2008-----------*/
		case 4:
		{	
						
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-5,5);
			}
		}
		break; 
		
		/*----------F5 CEC2008-----------*/
		case 5:
		{	
			
			
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-600,600);
			}
		}
		break; 
		
		/*----------F6 CEC2008-----------*/
		case 6:
		{	
			
			
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-32,32);
			}
		}
		break; 
		/*----------F7 CEC2008-----------*/
		case 7:
		{	
			
			
			for(int i=0; i< dimf; i++){
				lista[i] = Random_Real(-1,1);
			}
		}
		break; 
		 default:
			printf("Error");
		break;
    }		

}


/*-----------------------------Funcion controlar limite-------------------------*/
void controlar_limite (Sol &lista, int funcion, int dimf){
	
	switch ( funcion )
	{
        /*----------F1 CEC2008-----------*/
		case 1:  
        { 		
					
				for(int i=0; i< dimf; i++)
				{
					if(lista[i]<-100.0)
					{
						lista[i] = -100.0;	
					}else{
						if(lista[i]>100.0){
							lista[i] = 100.0;
						}
					}
				}
		}
		break;
        
		/*----------F2 CEC2008-----------*/
		case 2:
		{
			
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-100.0)
				{
					lista[i] = -100.0;	
				}else{
					if(lista[i]>100.0){
						lista[i] = 100.0;
					}
				}
			}
		}
		break; 
		
		/*----------F3 CEC2008-----------*/
		case 3:
		{
			
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-100.0)
				{
					lista[i] = -100.0;	
				}else{
					if(lista[i]>100.0){
						lista[i] = 100.0;
					}
				}
			}
		}
		break;
		
		/*----------F4 CEC2008-----------*/
		case 4:
		{
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-5.0)
				{
					lista[i] = -5.0;	
				}else{
					if(lista[i]>5.0){
						lista[i] = 5.0;
					}
				}
			}
		}
		break;
		
		/*----------F5 CEC2008-----------*/
		case 5:
		{
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-600.0)
				{
					lista[i] = -600.0;	
				}else{
					if(lista[i]>600.0){
						lista[i] = 600.0;
					}
				}
			}
		}
		break;
		
		/*----------F6 CEC2008-----------*/
		case 6:
		{
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-32.0)
				{
					lista[i] = -32.0;	
				}else{
					if(lista[i]>32.0){
						lista[i] = 32.0;
					}
				}
			}
		}
		break;
		/*----------F7 CEC2008-----------*/
		case 7:
		{
			
			for(int i=0; i< dimf; i++)
			{
				if(lista[i]<-1.0)
				{
					lista[i] = -1.0;	
				}else{
					if(lista[i]>1.0){
						lista[i] = 1.0;
					}
				}
			}
		}
		break;
		
		default:
			printf("Error");
		break;
    }		
	
}

/*------------------------Funcion de maximizacion o minimazacion--------------------------*/
// min = 0 o max = 1 //
int max_min (int funcion)
{
	switch ( funcion )
	{
		case 1:  
         	return 0;
		break;
        case 2:
			return 0;
		break;
		case 3:
			return 0;
		break;
		case 4:
			return 0;	
        break;
		case 5:
			return 0;	
        break;
		case 6:
			return 0;	
        break;
		case 7:
			return 0;	
        break;
		default:
			printf("Error");
		break;
    }	

}

/*--------------------------Funcion para calcular el tiempo---------------------------*/
double tiempo(time_t inicial){
    time_t final;
    time(&final);
    return difftime(final, inicial);

}

/*----------------------------Funcion random entero------------------------------------*/
int Random_Entero(int min, int max)
{
	double r = (double)rand()/(double(RAND_MAX)+1);
	
	r = r*(max-min) + min;

	return (int) r;

}

/*-------------------------------Funcion random real------------------------------------*/
double Random_Real(double min, double max)
{
	double r ;//= (double)rand()/(double(RAND_MAX)+1);
	
	r = drand48()*(max-min) + min;

	return r;

}


void reinicia (Sol &lista, int funcion, int dimf,int itera_total, int itera)
{
	int tau,b=2;
	double fit_ini = evaluar_fitness (lista, funcion,dimf),fitness;
	cout<<"*********************************************entra reinicia"<<endl;
	for(int i=0; i< dimf; i++){
		if(drand48()<0.2){
		tau = (drand48()>0.5 ?   0  : 1); //-1:1
		//lista[i] = lista[i] + tau*(100-(-100))*(1-pow(drand48(),pow(1-(double)itera/(double)itera_total,b)));
		if(tau==0)
			lista[i] = lista[i] + Random_Entero(0,100-lista[i])*(1-pow(drand48(),pow(1-(double)itera/(double)itera_total,b)));
		else
			lista[i] = lista[i] - Random_Entero(0,lista[i]-(-100))*(1-pow(drand48(),pow(1-(double)itera/(double)itera_total,b)));
		cout<<"----------------- "<<lista[i]<<"-----------------"<<endl;
		fitness = evaluar_fitness (lista, funcion,dimf);
		if(fit_ini > fitness){
			cout<<fit_ini <<" "<< fitness <<"--------------------------------------------------------------"<<endl;
			return;
		}
		}
	}


}