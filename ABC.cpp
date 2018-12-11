
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <new>
#include <string.h>
#include <limits>
#include <time.h>

#include "TypeDef.h"
#include "leer.h"
#include "benchmark.h"
#include "abejas.h"
#include "ABC.h"

using namespace std;

void almacenar(int iteraccion, Fuente mejor_fuente, Cfg parametros, int generacion, int numEval, double fuenteProm, double limiteProm);
void algoritmo_ABC (Cfg parametros);

int main (int argc, char **argv){
	Cfg param;
	param=leeArchivo(argv[1]);
	//param.funcion_o=atoi(argv[2]);
	//param.funcion_dim=atoi(argv[3]);
	//param.tipo_vecindario=atoi(argv[4]);
	algoritmo_ABC(param);
	return 1;
}
extern double numEval;
extern Fuente mejor_fuente;
extern int **vec_vecindario;
extern Fuente **matrizColonia;
extern Fuente **matrizColoniaAux;
extern double top_fitness;

/*----------Definimos variables necesarias-----------*/
int parcial = 0;
time_t t_inicial;
double t_mejor, t_total, bestEval;
/*----------Almacenar-----------*/
void almacenar(int iteraccion, double mejor_fitness, Cfg parametros, int generacion, int numEval, double fuenteProm, double limiteProm, double t_mejor,double t_total){
	FILE *fp;
	int j;
	short i;
	char nombre[50], linea[100];
	fp = fopen ("./res/parciales.txt", "a+");
	fprintf(fp, "%d	\t", parametros.funcion_o);
	fprintf(fp, "%d	\t", parametros.funcion_dim);
	fprintf(fp, "%d	\t", parametros.tipo_vecindario);
	//fprintf(fp, "%d	\t", parametros.recombina);
	fprintf(fp, "%6d	\t", parcial);
	fprintf(fp, "%6d	\t", generacion);
	fprintf(fp, "%7.0d	\t", numEval);
	fprintf(fp, "%3.3f	\t", mejor_fitness);
	fprintf(fp, "%3.3f	\t", fuenteProm);
	fprintf(fp, "%3.3f	\t", top_fitness);
	fprintf(fp, "%3.3f	\t", limiteProm);
	fprintf(fp, "%8d	\t", iteraccion);
	fprintf(fp, "%.3f\t", t_mejor);
	fprintf(fp, "%.3f\t", t_total);
	fprintf(fp, "\n");
	fclose(fp);
	parcial++;
}
/*----------Algoritmo ABC-----------*/
void algoritmo_ABC(Cfg parametros){
	int a;
	int tot;
	numEval = 0.0;
	bestEval= 0.0;
	t_mejor = 0.0;
	t_total = 0.0;
	int iteraccion = 0;
	int generacion = 1;
	Fuente mejor_fuente;
	double mejor_fitness = 0.0;
	int dimf = parametros.funcion_dim;
	double fuenteProm = 0.0, limiteProm=0.0;
	int dimFilas = parametros.dimension_filas;
	int dimColus = parametros.dimension_columnas;
	parametros.limite_permitido = parametros.numero_fuentes * dimf;

	srand(time(NULL)*10000);
	srand48(time(NULL)*10000);
	double totalEval = dimf*5000;

	mejor_fuente.solucion = (double*) calloc (dimf, sizeof(double));
	mejor_fitness = (max_min(parametros.funcion_o) == 1)? std::numeric_limits<double>::min() : std::numeric_limits<double>::max();
	matrizColonia = (Fuente **)malloc(dimFilas * sizeof(Fuente *));
	matrizColoniaAux = (Fuente **)malloc(dimFilas * sizeof(Fuente *));

	switch (parametros.tipo_vecindario) {
		case 1:
		tot = 12;
		break;
		case 2:
		tot = 18;
		break;
		case 3:
		tot = 16;
		default:
		break;
	}
	vec_vecindario = (int **)malloc(tot * sizeof(int));
	for(int k=0; k<tot; k++){
		vec_vecindario[k] = (int *)malloc(1 * sizeof(int));
	}
	for(int f=0; f<tot/2; f++){
		for(int c=0; c<2; c++){
			vec_vecindario[f][c] = 0;
		}
	}
	for(int f=0; f<dimFilas; f++){
		matrizColonia[f] = (Fuente *)malloc(dimColus * sizeof(Fuente));
		matrizColoniaAux[f] = (Fuente *)malloc(dimColus * sizeof(Fuente));
	}
	for(int f=0; f<dimFilas; f++){
		for(int c=0; c<dimColus; c++){
			matrizColonia[f][c].solucion = (double*) calloc (dimf, sizeof(double));
			matrizColonia[f][c].limite = 0;
			matrizColoniaAux[f][c].solucion = (double*) calloc (dimf, sizeof(double));
			matrizColoniaAux[f][c].limite = 0;
		}
	}
	for(int f=0; f<dimFilas; f++){
		for(int c=0; c<dimColus; c++){
			inicializacion(matrizColonia[f][c].solucion, parametros.funcion_o, dimf);
			matrizColonia[f][c].fitness = evaluar_fitness(matrizColonia[f][c].solucion, parametros.funcion_o,dimf);
		}
	}
	time(&t_inicial);
	while(numEval < totalEval){
		crear_empleada_matriz(parametros, dimf);
		observadora_matriz(parametros, dimf);
		exploradora_matriz(parametros, dimf);
		for(int f = 0; f < parametros.dimension_filas; f++){
			for(int c = 0; c < parametros.dimension_columnas; c++){
				if(max_min(parametros.funcion_o)){
					if (matrizColonia[f][c].fitness > mejor_fitness){
						mejor_fitness = matrizColonia[f][c].fitness;
						iteraccion = generacion;
						t_mejor = tiempo(t_inicial);
						bestEval = numEval;
						for (int j=0; j<dimf; j++){
							mejor_fuente.solucion[j] = matrizColonia[f][c].solucion[j];
						}
					}
				}else{
					if (matrizColonia[f][c].fitness < mejor_fitness){
						mejor_fitness = matrizColonia[f][c].fitness;
						iteraccion = generacion;
						t_mejor = tiempo(t_inicial);
						bestEval = numEval;
						for (int j=0; j<dimf; j++){
							mejor_fuente.solucion[j] = matrizColonia[f][c].solucion[j];
						}
					}
				}
			}
		}
		fuenteProm = 0.0;
		for(int f=0; f<dimFilas; f++){
			for(int c=0; c<dimColus; c++){
				fuenteProm += matrizColonia[f][c].fitness;
				limiteProm += matrizColonia[f][c].limite;
			}
		}
		fuenteProm = fuenteProm / parametros.numero_fuentes;
		limiteProm = limiteProm / parametros.numero_fuentes;
		t_total = tiempo(t_inicial);
		if(generacion % 500 == 0)
			almacenar(iteraccion, mejor_fitness, parametros, generacion, bestEval, fuenteProm, limiteProm, t_mejor, t_total);
		generacion = generacion + 1;
	}
	/*fuenteProm = 0.0;
	for(int f=0; f<dimFilas; f++){
		for(int c=0; c<dimColus; c++){
			fuenteProm += matrizColonia[f][c].fitness;
			limiteProm += matrizColonia[f][c].limite;
		}
	}*/
	//fuenteProm = fuenteProm / parametros.numero_fuentes;
	//limiteProm = limiteProm / parametros.numero_fuentes;
	//t_total = tiempo(t_inicial);
	//almacenar(iteraccion, mejor_fitness, parametros, generacion, bestEval, fuenteProm, limiteProm, t_mejor, t_total);
	/*
		cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * Resultados finales * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
		cout << " Funcion \t Dimension \t tipo_vecindario    Rec_Par \t Generacion \t Num eval \t Mejor    Fuente promedio    Limite promedio    N° iteracion" << endl;
		cout << "    " << parametros.funcion_o << "\t\t   " << parametros.funcion_dim << "\t\t    " << parametros.tipo_vecindario << "\t\t\t" << parametros.rec_par << "\t   " << generacion << "\t\t  " << numEval << " \t" << mejor_fuente.fitness << "\t" << fuenteProm/parametros.numero_fuentes << "\t\t   " << limiteProm/parametros.numero_fuentes << "\t\t   " << iteraccion <<  endl;
	*/
}
