#include <stdio.h>
#include <stdlib.h>
#include "abejas.h"
#include "benchmark.h"

int limite_fit;
Fuente mejor_fuente;
int **vec_vecindario;
double fitness_previo;
double top_fitness;
Fuente **matrizColonia, **matrizColoniaAux;

#define POSITIVO (1)
#define NEGATIVO (-(1))
#define ZERO (0)
void crear_empleada_matriz(Cfg parametros, int dimf){
	top_fitness = mejor_fitness(parametros.funcion_o-1);
	int ntodo, ncolu, nfila, ncolu1, nfila1;
	int dimFilas = parametros.dimension_filas;
	int dimColus = parametros.dimension_columnas;
	for (int f = 0; f < dimFilas; f++){
		for(int c = 0; c < dimColus; c++){
			conformar_vecindario(parametros.dimension_filas, parametros.dimension_columnas, f, c, parametros.tipo_vecindario);
			switch (parametros.tipo_vecindario) {
				case 1:
				ntodo = Random_Entero(0,99999999)%4;
				break;
				case 2:
				ntodo = Random_Entero(0,99999999)%7;
				break;
				case 3:
				ntodo = Random_Entero(0,99999999)%7;
				break;
			}
			nfila = vec_vecindario[ntodo][0];
			ncolu = vec_vecindario[ntodo][1];
			recombina_matriz(parametros, f, c, nfila, ncolu, parametros.recombina);
		}
	}
	actualizar_colonia_matriz(parametros, dimf);
}
void recombina_matriz(Cfg parametros, int f, int c, int filaux, int colaux, int par){
	int pos, j;
	Sol aux_sol;
	double l, aux_fit;
	int dimf = parametros.funcion_dim;
	switch(par){
		case 1:{
			/* Recombinacion tradicional */
			pos = Random_Entero(0,dimf);
			for (int j=0; j<dimf; j++) {
				if(j==pos) matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j] + (-1 + drand48()*(1+1)) * (matrizColonia[f][c].solucion[j] - matrizColonia[filaux][colaux].solucion[j]);
				else matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
			}
			controlar_limite (matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
		}
		break;
		case 2:{
		/* Recombinacion aritmetica */
			pos = Random_Entero(0,dimf);
			l=drand48();
			for (int j=0; j<dimf; j++)
				if(j==pos) matrizColoniaAux[f][c].solucion[j] = (1-l)*matrizColonia[f][c].solucion[j] + l * matrizColonia[filaux][colaux].solucion[j];
			else matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
			controlar_limite (matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
		}
		break;

		case 3:{
		/* Recombinacion binomial  */
			int ki;
			double ro;
			ro = parametros.rec_par;
			for (int j=0; j<dimf; j++){
				if(drand48()<=ro) matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
				else matrizColoniaAux[f][c].solucion[j] = matrizColonia[filaux][colaux].solucion[j];
			}
			ki = Random_Entero(0,dimf);
			if(matrizColoniaAux[f][c].solucion[ki] == matrizColonia[f][c].solucion[ki])
				matrizColoniaAux[f][c].solucion[ki] = matrizColonia[filaux][colaux].solucion[ki];
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
		}
		break;

		case 4:{
		/* Recombinacion linear */
			aux_sol = (double*) calloc (dimf, sizeof(double));
			pos = Random_Entero(0,dimf);
			for (int j=0; j<dimf; j++){
				if(j==pos){
					matrizColoniaAux[f][c].solucion[j] = 0.5*(matrizColonia[f][c].solucion[j]+matrizColonia[filaux][colaux].solucion[j]);
					aux_sol[j] = 1.5*matrizColonia[f][c].solucion[j]-0.5*matrizColonia[filaux][colaux].solucion[j];
				}else matrizColoniaAux[f][c].solucion[j] = aux_sol[j] = matrizColonia[f][c].solucion[j];
			}
			controlar_limite (matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			controlar_limite (aux_sol, parametros.funcion_o,parametros.funcion_dim);
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			aux_fit = evaluar_fitness(aux_sol, parametros.funcion_o,parametros.funcion_dim);

			if(matrizColoniaAux[f][c].fitness > aux_fit){
				for (int j=0; j<dimf; j++){
					if(j==pos) matrizColoniaAux[f][c].solucion[j] = (-1)*(0.5)*matrizColonia[f][c].solucion[j]+1.5*matrizColonia[filaux][colaux].solucion[j];
					else matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
				}
				controlar_limite (matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
				matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			}else{
				for (int j=0; j<dimf; j++){
					if(pos==j)	aux_sol[j] = (-1)*(0.5)*matrizColonia[f][c].solucion[j]+1.5*matrizColonia[filaux][colaux].solucion[j];
					else aux_sol[j] = matrizColonia[f][c].solucion[j];
				}
				controlar_limite (aux_sol, parametros.funcion_o,parametros.funcion_dim);
				aux_fit = evaluar_fitness(aux_sol, parametros.funcion_o,parametros.funcion_dim);
			}
			if(matrizColoniaAux[f][c].fitness > aux_fit){
				for (int j=0; j<dimf; j++){
					matrizColoniaAux[f][c].solucion[j] = aux_sol[j];
				}
				matrizColoniaAux[f][c].fitness = aux_fit;
			}
			free(aux_sol);
		}
		break;
		/* Recombinacion de un punto */
		case 5:{
			int point;
			point=Random_Entero(0,dimf);
			for(j=0;j<point;j++)
				matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
			for(;j<dimf;j++)
				matrizColoniaAux[f][c].solucion[j] = matrizColonia[filaux][colaux].solucion[j];
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
		}
		break;
		/* Recombinacion de n puntos */
		case 6:{
			int nP = 0;
			int nP1 = 0;
			int np_i = 0;
			nP = Random_Entero(0,dimf/4);
			while(np_i<nP){
				nP1=Random_Entero(nP1,dimf);
				for(;j<nP1;j++)
					matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
				nP1=Random_Entero(nP1,dimf);
				for(;j<nP1;j++)
					matrizColoniaAux[f][c].solucion[j] = matrizColonia[filaux][colaux].solucion[j];
				np_i+=2;
			}
			if(nP % 2 == 0){
				for(;j<dimf;j++)
					matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];

			}else for(;j<dimf;j++) matrizColoniaAux[f][c].solucion[j] = matrizColonia[filaux][colaux].solucion[j];

			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
		}
		break;
		/* Caso 7 */
		case 7:{
			l = drand48();
			pos = Random_Entero(0,dimf);
			aux_sol = (double*) calloc (dimf, sizeof(double));
			for (int j=0; j<dimf; j++){
				if(j==pos){
					matrizColoniaAux[f][c].solucion[j] = l * matrizColonia[f][c].solucion[j] + (1-l) * matrizColonia[filaux][colaux].solucion[j];
					aux_sol[j] =(1-l) * matrizColonia[f][c].solucion[j] + l * matrizColonia[filaux][colaux].solucion[j];
				}
				else matrizColoniaAux[f][c].solucion[j] = aux_sol[j] = matrizColonia[f][c].solucion[j];
			}
			controlar_limite (matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			controlar_limite (aux_sol, parametros.funcion_o,parametros.funcion_dim);
			matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			aux_fit = evaluar_fitness(aux_sol, parametros.funcion_o,parametros.funcion_dim);

			if(matrizColoniaAux[f][c].fitness > aux_fit){
				for (int j=0; j<dimf; j++){
					if(j==pos) matrizColoniaAux[f][c].solucion[j] = (matrizColonia[f][c].solucion[j] < matrizColonia[filaux][colaux].solucion[j]  ?  matrizColonia[f][c].solucion[j] : matrizColonia[filaux][colaux].solucion[j]);
					else matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
				}
				matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			}
			else{
				for (int j=0; j<dimf; j++){
					if(pos==j) aux_sol[j] = (matrizColonia[f][c].solucion[j] < matrizColonia[filaux][colaux].solucion[j]  ?  matrizColonia[f][c].solucion[j] : matrizColonia[filaux][colaux].solucion[j]);
					else  aux_sol[j] = matrizColonia[f][c].solucion[j];
				}
				aux_fit = evaluar_fitness(aux_sol, parametros.funcion_o,parametros.funcion_dim);
			}
			if(matrizColoniaAux[f][c].fitness > aux_fit){
				for (int j=0; j<dimf; j++){
					if(j==pos) matrizColoniaAux[f][c].solucion[j] = (matrizColonia[f][c].solucion[j] > matrizColonia[filaux][colaux].solucion[j]  ?  matrizColonia[f][c].solucion[j] : matrizColonia[filaux][colaux].solucion[j]);
					else matrizColoniaAux[f][c].solucion[j] = matrizColonia[f][c].solucion[j];
				}
				matrizColoniaAux[f][c].fitness = evaluar_fitness(matrizColoniaAux[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
			}
			else{
				for (int j=0; j<dimf; j++){
					if(pos==j)	aux_sol[j] = (matrizColonia[f][c].solucion[j] > matrizColonia[filaux][colaux].solucion[j]  ?  matrizColonia[f][c].solucion[j] : matrizColonia[filaux][colaux].solucion[j]);
					else  aux_sol[j] = matrizColonia[f][c].solucion[j];
				}
				aux_fit = evaluar_fitness(aux_sol, parametros.funcion_o,parametros.funcion_dim);
			}
			if(matrizColoniaAux[f][c].fitness > aux_fit){
				for (int j=0; j<dimf; j++)
					matrizColoniaAux[f][c].solucion[j] = aux_sol[j];
				matrizColoniaAux[f][c].fitness = aux_fit;
			}
			free(aux_sol);
		}
		break;
		default: printf("no recombina");
		break;
	}
}
void exploradora_matriz(Cfg parametros, int dimf){
	for(int f=0; f<parametros.dimension_filas ; f++){
		for(int c=0; c<parametros.dimension_columnas ; c++){
			if (matrizColonia[f][c].limite >= parametros.limite_permitido){
				inicializacion(matrizColonia[f][c].solucion, parametros.funcion_o,parametros.funcion_dim);
				matrizColonia[f][c].fitness = evaluar_fitness (matrizColonia[f][c].solucion, parametros.funcion_o,parametros.funcion_dim );
				matrizColonia[f][c].limite = 0;
			}
		}
	}
}

void torneo_matriz(Cfg parametros, int& f, int& c, int dimf){
	int pos, pos1;
	pos = Random_Entero(0,999999) % 7;
	pos1 = Random_Entero(0,999999) % 7;
	if(max_min(parametros.funcion_o)== 0){
		if ((matrizColonia[vec_vecindario[pos][0]][vec_vecindario[pos][1]].fitness) < (matrizColonia[vec_vecindario[pos1][0]][vec_vecindario[pos1][1]].fitness)){
			f = vec_vecindario[pos][0];
			c = vec_vecindario[pos][1];
		}else{
			f = vec_vecindario[pos1][0];
			c = vec_vecindario[pos1][1];
		}
	}else{
		if (matrizColonia[vec_vecindario[pos][0]][vec_vecindario[pos][1]].fitness > matrizColonia[vec_vecindario[pos1][0]][vec_vecindario[pos1][1]].fitness){
			f = vec_vecindario[pos][0];
			c = vec_vecindario[pos][1];
		}else{
			f = vec_vecindario[pos1][0];
			c = vec_vecindario[pos1][1];
		}

	}
}
void prepara_ruleta(Cfg parametros, double &p_promedio){
	for(int f=0; f<parametros.dimension_filas; f++){
		for (int c = 0; c < parametros.dimension_columnas; ++c){
			p_promedio += matrizColonia[f][c].fitness;
		}
	}
}
int ruleta_matriz(Cfg parametros, int& f, int& c, int dimf, double p_promedio){
	double r;
	r = drand48();
	float pi = 0.0;
	bool sigo = true;
	for(int f=0; f<parametros.dimension_filas && sigo; f++){
		for(int c=0; c<parametros.dimension_columnas && sigo; c++){
			pi += matrizColonia[f][c].fitness / p_promedio;
			if(pi >= r) sigo = false;
		}
	}
	return 0;
}

void observadora_matriz(Cfg parametros, int dimf){
	int k, i;
	double prom_col = 0;
	int faux = 0;
	int caux = 0;
	if(parametros.seleccion==2) prepara_ruleta(parametros,prom_col);
	for(int f=0; f<parametros.dimension_filas; f++){
		for(int c=0; c < parametros.dimension_columnas; c++){
			switch ( parametros.seleccion ){
				case 1:
				conformar_vecindario(parametros.dimension_filas, parametros.dimension_columnas, f, c, parametros.tipo_vecindario);
				faux = f;
				caux = c;
				torneo_matriz(parametros, faux, caux, dimf);
				break;
				case 2:
				ruleta_matriz(parametros, faux, caux, dimf, prom_col);
				break;
				default:
				break;
			}
			recombina_matriz(parametros, f, c, faux, caux, 1);
		}
	}
	actualizar_colonia_matriz(parametros, dimf);
}


void actualizar_colonia_matriz(Cfg parametros, int dimf){
	for(int f = 0; f < parametros.dimension_filas; f++){
		for(int c = 0; c < parametros.dimension_columnas; c++){
			if(max_min(parametros.funcion_o)== 0){
				if (matrizColoniaAux[f][c].fitness < matrizColonia[f][c].fitness){
					matrizColonia[f][c].fitness = matrizColoniaAux[f][c].fitness;
					matrizColonia[f][c].limite = 0;
					for (int j=0; j<dimf; j++)
						matrizColonia[f][c].solucion[j] = matrizColoniaAux[f][c].solucion[j];
				}else
				matrizColonia[f][c].limite = matrizColonia[f][c].limite + 1;
			}else{
				if (matrizColoniaAux[f][c].fitness > matrizColonia[f][c].fitness){
					matrizColonia[f][c].fitness = matrizColoniaAux[f][c].fitness;
					matrizColonia[f][c].limite = 0;
					for (int j=0; j<dimf; j++)
						matrizColonia[f][c].solucion[j] = matrizColoniaAux[f][c].solucion[j];
				}else
				matrizColonia[f][c].limite = matrizColonia[f][c].limite + 1;
			}
		}
	}
}
void conformar_vecindario(int maxF, int maxC, int f, int c, int metodo){
	switch(metodo){
		//Abajo, Izquierda, Arriba, Derecha, Centro
		case 1:{
			default:
			if(f == 0){
				if(c == 0){
					vec_vecindario[0][0] = f+1;
					vec_vecindario[0][1] = c;
					vec_vecindario[1][0] = f;
					vec_vecindario[1][1] = maxC-1;
					vec_vecindario[2][0] = maxF-1;
					vec_vecindario[2][1] = c;
					vec_vecindario[3][0] = f;
					vec_vecindario[3][1] = c+1;
					vec_vecindario[4][0] = f;
					vec_vecindario[4][1] = c;
				}else{
					if(c == maxC - 1){
						vec_vecindario[0][0] = f+1;
						vec_vecindario[0][1] = c;
						vec_vecindario[1][0] = f;
						vec_vecindario[1][1] = c-1;
						vec_vecindario[2][0] = maxF-1;
						vec_vecindario[2][1] = c;
						vec_vecindario[3][0] = 0;
						vec_vecindario[3][1] = 0;
						vec_vecindario[4][0] = f;
						vec_vecindario[4][1] = c;
					}else{
						vec_vecindario[0][0] = f+1;
						vec_vecindario[0][1] = c;
						vec_vecindario[1][0] = f;
						vec_vecindario[1][1] = c-1;
						vec_vecindario[2][0] = maxF-1;
						vec_vecindario[2][1] = c;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = c+1;
						vec_vecindario[4][0] = f;
						vec_vecindario[4][1] = c;
					}
				}
			}else{
				if(f == maxF - 1){
					if(c == 0){
						vec_vecindario[0][0] = 0;
						vec_vecindario[0][1] = c;
						vec_vecindario[1][0] = f;
						vec_vecindario[1][1] = maxC-1;
						vec_vecindario[2][0] = f-1;
						vec_vecindario[2][1] = c;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = 1;
						vec_vecindario[4][0] = f;
						vec_vecindario[4][1] = c;
					}else{
						if(c == maxC - 1){
							vec_vecindario[0][0] = 0;
							vec_vecindario[0][1] = c;
							vec_vecindario[1][0] = f;
							vec_vecindario[1][1] = c-1;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = 0;
							vec_vecindario[4][0] = f;
							vec_vecindario[4][1] = c;
						}else{
							vec_vecindario[0][0] = 0;
							vec_vecindario[0][1] = c;
							vec_vecindario[1][0] = f;
							vec_vecindario[1][1] = c-1;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = c+1;
							vec_vecindario[4][0] = f;
							vec_vecindario[4][1] = c;
						}
					}
				}else{
					if(c == 0){
						vec_vecindario[0][0] = f+1;
						vec_vecindario[0][1] = c;
						vec_vecindario[1][0] = f;
						vec_vecindario[1][1] = maxC-1;
						vec_vecindario[2][0] = f-1;
						vec_vecindario[2][1] = c;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = c+1;
						vec_vecindario[4][0] = f;
						vec_vecindario[4][1] = c;
					}else{
						if(c == maxC - 1){
							vec_vecindario[0][0] = f+1;
							vec_vecindario[0][1] = c;
							vec_vecindario[1][0] = f;
							vec_vecindario[1][1] = c-1;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = 0;
							vec_vecindario[4][0] = f;
							vec_vecindario[4][1] = c;
						}else{
							vec_vecindario[0][0] = f+1;
							vec_vecindario[0][1] = c;
							vec_vecindario[1][0] = f;
							vec_vecindario[1][1] = c-1;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = c+1;
							vec_vecindario[4][0] = f;
							vec_vecindario[4][1] = c;
						}
					}
				}
			}
		}
		break;
		//Izquierda, Arriba, Derecha, Abajo, Centro
		case 2:{
			if(f == 0){
				if(c == 0){
					vec_vecindario[0][0] = maxF-1;
					vec_vecindario[0][1] = maxC-1;
					vec_vecindario[1][0] = maxF-1;
					vec_vecindario[1][1] = c;
					vec_vecindario[2][0] = maxF-1;
					vec_vecindario[2][1] = c+1;
					vec_vecindario[3][0] = f;
					vec_vecindario[3][1] = c+1;
					vec_vecindario[4][0] = f+1;
					vec_vecindario[4][1] = c+1;
					vec_vecindario[5][0] = f+1;
					vec_vecindario[5][1] = c;
					vec_vecindario[6][0] = f+1;
					vec_vecindario[6][1] = maxC-1;
					vec_vecindario[7][0] = f;
					vec_vecindario[7][1] = maxC-1;
					vec_vecindario[8][0] = f;
					vec_vecindario[8][1] = c;
				}else{
					if(c == maxC - 1){
						vec_vecindario[0][0] = maxF-1;
						vec_vecindario[0][1] = c-1;
						vec_vecindario[1][0] = maxF-1;
						vec_vecindario[1][1] = c;
						vec_vecindario[2][0] = maxF-1;
						vec_vecindario[2][1] = 0;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = 1;
						vec_vecindario[4][0] = f+1;
						vec_vecindario[4][1] = 0;
						vec_vecindario[5][0] = f+1;
						vec_vecindario[5][1] = c;
						vec_vecindario[6][0] = f+1;
						vec_vecindario[6][1] = c-1;
						vec_vecindario[7][0] = f;
						vec_vecindario[7][1] = c-1;
						vec_vecindario[8][0] = f;
						vec_vecindario[8][1] = c;
					}else{
						vec_vecindario[0][0] = maxF-1;
						vec_vecindario[0][1] = c-1;
						vec_vecindario[1][0] = maxF-1;
						vec_vecindario[1][1] = c;
						vec_vecindario[2][0] = maxF-1;
						vec_vecindario[2][1] = c+1;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = c+1;
						vec_vecindario[4][0] = f+1;
						vec_vecindario[4][1] = c+1;
						vec_vecindario[5][0] = f+1;
						vec_vecindario[5][1] = c;
						vec_vecindario[6][0] = f+1;
						vec_vecindario[6][1] = c-1;
						vec_vecindario[7][0] = f;
						vec_vecindario[7][1] = c-1;
						vec_vecindario[8][0] = f;
						vec_vecindario[8][1] = c;
					}
				}
			}else{
				if(f == maxF - 1){
					if(c == 0){
						vec_vecindario[0][0] = f-1;
						vec_vecindario[0][1] = maxC-1;
						vec_vecindario[1][0] = f-1;
						vec_vecindario[1][1] = c;
						vec_vecindario[2][0] = f-1;
						vec_vecindario[2][1] = c+1;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = c+1;
						vec_vecindario[4][0] = 0;
						vec_vecindario[4][1] = c+1;
						vec_vecindario[5][0] = 0;
						vec_vecindario[5][1] = c;
						vec_vecindario[6][0] = 0;
						vec_vecindario[6][1] = maxC-1;
						vec_vecindario[7][0] = f;
						vec_vecindario[7][1] = maxC-1;
						vec_vecindario[8][0] = f;
						vec_vecindario[8][1] = c;
					}else{
						if(c == maxC - 1){
							vec_vecindario[0][0] = f-1;
							vec_vecindario[0][1] = c-1;
							vec_vecindario[1][0] = f-1;
							vec_vecindario[1][1] = c;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = 0;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = 0;
							vec_vecindario[4][0] = 0;
							vec_vecindario[4][1] = 0;
							vec_vecindario[5][0] = 0;
							vec_vecindario[5][1] = c;
							vec_vecindario[6][0] = 0;
							vec_vecindario[6][1] = c-1;
							vec_vecindario[7][0] = f;
							vec_vecindario[7][1] = c-1;
							vec_vecindario[8][0] = f;
							vec_vecindario[8][1] = c;
						}else{
							vec_vecindario[0][0] = f-1;
							vec_vecindario[0][1] = c-1;
							vec_vecindario[1][0] = f-1;
							vec_vecindario[1][1] = c;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c+1;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = c+1;
							vec_vecindario[4][0] = 0;
							vec_vecindario[4][1] = c+1;
							vec_vecindario[5][0] = 0;
							vec_vecindario[5][1] = c;
							vec_vecindario[6][0] = 0;
							vec_vecindario[6][1] = c-1;
							vec_vecindario[7][0] = f;
							vec_vecindario[7][1] = c-1;
							vec_vecindario[8][0] = f;
							vec_vecindario[8][1] = c;
						}
					}
				}else{
					if(c == 0){
						vec_vecindario[0][0] = f-1;
						vec_vecindario[0][1] = maxC-1;
						vec_vecindario[1][0] = f-1;
						vec_vecindario[1][1] = c;
						vec_vecindario[2][0] = f-1;
						vec_vecindario[2][1] = c+1;
						vec_vecindario[3][0] = f;
						vec_vecindario[3][1] = c+1;
						vec_vecindario[4][0] = f+1;
						vec_vecindario[4][1] = c+1;
						vec_vecindario[5][0] = f+1;
						vec_vecindario[5][1] = c;
						vec_vecindario[6][0] = f+1;
						vec_vecindario[6][1] = maxC-1;
						vec_vecindario[7][0] = f;
						vec_vecindario[7][1] = maxC-1;
						vec_vecindario[8][0] = f;
						vec_vecindario[8][1] = c;
					}else{
						if(c == maxC - 1){
							vec_vecindario[0][0] = f-1;
							vec_vecindario[0][1] = c-1;
							vec_vecindario[1][0] = f-1;
							vec_vecindario[1][1] = c;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = 0;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = 0;
							vec_vecindario[4][0] = f+1;
							vec_vecindario[4][1] = 0;
							vec_vecindario[5][0] = f+1;
							vec_vecindario[5][1] = c;
							vec_vecindario[6][0] = f+1;
							vec_vecindario[6][1] = c-1;
							vec_vecindario[7][0] = f;
							vec_vecindario[7][1] = c-1;
							vec_vecindario[8][0] = f;
							vec_vecindario[8][1] = c;
						}else{
							vec_vecindario[0][0] = f-1;
							vec_vecindario[0][1] = c-1;
							vec_vecindario[1][0] = f-1;
							vec_vecindario[1][1] = c;
							vec_vecindario[2][0] = f-1;
							vec_vecindario[2][1] = c+1;
							vec_vecindario[3][0] = f;
							vec_vecindario[3][1] = c+1;
							vec_vecindario[4][0] = f+1;
							vec_vecindario[4][1] = c+1;
							vec_vecindario[5][0] = f+1;
							vec_vecindario[5][1] = c;
							vec_vecindario[6][0] = f+1;
							vec_vecindario[6][1] = c-1;
							vec_vecindario[7][0] = f;
							vec_vecindario[7][1] = c-1;
							vec_vecindario[8][0] = f;
							vec_vecindario[8][1] = c;
						}
					}
				}
			}
		}
		break;
		//Izquierda, Arriba, Derecha, Abajo, Centro
		case 3:{
			// L9
			/*
			0 0 0 0 0 0 0
			0 0 0 x 0 0 0
			0 0 0 x 0 0 0
			0 x x X x x
			0 0 0 x 0 0 0
			0 0 0 x 0 0 0
			0 0 0 0 0 0 0
			*/
			//Izquierda
			if(c-1 < 0){
				vec_vecindario[0][0] = f;
				vec_vecindario[0][1] = maxC-1;
				vec_vecindario[1][0] = f;
				vec_vecindario[1][1] = maxC-2;
			}else{
				if(c-2 < 0){
					vec_vecindario[0][0] = f;
					vec_vecindario[0][1] = c-1;
					vec_vecindario[1][0] = f;
					vec_vecindario[1][1] = maxC-1;
				}else{
				vec_vecindario[0][0] = f;
				vec_vecindario[0][1] = c-1;
				vec_vecindario[1][0] = f;
				vec_vecindario[1][1] = c-2;
				}
			}
			//Arriba
			if(f-1 < 0){
				vec_vecindario[2][0] =  maxF-1;
				vec_vecindario[2][1] =  c;
				vec_vecindario[3][0] =  maxF-2;
				vec_vecindario[3][1] =  c;
			}else{
				if(f-2 < 0){
					vec_vecindario[2][0] =  f-1;
					vec_vecindario[2][1] =  c;
					vec_vecindario[3][0] =  maxF-1;
					vec_vecindario[3][1] =  c;
				}else{
				vec_vecindario[2][0] =  f-1;
				vec_vecindario[2][1] =  c;
				vec_vecindario[3][0] =  f-2;
				vec_vecindario[3][1] =  c;
				}
			}
			//Derecha
			if(f+1 > maxF-1){
				vec_vecindario[4][0] =  0;
				vec_vecindario[4][1] =  c;
				vec_vecindario[5][0] =  1;
				vec_vecindario[5][1] =  c;
			}else{
				if(f+2 > maxF-1){
					vec_vecindario[4][0] =  f+1;
					vec_vecindario[4][1] =  c;
					vec_vecindario[5][0] =  0;
					vec_vecindario[5][1] =  c;
				}else{
					vec_vecindario[4][0] =  f+1;
					vec_vecindario[4][1] =  c;
					vec_vecindario[5][0] =  f+2;
					vec_vecindario[5][1] =  c;
				}
			}
			//Abajo
			if(c+1 >= maxC-1){
				vec_vecindario[6][0] =  f;
				vec_vecindario[6][1] =  0;
				vec_vecindario[7][0] =  f;
				vec_vecindario[7][1] =  1;

			}else{
				if(c+2 >= maxC-1){
						vec_vecindario[6][0] =  f;
						vec_vecindario[6][1] =  c+1;
						vec_vecindario[7][0] =  f;
						vec_vecindario[7][1] =  0;
					}else{
					vec_vecindario[6][0] =  f;
					vec_vecindario[6][1] =  c+1;
					vec_vecindario[7][0] =  f;
					vec_vecindario[7][1] =  c+2;
				}
			}
			vec_vecindario[8][0] = f;
			vec_vecindario[8][1] = c;
			break;
		}
	}
}

void muestra_vecindario(int dimF, int dimC){
	for(int i=0; i<dimF; i++){
		for(int j=0; j<dimC; j++){
			cout << " " << matrizColonia[i][j].fitness;
		}
		cout << endl;
	}
}
void muestra_vecindario_sol(int dimF, int dimC, int dimf){
	for(int i=0; i<dimF; i++){
		for(int j=0; j<dimC; j++){
			for(int k=0; k<dimf; k++)
				cout << " " << matrizColonia[i][j].solucion[k] <<endl;
		}
		cout << endl;
	}
}
void muestra_vec_vecindario(int dimF, int dimC){
	for(int i=0; i<dimF; i++){
		for(int j=0; j<dimC; j++){
			cout << " " << vec_vecindario[i][j];
		}
		cout << endl;
	}
}
