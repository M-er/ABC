#ifndef leerArchivo
#define leerArchivo
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
struct Cfg{
	int   tamano;
	float rec_par;
	int   seleccion;
	int   funcion_o;
	int   recombina;
	int   funcion_dim;
	int   iteraciones;
	int   numero_fuentes;
	int   dimension_filas;
	int   tipo_vecindario;
	int   limite_permitido;
	int   dimension_columnas;
};
Cfg leeArchivo(string nombre);
#endif
