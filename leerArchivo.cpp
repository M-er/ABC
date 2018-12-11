#include "leer.h"
Cfg leeArchivo(string nombre)
{
	char linea[2000];
	Cfg resultado;
	ifstream f(nombre.c_str());
	if(!f){
		cout << "\n\t\t*** El archivo " <<  (nombre).c_str() << " no existe! ***" << endl << endl;
		exit(0);
	}
	int a=0;
	float pc=0.0;
	int parametro=0;
	unsigned long int b=0;
	while(f.getline(linea,2000,'\n')){
		a=-1;
		switch(parametro){
			case 0: sscanf(linea,"%d        %ld      ", &a, &b);
			resultado.dimension_filas = a;
			resultado.dimension_columnas = b;
			resultado.numero_fuentes = a*b;
			break;
			case 1: sscanf(linea,"%d                ", &a);
			resultado.limite_permitido = a;
			break;
			case 2: sscanf(linea,"%d                ", &a);
			resultado.iteraciones = a;
			break;
			case 3: sscanf(linea,"%d        %ld      ", &a, &b);
			resultado.seleccion = a;
			resultado.tamano = b;
			break;
			case 4: sscanf(linea,"%d        %ld      ", &a, &b);
			resultado.recombina = a;
			resultado.rec_par = b/100.0;
			break;
			/*case 5: sscanf(linea,"%d                ", &a);
			resultado.tipo_vecindario = a;
			break;*/
		}
		parametro++;
	}
	return resultado;
}
