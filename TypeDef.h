#ifndef TypeDef
#define TypeDef


typedef double * Sol;

/*----------Definimos una estructura para las abejas-----------*/
struct Fuente{
	Sol solucion;  
	double fitness;
	int limite;
};

#endif
