all: m
	gcj -C FastFractal.java
	gcj -c FastFractal.java UnitFunction1D.java FractalFunction1D.java RanTable.java DoubleDip.java RanQD1.java -o fractal.o
	gcjh --classpath . FastFractal
	g++ fastfractal_doubledip.cpp -c -o fastfractal_doubledip.o
	g++ fastfractal_doubledip.o fractal.o m.o -lgcj -o exe

m:
	#g++ -c ABC.cpp -o m.o
	g++  -Wno-deprecated -L../lib -lm ABC.cpp leerArchivo.cpp benchmark.cpp abejas.cpp  -o abejas

clean:
	rm -f *.class *.o FastFractal.h *~ exe
