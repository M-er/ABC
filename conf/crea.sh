#!/bin/bash
d=100
for j in {1..6}
do
  for v in {1..3}
  do
    for k in {1..3}
    do
      case "$k" in
        "1")
        d=100
        ;;
        "2")
        d=500
        ;;
        "3")
        d=1000
        ;;
      esac
      echo "
7 7    //Dimension (Fila, Columna)
10     //limite_permitido seteado por programa
100000 //iteraciones(generaciones)
1	2    //seleccion - 1) Torneo y 2) Ruleta
1	0    //recombinacion  (operador, parametro)
"$v"      //tipo_vecindario - 1) L5, 2) C9, 3) L9
"$d"    //dimension 100, 500, 1000
"$j"      //funcion, del 1 al 6
      " > 'abcF'$j'D'$k'V'$v'.cfg'
    done
  done
done
