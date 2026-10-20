clear
clear
clear
clear
clear
clear
clear
clear

echo "--> Compilando -"
gcc -fopenmp structures.h utilities.h searches.h GRASP.h main.c  -I/usr/include/igraph -L/usr/local/lib -ligraph -o TCC2

echo "Grafo teste"
G_SLICE=always-malloc G_DEBUG=gc-friendly  valgrind -v --track-origins=yes --tool=memcheck --leak-check=full --num-callers=40 --log-file=valgrind.log ./TCC2 "teste"

echo "Grafo b02"
G_SLICE=always-malloc G_DEBUG=gc-friendly  valgrind -v --track-origins=yes --tool=memcheck --leak-check=full --num-callers=40 --log-file=valgrind1.log ./TCC2 "b02"
