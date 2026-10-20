#include <iostream>
#include "raiz.h"
#include "scf_fl.h"

void SolveInstances();

int main()
{
    cout << "*********----Code Initialized----********* \n\n";
    SolveInstances();

    cout << "\n\n***********----Code Ended----********* \n\n";
    return 0;
}

void SolveInstances()
{
    /****Instance class*****/
    Instance ins;
    FILE *saida;
    double elapsedTime;
    struct timespec start,endIt;
    int qt_period,index_items,qt,id,id2,id3,id4,id5;
    int qt_items[4];
    int p[2];
    int cd[2],prob1[2],prob2[2];

    double cl1[2];
    double cl2 = 0.5;

    qt_items[0] = 15; qt_items[1] = 20; qt_items[2] = 15; qt_items[3] = 400;
    cl1[0] = 0.6; cl1[1] = 0.8;
    cd[0] = 50; cd[1] = 100;
    prob1[0] = 60; prob1[1] = 80;
    prob2[0] = 10; prob2[1] = 20;
    p[0] = 2; p[1] = 3;

    /******************************************************************************************/
    /*******************************Leitura dos dados da instancia*****************************/
    /******************************************************************************************/

    sprintf(ins.fileinfo,"../../Instancias/instancias_eli/mes6.dat");
    printf("%s\n",ins.fileinfo);

    clock_gettime(CLOCK_MONOTONIC, &start);
    elapsedTime = start.tv_sec;

    ins.le_dados_SM(ins.fileinfo);

    //ins.show_instance();
    //getchar();
    /******************************************************************************************/
    /*******************************Resolvendo com RelaxAndFix*********************************/
    /******************************************************************************************/
    Models * md21 = new Models(&ins,elapsedTime);
    sprintf(ins.pathname,"resultados/SCF_FL1.txt");
    md21->RodaGuimaraes();
    ins.desaloca_variaveis();
    clock_gettime(CLOCK_MONOTONIC, &endIt);
    elapsedTime = (endIt.tv_sec - start.tv_sec);
    elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

    ///Exit file with cplex results
    md21->saida = fopen(ins.pathname,"a");
    if(md21->saida == NULL) printf("Error in saida file!!!\n");
    fprintf(md21->saida,"eli_mesA\t%.4lf\t%.4lf\n",md21->BestSolution,elapsedTime);
    fclose(md21->saida);
    delete md21;
    /*******************************************************************************************/
}
