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




    printf("\nDa:\n");

    /*******************************Leitura dos dados da instancia*********************************************/

    //sprintf(ins.fileinfo,"../../Instancias/Ches_problems/ches05.txt");
    sprintf(ins.fileinfo,"../../Instancias/Ches_problems/ches05.txt");

    clock_gettime(CLOCK_MONOTONIC, &start);
    elapsedTime = start.tv_sec;
    ins.le_dados_SM(ins.fileinfo);
    //ins.show_instance();
    //getchar();
    /******************************************************************************************/

    /******************************Resolvendo com CPLEX*********************************/


    Models * md21 = new Models(&ins,elapsedTime);
    //sprintf(ins.pathname,"resultados/ches05.txt");
    sprintf(ins.pathname,"resultados/ches05.txt");
    md21->RodaGuimaraes();
    ins.desaloca_variaveis();
    clock_gettime(CLOCK_MONOTONIC, &endIt);
    elapsedTime = (endIt.tv_sec - start.tv_sec);
    elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

    ///Exit file with cplex results
    md21->saida = fopen(ins.pathname,"a");
    if(md21->saida == NULL) printf("Error in saida file!!!\n");
    //fprintf(md21->saida,"De\t%.4lf %.4lf %.4lf\n",md21->BestSolution,md21->best_bound,elapsedTime);
    fprintf(md21->saida,"Df\t%.4lf %.4lf %.4lf\n",md21->BestSolution,md21->best_bound,elapsedTime);
    fclose(md21->saida);
    delete md21;

    /***********************************************************************************/

}
