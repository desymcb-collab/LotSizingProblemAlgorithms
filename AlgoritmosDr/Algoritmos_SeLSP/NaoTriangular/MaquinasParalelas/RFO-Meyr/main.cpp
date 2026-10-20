#include <iostream>
#include "../Raiz/raiz.h"
#include "Base/rf.h"

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

    for(int i = 0; i <= 1; i++)
        for(int j = 1; j <= 9; j++)
        {
            printf("\nT%d%d.txt\n",i,j);

            /*******************************Leitura dos dados da instancia*********************************************/

            sprintf(ins.fileinfo,"../../../Instancias/T01_T12/T%d%d.txt",i,j);

            clock_gettime(CLOCK_MONOTONIC, &start);
            elapsedTime = start.tv_sec;
            ins.le_dados_M(ins.fileinfo);
            //ins.show_instance();
            //getchar();

            /******************************************************************************************/

            /******************************Resolvendo com CPLEX*********************************/

            RF * rf21 = new RF(&ins,elapsedTime);

            sprintf(ins.pathname,"resultados/T%d%d.txt",i,j);

            rf21->RelaxAndFix();

            ins.desaloca_variaveis();
            clock_gettime(CLOCK_MONOTONIC, &endIt);
            elapsedTime = (endIt.tv_sec - start.tv_sec);
            elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

            ///Exit file with cplex results
            rf21->saida = fopen(ins.pathname,"a");
            if(rf21->saida == NULL) printf("Error in saida file!!!\n");
            fprintf(rf21->saida,"T%d%d\t%.4lf %.4lf\n",i,j,rf21->BestSolution,elapsedTime);
            fclose(rf21->saida);
            delete rf21;

            /***********************************************************************************/
            if(i==1 && j==2)
                j=10;
        }
}
