#include <iostream>
#include "../Raiz/raiz.h"
#include "Base/rfks.h"

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

    int qt_items[4];
    int m[2];
    int cd[2],prob1[2],prob2[2];

    double cl1[2];
    double cl2 = 0.5;
    double elapsedTime;
    struct timespec start,endIt;

    qt_items[0] = 15; qt_items[1] = 20; qt_items[2] = 15; qt_items[3] = 400;
    cl1[0] = 0.6; cl1[1] = 0.8;
    cd[0] = 50; cd[1] = 100;
    prob1[0] = 60; prob1[1] = 80;
    prob2[0] = 10; prob2[1] = 20;
    m[0] = 2; m[1] = 3;

   for(int index_m = 0; index_m <= 1; index_m++) // 0 a 1 p
        for(int index_cl1 = 1; index_cl1 <= 1; index_cl1++) // 0 a 1 cl1
            for(int index_cd = 1; index_cd <= 1; index_cd++) // 0 a 1 cd
                for(int index_prob1 = 0; index_prob1 <= 0; index_prob1++) // 0 a 1 prob1
                    for(int index_prob2 = 1; index_prob2 <= 1; index_prob2++) // 0 a 1 prob2
                        for(int qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(int index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(int index = 1; index <= 5; index++)  // 0 a 9
                                {
                                    printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",m[index_m],qt_items[index_items],qt_period,cl1[index_cl1],cl2,cd[index_cd],prob1[index_prob1],prob2[index_prob2],index);
                                    /******************************************************************************************/
                                    /*******************************Leitura dos dados da instancia*****************************/
                                    /******************************************************************************************/

                                    sprintf(ins.fileinfo,"../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",m[index_m],qt_items[index_items],qt_period,cl1[index_cl1],cl2,cd[index_cd],prob1[index_prob1],prob2[index_prob2],index);

                                    clock_gettime(CLOCK_MONOTONIC, &start);
                                    elapsedTime = start.tv_sec;
                                    ins.le_dados_SM(ins.fileinfo);
                                    //ins.show_instance();
                                    //getchar();
                                    /******************************************************************************************/
                                    /*******************************Resolvendo com RelaxAndFix*********************************/
                                    /******************************************************************************************/
                                    RF * rf21 = new RF(&ins,elapsedTime);
                                    rf21->RelaxAndFix();
                                    ins.desaloca_variaveis();
                                    clock_gettime(CLOCK_MONOTONIC, &endIt);
                                    elapsedTime = (endIt.tv_sec - start.tv_sec);
                                    elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;
                                    sprintf(ins.pathname,"resultados/RFKS-p1_b2.txt");
                                    rf21->saida = fopen(ins.pathname,"a");
                                    if(rf21->saida == NULL) printf("Error in saida file!!!\n");
                                    fprintf(rf21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d\t%.4lf\t%.4lf\t%.4lf\n",m[index_m],qt_items[index_items],qt_period,cl1[index_cl1],cl2,cd[index_cd],prob1[index_prob1],prob2[index_prob2],index,rf21->solutionRF,rf21->BestSolution,elapsedTime);
                                    fclose(rf21->saida);
                                    delete rf21;
                                    /******************************************************************************************/

                            }

}
