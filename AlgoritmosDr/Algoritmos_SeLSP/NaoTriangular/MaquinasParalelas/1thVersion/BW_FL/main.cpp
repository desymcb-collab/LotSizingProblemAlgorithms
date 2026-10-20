#include <iostream>
#include "raiz.h"
#include "bw_fl.h"

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

    for(id5 = 0; id5 <= 0; id5++) // 0 a 1 p
        for(id = 1; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 0; id2 <= 0; id2++) // 0 a 1 cd
                for(id3 = 1; id3 <= 1; id3++) // 0 a 1 prob1
                    for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 5; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }

    for(id5 = 0; id5 <= 0; id5++) // 0 a 1 p
        for(id = 1; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
                for(id3 = 0; id3 <= 0; id3++) // 0 a 1 prob1
                    for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }

    for(id5 = 0; id5 <= 0; id5++) // 0 a 1 p
        for(id = 1; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
                for(id3 = 1; id3 <= 1; id3++) // 0 a 1 prob1
                    for(id4 = 0; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }

    for(id5 = 0; id5 <= 0; id5++) // 0 a 1 p
        for(id = 1; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
                for(id3 = 1; id3 <= 1; id3++) // 0 a 1 prob1
                    for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 1; index_items <= 1; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }

    for(id5 = 1; id5 <= 1; id5++) // 0 a 1 p
        for(id = 1; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 0; id2 <= 0; id2++) // 0 a 1 cd
                for(id3 = 1; id3 <= 1; id3++) // 0 a 1 prob1
                    for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 5; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }

    for(id5 = 1; id5 <= 1; id5++) // 0 a 1 p
        for(id = 0; id <= 1; id++) // 0 a 1 cl1
            for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
                for(id3 = 0; id3 <= 0; id3++) // 0 a 1 prob1
                    for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                        for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                            for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                for(qt = 1; qt <= 10; qt++)  // 0 a 9
                        {
                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                        clock_gettime(CLOCK_MONOTONIC, &start);
                        elapsedTime = start.tv_sec;
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();
                        //getchar();
                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/


                        Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/BW_FL1.txt");
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        clock_gettime(CLOCK_MONOTONIC, &endIt);
                        elapsedTime = (endIt.tv_sec - start.tv_sec);
                        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

                        ///Exit file with cplex results
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t%.4lf %.4lf %.4lf\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt,md21->BestSolution,md21->best_bound,elapsedTime);
                        fclose(md21->saida);
                        delete md21;

                        /***********************************************************************************/

                    }
}
