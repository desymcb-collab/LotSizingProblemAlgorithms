#include <iostream>
#include "raiz.h"
#include "../Base/rf.h"

int main(int argc, char *argv[])
{
    cout << "*********----Code Initialized----********* \n\n";
    //SolveInstances();
    /****Instance class*****/
    Instance ins;
    FILE *saida;

    double elapsedTime;
    struct timespec start,endIt;

    printf("argv=%d\n",argc);
    for(int i=0;i<argc;i++)
        printf("argv=[%d] %s\n",i,argv[i]);
    getchar();

    /******************************************************************************************/
    /*******************************Leitura dos dados da instancia*****************************/
    /******************************************************************************************/

    //sprintf(ins.fileinfo,"../../../../Instancias/instances_TNT/Data1-15-10-0.6-0.5-50-100-100-0-mod.txt");

    clock_gettime(CLOCK_MONOTONIC, &start);
    elapsedTime = start.tv_sec;
    ins.init_program(argc, argv);
    ins.le_dados_SM(ins.fileinfo);
    printf("%s\n",ins.fileinfo);
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
    sprintf(ins.pathname,"Ajuste/resultados/RFO-Ajuste11-02-19.txt");
    rf21->saida = fopen(ins.pathname,"a");
    if(rf21->saida == NULL) printf("Error in saida file!!!\n");
    fprintf(rf21->saida,"%s\t%.4lf\t%.4lf\n",ins.fileinfo,rf21->BestSolution,elapsedTime);
    fclose(rf21->saida);
    delete rf21;
    /******************************************************************************************/

    cout << "\n\n***********----Code Ended----********* \n\n";
    return 0;
}
