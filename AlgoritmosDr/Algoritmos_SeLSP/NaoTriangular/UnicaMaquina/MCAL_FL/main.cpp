#include <iostream>

#include "../Raiz/raiz.h"
#include "mcal_fl.h"

void SolveInstances();

int main()
{
    cout << "*********----Code Initialized 1----********* \n\n";
    SolveInstances();
    /*SolveInstances((char*)"AAA");
    SolveInstances((char*)"AAB");
    SolveInstances((char*)"ABA");
    SolveInstances((char*)"ABB");
    SolveInstances((char*)"NAA");
    SolveInstances((char*)"NAB");
    SolveInstances((char*)"NBA");
    SolveInstances((char*)"NBB");*/

    cout << "\n\n***********----Code Ended----********* \n\n";
    return 0;
}

void SolveInstances()
{
    /****Instance class*****/
    Instance ins;
    FILE *saida;
    int qt_period,index_items,qt,id,id2;
    int qt_items[4];
    int p = 1;
    int cd[2];
    double cl1[2];
    double cl2 = 0.5;

    qt_items[0] = 15; qt_items[1] = 25; qt_items[2] = 15; qt_items[3] = 400;
    cl1[0] = 0.6; cl1[1] = 0.8;
    cd[0] = 50; cd[1] = 100;

    for(id = 0; id <= 1; id++) // 0 a 1
        for(id2 = 0; id2 <= 1; id2++) // 0 a 1
            for(qt_period = 10; qt_period <= 15; qt_period = qt_period +5) // 5 a 15
                for(index_items = 0; index_items <= 1; index_items ++)  // 0 a 1
                    for(qt = 0; qt <= 9; qt++)  // 0 a 9
                    {
                        printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d.txt\n",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],100,100,qt);

                        /*******************************Leitura dos dados da instancia*********************************************/

                        sprintf(ins.fileinfo,"../../Instancias/instances/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.txt",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],100,100,qt);
                        ins.le_dados_SM(ins.fileinfo);
                        //ins.show_instance();

                        /******************************************************************************************/

                        /******************************Resolvendo com CPLEX*********************************/
                        /*Models * md11 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/SolucaoSamba.txt");
                        md11->saida = fopen(ins.pathname,"a");
                        if(md11->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md11->saida,"T%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md11->saida);
                        md11->RodaSamba();
                        delete md11;*/

                        /*Models * md21 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/Guimaraes_SingleMachine.txt");
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"T%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md21->saida);
                        md21->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        delete md21;*/

                        MCAL * md21 = new MCAL(&ins);
                        sprintf(ins.pathname,"resultados/MCAL_FL.txt");
                        md21->saida = fopen(ins.pathname,"a");
                        if(md21->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md21->saida,"Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d\t",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],100,100,qt);
                        fclose(md21->saida);
                        md21->RunMCAL();
                        ins.desaloca_variaveis();
                        delete md21;
                        //getchar();

                        /*Models * md31 = new Models(&ins);
                        sprintf(ins.pathname,"resultados/SolucaoFluxo.txt");
                        md31->saida = fopen(ins.pathname,"a");
                        if(md31->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md31->saida,"%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md31->saida);
                        md31->RodaFluxo();
                        delete md31;*/

                        /*Models_linear * md1 = new Models_linear(&ins);
                        sprintf(ins.pathname,"resultados/Guimaraes_SingleMachine.txt");
                        md1->saida = fopen(ins.pathname,"a");
                        if(md1->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md1->saida,"T%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md1->saida);
                        md1->RodaGuimaraes();
                        ins.desaloca_variaveis();
                        delete md1;*/

                        /*Models_linear * md2 = new Models_linear(&ins);
                        sprintf(ins.pathname,"resultados/SolucaoFLP_linear.txt");
                        md2->saida = fopen(ins.pathname,"a");
                        if(md2->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md2->saida,"%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md2->saida);
                        md2->RodaFLP();
                        delete md2;*/

                        /*Models_linear * md3 = new Models_linear(&ins);
                        sprintf(ins.pathname,"resultados/SolucaoFluxo.txt");
                        md3->saida = fopen(ins.pathname,"a");
                        if(md3->saida == NULL) printf("Error in saida file!!!\n");
                        fprintf(md3->saida,"%s0%d%d%d_%d\t",type,t,qt_period,qt_items[index_items],qt);
                        fclose(md3->saida);
                        md3->RodaFluxo();
                        delete md3;*/
                        /***********************************************************************************/

                    }
}
