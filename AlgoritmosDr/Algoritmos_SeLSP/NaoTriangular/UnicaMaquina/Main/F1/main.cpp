#include <iostream>

#include "../../Raiz/raiz.h"
#include "../../ModelosOriginais/F1/F1-Bernardo.h"

void SolveInstances(char *type);

int main()
{
    cout << "*********----Code Initialized 1----********* \n\n";
    SolveInstances((char*)"AAA");
    /*SolveInstances((char*)"AAB");
    SolveInstances((char*)"ABA");
    SolveInstances((char*)"ABB");
    SolveInstances((char*)"NAA");
    SolveInstances((char*)"NAB");
    SolveInstances((char*)"NBA");
    SolveInstances((char*)"NBB");*/

    cout << "\n\n***********----Code Ended----********* \n\n";
    return 0;
}

void SolveInstances(char *type)
{
    /****Instance class*****/
    cout << "Classe -- " << type << "\n\n";
    Instance ins;
    FILE *saida;
    int qt_plant,index_items,qt;
    int qt_items[4];
    int t = 12;

    qt_items[0] = 100; qt_items[1] = 200; qt_items[2] = 300; qt_items[3] = 400;

    for(qt_plant =2; qt_plant <= 2; qt_plant = qt_plant +2)
        for(index_items = 0; index_items <= 0; index_items ++)
            for(qt = 5; qt <= 5; qt++)
            {
                cout << "Instancia --- T" << type << "0" << t << qt_plant << qt_items[index_items] << "_" << qt << "\n\n";

                /*******************************Leitura dos dados da instancia*********************************************/

                sprintf(ins.fileinfo,"../../../CLSD-PM-Data-Files-OPLformat_James/SingleMachineCapVar/T%s0%d%d%d_%d.dat",type,t,qt_plant,qt_items[index_items],qt);
                ins.le_dados_SM(ins.fileinfo);
                //ins.show_instance();

                /******************************************************************************************/

                /******************************Resolvendo com CPLEX*********************************/
                /*Models * md11 = new Models(&ins);
                sprintf(ins.pathname,"resultados/SolucaoSamba.txt");
                md11->saida = fopen(ins.pathname,"a");
                if(md11->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md11->saida,"T%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md11->saida);
                md11->RodaSamba();
                delete md11;*/

                /*Models * md21 = new Models(&ins);
                sprintf(ins.pathname,"resultados/Guimaraes_SingleMachine.txt");
                md21->saida = fopen(ins.pathname,"a");
                if(md21->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md21->saida,"T%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md21->saida);
                md21->RodaGuimaraes();
                ins.desaloca_variaveis();
                delete md21;*/

                F1 * md21 = new F1(&ins);
                /*sprintf(ins.pathname,"resultados/Classic_SingleMachine.txt");
                md21->saida = fopen(ins.pathname,"a");
                if(md21->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md21->saida,"T%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md21->saida);*/
                md21->RunF1();
                ins.desaloca_variaveis();
                delete md21;

                /*Models * md31 = new Models(&ins);
                sprintf(ins.pathname,"resultados/SolucaoFluxo.txt");
                md31->saida = fopen(ins.pathname,"a");
                if(md31->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md31->saida,"%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md31->saida);
                md31->RodaFluxo();
                delete md31;*/

                /*Models_linear * md1 = new Models_linear(&ins);
                sprintf(ins.pathname,"resultados/Guimaraes_SingleMachine.txt");
                md1->saida = fopen(ins.pathname,"a");
                if(md1->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md1->saida,"T%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md1->saida);
                md1->RodaGuimaraes();
                ins.desaloca_variaveis();
                delete md1;*/

                /*Models_linear * md2 = new Models_linear(&ins);
                sprintf(ins.pathname,"resultados/SolucaoFLP_linear.txt");
                md2->saida = fopen(ins.pathname,"a");
                if(md2->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md2->saida,"%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md2->saida);
                md2->RodaFLP();
                delete md2;*/

                /*Models_linear * md3 = new Models_linear(&ins);
                sprintf(ins.pathname,"resultados/SolucaoFluxo.txt");
                md3->saida = fopen(ins.pathname,"a");
                if(md3->saida == NULL) printf("Error in saida file!!!\n");
                fprintf(md3->saida,"%s0%d%d%d_%d\t",type,t,qt_plant,qt_items[index_items],qt);
                fclose(md3->saida);
                md3->RodaFluxo();
                delete md3;*/
                /***********************************************************************************/

            }
}
