#include <iostream>
#include "raiz.h"
#include "scf_fl.h"

void SolveInstances(const char* instance_name);

int main(int argc, char ** argv)
{
    cout << "*********----Code Initialized----********* \n\n";
    //SolveInstances(argv[1]);
    SolveInstances("Data3-15-10-0.8-0.5-100-60-20-1-mod");

    cout << "\n***********----Code Ended----**********\n";
    return 0;
}

void SolveInstances(const char* instance_name)
{
    /****Instance class*****/
    Instance ins;
    FILE *saida;
    double elapsedTime;
    struct timespec start,endIt;
    
    printf("\n%s\n",instance_name);
    try
    {

        /*******************************Leitura dos dados da instancia*********************************************/

        sprintf(ins.fileinfo,"../Instancias/MultiMachineCapVar/%s.dat",instance_name);

        clock_gettime(CLOCK_MONOTONIC, &start);
        elapsedTime = start.tv_sec;
        ins.le_dados_SM(ins.fileinfo);
        //ins.show_instance();
        //getchar();
        /******************************************************************************************/

        /******************************Resolvendo com CPLEX*********************************/

        Models * md21 = new Models(&ins,elapsedTime);
        md21->RodaGuimaraes();
        ins.desaloca_variaveis();
        clock_gettime(CLOCK_MONOTONIC, &endIt);
        elapsedTime = (endIt.tv_sec - start.tv_sec);
        elapsedTime += (endIt.tv_nsec - start.tv_nsec) / 1000000000.0;

        ///Exit file with cplex results
        saida = fopen("resultados/cplex.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\n",instance_name,md21->best_bound,md21->BestSolution,elapsedTime);
        fclose(saida);
        delete md21;

        /***********************************************************************************/
    }              
    catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
    }
    catch (...) {
        cerr << "Error Cpp" << endl;
    }                  
}
