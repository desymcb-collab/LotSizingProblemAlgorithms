#include <iostream>
#include "raiz.h"
#include "rf.h"
#include "fo.h"
#include "rfoks.h"
#include "rfopr.h"

void SolveInstances(const char* instance_name);

int main(int argc, char ** argv)
{
    cout << "*********----Code Initialized----********* \n\n";
    SolveInstances(argv[1]);
    //SolveInstances("Data2-15-10-0.8-0.5-100-60-20-1-mod");

    cout << "\n***********----Code Ended----**********\n";
    return 0;
}

void SolveInstances(const char* instance_name)
{
    /****Instance class*****/
    Instance ins;
    FILE *saida;

    double elapsedTimeRF, elapsedTimeFO, elapsedTimeRFOKS, elapsedTimeRFOPR;
    double additionalTimeFO, additionalTimeRFOKS, additionalTimeRFOPR;
    struct timespec startRF,endItRF;
    struct timespec startFO,endItFO;
    struct timespec startRFOKS,endItRFOKS;
    struct timespec startRFOPR,endItRFOPR;

    printf("\n%s\n",instance_name);
    try
    {
        /******************************************************************************************/
        /*******************************Leitura dos dados da instancia*****************************/
        /******************************************************************************************/

        sprintf(ins.fileinfo,"../Instancias/MultiMachineCapVar/%s.dat",instance_name);      //path to the instance

        clock_gettime(CLOCK_MONOTONIC, &startRF);
        elapsedTimeRF = startRF.tv_sec;
        ins.le_dados_SM(ins.fileinfo);
        //ins.show_instance();
        //getchar();
        
        /******************************************************************************************
                                        Relax-and-fix heuristic step
        /******************************************************************************************/

        RF * rf = new RF(&ins,elapsedTimeRF);
        rf->RelaxAndFix();
        ins.LowerBound = rf->LowerBound;
        ins.solutionRF = rf->solutionRF;
        
        clock_gettime(CLOCK_MONOTONIC, &endItRF);
        elapsedTimeRF = (endItRF.tv_sec - startRF.tv_sec);
        elapsedTimeRF += (endItRF.tv_nsec - startRF.tv_nsec) / 1000000000.0;
        delete rf;
        //ins.show_solution();
        /******************************************************************************************
                                        Fix-and-optimize heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startFO);
        elapsedTimeFO = startFO.tv_sec;
        
        if(elapsedTimeRF < 900){additionalTimeFO = 900 - elapsedTimeRF;}
        else{additionalTimeFO = 0;}

        FO * fo = new FO(&ins,elapsedTimeFO,additionalTimeFO);
        fo->FixAndOptimize();
        ins.solutionRFO = fo->solutionFO;

        clock_gettime(CLOCK_MONOTONIC, &endItFO);
        elapsedTimeFO = (endItFO.tv_sec - startFO.tv_sec);
        elapsedTimeFO += (endItFO.tv_nsec - startFO.tv_nsec) / 1000000000.0;
        elapsedTimeFO += elapsedTimeRF;
        delete fo;

        /******************************************************************************************
                                Kernel search heuristic step after RFO
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFOKS);
        elapsedTimeRFOKS = startRFOKS.tv_sec;
        
        if(elapsedTimeFO < 1800){additionalTimeRFOKS = 1800 - elapsedTimeFO;}
        else{additionalTimeRFOKS = 0;}

        RFOKS * rfoks = new RFOKS(&ins,elapsedTimeRFOKS,additionalTimeRFOKS);
        rfoks->KernelSearch();
        ins.solutionRFOKS = rfoks->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItRFOKS);
        elapsedTimeRFOKS = (endItRFOKS.tv_sec - startRFOKS.tv_sec);
        elapsedTimeRFOKS += (endItRFOKS.tv_nsec - startRFOKS.tv_nsec) / 1000000000.0;
        elapsedTimeRFOKS += elapsedTimeFO;
        delete rfoks;

        /******************************************************************************************
                                 RFO + Path relinking heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFOPR);
        elapsedTimeRFOPR = startRFOPR.tv_sec;
        
        if(elapsedTimeFO < 1800){additionalTimeRFOPR = 1800 - elapsedTimeFO;}
        else{additionalTimeRFOPR = 0;}

        RFOPR * rfopr = new RFOPR(&ins,elapsedTimeRFOPR,additionalTimeRFOPR);
        rfopr->PathRelinking();
        ins.solutionRFOPR = rfopr->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItRFOPR);
        elapsedTimeRFOPR = (endItRFOPR.tv_sec - startRFOPR.tv_sec);
        elapsedTimeRFOPR += (endItRFOPR.tv_nsec - startRFOPR.tv_nsec) / 1000000000.0;
        elapsedTimeRFOPR += elapsedTimeFO;
        delete rfopr;

        /******************************************************************************************
                                            Instances outputs
        ******************************************************************************************/
        ins.desaloca_variaveis();
        saida = fopen("solutions/p1_b2-3600s.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,ins.LowerBound,ins.solutionRF,ins.solutionRFO,ins.solutionRFOKS,ins.solutionRFOPR);
        fclose(saida);

        saida = fopen("times/p1_b2-3600s.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,elapsedTimeRF,elapsedTimeFO,elapsedTimeRFOKS,elapsedTimeRFOPR);
        fclose(saida);
        
        cout << "ElapsedTimeRF: " << elapsedTimeRF << "\n" << "ElapsedTimeRFO: " << elapsedTimeFO << "\n";
        cout << "ElapsedTimeRFOKS: " << elapsedTimeRFOKS << "\n";
        cout << "ElapsedTimeRFOPR: " << elapsedTimeRFOPR << "\n";
        /******************************************************************************************/
    }
    catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
    }
    catch (...) {
        cerr << "Error Cpp" << endl;
    }
}
