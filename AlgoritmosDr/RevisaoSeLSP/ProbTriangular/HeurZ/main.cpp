#include <iostream>
#include "raiz.h"
#include "rf.h"
#include "rfII.h"
#include "fo.h"
#include "foII.h"
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

    double elapsedTimeRF, elapsedTimeFO, elapsedTimeRF2, elapsedTimeFO2, elapsedTimeRFOKS, elapsedTimeRFOPR;
    double additionalTimeFO, additionalTimeFO2, additionalTimeRFOKS, additionalTimeRFOPR;
    struct timespec startRF,endItRF;
    struct timespec startRF2,endItRF2;
    struct timespec startFO,endItFO;
    struct timespec startFO2,endItFO2;
    struct timespec startRFOKS,endItRFOKS;
    struct timespec startRFOPR,endItRFOPR;

    printf("\n%s\n",instance_name);
    try
    {
        /******************************************************************************************/
        /*******************************Leitura dos dados da instancia*****************************/
        /******************************************************************************************/

        sprintf(ins.fileinfo,"../../Instancias/MultiMachineCapVar/%s.dat",instance_name);      //path to the instance

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
        
        if(elapsedTimeRF < 600){additionalTimeFO = 600 - elapsedTimeRF;}
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
                                        Relax-and-fix heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRF2);
        elapsedTimeRF2 = startRF2.tv_sec;
       
        RFII * rf2 = new RFII(&ins,elapsedTimeRF2);
        rf2->RelaxAndFix();
        ins.LowerBound = rf2->LowerBound;
        ins.solutionRF2 = rf2->solutionRF;
        
        clock_gettime(CLOCK_MONOTONIC, &endItRF2);
        elapsedTimeRF2 = (endItRF2.tv_sec - startRF2.tv_sec);
        elapsedTimeRF2 += (endItRF2.tv_nsec - startRF2.tv_nsec) / 1000000000.0;
        delete rf2;
        //ins.show_solution();
        /******************************************************************************************
                                        Fix-and-optimize heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startFO2);
        elapsedTimeFO2 = startFO2.tv_sec;
        
        if(elapsedTimeRF2 < 600){additionalTimeFO2 = 600 - elapsedTimeRF2;}
        else{additionalTimeFO2 = 0;}

        FOII * fo2 = new FOII(&ins,elapsedTimeFO2,additionalTimeFO2);
        fo2->FixAndOptimize();
        ins.solutionRFO2 = fo2->solutionFO;

        clock_gettime(CLOCK_MONOTONIC, &endItFO2);
        elapsedTimeFO2 = (endItFO2.tv_sec - startFO2.tv_sec);
        elapsedTimeFO2 += (endItFO2.tv_nsec - startFO2.tv_nsec) / 1000000000.0;
        elapsedTimeFO2 += elapsedTimeRF2;
        delete fo2;
        /******************************************************************************************
                                Kernel search heuristic step after RFO
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFOKS);
        elapsedTimeRFOKS = startRFOKS.tv_sec;
        
        if(elapsedTimeFO < 1200){additionalTimeRFOKS = 1200 - elapsedTimeFO;}
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
        
        if((elapsedTimeFO+elapsedTimeFO2) < 2400){additionalTimeRFOPR = 2400 - (elapsedTimeFO+elapsedTimeFO2);}
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
        saida = fopen("resultados/p1_b4-3600s.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
            instance_name,
            ins.solutionRF,ins.solutionRFO, elapsedTimeFO,
            ins.solutionRF2,ins.solutionRFO2, elapsedTimeFO2,
            ins.solutionRFOKS, elapsedTimeRFOKS, ins.solutionRFOPR, elapsedTimeRFOPR);
        fclose(saida);

        saida = fopen("resultados/solutions.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
            instance_name, ins.solutionRF,ins.solutionRFO, ins.solutionRF2,ins.solutionRFO2,
            ins.solutionRFOKS, ins.solutionRFOPR);
        fclose(saida);

        saida = fopen("resultados/times.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,elapsedTimeRF,elapsedTimeFO,
            elapsedTimeRF2, elapsedTimeFO2, elapsedTimeRFOKS,elapsedTimeRFOPR);
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
