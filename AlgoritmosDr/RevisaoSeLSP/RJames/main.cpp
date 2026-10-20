#include <iostream>
#include "raiz.h"
#include "rf.h"
#include "fo.h"
#include "ks.h"
#include "rfks.h"
#include "kspr.h"

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

    double elapsedTimeRF, elapsedTimeFO, elapsedTimeKS, elapsedTimeRFKS, elapsedTimeKSPR;
    double additionalTimeFO, additionalTimeKS, additionalTimeRFKS, additionalTimeKSPR;
    struct timespec startRF,endItRF;
    struct timespec startFO,endItFO;
    struct timespec startKS,endItKS;
    struct timespec startRFKS,endItRFKS;
    struct timespec startKSPR,endItKSPR;

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
        
        if(elapsedTimeRF < 1800){additionalTimeFO = 1800 - elapsedTimeRF;}
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
                                    independent Kernel search heuristic
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startKS);
        elapsedTimeKS = startKS.tv_sec;

        KS * ks = new KS(&ins,elapsedTimeKS);
        ks->KernelSearch();
        ins.solutionKS = ks->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItKS);
        elapsedTimeKS = (endItKS.tv_sec - startKS.tv_sec);
        elapsedTimeKS += (endItKS.tv_nsec - startKS.tv_nsec) / 1000000000.0;
        delete ks;

        /******************************************************************************************
                                Kernel search heuristic step after RF
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFKS);
        elapsedTimeRFKS = startRFKS.tv_sec;
        
        if(elapsedTimeRF < 1800){additionalTimeRFKS = 1800 - elapsedTimeRF;}
        else{additionalTimeRFKS = 0;}

        RFKS * rfks = new RFKS(&ins,elapsedTimeRFKS,additionalTimeRFKS);
        rfks->KernelSearch();
        ins.solutionRFKS = rfks->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItRFKS);
        elapsedTimeRFKS = (endItRFKS.tv_sec - startRFKS.tv_sec);
        elapsedTimeRFKS += (endItRFKS.tv_nsec - startRFKS.tv_nsec) / 1000000000.0;
        elapsedTimeRFKS += elapsedTimeRF;
        delete rfks;

        /******************************************************************************************
                                 KS + Path relinking heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startKSPR);
        elapsedTimeKSPR = startKSPR.tv_sec;
        
        if(elapsedTimeKS < 1800){additionalTimeKSPR = 1800 - elapsedTimeKS;}
        else{additionalTimeKSPR = 0;}

        KSPR * kspr = new KSPR(&ins,elapsedTimeKSPR,additionalTimeKSPR);
        kspr->PathRelinking();
        ins.solutionKSPR = kspr->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItKSPR);
        elapsedTimeKSPR = (endItKSPR.tv_sec - startKSPR.tv_sec);
        elapsedTimeKSPR += (endItKSPR.tv_nsec - startKSPR.tv_nsec) / 1000000000.0;
        elapsedTimeKSPR += elapsedTimeRF;
        elapsedTimeKSPR += elapsedTimeKS;
        delete kspr;

        /******************************************************************************************
                                            Instances outputs
        ******************************************************************************************/
        ins.desaloca_variaveis();
        saida = fopen("solutions/p1_b4-3600s.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,ins.LowerBound,ins.solutionRF,ins.solutionRFO,ins.solutionKS,ins.solutionRFKS,ins.solutionKSPR);
        fclose(saida);

        saida = fopen("times/p1_b4-3600s.txt","a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,elapsedTimeRF,elapsedTimeFO,elapsedTimeKS,elapsedTimeRFKS,elapsedTimeKSPR);
        fclose(saida);
        
        cout << "ElapsedTimeRF: " << elapsedTimeRF << "\n" << "ElapsedTimeRFO: " << elapsedTimeFO << "\n";
        cout << "ElapsedTimeKS: " << elapsedTimeKS << "\n"<< "ElapsedTimeRFKS: " << elapsedTimeRFKS << "\n";
        cout << "ElapsedTimeKSPR: " << elapsedTimeKSPR << "\n";
        /******************************************************************************************/
    }
    catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
    }
    catch (...) {
        cerr << "Error Cpp" << endl;
    }
}
