#include <iostream>
#include "raiz.h"
#include "rf.h"
#include "fo.h"
#include "rfoks.h"
#include "rfopr.h"

void SolveInstances(const char* instance_name,const char* instance_type);

int main(int argc, char ** argv)
{
    cout << "*********----Code Initialized----********* \n\n";
    SolveInstances(argv[1],argv[2]);
    //SolveInstances("Data2-15-10-0.8-0.5-100-60-20-1-mod");

    cout << "\n***********----Code Ended----**********\n";
    return 0;
}

void SolveInstances(const char* instance_name,const char* instance_type)
{
    /****Instance class*****/
    Instance ins;
    FILE *saida;
    bool new_sol;

    double elapsedTimeRF, elapsedTimeFO, elapsedTimeRFKS, elapsedTimeRFPR;
    double additionalTimeFO, additionalTimeRFKS, additionalTimeRFPR;
    struct timespec startRF,endItRF;
    struct timespec startFO,endItFO;
    struct timespec startRFKS,endItRFKS;
    struct timespec startRFPR,endItRFPR;

    printf("\n%s\n",instance_name);
    try
    {
        /******************************************************************************************/
        /*******************************Leitura dos dados da instancia*****************************/
        /******************************************************************************************/

        sprintf(ins.fileinfo,"../../../Instancias/SingleMachineMod/%s.txt",instance_name);      //path to the instance

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
        
        if(elapsedTimeRF < 1200){additionalTimeFO = 2400 - elapsedTimeRF;}
        else{additionalTimeFO = 1200;}

        FO * fo = new FO(&ins,elapsedTimeFO,additionalTimeFO);
        fo->FixAndOptimize();
        ins.solutionRFO = ins.vet[3].sol;

        clock_gettime(CLOCK_MONOTONIC, &endItFO);
        elapsedTimeFO = (endItFO.tv_sec - startFO.tv_sec);
        elapsedTimeFO += (endItFO.tv_nsec - startFO.tv_nsec) / 1000000000.0;
        elapsedTimeFO += elapsedTimeRF;
        delete fo;
        
        /******************************************************************************************
                                Kernel search heuristic step after RFO
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFKS);
        elapsedTimeRFKS = startRFKS.tv_sec;
        
        if(elapsedTimeFO < 2400){additionalTimeRFKS = 2400 - elapsedTimeFO;}
        else{additionalTimeRFKS = 0;}

        RFOKS * rfks = new RFOKS(&ins,elapsedTimeRFKS,additionalTimeRFKS);
        rfks->KernelSearch();
        ins.solutionRFKS = rfks->BestSolution;

        clock_gettime(CLOCK_MONOTONIC, &endItRFKS);
        elapsedTimeRFKS = (endItRFKS.tv_sec - startRFKS.tv_sec);
        elapsedTimeRFKS += (endItRFKS.tv_nsec - startRFKS.tv_nsec) / 1000000000.0;
        elapsedTimeRFKS += elapsedTimeFO;
        delete rfks;

        /******************************************************************************************
                                 RFO + Path relinking heuristic step
        /******************************************************************************************/
        clock_gettime(CLOCK_MONOTONIC, &startRFPR);
        elapsedTimeRFPR = startRFPR.tv_sec;

        if(elapsedTimeFO < 2400){additionalTimeRFPR = 3600 - elapsedTimeFO;}
        else{additionalTimeRFPR = 1200;}

        RFOPR * rfpr = new RFOPR(&ins,elapsedTimeRFPR,additionalTimeRFPR);
        rfpr->PathRelinking();
        new_sol = rfpr->new_sol;
        ins.solutionRFPR = ins.vet[2].sol;
        clock_gettime(CLOCK_MONOTONIC, &endItRFPR);
        elapsedTimeRFPR = (endItRFPR.tv_sec - startRFPR.tv_sec);
        elapsedTimeRFPR += (endItRFPR.tv_nsec - startRFPR.tv_nsec) / 1000000000.0;
        elapsedTimeRFPR += elapsedTimeFO;
        delete rfpr;

        while(new_sol)
        {
            struct timespec start,endIt;
            double elapsedTime;
            new_sol = false;
            if(elapsedTimeRFPR < 3600)
            {
                clock_gettime(CLOCK_MONOTONIC, &start);
                elapsedTime = start.tv_sec;
                additionalTimeFO = 3600 - elapsedTimeRFPR;
                FO * fo = new FO(&ins,elapsedTime,additionalTimeFO);
                fo->FixAndOptimize();
                if(ins.solutionRFPR > ins.vet[3].sol)
                {
                    ins.solutionRFPR = ins.vet[3].sol;
                    new_sol = true;
                }
                delete fo;
            }

            clock_gettime(CLOCK_MONOTONIC, &endItRFPR);
            elapsedTimeRFPR = (endItRFPR.tv_sec - startRFPR.tv_sec);
            elapsedTimeRFPR += (endItRFPR.tv_nsec - startRFPR.tv_nsec) / 1000000000.0;
            elapsedTimeRFPR += elapsedTimeFO;

            if(elapsedTimeRFPR < 3600 && new_sol == true)
            {
                clock_gettime(CLOCK_MONOTONIC, &start);
                elapsedTime = start.tv_sec;
                additionalTimeRFPR = 3600 - elapsedTimeRFPR;
                RFOPR * rfpr = new RFOPR(&ins,elapsedTime,additionalTimeRFPR);
                rfpr->PathRelinking();
                if(ins.solutionRFPR > ins.vet[2].sol)
                {
                    ins.solutionRFPR = ins.vet[2].sol;
                    new_sol = true;
                }
                delete rfpr;
            }

            clock_gettime(CLOCK_MONOTONIC, &endItRFPR);
            elapsedTimeRFPR = (endItRFPR.tv_sec - startRFPR.tv_sec);
            elapsedTimeRFPR += (endItRFPR.tv_nsec - startRFPR.tv_nsec) / 1000000000.0;
            elapsedTimeRFPR += elapsedTimeFO;            
        }

        /******************************************************************************************
                                            Instances outputs
        ******************************************************************************************/
        ins.desaloca_variaveis();
        sprintf(ins.pathname,"resultados/%s.txt",instance_type);
        saida = fopen(ins.pathname,"a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,ins.solutionRF,elapsedTimeRF,ins.solutionRFO,elapsedTimeFO,ins.solutionRFKS,elapsedTimeRFKS,ins.solutionRFPR,elapsedTimeRFPR);
        fclose(saida);

        sprintf(ins.pathname,"resultados/Solutions-%s.txt",instance_type);
        saida = fopen(ins.pathname,"a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,ins.LowerBound,ins.solutionRF,ins.solutionRFO,ins.solutionRFKS,ins.solutionRFPR);
        fclose(saida);

        sprintf(ins.pathname,"resultados/Times-%s.txt",instance_type);
        saida = fopen(ins.pathname,"a");
        if(saida == NULL) printf("Error in saida file!!!\n");
        fprintf(saida,"%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",instance_name,elapsedTimeRF,elapsedTimeFO,elapsedTimeRFKS,elapsedTimeRFPR);
        fclose(saida);
        
        cout << "ElapsedTimeRF: " << elapsedTimeRF << "\n";
        cout << "ElapsedTimeFO: " << elapsedTimeFO << "\n";
        cout << "ElapsedTimeRFOKS: " << elapsedTimeRFKS << "\n";
        cout << "ElapsedTimeRFOPR: " << elapsedTimeRFPR << "\n";
        /******************************************************************************************/
    }
    catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
    }
    catch (...) {
        cerr << "Error Cpp" << endl;
    }
}
