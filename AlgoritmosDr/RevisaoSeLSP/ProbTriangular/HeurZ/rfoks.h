#ifndef RFOKS_H_INCLUDED
#define RFOKS_H_INCLUDED

#include <iostream>
#include "raiz.h"
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloexpression.h>
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
#include <map>
#include <utility>
#include <sys/time.h>

typedef IloArray<IloNumVarArray> Num2D;
typedef IloArray<Num2D> Num3D;
typedef IloArray<Num3D> Num4D;
typedef IloArray<Num4D> Num5D;

typedef IloArray<IloExprArray> Exp2D;
typedef IloArray<Exp2D> Exp3D;
typedef IloArray<Exp3D> Exp4D;
typedef IloArray<Exp4D> Exp5D;

class RFOKS
{
    public:
        RFOKS(Instance * instance, double st, double addT){ins=instance;N=ins->N;M=ins->M;T=ins->T;start=st;LowerBound=ins->LowerBound;additionalTimeKS=addT;}
        ~RFOKS(){}

        Instance * ins;
        FILE * saida;
        int N,M,T,nm,binary,bk;
        double BestSolution,
               LowerBound,
               tempo,start,
               solutionKS,
               elapsedTime,
               total_time,
               subKS_time,
               additionalTimeKS;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num4D Fx;
        Num3D Gx;
        Num3D X_ini;
        Num3D X_co;
        Num4D Xx;
        Num4D Yx;
        Num3D Zx;
        Num2D R;

        IloRangeArray constraints;

        void KernelSearch();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void getFOptValue();

        void CleanModel();

        void SetParameters();

        void store_solution();

        /************************************************************************
                            Kernel Search - Auxiliary functions
        ************************************************************************/
        int qtd_buckets, qtd_var_bucket;
        int busca_binaria_ks(double LS, LRC_ks vetor);
        int busca_elemento(LRC_ks vetor,int t,int j,int i, int k);
        void Insere_ks(LRC_ks *lrc, double campo1, int campo2, int campo3, int campo4, int campo5, int kernel);
        void retira_elemento(LRC_ks *lrc, int campo2, int campo3, int campo4, int campo5);

         /************************************************************************
                          Kernel Search - Local Search
        ************************************************************************/
        void Update_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket);
        void Define_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket);

        void Phase1(LRC_ks *kernel,LRC_ks *bucket);
        bool Prob_Reduzido_kernel(LRC_ks bucket);
        bool Prob_Reduzido(LRC_ks bucket);
        void FLP_Reduzido(IloModel model, LRC_ks bucket);    //inspired in the facility location problem

        void Phase2(LRC_ks *kernel,LRC_ks *bucket);        
        bool Prob_Reduzido_bucket(LRC_ks kernel, LRC_ks bucket, int id);
        bool Prob_bucket(LRC_ks kernel, LRC_ks bucket,int id);
        void FLP_Reduzido_bucket(IloModel model, LRC_ks kernel, LRC_ks bucket, int id);    //inspired in the facility location problem

};

#endif // RFOKS_H_INCLUDED
