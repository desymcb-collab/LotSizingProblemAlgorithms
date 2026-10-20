#ifndef RF_H_INCLUDED
#define RF_H_INCLUDED

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


class RF
{
    public:
        RF(Instance * instance, double st){ins=instance;N=ins->N;M=ins->M;T=ins->T;start=st;}
        ~RF(){}

        Instance * ins;
        FILE * saida;
        
        int N,M,T,binary;
        int jump;
        double LowerBound,
               tempo,start,
               solutionRF,
               elapsedTime,
               subRF_time;

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

        bool RelaxAndFix();

        void SolveRelaxedModel();                  //Start solution of the relaxed model

        bool SolveRFModel(int t, int id);          //Solution of the t-th relax-and-fix model

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/

        void RelaxedModel(IloModel model);

        void RFModel(IloModel model, int l, int id);

        void getRelaxedValue();

        bool getRFSolution(int t, int id);

        void getReducedCost();

        void CleanModel();

        void SetParameters();

        void store_solution();
};

#endif // RF_H_INCLUDED