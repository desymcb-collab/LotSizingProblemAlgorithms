#ifndef RFOPR_H_INCLUDED
#define RFOPR_H_INCLUDED

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

class RFOPR
{
    public:
        RFOPR(Instance * instance, double st, double addT){ins=instance;N=ins->N;T=ins->T;start=st;additionalTimePR=addT;}
        ~RFOPR(){}

        Instance * ins;
        FILE * saida;
        bool new_sol;
        int N,T;
        double BestSolution,
               best_bound,
               tempo,start,
               elapsedTime,
               subPR_time,
               total_time,
               additionalTimePR,
               solutionPR;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num3D Fx;
        Num2D Gx;
        Num2D X_ini;
        Num2D X_co;
        Num3D Xx;
        Num3D Yx;
        Num2D Zx;
        IloNumVarArray R;

        IloRangeArray constraints;

        void PathRelinking();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/

        void PRModel(IloModel model);

        void RFModel(IloModel model, int l, int jump, int id);

        void getFOptValue();

        void CleanModel();

        void SetParameters();

        void store_solution();
};

#endif // RFOPR_H_INCLUDED
