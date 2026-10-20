#ifndef FO_H_INCLUDED
#define FO_H_INCLUDED

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

class FO
{
    public:
        FO(Instance * instance, double st, double addFO){ins=instance;N=ins->N;M=ins->M;T=ins->T;start=st;additionalTimeFO=addFO;}
        ~FO(){}

        Instance * ins;
        FILE * saida;
        
        int N,M,T,binary,jump;
        double additionalTimeFO;
        double tempo,start,
               solutionFO,
               elapsedTime,
               subFO_time;

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

        void FixAndOptimize();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/

        bool SolveFO(bool *re, int t, int id);

        void FOModel(IloModel model, int l, int id);

        void getFOptValue();

        void CleanModel();

        void SetParameters();

        void store_solution();

        bool verify_duplication();

};
#endif // FO_H_INCLUDED