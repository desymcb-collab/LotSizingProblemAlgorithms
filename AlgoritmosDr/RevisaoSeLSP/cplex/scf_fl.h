#ifndef SCF_FL_H_INCLUDED
#define SCF_FL_H_INCLUDED

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

class Models
{
    public:
        Models(Instance *instance,double st){ins = instance;N=ins->N;M=ins->M;T=ins->T;start=st;}
        ~Models(){}

        Instance * ins;
        FILE *saida;
        double BestSolution,
               best_bound,
               tempo,start;
        int N,M,T;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num4D Fx;
        Num3D Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num3D X_ini;        //variable from Sambasivan's model
        Num3D X_co;        //variable from Sambasivan's model
        Num4D Xx;        //variable from Sambasivan's model
        Num4D Yx;
        Num3D Zx;
        Num2D R;

        IloRangeArray constraints;

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RodaGuimaraes();

        bool SolveGuimaraes();

        void set_parameters();

        /************************************************************************
                                    Auxiliary functions
        ************************************************************************/
        void printSolution();

        void getValue();

        private:

		void Guimaraes_Model(IloModel model);

		void CleanGM();
};

#endif // SCF_FL_H_INCLUDED
