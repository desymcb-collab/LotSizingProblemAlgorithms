#ifndef F1_RELAX_H_INCLUDED
#define F1_RELAX_H_INCLUDED

#include <iostream>
#include "../../Raiz/raiz.h"
#include <ilcplex/ilocplex.h>
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

class F1
{
    public:
        F1(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~F1(){}

        Instance * ins;
        FILE *saida;

        double BestSolution,
               best_bound,
               tempo;
        int N,T;

        std::vector< std::vector<double> > X;
        std::vector< std::vector<double> > V;
        std::vector< std::vector< std::vector<double> > > Y;
        std::vector< std::vector<double> > Z;
        std::vector< std::vector<double> > I;

        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloEnv env2;
        IloObjective obj;
        IloCplex cplex;
        IloCplex UPcplex;
        Num2D Ix;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D Xx;        //variable from Sambasivan's model
        Num2D Vx;
        Num3D Yx;
        Num2D Zx;
        Num2D UPI;        //variable from Sambasivan's model
        Num2D UPX;        //variable from Sambasivan's model
        Num2D UPV;
        Num3D UPY;
        Num2D UPZ;

        Exp2D expr1;
        IloExprArray expr2;
        Exp2D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        Exp3D expr6;
        IloExpr expr7;
        IloExprArray expr8;

        IloRangeArray constraints;

        //std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/

        void obj_function();

        void print_solution();

        void alocate_F1();

        void free_F1();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RunF1();

        //void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        private:

        void SolveF1();

        void set_parameters();

		void F1_Model(IloModel model);

		void F1_UncapacitatedModel(IloModel UncapP);

		bool Algoritmo_Barany(IloModel model);

		void insereRestricoesBarany(IloModel model, int conj_S[], int i, int l);

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // CLASSIC_MODELS_H_INCLUDED
