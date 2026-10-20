#ifndef MCAL_FL_H_INCLUDED
#define MCAL_FL_H_INCLUDED

#include <iostream>
#include "../Raiz/raiz.h"
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

class MCAL
{
    public:
        MCAL(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~MCAL(){}

        Instance * ins;
        FILE *saida;

        double BestSolution,
               best_bound,
               tempo,inicio;
        int N,T;

        std::vector< std::vector< std::vector<double> > > X;
        std::vector< std::vector<double> > V;
        std::vector< std::vector< std::vector<int> > > Y;
        std::vector< std::vector< std::vector<int> > > Q;
        std::vector< std::vector<int> > Z;
        std::vector< std::vector<int> > G;
        std::vector< std::vector<double> > NetG;

        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        IloCplex cplexRelaxed;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D X_ini;        //variable from Sambasivan's model
        Num2D X_co;        //variable from Sambasivan's model
        Num3D Xx;        //variable from Sambasivan's model
        Num2D Vx;
        Num3D Qx;
        Num3D Yx;
        Num2D Zx;
        Num2D Gx;
        IloNumVarArray R;

        Exp2D expr1;
        IloExprArray expr2;
        Exp3D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        Exp3D expr6;
        Exp3D expr7;
        Exp3D expr8;
        Exp2D expr9;
        Exp2D expr10;
        Exp2D expr11;
        Exp2D expr12;
        Exp2D expr13;
        Exp3D expr14;
        IloExpr expr15;

        IloRangeArray constraints;

        //std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/

        void obj_function();

        void print_solution();

        void printErro(int t);

        void alocate_MCAL();

        void free_MCAL();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RunMCAL();

        void getValue();

        bool Algoritmo_FordFulkerson(IloModel model);

		void addMinCutConstraints(IloModel model, int conj_S[], int i, double mincut);

		void NetworkG( int t);

		bool bfs(int st, int last, int parent[]);

        void dfs(int st, bool visitedS[], bool visited[]);

        private:

        void SolveMCAL();

        void set_parameters();

		void MCAL_Model(IloModel model);

};

#endif // MCAL_FL_H_INCLUDED
