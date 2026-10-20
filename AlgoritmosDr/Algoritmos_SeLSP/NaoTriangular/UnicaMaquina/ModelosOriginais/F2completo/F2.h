#ifndef F2_H_INCLUDED
#define F2_H_INCLUDED

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

class F2
{
    public:
        F2(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~F2(){}

        Instance * ins;
        FILE *saida;

        double BestSolution,
               best_bound,
               tempo;
        int N,T;

        std::vector< std::vector<double> > X;
        //std::vector< std::vector<double> > V;
        std::vector< std::vector< std::vector<double> > > Y;
        std::vector< std::vector<double> > Z;
        std::vector< std::vector<double> > I;
        std::vector< std::vector<double> > NetG;
        std::vector<double> G;

        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num2D Ix;
        IloNumVarArray Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D Xx;        //variable from Sambasivan's model
        Num2D Vx;
        Num3D Yx;
        Num2D Zx;

        Exp2D expr;
        Exp2D expr1;
        IloExprArray expr2;
        Exp2D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        Exp3D expr6;
        IloExpr expr7;
        IloExprArray expr8;
        IloExprArray expr9;
        Exp2D expr10;
        Exp2D expr11;
        Exp2D expr12;
        IloExpr expr13;

        IloRangeArray constraints;
        IloRangeArray CapacityConstraints;

        //std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/

        void obj_function();

        void print_solution();

        void getValue();

        void alocate_F2();

        void free_F2();

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RunF2();

        //void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        private:

        void SolveF2();

        void set_parameters();

		void F2_Model(IloModel model);

		//void F2_CapacityConstraints(IloModel model);

		bool Algoritmo_Barany(IloModel model);

		void addBaranysConstraints(IloModel model, int conj_S[], int i, int l);

		bool Algoritmo_FordFulkerson(IloModel model);

		void addMinCutConstraints(IloModel model, int conj_S[], int l, int i);

		void NetworkG(int z, int t);

		bool bfs(int st, int last, int parent[]);

        void dfs(int st, bool visitedS[], bool visited[]);

        int dfsY(int st, bool visited[]);

        int stZ(int t);

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // CLASSIC_MODELS_H_INCLUDED
