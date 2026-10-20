#ifndef MCF_H_INCLUDED
#define MCF_H_INCLUDED

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

class Models
{
    public:
        Models(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~Models(){}

        Instance * ins;
        FILE *saida;
        double BestSolution,
               best_bound,
               tempo;
        int N,T;

        std::vector< std::vector<double> > X;
        std::vector< std::vector< std::vector<int> > > Y;
        std::vector< std::vector<int> > Z;
        std::vector< std::vector<double> > I;
        std::vector< std::vector<int> > G;
        std::vector< std::vector< std::vector< std::vector<int> > > > F;
        std::vector< std::vector<double> > NetG;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num4D Fx;
        Num2D Ix;
        Num2D Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D X_ini;        //variable from Sambasivan's model
        Num2D X_co;        //variable from Sambasivan's model
        Num2D Xx;        //variable from Sambasivan's model
        Num3D Yx;
        Num2D Zx;
        IloNumVarArray R;

        Exp2D expr1;
        IloExprArray expr2;
        Exp2D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        Exp2D expr6;
        Exp3D expr7;
        Exp3D expr8;
        Exp4D expr9;
        Exp2D expr10;
        Exp2D expr11;
        Exp2D expr12;
        Exp2D expr13;

        Exp2D expr14;
        Exp2D expr15;
        Exp2D expr16;
        Exp3D expr17;
        IloExpr expr18;

        IloRangeArray constraints;

        std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;


        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RodaMCF();

        //void RodaFluxo();

        bool SolveMCF();

        //void SolveFluxo(LS_Solution_linear * s);

        void set_parameters();

        //void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/
        void aloca_MCF();

        void funcao_objetivo();

        void libera_MCF();

        void printSolution();

        void printSeq();

        void getValue();

        bool Algoritmo_FordFulkerson(IloModel model);

		void addMinCutConstraints(IloModel model, int conj_S[], int k, int i, double mincut);

		void NetworkG( int t);

		bool bfs(int st, int last, int parent[]);

        void dfs(int st, bool visitedS[], bool visited[]);

        private:

		void MCF_Model(IloModel model);

		void CleanGM();

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // SCF_H_INCLUDED
