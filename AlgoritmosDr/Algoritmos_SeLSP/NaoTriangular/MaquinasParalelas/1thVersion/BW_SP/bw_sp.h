#ifndef BW_SP_H_INCLUDED
#define BW_SP_H_INCLUDED

#include <iostream>
#include "raiz.h"
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
        Models(Instance *instance){ins = instance;N=ins->N;M=ins->M;T=ins->T;}
        ~Models(){}

        Instance * ins;
        FILE *saida;
        double BestSolution,
               best_bound,
               tempo,inicio;
        int N,M,T;

        std::vector< std::vector< std::vector< std::vector<double> > > > X;
        std::vector< std::vector< std::vector < std::vector<int> > > > Y;
        std::vector< std::vector< std::vector<int> > > Z;
        std::vector< std::vector< std::vector<int> > > G;
        std::vector< std::vector<double > > NetG;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num3D Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num3D X_ini;        //variable from Sambasivan's model
        Num3D X_co;        //variable from Sambasivan's model
        Num4D Xx;        //variable from Sambasivan's model
        Num4D Yx;
        Num3D Zx;
        Num2D R;

        IloExprArray expr1;
        Exp2D expr2;
        IloExprArray expr3;
        Exp2D expr4;
        Exp3D expr5;
        Exp2D expr6;
        Exp3D expr7;
        Exp2D expr8;
        Exp3D expr9;
        Exp3D expr10;
        Exp4D expr11;
        Exp3D expr12;
        Exp3D expr13;

        Exp3D expr14;
        Exp3D expr15;
        Exp3D expr16;
        Exp4D expr17;
        IloExpr expr18;

        IloRangeArray constraints;

        std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;


        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RodaGuimaraes();

        //void RodaFluxo();

        bool SolveGuimaraes();

        //void SolveFluxo(LS_Solution_linear * s);

        void set_parameters();

        //void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/
        void aloca_Guimaraes();

        void funcao_objetivo();

        void libera_Guimaraes();

        void printSolution();

        void printSeq();

        void getValue();

        bool Algoritmo_FordFulkerson(IloModel model);

		void addMinCutConstraints(IloModel model, int conj_S[], int l, int k, int i, double mincut);

		void NetworkG(int k, int t);

		bool bfs(int st, int last, int parent[]);

        void dfs(int st, bool visitedS[], bool visited[]);


        private:

		void Guimaraes_Model(IloModel model);

		void CleanGM();

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // BW_SP_H_INCLUDED
