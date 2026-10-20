#ifndef CLASSIC_MODELS_H_INCLUDED
#define CLASSIC_MODELS_H_INCLUDED

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

typedef struct item
{
    int itm;
    struct item *I_prox;
}Item;
typedef Item *tipoItem;

typedef struct setI
{
    int n_items;
    Item *IS;
}ISet;
typedef ISet *Set_S;

typedef struct sl
{
    int n_sets;
    ISet *S;
    ISet *NS;
}s_list;
typedef s_list *SetList;


class CM
{
    public:
        CM(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~CM(){}

        Instance * ins;
        FILE *saida;

        double BestSolution,
               best_bound,
               tempo;
        int N,T;

        std::vector< std::vector<double> > X;
        std::vector< std::vector<double> > V;
        std::vector< std::vector< std::vector<int> > > Y;
        std::vector< std::vector<int> > Z;
        std::vector< std::vector<double> > I;
        std::vector<int>  G;
        std::vector< std::vector<int> > cY;
        std::vector< std::vector<int> > WG;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        IloCplex cplexRelaxed;
        Num2D Ix;
        IloNumVarArray Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D Xx;        //variable from Sambasivan's model
        Num2D Vx;
        Num3D Yx;
        Num2D Zx;

        Exp2D expr1;
        IloExprArray expr2;
        Exp2D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        Exp2D expr6;
        Exp2D expr7;
        IloExprArray expr8;
        IloExprArray expr9;
        IloExpr expr10;
        IloExprArray expr11;
        IloExprArray expr12;
        //IloExprArray expr13;
        Exp3D expr13;
        Exp2D expr14;

        IloRangeArray constraints;

        std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/

        void funcao_objetivo();

        void aloca_Classic();

        void libera_Classic();

        SetList libera_conjuntos(SetList conj);

        void libera_conjunto(Set_S conj);

        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RodaClassic();

        //void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        private:

        void copiaY(int t);

        void WeightedGraph(int t);

        void insereItem(Set_S st, int a);

        void insereSet(SetList conj);

        void insereRestricoes1(IloModel model,int conj_S[], int k, int t);

        void insereRestricoes2(IloModel model, int conj_S[], int l, int i);

        SetList inicializaList(SetList conj);

        bool bfs(int st, int last, int parent[]);

        void dfs(int st, bool visited[]);

        int dfsY(int st);

        int stZ(int t);

        bool Algoritmo_FordFulkerson(IloModel model);

        bool Algoritmo_Barany(IloModel model);

        void SolveClassic();

        void set_parameters();

		void Classic_Model(IloModel model);

		void CleanCM();



		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // CLASSIC_MODELS_H_INCLUDED
