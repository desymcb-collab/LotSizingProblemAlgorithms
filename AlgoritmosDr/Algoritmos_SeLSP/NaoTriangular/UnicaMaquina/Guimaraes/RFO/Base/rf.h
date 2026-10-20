#ifndef RF_H_INCLUDED
#define RF_H_INCLUDED

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

#define MAX_N 400
#define MAX_T 12
#define MAX_M 2


typedef IloArray<IloNumVarArray> Num2D;
typedef IloArray<Num2D> Num3D;
typedef IloArray<Num3D> Num4D;
typedef IloArray<Num4D> Num5D;

typedef IloArray<IloExprArray> Exp2D;
typedef IloArray<Exp2D> Exp3D;
typedef IloArray<Exp3D> Exp4D;
typedef IloArray<Exp4D> Exp5D;

/************************************************************************
                    Estutura de dados kernel search
************************************************************************/

struct dado_ks {
    int i1;
    int t1; //indica qual caminho que foi o escolhido (periodo)
    int k1; //indica qual caminho que foi o escolhido (maquina)
    int j1; //indica qual caminho que foi o escolhido (item)
    double cost; //fornece o valor deste caminho
    int kernel;
};

struct LRC_ks {
    int nelementos; //indica o número de elementos da lrc
    dado_ks campo[MAX_N * MAX_T * MAX_M]; //dados da estrutura
};

class RF
{
    public:
        RF(Instance * instance, double st){ins=instance;N=ins->N;M=ins->M;T=ins->T;start=st;jump=ins->jump;jump_fo=ins->jump_fo;binary=ins->binary;}
        ~RF(){}

        Instance * ins;
        FILE * saida;
        guarda vet[vetKS];
        int N,M,T,nm,binary,bk;
        int jump,jump_fo;
        double BestSolution,
               best_bound,
               tempo,start,
               solutionRF,
               elapsed_time,pt;

        std::vector< std::vector< std::vector < double > > > X;
        std::vector< std::vector< std::vector<double> > > V;
        std::vector< std::vector< std::vector < double > > > Y;
        std::vector< std::vector< std::vector < double > > > reducedCostY;
        std::vector< std::vector< double > > Z;
        std::vector< std::vector< int > > G;
        std::vector< std::vector< std::vector < double > > > F;
        std::vector< double > NetG;

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

        Exp2D expr1;
        IloExprArray expr2;
        Exp3D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        IloExprArray expr6;
        Exp2D expr7;
        Exp2D expr8;
        Exp3D expr9;
        Exp2D expr10;
        Exp2D expr11;
        Exp3D expr12;
        Exp3D expr13;
        Exp2D expr14;
        Exp2D expr15;
        Exp2D expr16;
        Exp3D expr17;
        IloExpr expr18;
        IloExpr expr_ks_18;

        IloRangeArray constraints;

        std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;


        void RelaxAndFix();

        void SolveRelaxedModel();          //Start solution of the relaxed model

        bool SolveRFModel(int t, int jump, int id);          //Start solution of the relaxed model

        bool FixAndOptimize(bool *re, int t, int id);



        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void AllocateVariables();               //Allcate variables of the problem

        void FreeVariables();                   //Deallocate variables of the problem

        void RelaxedModel(IloModel model);

        void FixedModel(IloModel model, int l, int id);

        void RFModel(IloModel model, int l, int jump, int id);

        void FixedVariables(int t);

        void BinaryVariables(int t);

        void getFirstItems();

        void getRelaxedValue();

        void getFOptValue();

        bool getRFValue(int t, int id);

        void getReducedCost();

        void CleanModel();

        void SetParameters();

        void store_solution();

        bool verify_duplication();

        void RodaFLP_linear();
        void SolveFLP_linear();
        void FLP_Model_linear(IloModel model);   //inspired in the facility location problem

};

#endif // RF_H_INCLUDED
