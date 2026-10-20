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
        RF(Instance * instance, double st){ins=instance;N=ins->N;M=ins->M;T=ins->T;start=st;}
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
               elapsed_time,pt;;

        std::vector< std::vector< std::vector < std::vector<double> > > > X;
        std::vector< std::vector< std::vector<double> > > V;
        std::vector< std::vector< std::vector < std::vector<double> > > > Y;
        std::vector< std::vector< std::vector < std::vector<double> > > > reducedCostY;
        std::vector< std::vector< std::vector<double> > > Z;
        std::vector< std::vector< std::vector<int> > > G;
        std::vector< std::vector< std::vector < std::vector<double> > > > F;
        std::vector< std::vector<double > > NetG;
        std::vector< int > firstItem;


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

        Exp2D expr1;
        Exp2D expr2;
        Exp4D expr3;
        Exp2D expr4;
        Exp3D expr5;
        Exp2D expr6;
        Exp3D expr7;
        Exp3D expr8;
        Exp4D expr9;
        Exp3D expr10;
        Exp3D expr11;
        Exp3D expr12;
        Exp3D expr13;
        Exp3D expr14;
        Exp3D expr15;
        Exp3D expr16;
        Exp4D expr17;
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

        void KernelSearch();

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

        /************************************************************************
                            Kernel Search - Auxiliary functions
        ************************************************************************/
        int qtd_buckets;
        int busca_binaria_ks(double LS, LRC_ks vetor);
        int busca_elemento(LRC_ks vetor,int t,int j,int i, int k);
        void Insere_ks(LRC_ks *lrc, double campo1, int campo2, int campo3, int campo4, int campo5, int kernel);
        void retira_elemento(LRC_ks *lrc, int campo2, int campo3, int campo4, int campo5);

        /************************************************************************
                          Kernel Search - Local Search
        ************************************************************************/
        void Define_Kernel_Bucket1(LRC_ks *kernel, LRC_ks *bucket, int rt);
        void Define_Kernel_Bucket2(LRC_ks *kernel, LRC_ks *bucket,int rt);
        void Phase1(LRC_ks *kernel,LRC_ks *bucket);
        void Phase2(LRC_ks *kernel,LRC_ks *bucket, int NB);
        void Update_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket);

        void RodaFLP_linear();
        void SolveFLP_linear();
        void FLP_Model_linear(IloModel model);   //inspired in the facility location problem

        bool Prob_Reduzido_kernel(LRC_ks bucket);
        bool Prob_Reduzido(LRC_ks bucket);
        void FLP_Reduzido(IloModel model, LRC_ks bucket);    //inspired in the facility location problem

        bool Prob_Reduzido_bucket(LRC_ks kernel, LRC_ks bucket, int tam, int id);
        bool Prob_bucket(LRC_ks kernel, LRC_ks bucket,int tam,int id);
        void FLP_Reduzido_bucket(IloModel model, LRC_ks kernel, LRC_ks bucket, int tam, int id);    //inspired in the facility location problem


};

#endif // RF_H_INCLUDED
