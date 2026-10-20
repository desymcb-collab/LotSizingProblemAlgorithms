#ifndef SCF_FL_H_INCLUDED
#define SCF_FL_H_INCLUDED

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
        Models(Instance *instance,double st){ins = instance;N=ins->N;M=ins->M;T=ins->T;start=st;}
        ~Models(){}

        Instance * ins;
        FILE *saida;
        double BestSolution,
               best_bound,
               tempo,start;
        int N,M,T;

        std::vector< std::vector< std::vector < std::vector<double> > > > X;
        std::vector< std::vector< std::vector<double> > > V;
        std::vector< std::vector< std::vector < std::vector< std::vector<int> > > > > Y;
        std::vector< std::vector< std::vector<int> > > Z;
        std::vector< std::vector< std::vector<int> > > L;
        std::vector< std::vector< double > > I;
        std::vector< std::vector< double > > SALE;
        std::vector< std::vector< std::vector<int> > > G;
        std::vector< std::vector< std::vector < std::vector<double> > > > F;
        std::vector< std::vector<double > > NetG;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num4D Fx;
        Num3D Gx;
        Num2D R;
        Num3D X_b;        //variable from Sambasivan's model
        Num3D X_a;        //variable from Sambasivan's model
        Num4D Xx;        //variable from Sambasivan's model
        Num5D Yx;
        Num3D Zx;
        Num3D Lx;
        Num2D Ix;
        Num2D SALEx;

        Exp2D expr1;
        Exp2D expr2;
        Exp2D expr3;
        Exp2D expr4;
        Exp2D expr5;
        Exp3D expr6;
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
        Exp3D expr17;
        Exp3D expr18;
        Exp3D expr19;
        Exp2D expr20;
        Exp3D expr21;
        Exp3D expr22;
        Exp2D expr23;
        Exp3D expr24;
        Exp2D expr25;
        Exp3D expr26;
        Exp3D expr27;

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

        private:

		void Guimaraes_Model(IloModel model);

		void CleanGM();

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();

};

#endif // SCF_FL_H_INCLUDED
