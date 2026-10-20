#ifndef MODELS_LINEAR_H_INCLUDED
#define MODELS_LINEAR_H_INCLUDED

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

class Models_linear
{
    public:
        Models_linear(Instance *instance){ins = instance;N=ins->N;T=ins->T;}
        ~Models_linear(){}

        Instance * ins;
        FILE *saida;
        double BestSolution,
               best_bound,
               tempo;
        int N,T;

        std::vector< std::vector<double> > X;
        std::vector< std::vector< std::vector<double> > > Y;
        std::vector< std::vector<double> > Z;
        std::vector< std::vector<double> > I;
        std::vector< std::vector<double> > G;
        std::vector< std::vector< std::vector<double> > > F;


        /************************************************************************
                    Variables to build the model on CPLEX solver
        ************************************************************************/

        IloEnv env;
        IloObjective obj;
        IloCplex cplex;
        Num3D Fx;
        Num2D Ix;
        Num2D Gx;
        //Num4D X_FLP;    //variable from facility location problem
        Num2D Xx;        //variable from Sambasivan's model
        Num3D Yx;
        Num2D Zx;

        Exp2D expr1;
        IloExprArray expr2;
        Exp2D expr3;
        IloExprArray expr4;
        Exp2D expr5;
        IloExprArray expr6;
        Exp2D expr7;
        Exp2D expr8;
        Exp3D expr9;
        IloExprArray expr10;
        IloExprArray expr11;
        
        IloRangeArray constraints;

        std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > cost;


        /************************************************************************
                         Functions to build and solve the model
        ************************************************************************/
        void RodaGuimaraes();

        //void RodaFLP();

        //void RodaFluxo();

        bool SolveGuimaraes();

        //void SolveFLP(LS_Solution_linear * s);

        //void SolveFluxo(LS_Solution_linear * s);

        void set_parameters();

        void InitCost();

        //void ShowSolution(LS_Solution_linear * n);

        /************************************************************************
                        Functions to alocate and dealocate memory
        ************************************************************************/
        void aloca_Guimaraes();

        void funcao_objetivo();

        void libera_Guimaraes();

        private:

		void Guimaraes_Model(IloModel model);

		void CleanGM();

		//void FLP_Model(IloModel model);   //inspired in the facility location problem

		//void CleanFLP();

		//void Fluxo_Model(IloModel model);  //inspired in flow problem

		//void CleanFM();       

};

#endif // MODELS_LINEAR_H_INCLUDED
