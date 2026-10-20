#include "F1-AddBarany.h"
ILOSTLBEGIN

void F1::RunF1()
{
/****Variables initialization****/

    double inicio,fim;
    bool re = false;
    best_bound = -INF;
    BestSolution = INF;
    alocate_F1();

/****Solution method (CPLEX) and processing time****/

    inicio = time(NULL);
    SolveF1();
    env.end();
    fim = time(NULL);
    tempo = difftime(fim,inicio);

/********************Final results and memory dealocation************************/

    cout << "\nBest Solution: " << BestSolution << "\tTempo: " << tempo << "\n";
    //obj_function();
    //print_solution();
    free_F1();

/***********************Exit file with cplex results*****************************/

    /*sprintf(ins->pathname,"resultados/F1_SingleMachine.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf\n",BestSolution,tempo);
    fclose(saida);*/

/********************************************************************************/
}

void F1::SolveF1()
{
    bool ro = false, re = false;

    IloModel model(env);
    F1_UncapacitatedModel(model);                                      //Constroi modelo classico
    set_parameters();

    re = cplex.solve();

    ///Get solution of the current problem
    if(re)
    {
        getValue();
        if(!Algoritmo_Barany(model))
        {
            F1_CapacityConstraints(model);
            cplex.exportModel("lp/F1-Relax.lp");
            re = cplex.solve();
            double aux;
            BestSolution = re?cplex.getObjValue():INF;
        }
        model.end();
    }
}

void F1::F1_UncapacitatedModel(IloModel UncapP)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    double menor = INF;
    for(int i=1; i<=N; i++)
        if(ins->ProdTime[i] < menor)
            menor = ins->ProdTime[i];

    Xx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            int cap = 0;
            cap = int(ins->Capacity[t]/menor);
            Xx[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"Xx%d_%d",i,t);
            Xx[i][t].setName(name);
        }
    }

    Vx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Vx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            int cap = 0;
            cap = int(ins->Capacity[t]/menor);
            Vx[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"Vx%d_%d",i,t);
            Vx[i][t].setName(name);
        }
    }

    Yx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Yx[i][j][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Yx%d_%d_%d",i,j,t);
                Yx[i][j][t].setName(name);
            }
        }
    }

    Zx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Zx[i] = IloNumVarArray(env,T+2);
        for(int t=0; t<=T+1; t++)
        {
            Zx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Zx%d_%d",i,t);
            Zx[i][t].setName(name);
        }
    }

    Ix = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Ix[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            Ix[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"Ix%d_%d",i,t);
            Ix[i][t].setName(name);
        }
    }

    ///Insert objective function in the model
    IloExpr obj1(env);

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
                obj1 += ins->SetupCost[i][j]*Yx[i][j][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            obj1 += ins->HoldingCost[i]*Ix[i][t];

    obj = IloMinimize(env, obj1);
    UncapP.add(obj);
    obj1.end();

    ///Insert constraints in the model

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr1[i][t] += Xx[i][t];
            expr1[i][t] += Ix[i][t-1];
            expr1[i][t] -= Ix[i][t];
            constraints.add(expr1[i][t] == ins->Demand[i][t]);
        }
    expr1.end();

    expr3 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr3[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            double cap;
            //cap = int((ins->Capacity[t])/menor);
            //if(cap > ins->MaxProd[i][t])
            cap = ins->MaxProd[i][t];
            expr3[i][t] += Xx[i][t];
            for(int j=1; j<=N; j++)
                expr3[i][t] -= cap*Yx[j][i][t];
            expr3[i][t] -= cap*Zx[i][t];
            constraints.add(expr3[i][t] <= 0);
        }
    expr3.end();

    expr4 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr4[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            expr4[t] += Zx[i][t];
        constraints.add(expr4[t] == 1);
    }
    expr4.end();

    expr5 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr5[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr5[i][t] += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr5[i][t] += Yx[j][i][t];
                expr5[i][t] -= Yx[i][j][t];
            }
            expr5[i][t] -= Zx[i][t+1];
            constraints.add(expr5[i][t] == 0);
        }
    expr5.end();

    expr6 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr6[i] = Exp2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr6[i][j] = IloExprArray(env,N+1);
            for(int t=0; t<=T; t++)
                expr6[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            if(i!=j)
                for(int t=1; t<=T; t++)
                {
                    expr6[i][j][t] += Vx[i][t];
                    expr6[i][j][t] += N*Yx[i][j][t];
                    expr6[i][j][t] -= N*Zx[i][t];
                    expr6[i][j][t] -= Vx[j][t];
                    constraints.add(expr6[i][j][t] <= (N-1));
                }
    expr6.end();

    UncapP.add(constraints);
    constraints.end();

    Zx[1][1].setBounds(1,1);
    for(int i=1; i<=N; i++)
    {
        Ix[i][0].setBounds(0,0);
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
    }

    cplex = IloCplex(UncapP);
//    UPcplex.exportModel("lp/SM_NTS_F1.lp");
}

void F1::F1_CapacityConstraints(IloModel model)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    expr2 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr2[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int k=1; k<=N; k++)
        {
            expr2[t] += ins->ProdTime[k]*Xx[k][t];
            for(int l=1; l<=N; l++)
                expr2[t] += ins->SetupTime[l][k]*Yx[l][k][t];
        }
        constraints.add(expr2[t] <= ins->Capacity[t]);
    }

    model.add(constraints);
    constraints.end();
//    cplex.exportModel("lp/SM_NTS_F1.lp");
}

bool F1::Algoritmo_Barany(IloModel model)
{
    bool ro = true;
    int conj_S[N+1][T+1][T+1];

    for(int i=1; i<=N; i++)
        for(int l=1; l<=T; l++)
            for(int t=1; t<=l; t++)
            {
                double soma=0.0,yz=0.0;
                for(int u=t; u<=l; u++)
                    soma += ins->Demand[i][u];

                for(int j=1; j<=N; j++)
                    yz += Y[j][i][t];
                yz += Z[i][t];
                soma = soma*yz;

                if(X[i][t] <= soma)
                    conj_S[i][l][t] = 1;
                else
                    conj_S[i][l][t] = 0;
            }


    ///Verify if constraint is satisfied
    for(int i=1; i<=N; i++)
        for(int l=1; l<=T; l++)
        {
            double somaT=0.0,Tdemand=0.0;
            double somaX=0.0,soma,yz;
            for(int t=1; t<=l; t++)
            {
                Tdemand += ins->Demand[i][t];
                if(conj_S[i][l][t] == 0)
                {
                    soma=0.0;
                    yz=0.0;
                    for(int u=t; u<=l; u++)
                        soma += ins->Demand[i][u];
                    for(int j=1; j<=N; j++)
                        yz += Y[j][i][t];
                    yz += Z[i][t];
                    somaT = soma*yz;
                }

            }
            soma = somaT - Tdemand;
            if(soma < 0.0 && somaT > zero)
            {
                ro = false;
                printf("Inequacao invalida. (metodo do Barani)\n");
                ///Insert new constraint to the model
                addBaranysConstraints(model,conj_S[i][l],i,l);
            }
        }

    return ro;
}

void F1::addBaranysConstraints(IloModel model, int conj_S[], int i, int l)
{
    double Tdemand = 0.0;
    expr7 = IloExpr(env);
    expr8 = IloExprArray(env,l+1);
        for(int t=0;t<=l;t++)
            expr8[t] = IloExpr(env);

    for(int t=1; t<=l; t++)
    {
        Tdemand += ins->Demand[i][t];
        if(conj_S[t] == 0)
        {
            double soma = 0.0;
            for(int u=t; u<=l; u++)
                soma += ins->Demand[i][u];
            //printf("%.2lf ( ",soma);
            for(int j=1; j<=N; j++)
            {
                expr8[t] += Yx[j][i][t];
                //printf("+ Y[%d][%d][%d] ",j,i,t);
            }
            expr8[t] += Zx[i][t];
            //printf("+ Z[%d][%d]) ",i,t);
            expr7 += soma*expr8[t];
        }else
            expr7 += Xx[i][t];
    }
    //printf(">= %.2lf\n",Tdemand);
    model.add(expr7 >= Tdemand);
    expr7.end();
    expr8.end();
}

void F1::getValue()
{
    for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                X[i][t] = cplex.getValue(Xx[i][t]);
                V[i][t] = cplex.getValue(Vx[i][t]);
                Z[i][t] = cplex.getValue(Zx[i][t]);
                I[i][t] = cplex.getValue(Ix[i][t]);
                for(int j=1; j<=N; j++)
                    Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
            }
        }
        for(int i=1; i<=N; i++)
            Z[i][T+1] = cplex.getValue(Zx[i][T+1]);
}

void F1::obj_function()
{
    double fo;
    fo = 0;

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                fo += ins->SetupCost[i][j]*Y[i][j][t];
                //printf("%.6lf %.2lf %.2lf \n",fo,ins->SetupCost[i][j],Y[i][j][t]);
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            fo += ins->HoldingCost[i]*I[i][t];
            //printf("%.6lf %.2lf %.2lf \n",fo,ins->HoldingCost[i],I[i][t]);
        }

    printf("Funcao Objetivo: %.2lf\n",fo);
}

void F1::print_solution()
{
    for(int t=1; t<=T; t++)
    {
        for(int i=0; i<=N; i++)
        {
            printf("\n(i,t)%d %d\n",i,t);
            printf("(X,Z,I)%.2lf %d %.2lf\n",X[i][t],Z[i][t],I[i][t]);
            for(int j=1; j<=N; j++)
                printf("\n(i,j,t) %d %d %d - Y: %d",i,j,t,Y[i][j][t]);
            printf("\n");
        }
    }
}

void F1::alocate_F1()
{
    X.resize(N+1);
    V.resize(N+1);
    Z.resize(N+1);
    I.resize(N+1);
    Y.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        X[i].resize(T+1,0);
        V[i].resize(T+1,0);
        Z[i].resize(T+2,0);
        I[i].resize(T+1,0);
        Y[i].resize(N+1);
        for(int j=0;j<=N;j++)
            Y[i][j].resize(T+1,0);
    }
}

void F1::free_F1()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
            Y[i][j].clear();
        V[i].clear();
        X[i].clear();
        Z[i].clear();
        I[i].clear();
        Y[i].clear();
    }
    V.clear();
    X.clear();
    Y.clear();
    Z.clear();
    I.clear();
}

void F1::set_parameters()
{
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.setParam(IloCplex::EpGap, 0.000001);
    cplex.setParam(IloCplex::TreLim, 4000.0);
    //cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::Threads,4);
    //cplex.setParam(IloCplex::WorkMem, 128.0);
    //cplex.setParam(IloCplex::NumericalEmphasis, 1);
    //cplex.setParam(IloCplex::EpRHS,1e-9);
    //cplex.setParam(IloCplex::EpOpt,1e-9);
    //cplex.setParam(IloCplex::EpInt,1e-9);
    //cplex.setParam(IloCplex::EpMrk,0.99999);
    //cplex.setParam(IloCplex::MIPEmphasis, 4);
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
    //cplex.setParam(IloCplex::NodeAlg, IloCplex::Dual);
    //cplex.setParam(IloCplex::IntSolLim, 1);
    //cplex.setParam(IloCplex::NodeLim, 1);
    //cplex.setParam(IloCplex::NodeSel, 1);
    //cplex.setParam(IloCplex::BBInterval, 1);
    //cplex.setParam(IloCplex::Reduce, 0); // 0 No primal and dual reductions, 1 Only primal reductions, 2 Only dual reductions, 3 Both primal and dual reductions  (default)
    //cplex.setParam(IloCplex::PreLinear, 1); // 0 Perform only linear reductions, 1 Perform full reductions (default)
    //cplex.setParam(IloCplex::RepeatPresolve, 0); // -1 Automatic: let CPLEX choose (default), 0 Turn off represolve, 1 Represolve without cuts, 2 Represolve with cuts, 3 Represolve with cuts and allow new root cuts
    //cplex.setParam(IloCplex::Probe, -1); // -1 No probing 0 Automatic: let CPLEX choose (default), 1 Moderate probing level, 2 Aggressive probing level, 3 Very aggressive probing level
    //cplex.setParam(IloCplex::ProbeTime, 10);
    //cplex.setParam(IloCplex::BndStrenInd, 0); // -1 Automatic: let CPLEX choose; default, 0 Do not apply bound strengthening, 1 Apply bound strengthening
    //cplex.setParam(IloCplex::ParallelMode, -1); // -1 opportunistic, 0 defaul let cplex decide, 1 deterministic
    //cplex.setParam(IloCplex::FPHeur,       -1); // -1 Do not apply the feasibility pump heuristic, 0 Automatic: let CPLEX choose; default, 1 Apply the feasibility pump heuristic with an emphasis on finding a feasible solution, 2 Apply the feasibility pump heuristic with an emphasis on finding a feasible solution with a good objective value
    //cplex.setParam(IloCplex::HeurFreq,     -1); //-1 off, 0 default let cplex decide, any positive invokes the heuristic with this frequency
    //cplex.setParam(IloCplex::StrongItLim,   0); // 0 default let cplex decide, any positive limits simplex iterations for strong branching
    //cplex.setParam(IloCplex::PreInd, 0); // 0 off, 1 apply presolve (default)
    //cplex.setParam(IloCplex::RelaxPreInd, 1); // -1 let cplex decide (default), 0 off, 1 apply presolve
    //cplex.setParam(IloCplex::PreslvNd, -1); // -1 No node presolve, 0 Automatic: let CPLEX choose; default, 1 Force presolve at nodes, 2 Perform probing on integer-infeasible variables, 3 Perform aggressive node probing
    //cplex.setParam(IloCplex::AggInd, -1); // -1 default automatic (1 for LP, infinite for MIP), 0 off, positive: max nbr of aggregations
    //cuts: -1 off, 0 default let cplex choose, 1 generate, 2 generate aggressively
    //cplex.setParam(IloCplex::Cliques,      0);
    //cplex.setParam(IloCplex::Covers,       0);
    //cplex.setParam(IloCplex::DisjCuts,     0);
    //cplex.setParam(IloCplex::FlowCovers,   0);
    //cplex.setParam(IloCplex::FlowPaths,    0);
    //cplex.setParam(IloCplex::FracCuts,     0);
    //cplex.setParam(IloCplex::GUBCovers,    0);
    //cplex.setParam(IloCplex::ImplBd,       0);
    //cplex.setParam(IloCplex::ZeroHalfCuts, 0);
    //cplex.setParam(IloCplex::MIRCuts,      0);
    //cplex.setParam(IloCplex::LiftProjCuts, 0);
    //cplex.setParam(IloCplex::MCFCuts,      0);
}
