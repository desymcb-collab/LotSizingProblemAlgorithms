#include "models.h"

void Models::RodaGuimaraes()
{
    double inicio,fim;
    bool re = false;
    best_bound = INF;
    BestSolution = INF;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    aloca_Guimaraes();
    inicio = time(NULL);
    env = IloEnv();
    re = SolveGuimaraes();
    env.end();
    fim = time(NULL);
    tempo = difftime(fim,inicio);
    //funcao_objetivo();
    
    for(int t=1; t<=T; t++)
        for(int i=0; i<=N; i++)
        {
            printf("\n%d %d\n",i,t);
            printf("%.2lf %d %d %.2lf\n",X[i][t],Z[i][t],G[i][t],I[i][t]);  

            for(int j=1; j<=N; j++)
                printf("\n%d %d %d - %.2lf %d",i,j,t,F[i][j][t],Y[i][j][t]);
            
            printf("\n");
        }
    libera_Guimaraes();

    if(re)
        cout << "\nBest Solution: " << BestSolution << "\tTempo: " << tempo << "\n";
        //cout << "\nBest Solution: " << BestSolution << "\tBest Bound: " << best_bound << "\tTempo: " << time << "\n";

    ///Exit file with cplex results

    sprintf(ins->pathname,"resultados/Guimaraes_SingleMachine.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf\n",BestSolution,tempo);
    fclose(saida);
}

bool Models::SolveGuimaraes()
{
    bool re = false;
    IloModel model(env);
    Guimaraes_Model(model);                                         //Constroi modelo classico

    set_parameters();

    re = cplex.solve();
    BestSolution = re?cplex.getObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF

    ///Get solution of the current problem
    if(re)
    {
        for(int i=1; i<=N; i++)
            for(int t=1; t<=T; t++)
            {
                X[i][t] = cplex.getValue(Xx[i][t]);
                Z[i][t] = IloRound(cplex.getValue(Zx[i][t]));
                G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
                I[i][t] = cplex.getValue(Ix[i][t]);
                for(int j=1; j<=N; j++)
                {
                    F[i][j][t] = cplex.getValue(Fx[i][j][t]);
                    Y[i][j][t] = IloRound(cplex.getValue(Yx[i][j][t]));
                }
            }
        for(int t=1; t<=T; t++)
            for(int j=1; j<=N; j++)
                F[0][j][t] = cplex.getValue(Fx[0][j][t]);
            
        CleanGM();
        model.end();
    }
    return re;

}

void Models::Guimaraes_Model(IloModel model)
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
            Xx[i][t] = IloNumVar(env,0,cap,ILOFLOAT);
            sprintf(name,"Xx%d_%d",i,t);
            Xx[i][t].setName(name);
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
                Yx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Yx%d_%d_%d",i,j,t);
                Yx[i][j][t].setName(name);
            }
        }
    }

    Zx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Zx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
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

    Fx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Fx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Fx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Fx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"Fx%d_%d_%d",i,j,t);
                Fx[i][j][t].setName(name);
            }
        }
    }

    Gx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Gx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            Gx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Gx%d_%d",i,t);
            Gx[i][t].setName(name);
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

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
    model.add(obj);
    obj1.end();

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            expr1[i][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr1[i][t] += Xx[i][t];
            expr1[i][t] += Ix[i][t-1];
            expr1[i][t] -= Ix[i][t];
            constraints.add(expr1[i][t] == ins->Demand[i][t]);
        }

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

    /*expr3 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr3[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
    {
        
        for(int t=1; t<=T; t++)
        {
            expr3[i][t] += Xx[i][t];
            for(int j=1; j<=N; j++)
                expr3[i][t] -= 25*Yx[j][i][t];
            expr3[i][t] -= 25*Zx[i][t];
            constraints.add(expr3[i][t] >= 0);
        }
    }
    expr3.end();*/

    expr3 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr3[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
    {
        
        for(int t=1; t<=T; t++)
        {
            int cap;
            cap = int((ins->Capacity[t])/menor);
            if(cap > ins->MaxProd[i][t])
                cap = ins->MaxProd[i][t];
            expr3[i][t] += Xx[i][t];
            for(int j=1; j<=N; j++)
                expr3[i][t] -= cap*Yx[j][i][t];
            expr3[i][t] -= cap*Zx[i][t];
            constraints.add(expr3[i][t] <= 0);
        }
    }

    expr4 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr4[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            expr4[t] += Zx[i][t];
        constraints.add(expr4[t] == 1);
    }

    expr5 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr5[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<T; t++)
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
        expr5[i][T] += Zx[i][T];
        for(int j=1; j<=N; j++)
        {
            expr5[i][T] += Yx[j][i][T];
            expr5[i][T] -= Yx[i][j][T];
        }
        constraints.add(expr5[i][T] >= 0);
    }

    expr6 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr6[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int j=1; j<=N; j++)
        {
            expr6[t] += Fx[0][j][t];
            expr6[t] -= Gx[j][t];
        }
        constraints.add(expr6[t] == 0);
    }

    expr7 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr7[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr7[i][t] -= Gx[i][t];
            expr7[i][t] += Fx[0][i][t];
            for(int j=1; j<=N; j++)
            {
                expr7[i][t] += Fx[j][i][t];
                expr7[i][t] -= Fx[i][j][t];
            }
            constraints.add(expr7[i][t] == 0);
        }

    expr8 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr8[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr8[i][t] += Fx[0][i][t];
            expr8[i][t] -= N*Zx[i][t];
            constraints.add(expr8[i][t] <= 0);
        }

    expr9 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr9[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                expr9[i][j][t] += Fx[i][j][t];
                expr9[i][j][t] -= N*Yx[i][j][t];
                constraints.add(expr9[i][j][t] <= 0);
            }


    expr10 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr10[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            expr10[i][t] += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr10[i][t] += Yx[j][i][t];
            }
            expr10[i][t] -= Gx[i][t];
            constraints.add(expr10[i][t] >= 0);
        }
    }

    expr11 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr11[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                expr11[i][t] += Yx[j][i][t];
            }
            expr11[i][t] += Zx[i][t];
            expr11[i][t] -= INF*Gx[i][t];
            constraints.add(expr11[i][t] <= 0);
        }
        
    }

    /*expr12 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr12[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            expr12[t] += Fx[0][i][t];
     
        constraints.add(expr12[t] >= 1);
    }
    */


    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Ix[i][0].setBounds(0,0);
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
        {
            Fx[i][i][t].setBounds(0,0);
            Yx[i][i][t].setBounds(0,0);
        }
    }

    //Zx[1][1].setBounds(1,1);

    cplex = IloCplex(model);
    cplex.exportModel("lp/SM_NTS.lp");
}

void Models::CleanGM()
{
    obj.end();
    constraints.end();

    Fx.end();
    Ix.end();
    Gx.end();
    Xx.end();
    Yx.end();
    Zx.end();

    expr1.end();
    expr2.end();
    expr3.end();
    expr4.end();
    expr5.end();
    expr6.end();
    expr7.end();
    expr8.end();
    expr9.end();
    expr10.end();
    expr11.end();
    expr12.end();
    cplex.end();
}

void Models::funcao_objetivo()
{
    double fo;
    fo = 0;

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                fo += ins->SetupCost[i][j]*double(Y[i][j][t]);
                //printf("%.6lf %.2lf %.2lf \n",fo,ins->SetupCost[i][j],Y[i][j][t]);
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            fo += ins->HoldingCost[i]*double(I[i][t]);
            //printf("%.6lf %.2lf %.2lf \n",fo,ins->HoldingCost[i],I[i][t]);
        }

    printf("Funcao Objetivo: %.2lf\n",fo);

}

void Models::aloca_Guimaraes()
{
    X.resize(N+1);
    Z.resize(N+1);
    I.resize(N+1);
    G.resize(N+1);
    F.resize(N+1);
    Y.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        X[i].resize(T+1,0);
        Z[i].resize(T+1,0);
        I[i].resize(T+1,0);
        G[i].resize(T+1,0);
        Y[i].resize(N+1);
        F[i].resize(N+1);
        for(int j=0;j<=N;j++)
        {
            Y[i][j].resize(T+1,0);
            F[i][j].resize(T+1,0);
        }
        
    }
}

void Models::libera_Guimaraes()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            
            Y[i][j].clear();
            F[i][j].clear();
        }
        X[i].clear();
        Z[i].clear();
        I[i].clear();
        G[i].clear();
        Y[i].clear();
        F[i].clear();
    }
    X.clear();
    Y.clear();
    Z.clear();
    I.clear();
    F.clear();
    G.clear();
}


void Models::set_parameters()
{
    //cplex.setParam(IloCplex::Threads,4);
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.setParam(IloCplex::EpGap, 0.000001);
    cplex.setParam(IloCplex::TreLim, 4000.0);
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


    /*cplex.setParam(IloCplex::Cliques,      0);
    cplex.setParam(IloCplex::Covers,       0);
    cplex.setParam(IloCplex::DisjCuts,     0);
    cplex.setParam(IloCplex::FlowCovers,   0);
    cplex.setParam(IloCplex::FlowPaths,    0);
    cplex.setParam(IloCplex::FracCuts,     0);
    cplex.setParam(IloCplex::GUBCovers,    0);
    cplex.setParam(IloCplex::ImplBd,       0);
    cplex.setParam(IloCplex::ZeroHalfCuts, 0);
    cplex.setParam(IloCplex::MIRCuts,      0);
    cplex.setParam(IloCplex::LiftProjCuts, 0);
    cplex.setParam(IloCplex::MCFCuts,      0);*/
    //cplex.setOut(env.getNullStream());
}

/*void Models::InitCost()
{
    cost.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        cost[i].resize(M+1);
        for(int j=0;j<=M;j++)
        {
            cost[i][j].resize(T+1);
            for(int t=0;t<=T;t++)
            {
                cost[i][j][t].resize(M+1);
                for(int k=0;k<=M;k++)
                    cost[i][j][t][k].resize(T+1,0);
            }
        }
    }

    double aux1,aux2,aux3,menor,pmenor;

    for(int ii=1; ii<=N; ii++)
    {
        for(int jj = 1; jj<=M; jj++)
        {
            for(int tt=1; tt<=T; tt++)
            {
                for(int kk=1; kk<=M; kk++)
                {
                    for(int uu=1; uu<=T; uu++)
                    {
                        if(uu>=tt)
                        {
                            aux1 = (uu-tt)*ins->hpdc[ii][jj] + ins->rpdc[jj][kk];
                            aux2 = (uu-tt)*ins->hpdc[ii][kk] + ins->rpdc[jj][kk];
                            menor = INF;
                            for(int mi=1; mi<=M; mi++)
                            {
                                aux3 = (uu-tt)*ins->hpdc[ii][mi] + ins->rpdc[jj][mi] + ins->rpdc[mi][kk];
                                if(aux3 < menor)
                                {
                                    menor = aux3;
                                    pmenor = mi;
                                }
                            }
                            if(aux1 < aux2 && aux1 < menor)
                                cost[ii][jj][tt][kk][uu] = ins->cpdc[ii][jj][tt] + aux1;
                            else
                            {
                                if(aux2 < aux1 && aux2 < menor)
                                    cost[ii][jj][tt][kk][uu] = ins->cpdc[ii][jj][tt] + aux2;
                                else
                                    cost[ii][jj][tt][kk][uu] = ins->cpdc[ii][jj][tt] + menor;
                            }

                        }
                        else
                            cost[ii][jj][tt][kk][uu] = 0;
                    }
                }
            }
        }
    }
}*/
