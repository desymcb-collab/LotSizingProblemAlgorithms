#include "rfopr.h"

void RFOPR::PathRelinking()
{
    struct timespec fim;
    bool re = true;
    total_time = additionalTimePR;
    new_sol = false;
    env = IloEnv();

    while(re)
    {
        re = false;
        clock_gettime(CLOCK_MONOTONIC, &fim);
        elapsedTime = (fim.tv_sec - start);
        subPR_time = floor((total_time-elapsedTime)/2);

        if((elapsedTime < total_time) && (subPR_time != 0))
        {
            IloModel model(env);
            PRModel(model);                                       //Constroi modelo classico

            SetParameters();

            re = cplex.solve();

            if(re)
            {
                getFOptValue();
                if(BestSolution < (ins->vet[2].sol + 0.01))
                {
                    store_solution();
                }else re = false;
                CleanModel();
                model.end();
                //cout << "PR: " << BestSolution << "\n";
            }

            double aux;
            aux = fabs(ins->vet[2].sol - ins->vet[3].sol);
            if(aux < zero)
            {
                re = false;
                //printf("Equal solutions\n");
            }
        }
    }

    BestSolution = ins->vet[2].sol;

    env.end();
    cout << "RFPR: " << BestSolution << "\n";
}

void RFOPR::PRModel(IloModel model)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    double menor = INF;
    for(int i=1; i<=N; i++)
        if(ins->ProdTime[i] < menor)
            menor = ins->ProdTime[i];

    Xx = Num3D(env,N+1);
    X_ini = Num2D(env,N+1);
    X_co = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = Num2D(env,T+1);
        X_ini[i] = IloNumVarArray(env,T+1);
        X_co[i] = IloNumVarArray(env,T+1);
        
        for(int t=0; t<=T; t++)
        {
            Xx[i][t] = IloNumVarArray(env,T+1);

            X_ini[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"X_ini_%hu_%hu",(unsigned short)i,(unsigned short)t);
            X_ini[i][t].setName(name);

            X_co[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"X_co_%hu_%hu",(unsigned short)i,(unsigned short)t);
            X_co[i][t].setName(name);

            for(int u=0; u<=T; u++)
            {
                Xx[i][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"Xx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)t,(unsigned short)u);
                Xx[i][t][u].setName(name);
            }
        }
        
    }

    R = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        R[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"R_%hu",(unsigned short)t);
        R[t].setName(name);
    }
    
    Zx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Zx[i] = IloNumVarArray(env,T+2);
        for(int t=0; t<=T+1; t++)
        {
            Zx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Zx_%hu_%hu",(unsigned short)i,(unsigned short)t);
            Zx[i][t].setName(name);
        }
    }

    Yx = Num3D(env,N+1);
    Fx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num2D(env,N+1);
        Fx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = IloNumVarArray(env,T+1);
            Fx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Yx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Yx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)t);
                Yx[i][j][t].setName(name);

                Fx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"Fx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)t);
                Fx[i][j][t].setName(name);
            }
            
        }
    }

    Gx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Gx[i] = IloNumVarArray(env,T+2);
        for(int t=0; t<=T+1; t++)
        {
            Gx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Gx_%hu_%hu",(unsigned short)i,(unsigned short)t);
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
            for(int u=t; u<=T; u++)
                obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();


    IloExpr expr(env);

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int t=1; t<=u; t++)
                expr += Xx[i][t][u];
            constraints.add(expr == ins->Demand[i][u]);
            expr.clear();
        }

    
    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int u=t; u<=T; u++)
                expr += ins->ProdTime[i]*Xx[i][t][u];
            for(int l=1; l<=N; l++)
                expr += ins->SetupTime[l][i]*Yx[l][i][t];
        }
        constraints.add(expr <= ins->Capacity[t]);
        expr.clear();
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int u=t; u<=T; u++)
            {
                int cap;
                cap = int((ins->Capacity[t])/menor);
                if(cap > ins->Demand[i][u])
                    cap = ins->Demand[i][u];
                expr += Xx[i][t][u];
                for(int j=1; j<=N; j++)
                    expr -= ins->MaxProd[i][t]*Yx[j][i][t];
                expr -= ins->MaxProd[i][t]*Zx[i][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    
    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            expr += Zx[i][t];
        constraints.add(expr == 1);
        expr.clear();
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr += Yx[j][i][t];
                expr -= Yx[i][j][t];
            }
            expr -= Zx[i][t+1];
            constraints.add(expr == 0);
            expr.clear();
        }
    
    for(int t=1; t<=T; t++)
    {
        for(int j=1; j<=N; j++)
        {
            expr += Fx[0][j][t];
            expr -= Gx[j][t];
        }
        constraints.add(expr == 0);
        expr.clear();
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr -= Gx[i][t];
            expr += Fx[0][i][t];
            for(int j=1; j<=N; j++)
            {
                expr += Fx[j][i][t];
                expr -= Fx[i][j][t];
            }
            constraints.add(expr == 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr += Fx[0][i][t];
            expr -= N*Zx[i][t];
            constraints.add(expr <= 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                expr += Fx[i][j][t];
                expr -= N*Yx[i][j][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    for(int t=1; t<=T; t++)
        for(int i=1; i<=N; i++)
        {
            expr += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr += Yx[j][i][t];
            }
            expr -= Gx[i][t];
            constraints.add(expr >= 0);
            expr.clear();
        }

    for(int t=1; t<=T; t++)
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                expr += Yx[j][i][t];
            }
            expr += Zx[i][t];
            expr -= BigM*Gx[i][t];
            constraints.add(expr <= 0);
            expr.clear();
        }

    /*for(int t=1; t<=T; t++)
        
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[j][i][t];
                constraints.add(expr <= 1);
                expr.clear();
            }

    for(int t=1; t<=T; t++)
        
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[i][j][t];
                constraints.add(expr <= 1);
                expr.clear();
            }*/

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                expr += Xx[i][t][u];
            expr -= X_co[i][t];
            expr -= X_ini[i][t];
            constraints.add(expr == 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr += X_co[i][t];
            expr -= BigM*Zx[i][t];
            constraints.add(expr <= 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr += X_ini[i][t];
            for(int j=1; j<=N; j++)
                expr -= 25*Yx[j][i][t];

            expr += 25*Zx[i][t+1];
            constraints.add(expr >= 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int s=t+1; s<=T; s++)
            {
                expr += X_ini[i][t];
                for(int u=t+1; u<=s; u++)
                    expr += X_co[i][u];
                for(int j=1; j<=N; j++)
                    expr -= 25*Yx[j][i][t];
                for(int u=t+1; u<=s-1; u++)
                    expr += BigM*R[u];
                expr -= BigM*R[s];
                constraints.add(expr >= -BigM);
                expr.clear();
            }
    
    expr.end();
    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Zx[i][0].setBounds(0,0);
        
        for(int t=1; t<=T; t++){
            Yx[i][i][t].setBounds(0,0);
            if(t < int(2*(T/3)))
                for(int j=1; j<=N; j++)
                    if(ins->vet[2].ystore[i][j][t] < 0.001 && ins->vet[3].ystore[i][j][t] < 0.001 && ins->vet[1].ystore[i][j][t] < 0.001 && i!=ins->nontri[0] && i!= ins->nontri[1] && i!=ins->nontri[2] && j!=ins->nontri[0] && j!= ins->nontri[1] && j!=ins->nontri[2])
                    {
                        Yx[i][j][t].setBounds(0,0);
                    }
            //if(ins->vet[2].zstore[i][t]==0 && ins->vet[1].zstore[i][t]==0)
                //Zx[i][t].setBounds(0,0);
        }
    }

    cplex = IloCplex(model);
    //cplex.exportModel("lp/RF.lp");
}

void RFOPR::CleanModel()
{
    obj.end();
    constraints.end();
    Fx.end();
    Gx.end();
    Xx.end();
    Yx.end();
    Zx.end();
    cplex.end();
}

void RFOPR::getFOptValue()
{
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                ins->X[i][t][u] = cplex.getValue(Xx[i][t][u]);
            ins->Z[i][t] = IloRound(cplex.getValue(Zx[i][t]));
            ins->G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
            ins->F[0][i][t] = cplex.getValue(Fx[0][i][t]);
            for(int j=1; j<=N; j++)
            {
                ins->F[i][j][t] = cplex.getValue(Fx[i][j][t]);
                ins->Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t]));
            }
        }
        ins->Z[i][T+1] = IloRound(cplex.getValue(Zx[i][T+1]));
    }
    BestSolution = cplex.getObjValue();
}

void RFOPR::store_solution()
{
    int id;

    id = 2;

    new_sol = true;
    ins->vet[id].sol = BestSolution;
    for(int i = 1; i <= N; i++)
    {
        for(int t = 1; t <= T; t++)
        {
            ins->vet[id].zstore[i][t] = ins->Z[i][t];
            for(int u = t; u <= T; u++)
                ins->vet[id].xstore[i][t][u] = ins->X[i][t][u];
            for(int j = 1; j <= N; j++)
                ins->vet[id].ystore[i][j][t] = ins->Y[i][j][t];
        }
        ins->vet[id].zstore[i][T+1] = ins->Z[i][T+1];
    }
}


void RFOPR::SetParameters()
{
    cplex.setParam(IloCplex::Param::Threads,1);
    cplex.setParam(IloCplex::Param::TimeLimit,subPR_time);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.000000001);
    //cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 8000.0);
    //Mip Display
    cplex.setParam(IloCplex::Param::MIP::Display, 3);
    //Mip Interval
    //cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
    //cplex.setOut(env.getNullStream());
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
}
