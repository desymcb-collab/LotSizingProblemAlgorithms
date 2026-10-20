#include "rf.h"

bool RF::RelaxAndFix()
{
    struct timespec fim;
    bool re = false, lt = true;
    int t,it;
    double total_time = 1800;

    env = IloEnv();
    subRF_time = 200;
    SolveRelaxedModel();
    jump = 1;           
    binary = 4;
    t = binary;

    it = ceil(double(T-binary)/double(jump))+1;

    while(t<=T)
    {
        clock_gettime(CLOCK_MONOTONIC, &fim);
        elapsedTime = (fim.tv_sec - start);
        subRF_time = floor((total_time-elapsedTime)/it);
        it -= 1;
        if(elapsedTime < total_time)
        {
            re = SolveRFModel(t,it);
            clock_gettime(CLOCK_MONOTONIC, &fim);
            elapsedTime = (fim.tv_sec - start);
            t += jump;
            if(t==T)
                lt = false;
            if(t>T && lt)
            {
                lt = false;
                t = T;
            }
        }else t = T+2;
    }

    env.end();
    //cout << "Number of feasible solutions: " << nm << "\n";
    cout << "RF: " << solutionRF << "\n";
    return re;
}

/************************************************************************
            Functions to solve and keep records of linear solution
************************************************************************/

void RF::SolveRelaxedModel()
{
    bool re = false;

    IloModel model(env);

    RelaxedModel(model);                                         //Constroi modelo classico
    SetParameters();
    //cplex.setOut(env.getNullStream());
    re = cplex.solve();

    ///Get solution of the linear relaxation of the original problem
    if(re)
    {
        getRelaxedValue();
        getReducedCost();
        CleanModel();
        model.end();
    }else LowerBound = -INF;
    cout << "LB: " << LowerBound << "\n";
}

void RF::RelaxedModel(IloModel model)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    double menor = INF;
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            if(ins->ProdTime[i][k] < menor)
                menor = ins->ProdTime[i][k];

    Xx = Num4D(env,N+1);
    X_ini = Num3D(env,N+1);
    X_co = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = Num3D(env,M+1);
        X_ini[i] = Num2D(env,M+1);
        X_co[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Xx[i][k] = Num2D(env,T+1);
            X_ini[i][k] = IloNumVarArray(env,T+1);
            X_co[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Xx[i][k][t] = IloNumVarArray(env,T+1);

                X_ini[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_ini_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t,(unsigned short)u);
                    Xx[i][k][t][u].setName(name);
                }
            }
        }
    }

    R = Num2D(env,M+1);
    for(int i=0; i<=M; i++)
    {
        R[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            R[i][t] = IloNumVar(env,0,1,ILOFLOAT);
            sprintf(name,"R_%hu_%hu",(unsigned short)i,(unsigned short)t);
            R[i][t].setName(name);
        }
    }

    Zx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Zx[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Zx[i][k] = IloNumVarArray(env,T+2);
            for(int t=0; t<=T+1; t++)
            {
                Zx[i][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                sprintf(name,"Zx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                Zx[i][k][t].setName(name);
            }
        }
    }

    Yx = Num4D(env,N+1);
    Fx = Num4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num3D(env,N+1);
        Fx[i] = Num3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = Num2D(env,M+1);
            Fx[i][j] = Num2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                Yx[i][j][k] = IloNumVarArray(env,T+1);
                Fx[i][j][k] = IloNumVarArray(env,T+1);
                for(int t=0; t<=T; t++)
                {
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Yx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)k,(unsigned short)t);
                    Yx[i][j][k][t].setName(name);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)k,(unsigned short)t);
                    Fx[i][j][k][t].setName(name);
                }
            }
        }
    }

    Gx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Gx[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Gx[i][k] = IloNumVarArray(env,T+2);
            for(int t=0; t<=T+1; t++)
            {
                Gx[i][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                sprintf(name,"Gx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                Gx[i][k][t].setName(name);
            }
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();


    IloExpr expr(env);

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr += Xx[i][k][t][u];
            constraints.add(expr == ins->Demand[i][u]);
            expr.clear();
        }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                    expr += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
            }
            constraints.add(expr <= ins->Capacity[k][t]);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                {
                    int cap;
                    cap = int((ins->Capacity[k][t])/menor);
                    if(cap > ins->Demand[i][u])
                        cap = ins->Demand[i][u];
                    expr += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr <= 0);
                    expr.clear();
                }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr += Zx[i][k][t];
            constraints.add(expr == 1);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                    expr -= Yx[i][j][k][t];
                }
                expr -= Zx[i][k][t+1];
                constraints.add(expr == 0);
                expr.clear();
            }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr += Fx[0][j][k][t];
                expr -= Gx[j][k][t];
            }
            constraints.add(expr == 0);
            expr.clear();
        }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr -= Gx[i][k][t];
                expr += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Fx[j][i][k][t];
                    expr -= Fx[i][j][k][t];
                }
                constraints.add(expr == 0);
                expr.clear();
            }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += Fx[0][i][k][t];
                expr -= N*Zx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr += Fx[i][j][k][t];
                    expr -= N*Yx[i][j][k][t];
                    constraints.add(expr <= 0);
                    expr.clear();
                }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                }
                expr -= Gx[i][k][t];
                constraints.add(expr >= 0);
                expr.clear();
            }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                }
                expr += Zx[i][k][t];
                expr -= BigM*Gx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    /*for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[j][i][k][t];
                constraints.add(expr <= 1);
                expr.clear();
            }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[i][j][k][t];
                constraints.add(expr <= 1);
                expr.clear();
            }*/

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr += Xx[i][k][t][u];
                expr -= X_co[i][k][t];
                expr -= X_ini[i][k][t];
                constraints.add(expr == 0);
                expr.clear();
            }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += X_co[i][k][t];
                expr -= BigM*Zx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr -= 25*Yx[j][i][k][t];

                expr += 25*Zx[i][k][t+1];
                constraints.add(expr >= 0);
                expr.clear();
            }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr -= 25*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr += BigM*R[k][u];
                    expr -= BigM*R[k][s];
                    constraints.add(expr >= -BigM);
                    expr.clear();
                }
    expr.end();
    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        for(int k=1; k<=M; k++)
        {
            Zx[i][k][0].setBounds(0,0);

            if(ins->MachineProd[i][k] == 0)
            {
                for(int t=1; t<=T; t++)
                {
                    Zx[i][k][t].setBounds(0,0);
                    Gx[i][k][t].setBounds(0,0);
                    for(int j=1; j<=N; j++)
                        Yx[j][i][k][t].setBounds(0,0);
                }
                Zx[i][k][T+1].setBounds(0,0);
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);
        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/Relax1.lp");
    //cplex.exportModel(name);
}

/************************************************************************
            Functions to solve and keep records of RF heuristic
************************************************************************/

bool RF::SolveRFModel(int t, int id)
{
    bool re = false;
    IloModel model(env);
    RFModel(model,t);                           //build subproblem model

    SetParameters();
    //cplex.setOut(env.getNullStream());
    re = cplex.solve();
    solutionRF = re?cplex.getObjValue():INF;       //If re = true then solutionRF = cplex.getObjValue(), otherwise, solutionRF = INF
    
    ///Get solution of the current problem
    if(re)
    {
        re = getRFSolution(t,id);
        CleanModel();
        model.end();
        //cout << "RF: " << solutionRF << "\n";
    }
    return re;
}

void RF::RFModel(IloModel model, int l, int id)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    double menor = INF;
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            if(ins->ProdTime[i][k] < menor)
                menor = ins->ProdTime[i][k];

    Xx = Num4D(env,N+1);
    X_ini = Num3D(env,N+1);
    X_co = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = Num3D(env,M+1);
        X_ini[i] = Num2D(env,M+1);
        X_co[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Xx[i][k] = Num2D(env,T+1);
            X_ini[i][k] = IloNumVarArray(env,T+1);
            X_co[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Xx[i][k][t] = IloNumVarArray(env,T+1);

                X_ini[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_ini_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t,(unsigned short)u);
                    Xx[i][k][t][u].setName(name);
                }
            }
        }
    }

    R = Num2D(env,M+1);
    for(int i=0; i<=M; i++)
    {
        R[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            if(t<=l)
            {
                R[i][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"R_%hu_%hu",(unsigned short)i,(unsigned short)t);
                R[i][t].setName(name); 
            }else
            {
                R[i][t] = IloNumVar(env,0,1,ILOFLOAT);
                sprintf(name,"R_%hu_%hu",(unsigned short)i,(unsigned short)t);
                R[i][t].setName(name);
            }
        }
    }

    Zx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Zx[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Zx[i][k] = IloNumVarArray(env,T+2);
            for(int t=0; t<=T+1; t++)
            {
                if(t<=l)
                {
                    Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                    sprintf(name,"Zx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                    Zx[i][k][t].setName(name);
                }
                else
                {
                    Zx[i][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                    sprintf(name,"Zx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                    Zx[i][k][t].setName(name);
                }
            }
        }
    }

    Yx = Num4D(env,N+1);
    Fx = Num4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num3D(env,N+1);
        Fx[i] = Num3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = Num2D(env,M+1);
            Fx[i][j] = Num2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                Yx[i][j][k] = IloNumVarArray(env,T+1);
                Fx[i][j][k] = IloNumVarArray(env,T+1);
                for(int t=0; t<=T; t++)
                {
                    if( (t<=l) )// || (reducedCostY[i][j][k][t] <= 0))
                    {
                        Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                        sprintf(name,"Yx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)k,(unsigned short)t);
                        Yx[i][j][k][t].setName(name);
                    }
                    else
                    {
                        Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                        sprintf(name,"Yx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)k,(unsigned short)t);
                        Yx[i][j][k][t].setName(name);
                    }

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx_%hu_%hu_%hu_%hu",(unsigned short)i,(unsigned short)j,(unsigned short)k,(unsigned short)t);
                    Fx[i][j][k][t].setName(name);
                }
            }
        }
    }

    Gx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Gx[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Gx[i][k] = IloNumVarArray(env,T+2);
            for(int t=0; t<=T+1; t++)
            {
                if(t <= l)
                {
                    Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Gx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                Gx[i][k][t].setName(name);    
                }else
                {
                    Gx[i][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                    sprintf(name,"Gx_%hu_%hu_%hu",(unsigned short)i,(unsigned short)k,(unsigned short)t);
                    Gx[i][k][t].setName(name);    
                }             
                
            }
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();

    IloExpr expr(env);

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr += Xx[i][k][t][u];
            constraints.add(expr == ins->Demand[i][u]);
            expr.clear();
        }


    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                    expr += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
            }
            constraints.add(expr <= ins->Capacity[k][t]);
            expr.clear();
        }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                {
                    int cap;
                    cap = int((ins->Capacity[k][t])/menor);
                    if(cap > ins->Demand[i][u])
                        cap = ins->Demand[i][u];
                    expr += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr <= 0);
                    expr.clear();
                }


    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr += Zx[i][k][t];
            constraints.add(expr == 1);
            expr.clear();
        }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                    expr -= Yx[i][j][k][t];
                }
                expr -= Zx[i][k][t+1];
                constraints.add(expr == 0);
                expr.clear();
            }


    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr += Fx[0][j][k][t];
                expr -= Gx[j][k][t];
            }
            constraints.add(expr == 0);
            expr.clear();
        }

 
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr -= Gx[i][k][t];
                expr += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Fx[j][i][k][t];
                    expr -= Fx[i][j][k][t];
                }
                constraints.add(expr == 0);
                expr.clear();
            }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += Fx[0][i][k][t];
                expr -= N*Zx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }


    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr += Fx[i][j][k][t];
                    expr -= N*Yx[i][j][k][t];
                    constraints.add(expr <= 0);
                    expr.clear();
                }


    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                }
                expr -= Gx[i][k][t];
                constraints.add(expr >= 0);
                expr.clear();
            }


    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr += Yx[j][i][k][t];
                }
                expr += Zx[i][k][t];
                expr -= BigM*Gx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }

    /*for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[j][i][k][t];
                constraints.add(expr <= 1);
                expr.clear();
            }


    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr += Yx[i][j][k][t];
                constraints.add(expr <= 1);
                expr.clear();
            }*/

 
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr += Xx[i][k][t][u];
                expr -= X_co[i][k][t];
                expr -= X_ini[i][k][t];
                constraints.add(expr == 0);
                expr.clear();
            }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += X_co[i][k][t];
                expr -= BigM*Zx[i][k][t];
                constraints.add(expr <= 0);
                expr.clear();
            }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr -= 25*Yx[j][i][k][t];

                expr += 25*Zx[i][k][t+1];
                constraints.add(expr >= 0);
                expr.clear();
            }


    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr -= 25*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr += BigM*R[k][u];
                    expr -= BigM*R[k][s];
                    constraints.add(expr >= -BigM);
                    expr.clear();
                }
    expr.end();
    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        for(int k=1; k<=M; k++)
        {
            Zx[i][k][0].setBounds(0,0);
            if(ins->MachineProd[i][k]==0)
            {
                for(int t=1; t<=T; t++)
                {
                    Zx[i][k][t].setBounds(0,0);
                    Gx[i][k][t].setBounds(0,0);
                    for(int j=1; j<=N; j++)
                        Yx[j][i][k][t].setBounds(0,0);
                }
                Zx[i][k][T+1].setBounds(0,0);
            }

            for(int t=1; t<=T; t++)
            {
                Yx[i][i][k][t].setBounds(0,0);
            }

            for(int j=1; j<=N; j++)
            {
                if(i!=j)
                {
                    for(int t=1; t<=T; t++)
                    {
                        if(t <= l-binary)
                        {
                            //printf("Variable %d_%d_%d_%d : %.4lf\n",i,j,k,t,Y[i][j][k][t]);
                            if(ins->Y[i][j][k][t] == 0.0)
                                Yx[i][j][k][t].setBounds(0,0);
                            else
                            {
                                int b;
                                b = int(ins->Y[i][j][k][t]);
                                Yx[i][j][k][t].setBounds(b,b);
                            }
                        }
                    }
                }
            }
        }
    }

    cplex = IloCplex(model);
    //cplex.exportModel("resultados/RF.lp");
    //getchar();
}

/*****************************************************************************
                                Auxiliary functions
*****************************************************************************/

void RF::CleanModel()
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

void RF::getRelaxedValue()
{
    bool re = true;

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    ins->vet[0].xstore[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                ins->vet[0].zstore[i][k][t] = cplex.getValue(Zx[i][k][t]);

                for(int j=1; j<=N; j++)
                {

                    ins->vet[0].ystore[i][j][k][t] = cplex.getValue(Yx[i][j][k][t]);
                    //printf("\nVariable %d_%d_%d_%d : %.4lf\n",i,j,k,t,Y[i][j][k][t]);
                    //getchar();
                }
            }
            ins->vet[0].zstore[i][k][T+1] = cplex.getValue(Zx[i][k][T+1]);
        }
    LowerBound = cplex.getObjValue();
}

void RF::getReducedCost()
{
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int j=1; j<=N; j++)
                for(int t=1; t<=T; t++)
                {
                    ins->reducedCostY[i][j][k][t] = cplex.getReducedCost(Yx[i][j][k][t]);
                    //printf("[%d][%d][%d][%d]: %.2lf!!!\n",i,j,k,t,ins->reducedCostY[i][j][k][t]);
                }
    //getchar();
}

bool RF::getRFSolution(int l, int id)
{
    bool re = true;

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    ins->X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                ins->Z[i][k][t] = 1.0*IloRound(cplex.getValue(Zx[i][k][t]));
                ins->G[i][k][t] = 1.0*IloRound(cplex.getValue(Gx[i][k][t]));
                ins->F[0][i][k][t] = cplex.getValue(Fx[0][i][k][t]);
                for(int j=1; j<=N; j++)
                {
                    ins->F[i][j][k][t] = cplex.getValue(Fx[i][j][k][t]);
                    
                    if(t > l)
                    {
                        ins->Y[i][j][k][t] = cplex.getValue(Yx[i][j][k][t]);
                    }else
                        ins->Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t]));
                        
                }
            }
            ins->Z[i][k][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][k][T+1]));
        }

    solutionRF = cplex.getObjValue();
    if(id == 3)
    {
        for(int i = 1; i <= N; i++)
            for(int k = 1; k <= M; k++)
            {
                for(int t = 1; t <= T; t++)
                {
                    ins->vet[1].zstore[i][k][t] = ins->Z[i][k][t];
                    for(int u = t; u <= T; u++)
                        ins->vet[1].xstore[i][k][t][u] = ins->X[i][k][t][u];
                    for(int j = 1; j <= N; j++)
                    {
                        ins->vet[1].ystore[i][j][k][t] = ins->Y[i][j][k][t];
                    }
                }
                ins->vet[1].zstore[i][k][T+1] = ins->Z[i][k][T+1];
            }
    }
    if(id == 0)
        store_solution();

    return re;
}

void RF::store_solution()
{
    if(ins->vet[2].sol > solutionRF)
    {
        ins->vet[2].sol = solutionRF;
        for(int i = 1; i <= N; i++)
            for(int k = 1; k <= M; k++)
            {
                for(int t = 1; t <= T; t++)
                {
                    ins->vet[2].zstore[i][k][t] = ins->Z[i][k][t];
                    for(int u = t; u <= T; u++)
                        ins->vet[2].xstore[i][k][t][u] = ins->X[i][k][t][u];
                    for(int j = 1; j <= N; j++)
                    {
                        ins->vet[2].ystore[i][j][k][t] = ins->Y[i][j][k][t];
                    }
                    
                }
                ins->vet[2].zstore[i][k][T+1] = ins->Z[i][k][T+1];

            }
    }
}

void RF::SetParameters()
{
    cplex.setParam(IloCplex::Param::Threads,1);
    cplex.setParam(IloCplex::Param::TimeLimit,subRF_time);
    //cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.000000001);
    cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 8000.0);
    //Mip Display
    cplex.setParam(IloCplex::Param::MIP::Display, 3);
    //Mip Interval
    cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
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
