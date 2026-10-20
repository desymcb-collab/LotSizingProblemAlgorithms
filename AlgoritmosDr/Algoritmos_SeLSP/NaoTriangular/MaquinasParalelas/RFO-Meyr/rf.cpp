#include "rf.h"

void RF::RelaxAndFix()
{
    struct timespec fim;
    double elapsedTime;
    bool fo, re = false, lt = true;
    best_bound = -INF;
    BestSolution = INF;

    AllocateVariables();
    env = IloEnv();
    nm = 0;
    pt = 0.7;
    bk = 3;

    //SolveRelaxedModel();
    int jump = 1;           //
    int jump_fo = 1;           //
    int t=1;
    binary = 4;
    while(t<=T)
    {
        clock_gettime(CLOCK_MONOTONIC, &fim);
        elapsedTime = (fim.tv_sec - start);
        if(elapsedTime < 2000)
        {
            re = SolveRFModel(t,jump,1);
            if(re)
                fo = true;
            clock_gettime(CLOCK_MONOTONIC, &fim);
            elapsedTime = (fim.tv_sec - start);
            t += jump;
            if(t>T && lt)
            {
                lt = false;
                t = T;
            }
        }//else t = T+2;
    }

    lt = false;
    while(fo)
    {
        lt = true;
        fo = false;
        for(int t=1; t<=T; t++)
        {
            clock_gettime(CLOCK_MONOTONIC, &fim);
            elapsedTime = (fim.tv_sec - start);
            if(elapsedTime < 3410)
            {
                bool bt;
                bt = false;
                re = FixAndOptimize(&bt,t,1);

                if(bt && fo!=true)
                    fo = true;

                clock_gettime(CLOCK_MONOTONIC, &fim);
                elapsedTime = (fim.tv_sec - start);
                //t += jump_fo;

            }
        }
        if(fo)
        {
            fo = false;
            for(int t=T; t>=1; t--)
            {
                clock_gettime(CLOCK_MONOTONIC, &fim);
                elapsedTime = (fim.tv_sec - start);
                if(elapsedTime < 3410)
                {
                    bool bt;
                    bt = false;
                    re = FixAndOptimize(&bt,t,1);

                    if(bt && fo!=true)
                        fo = true;

                    clock_gettime(CLOCK_MONOTONIC, &fim);
                    elapsedTime = (fim.tv_sec - start);
                    //t -= jump_fo;

                }
            }
        }
    }

//    if(lt)
//        KernelSearch();


    env.end();
    FreeVariables();
    //cout << "Number of feasible solutions: " << nm << "\n";
    cout << "Best solution: " << BestSolution;

}

void RF::SolveRelaxedModel()
{
    bool re = false;
    char name[20];
    IloModel model(env);

    RelaxedModel(model);                                         //Constroi modelo classico
    SetParameters();
    re = cplex.solve();

    ///Get solution of the linear relaxation of the original problem
    if(re)
    {
        /*getFirstItems();
        for(int k=1;k<=M;k++)
            Zx[firstItem[k]][k][1].setBounds(1,1);*/

        //sprintf(name,"lp/Relax.lp");
        //cplex.exportModel(name);
        //re = cplex.solve();
        //best_bound = re?cplex.getObjValue():INF;

            //getRelaxedValue();
        getReducedCost();
        CleanModel();
        model.end();
    }
    //cout << "Relaxed solution: " << best_bound << "\n\n";
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            sprintf(name,"R%d_%d",i,t);
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
                sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                    sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                    Yx[i][j][k][t].setName(name);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                sprintf(name,"Gx%d_%d_%d",i,k,t);
                Gx[i][k][t].setName(name);
            }
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
            {
                if(i==1)
                    obj1 += ins->SetupCost[0][j][k]*Zx[j][k][1];
                for(int t=1; t<=T; t++)
                {
                    obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];
                }
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                {
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
                }
                if(t==1)
                    expr2[k][t] += ins->SetupTime[0][i][k]*Zx[i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= ins->minl[i][k]*Yx[j][i][k][t];

                expr16[i][k][t] += ins->minl[i][k]*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= ins->minl[i][k]*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

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

bool RF::FixAndOptimize(bool *bt, int t, int id)
{
    bool re = false;
    IloModel model(env);

    FixedModel(model,t,id);                                         //Constroi modelo classico

    SetParameters();
    cplex.setOut(env.getNullStream());
    IloNumArray val(env);
    IloNumVarArray start(env);
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int j=1; j<=N; j++)
                {
                    start.add(Yx[i][j][k][t]);
                    val.add(vet[0].ystore[i][j][k][t]);
                }

                for(int u=t; u<=T; u++)
                {
                    start.add(Xx[i][k][t][u]);
                    val.add(vet[0].xstore[i][k][t][u]);
                }
                start.add(Zx[i][k][t]);
                val.add(vet[0].zstore[i][k][t]);
            }

    cplex.addMIPStart(start,val);
    start.end();
    val.end();
    re = cplex.solve();
    solutionRF = re?cplex.getObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //best_bound = re?cplex.getBestObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //cout << "FixAndOpt solution item " << t << ": " << solutionRF << " -- " << cplex.getSolnPoolNsolns() << "\n";
    ///Get solution of the current problem
    if(re)
    {
        if(BestSolution > solutionRF)
            *bt = true;
        getFOptValue();
        CleanModel();
        model.end();
        cout << "FO: " << BestSolution << "\n";
    }
    return re;
}

void RF::FixedModel(IloModel model, int l, int id)
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"R%d_%d",i,t);
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                    Yx[i][j][k][t].setName(name);
                    if(t<l)
                        Yx[i][j][k][t].setBounds(Y[i][j][k][t],Y[i][j][k][t]);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Gx%d_%d_%d",i,k,t);
                Gx[i][k][t].setName(name);
            }
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
            {
                if(i==1)
                    obj1 += ins->SetupCost[0][j][k]*Zx[j][k][1];
                for(int t=1; t<=T; t++)
                {
                    obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];
                }
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                {
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
                }
                if(t==1)
                    expr2[k][t] += ins->SetupTime[0][i][k]*Zx[i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= ins->minl[i][k]*Yx[j][i][k][t];

                expr16[i][k][t] += ins->minl[i][k]*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= ins->minl[i][k]*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

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
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);

        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF_%d.lp",l);
    //cplex.exportModel(name);
}

bool RF::SolveRFModel(int t, int jump, int id)
{
    bool re = false;
    IloModel model(env);
    RFModel(model,t,jump,id);                                         //Constroi modelo classico

    SetParameters();
    cplex.setOut(env.getNullStream());
    re = cplex.solve();
    solutionRF = re?cplex.getObjValue():INF;                        //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //best_bound = re?cplex.getBestObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //cout << "RF solution period " << t << ": " << solutionRF << " -- Quantity of fea sol: " << cplex.getSolnPoolNsolns() << "\n";
    ///Get solution of the current problem
    if(re)
    {
        re = getRFValue(t,id);
        CleanModel();
        model.end();
        cout << "RF: " << BestSolution << "\n";
    }
    return re;

}

void RF::RFModel(IloModel model, int l, int jump, int id)
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"R%d_%d",i,t);
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
                if(t<=l)
                {
                    Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                    sprintf(name,"Zx%d_%d_%d",i,k,t);
                    Zx[i][k][t].setName(name);
                }
                else
                {
                    Zx[i][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                    sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                        sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                        Yx[i][j][k][t].setName(name);
                    }
                    else
                    {
                        Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                        sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                        Yx[i][j][k][t].setName(name);
                    }

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Gx%d_%d_%d",i,k,t);
                Gx[i][k][t].setName(name);
            }
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
            {
                if(i==1)
                    obj1 += ins->SetupCost[0][j][k]*Zx[j][k][1];
                for(int t=1; t<=T; t++)
                {
                    obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];
                }
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    obj = IloMinimize(env, obj1);
    model.add(obj);
    obj1.end();

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                {
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
                }
                if(t==1)
                    expr2[k][t] += ins->SetupTime[0][i][k]*Zx[i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= ins->minl[i][k]*Yx[j][i][k][t];

                expr16[i][k][t] += ins->minl[i][k]*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= ins->minl[i][k]*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

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
                /*if(t <= l-2)
                {
                    if(Z[i][k][t] == 0.0)
                        Zx[i][k][t].setBounds(0,0);
                    else
                        Zx[i][k][t].setBounds(1,1);
                }*/
            }

            for(int j=1; j<=N; j++)
            {
                if(i!=j)
                for(int t=1; t<=T; t++)
                {
                    if(id == 1)
                    {
                        if(t <= l-binary)
                        {
                            //printf("Variable %d_%d_%d_%d : %.4lf\n",i,j,k,t,Y[i][j][k][t]);
                            if(Y[i][j][k][t] == 0.0)
                                Yx[i][j][k][t].setBounds(0,0);
                            else
                            {
                                int b;
                                b = int(Y[i][j][k][t]);
                                Yx[i][j][k][t].setBounds(b,b);
                            }
                        }
                    }
                }
            }
        }
    }

    cplex = IloCplex(model);
    //cplex.exportModel("lp/SM_NTS.lp");
    //getchar();
}

void RF::Update_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket)
{

    int l;
    for(int i = 1; i<=N; i++) {
        for(int j = 1; j<=N; j++) {
            for(int k = 1; k<=M; k++) {
                for(int t=1; t<=T; t++) {
                    if(Y[i][j][k][t] == 1){
                        l=busca_elemento(*kernel,t,j,i,k);
                        if(l<0)
                        {
                            Insere_ks(kernel,Y[i][j][k][t],t,j,i,k,1);
                            retira_elemento(bucket,t,j,i,k);
                        }
                    }
                /*else
                {
                    l=busca_elemento(*kernel,t,j,i);
                    if(l > 0){
                        kernel->campo[l].kernel -= 1;
                        if(kernel->campo[l].kernel == 0){
                            retira_elemento(kernel,t,j,i);
                            double custo;
                            custo = ins->spdc[i][j]*double(s->y[i][j][t]) + ins->cpdc[i][j][t]*double(s->x[i][j][t]);
                            Insere_ks(bucket,custo,t,j,i,0);
                        }
                    }
                    else{
                          double custo;
                          custo = ins->spdc[j][i]*double(s->y[i][j][t]) + ins->cpdc[i][j][t]*double(s->x[i][j][t]);
                          Insere_ks(bucket,custo,t,j,i,0);
                    }
                }*/
                }
            }
        }
    }

}

/*******************************************************************************/
/*******************************************************************************/

int RF::busca_binaria_ks(double LS, LRC_ks vetor)
{
    int e,m,d;

    e=0;
    d= vetor.nelementos-1;
    if (d==-1)
        return(0);
    while(e <= d ) {
        m = (e+d)/2;
        if(vetor.campo[m].cost <= LS) {
            if(m+1 < vetor.nelementos) {
                if(vetor.campo[m+1].cost > LS)
                    return(m+1);
                else e = m+1;
            } //if
            else if (m+1 == vetor.nelementos)
                return(m+1);
            else return(m);
        } //if
        else d = m-1;
    }//while
    return(m);
}

void RF::Insere_ks(LRC_ks *lrc, double campo1, int campo2, int campo3, int campo4, int campo5, int kernel)
{ //campos da estrutura
    int ii, aux;

    aux = busca_binaria_ks(campo1,*lrc);
    for(ii=lrc->nelementos-1; ii>=aux; ii--)
        lrc->campo[ii+1] = lrc->campo[ii];
    lrc->campo[aux].cost = campo1;
    lrc->campo[aux].t1 = campo2;
    lrc->campo[aux].j1 = campo3;
    lrc->campo[aux].i1 = campo4;
    lrc->campo[aux].k1 = campo5;
    lrc->campo[aux].kernel = kernel;
    lrc->nelementos++;
}

int RF::busca_elemento(LRC_ks vetor,int t,int j,int i,int k)
{
    int e,m,d;
    m=-1;
    e=0;
    d= vetor.nelementos-1;
    if (d==-1)
        return(0);
    while(e <= d ) {
        if((vetor.campo[e].t1==t) && (vetor.campo[e].k1==k) && (vetor.campo[e].i1==i) && (vetor.campo[e].j1==j)){
            return (vetor.nelementos);
        }
        e+=1;
    }//while
    return(m);
}

void RF::retira_elemento(LRC_ks *lrc, int campo2, int campo3,int campo4, int campo5)
{ //campos da estrutura
    int ii, aux;

    aux = busca_elemento(*lrc,campo2,campo3,campo4,campo5);
    lrc->campo[aux].kernel = 1;
//    for(ii=aux; ii<=lrc->nelementos-1; ii++)
//        lrc->campo[ii] = lrc->campo[ii+1];
//    lrc->nelementos--;
}

/*****************************************************************************/
/*****************************************************************************/

void RF::KernelSearch()
{
    double tempo2;
    int NB;

    LRC_ks kernel,bucket;

    Phase1(&kernel,&bucket);

    qtd_buckets = bk;
    NB = ceil(bucket.nelementos/bk);

    Phase2(&kernel,&bucket,NB);
}

void RF::Phase1(LRC_ks *kernel, LRC_ks *bucket)
{
    struct timespec fim;
    double tempo;

    RodaFLP_linear();                   ///Solve linear relaxation of model 2

    kernel->nelementos = 0;
    bucket->nelementos = 0;

    int rt=1;
    while(kernel->nelementos == 0)
    {
        if(solutionRF!=INF)
            Define_Kernel_Bucket1(kernel,bucket,rt);
        else
            Define_Kernel_Bucket2(kernel,bucket,rt);
        rt++;
    }

    //printf("Quantidade de elementos no kernel: %d\n",kernel->nelementos);
    //printf("Quantidade de elementos no bucket: %d\n",bucket->nelementos);
    clock_gettime(CLOCK_MONOTONIC, &fim);
    tempo = (fim.tv_sec - start);
    //tempo = difftime(tempo2,start);
    if(tempo < 3510)
    {
        Prob_Reduzido_kernel(*bucket);
        printf("Kernel --- Sol: %lf\n",solutionRF);
    }
}

void RF::Phase2(LRC_ks *kernel,LRC_ks *bucket, int NB)
{
    struct timespec fim;
    double tempo;

    for(int ii=1;ii<=qtd_buckets;ii++)
    {
        clock_gettime(CLOCK_MONOTONIC, &fim);
        tempo = (fim.tv_sec - start);
        //tempo = difftime(tempo2,start);
        if(tempo < 3510)
        {
            Prob_Reduzido_bucket(*kernel,*bucket,NB,ii);
            printf("Bucket: %d --- Sol: %lf\n",ii,solutionRF);
            if(solutionRF < INF)
            {
                Update_Kernel_Bucket(kernel,bucket);
            }
        }else ii = qtd_buckets+1;
    }
}

void RF::Define_Kernel_Bucket1(LRC_ks *kernel, LRC_ks *bucket, int rt)
{
    /*int aux[N+1][N+1][M+1][T+1];
    int taxa = floor((pt*vetKS)/rt);

    for(int i = 0; i<=N; i++)
        for(int j = 0; j<=N; j++)
            for(int k = 0; k<=M; k++)
                for(int t=0; t<=T; t++)
                    aux[i][j][k][t] = 0;

    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    for(int l=0; l < nm; l++)
                        if(vet[l].ystore[i][j][k][t] == 1 || vet[l].zstore[j][k][t] == 1)
                            aux[i][j][k][t] += 1;*/
    int cont = 0;
    int cont2 = 0;
    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    if(reducedCostY[i][j][k][t] <= 0 || vet[0].ystore[i][j][k][t] == 1 || vet[0].zstore[j][k][t] == 1)    //criterio para criacao do kernel (a variavel y tem valor zero em mais de taxa% das solucoes)
                    {
                        Insere_ks(kernel,reducedCostY[i][j][k][t],t,j,i,k,1);  //as variaveis que obtiverem solucao menor que o criterio durante mais de 5 solucao serao excluidas do kernel
                        cont++;
                    }else{
                        Insere_ks(bucket,reducedCostY[i][j][k][t],t,j,i,k,0);   // 0 indica se a variavel do bucket foi inserida no kernel durante alguma iteracao
                        cont2++;
                    }
    //printf("Cont: %d --- Cont2: %d\n",cont,cont2);
}

void RF::Define_Kernel_Bucket2(LRC_ks *kernel, LRC_ks *bucket, int rt)
{
    int aux[N+1][N+1][M+1][T+1];
    int taxa = floor((pt*vetKS)/rt);

    for(int i = 0; i<=N; i++)
        for(int j = 0; j<=N; j++)
            for(int k = 0; k<=M; k++)
                for(int t=0; t<=T; t++)
                    aux[i][j][k][t] = 0;

    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    for(int l=0; l < nm; l++)
                        if(vet[l].ystore[i][j][k][t] == 1 || vet[k].zstore[i][k][t] == 1)
                            aux[i][j][k][t] += 1;

    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    if(aux[i][j][k][t] >= taxa || vet[0].ystore[i][j][k][t] == 1)    //criterio para criacao do kernel (a variavel y tem valor zero em mais de taxa% das solucoes)
                    {
                        Insere_ks(kernel,(vetKS-aux[i][j][k][t]),t,j,i,k,1);  //as variaveis que obtiverem solucao menor que o criterio durante mais de 5 solucao serao excluidas do kernel
                    }else{
                        Insere_ks(bucket,(vetKS-aux[i][j][k][t]),t,j,i,k,0);   // 0 indica se a variavel do bucket foi inserida no kernel durante alguma iteracao
                    }
}

bool RF::Prob_Reduzido_kernel(LRC_ks ks_data)
{
    double end;
    bool re = false;
    LRC_ks bucket;
    bucket = ks_data;
    tempo = time(NULL);                                 //Inicia contagem de tempo para resolver problema

    env = IloEnv();
    re = Prob_Reduzido(bucket);                       //Variante do modelo de Sambasivan para a utilizacao de carry-over
    env.end();

    end = time(NULL);                                   //Encerra contagem de tempo para resolver o problema
    elapsed_time = difftime(end,tempo);

    //cout << "SolucaoFLP: " << LS->BestSolution << "\tLimitante Inferior: " << LS->LowerBound << "\n";

    return re;
}

bool RF::Prob_Reduzido(LRC_ks bucket)
{
    bool re = false;

    IloModel model(env);

    FLP_Reduzido(model,bucket);   //Constroi modelo localizacao de facilidades

    SetParameters();

    //cplex.setParam(IloCplex::TiLim, 800);
    switch(bk)
    {
    case 3:
        cplex.setParam(IloCplex::TiLim, 800);
        break;
	case 5:
	    cplex.setParam(IloCplex::TiLim, 550);
	    break;
    case 10:
        cplex.setParam(IloCplex::TiLim, 250);
        break;
    case 15:
        cplex.setParam(IloCplex::TiLim, 150);
        break;
    case 20:
        cplex.setParam(IloCplex::TiLim, 120);
        break;
    }

    IloNumArray val(env);
    IloNumVarArray start(env);

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int j=1; j<=N; j++)
                {
                    start.add(Yx[i][j][k][t]);
                    val.add(vet[0].ystore[i][j][k][t]);
                }

                for(int u=t; u<=T; u++)
                {
                    start.add(Xx[i][k][t][u]);
                    val.add(vet[0].xstore[i][k][t][u]);
                }
                start.add(Zx[i][k][t]);
                val.add(vet[0].zstore[i][k][t]);
            }

    cplex.addMIPStart(start,val);
    start.end();
    val.end();

    re = cplex.solve();
    ///Get solution of the current problem

    solutionRF = re?cplex.getObjValue():INF;

    if(re)
    {
        getFOptValue();
        CleanModel();
        model.end();
    }

    return re;
}

void RF::FLP_Reduzido(IloModel model, LRC_ks bucket)
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"R%d_%d",i,t);
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                    Yx[i][j][k][t].setName(name);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Gx%d_%d_%d",i,k,t);
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

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= 25*Yx[j][i][k][t];

                expr16[i][k][t] += 25*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= 25*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

    model.add(constraints);

    if(bucket.nelementos >= 1)
        for(int i=0; i<bucket.nelementos; i++)
            Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1].setBounds(0,0);

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
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);

        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF_%d.lp",l);
    //cplex.exportModel(name);
}

bool RF::Prob_Reduzido_bucket(LRC_ks ks_kernel, LRC_ks ks_data, int tam, int id)
{
    bool re = false;
    LRC_ks bucket,kernel;
    bucket = ks_data;
    kernel = ks_kernel;

    env = IloEnv();
    re = Prob_bucket(kernel,bucket,tam,id);                                    //Variante do modelo de Sambasivan para a utilizacao de carry-over
    env.end();
    //cout << "SolucaoFLP: " << LS->BestSolution << "\tLimitante Inferior: " << LS->LowerBound << "\n";

    return re;
}

bool RF::Prob_bucket(LRC_ks kernel, LRC_ks bucket, int tam, int id)
{
    bool re = false;

    IloModel model(env);

    FLP_Reduzido_bucket(model,kernel,bucket,tam,id);   //Constroi modelo localizacao de facilidades

    SetParameters();

    switch(bk)
    {
    case 3:
        cplex.setParam(IloCplex::TiLim, 800);
        break;
	case 5:
	    cplex.setParam(IloCplex::TiLim, 550);
	    break;
    case 10:
        cplex.setParam(IloCplex::TiLim, 250);
        break;
    case 15:
        cplex.setParam(IloCplex::TiLim, 150);
        break;
    case 20:
        cplex.setParam(IloCplex::TiLim, 120);
        break;
    }

    IloNumArray val(env);
    IloNumVarArray start(env);
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int j=1; j<=N; j++)
                {
                    start.add(Yx[i][j][k][t]);
                    val.add(vet[0].ystore[i][j][k][t]);
                }

                for(int u=t; u<=T; u++)
                {
                    start.add(Xx[i][k][t][u]);
                    val.add(vet[0].xstore[i][k][t][u]);
                }
                start.add(Zx[i][k][t]);
                val.add(vet[0].zstore[i][k][t]);
            }

    cplex.addMIPStart(start,val);
    start.end();
    val.end();

    re = cplex.solve();
    ///Get solution of the current problem

    solutionRF = re?cplex.getObjValue():INF;

    if(re)
    {
        getFOptValue();
        CleanModel();
        model.end();
    }

    return re;
}

void RF::FLP_Reduzido_bucket(IloModel model, LRC_ks kernel, LRC_ks bucket, int tam, int id)
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"R%d_%d",i,t);
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                    Yx[i][j][k][t].setName(name);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Gx%d_%d_%d",i,k,t);
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

    IloExpr expr_ks_2(env);

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    expr_ks_2 += ins->SetupCost[i][j][k]*Yx[i][j][k][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    expr_ks_2 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    constraints.add(expr_ks_2 <= (BestSolution -0.0001));

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= 25*Yx[j][i][k][t];

                expr16[i][k][t] += 25*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= 25*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

    expr_ks_18 = IloExpr(env);
    int aux;
    aux = id*tam;
    if(aux > bucket.nelementos)
        aux = bucket.nelementos;
    for(int i=(id-1)*tam; i < aux; i++)
            expr_ks_18 += Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1];

    //for(int i=0; i < kernel.nelementos; i++)
      //  expr_ks_8 += z[kernel.campo[i].item][kernel.campo[i].j1][kernel.campo[i].t1];

    constraints.add(expr_ks_18 >= 1);

    model.add(constraints);

    for(int i=0; i < (id-1)*tam; i++)
        if(bucket.campo[i].kernel == 0)
            Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1].setBounds(0,0);
    for(int i=id*tam; i < bucket.nelementos; i++)
        Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1].setBounds(0,0);

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
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);

        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF_%d.lp",l);
    //cplex.exportModel(name);
}

void RF::RodaFLP_linear()
{
    env = IloEnv();
    SolveFLP_linear();                                      //Variante do modelo de Sambasivan para a utilizacao de carry-over
    env.end();
}

void RF::SolveFLP_linear()
{
    bool re = false;

    IloModel model(env);

    FLP_Model_linear(model);   //Constroi modelo localizacao de facilidades

    SetParameters();

    cplex.setParam(IloCplex::TiLim, 100);

    re = cplex.solve();

    ///Get solution of the current problem
    if(re)
    {
        solutionRF = cplex.getObjValue();
        //getReducedCost();
        for(int i=1; i<=N; i++)
            for(int k=1; k<=M; k++)
                for(int j=1; j<=N; j++)
                    for(int t=1; t<=T; t++)
                        reducedCostY[i][j][k][t] = cplex.getReducedCost(Yx[i][j][k][t]);
        CleanModel();
        model.end();
    }else solutionRF = INF;
}

void RF::FLP_Model_linear(IloModel model)
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
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);

                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);

                for(int u=0; u<=T; u++)
                {
                    Xx[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Xx%d_%d_%d_%d",i,k,t,u);
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
            sprintf(name,"R%d_%d",i,t);
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
                sprintf(name,"Zx%d_%d_%d",i,k,t);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                    sprintf(name,"Yx%d_%d_%d_%d",i,j,k,t);
                    Yx[i][j][k][t].setName(name);

                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
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
                sprintf(name,"Gx%d_%d_%d",i,k,t);
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

    expr1 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr1[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int u=1; u<=T; u++)
        {
            for(int k=1; k<=M; k++)
                for(int t=1; t<=u; t++)
                    expr1[i][u] += Xx[i][k][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr2[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                for(int u=t; u<=T; u++)
                    expr2[k][t] += ins->ProdTime[i][k]*Xx[i][k][t][u];
                for(int l=1; l<=N; l++)
                    expr2[k][t] += ins->SetupTime[l][i][k]*Yx[l][i][k][t];
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr3[i][k] = Exp2D(env,T+1);
            for(int t=0; t<=T; t++)
            {
                expr3[i][k][t] = IloExprArray(env,T+1);
                for(int u=0; u<=T; u++)
                    expr3[i][k][t][u] = IloExpr(env);
            }
        }
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
                    expr3[i][k][t][u] += Xx[i][k][t][u];
                    for(int j=1; j<=N; j++)
                        expr3[i][k][t][u] -= ins->MaxProd[i][t]*Yx[j][i][k][t];
                    expr3[i][k][t][u] -= ins->MaxProd[i][t]*Zx[i][k][t];
                    constraints.add(expr3[i][k][t][u] <= 0);
                }

    expr4 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr4[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr4[k][t] += Zx[i][k][t];
            constraints.add(expr4[k][t] == 1);
        }

    expr5 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr5[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr5[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr5[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr5[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr5[i][k][t] += Yx[j][i][k][t];
                    expr5[i][k][t] -= Yx[i][j][k][t];
                }
                expr5[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr5[i][k][t] == 0);
            }

    expr6 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr6[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr6[k][t] += Fx[0][j][k][t];
                expr6[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr7[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][k][t] -= Gx[i][k][t];
                expr7[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += Fx[j][i][k][t];
                    expr7[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr7[i][k][t] == 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr8[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][k][t] += Fx[0][i][k][t];
                expr8[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr8[i][k][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr9[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr9[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr9[i][j][k][t] += Fx[i][j][k][t];
                    expr9[i][j][k][t] -= N*Yx[i][j][k][t];
                    constraints.add(expr9[i][j][k][t] <= 0);
                }


    expr10 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr10[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr10[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                expr10[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr10[i][k][t] += Yx[j][i][k][t];
                }
                expr10[i][k][t] -= Gx[i][k][t];
                constraints.add(expr10[i][k][t] >= 0);
            }

    expr11 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr11[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr11[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                {
                    expr11[i][k][t] += Yx[j][i][k][t];
                }
                expr11[i][k][t] += Zx[i][k][t];
                expr11[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr11[i][k][t] <= 0);
            }

    /*expr12 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr12[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr12[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr12[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr12[i][k][t] <= 1);
            }

    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr13[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr13[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr13[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr13[i][k][t] <= 1);
            }*/

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr14[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    expr14[i][k][t] += Xx[i][k][t][u];
                expr14[i][k][t] -= X_co[i][k][t];
                expr14[i][k][t] -= X_ini[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
            }

    expr15 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr15[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr15[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr15[i][k][t] += X_co[i][k][t];
                expr15[i][k][t] -= BigM*Z[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    expr16 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr16[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr16[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr16[i][k][t] += X_ini[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= 25*Yx[j][i][k][t];

                expr16[i][k][t] += 25*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp3D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = Exp2D(env,T+1);
            for(int j=0; j<=T; j++)
            {
                expr17[i][k][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr17[i][k][j][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int s=t+1; s<=T; s++)
                {
                    expr17[i][k][t][s] += X_ini[i][k][t];
                    for(int u=t+1; u<=s; u++)
                        expr17[i][k][t][s] += X_co[i][k][u];
                    for(int j=1; j<=N; j++)
                        expr17[i][k][t][s] -= 25*Yx[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }

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
            }
        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF.lp");
    //cplex.exportModel(name);
}

void RF::CleanModel()
{
    obj.end();
    constraints.end();

    Fx.end();
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
    expr13.end();
    expr14.end();
    expr15.end();
    expr16.end();
    expr17.end();

    cplex.end();
}

void RF::FixedVariables(int t)
{
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
            {
                if(reducedCostY[i][j][k][t] <= 0)
                    Y[i][j][k][t] = -1;
            }
}

void RF::getFirstItems()
{
    double initItem;
    //printf("Qtd de sol encontradas %d!\n",pool);

    for(int k=1; k<=M; k++)
    {
        initItem = 0.0;
        for(int i=1; i<=N; i++)
        {
            Z[i][k][1] = cplex.getValue(Zx[i][k][1]);
            if(initItem <= Z[i][k][1])
            {
                initItem = Z[i][k][1];
                firstItem[k] = i;
            }
        }
        printf("Primeiro item maquina %d: %d\n",k,firstItem[k]);
    }
}

void RF::getRelaxedValue()
{
    bool re = true;
    double initItem[M+1];
    //printf("Qtd de sol encontradas %d!\n",pool);

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                Z[i][k][t] = cplex.getValue(Zx[i][k][t]);

                for(int j=1; j<=N; j++)
                {

                    Y[i][j][k][t] = cplex.getValue(Yx[i][j][k][t]);
                    //printf("\nVariable %d_%d_%d_%d : %.4lf\n",i,j,k,t,Y[i][j][k][t]);
                    //getchar();
                    double aux;
                    aux = cplex.getValue(Yx[i][j][k][t]);
                    aux = aux - floor(aux+0.1);
                    if(fabs(aux) > 0.01 && re==true)
                    {
                        re = false;
                        //printf("\nContinuous value %.4lf %.4lf %.4lf\n",Y[i][j][k][t],aux,fabs(aux));
                    }
                }
            }
            Z[i][k][T+1] = cplex.getValue(Zx[i][k][T+1]);
        }
    solutionRF = cplex.getObjValue();
    if(re)
    {
        store_solution();
        printf("Optimal solution...you are not doing this right!!\n");
        if(BestSolution > solutionRF)
            BestSolution = solutionRF;
    }

}

void RF::getReducedCost()
{
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int j=1; j<=N; j++)
                for(int t=1; t<=T; t++)
                    reducedCostY[i][j][k][t] = cplex.getReducedCost(Yx[i][j][k][t]);
}

bool RF::getRFValue(int l, int id)
{
    int pool;
    bool re = true;
    pool = cplex.getSolnPoolNsolns();
    //printf("Qtd de sol encontradas %d!\n",pool);
    for(int p=0;p<pool;p++)
    {
        for(int i=1; i<=N; i++)
            for(int k=1; k<=M; k++)
            {
                for(int t=1; t<=T; t++)
                {
                    for(int u=t; u<=T; u++)
                        X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u],p);

                    Z[i][k][t] = 1.0*IloRound(cplex.getValue(Zx[i][k][t],p));

                    for(int j=1; j<=N; j++)
                    {
                        if(id==1)
                        {
                            if(t > l)
                            {
                                if(reducedCostY[i][j][k][t]>0)
                                    Y[i][j][k][t] = cplex.getValue(Yx[i][j][k][t],p);
                                else
                                    Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t],p));
                                double aux;
                                aux = cplex.getValue(Yx[i][j][k][t],p);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    //printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][k][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t],p));
                        }
                    }
                }
                Z[i][k][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][k][T+1],p));
            }
        solutionRF = cplex.getObjValue(p);
        if(re)
        {
            store_solution();
            //printf("Solution stored %.4lf!\n",solutionRF);
            if(BestSolution > solutionRF)
                BestSolution = solutionRF;
        }
        //re = true;
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                Z[i][k][t] = 1.0*IloRound(cplex.getValue(Zx[i][k][t]));
                G[i][k][t] = IloRound(cplex.getValue(Gx[i][k][t]));
                F[0][i][k][t] = cplex.getValue(Fx[0][i][k][t]);
                for(int j=1; j<=N; j++)
                {
                    F[i][j][k][t] = cplex.getValue(Fx[i][j][k][t]);
                    if(id==1)
                        {
                            if(t > l)
                            {
                                if(reducedCostY[i][j][k][t]>0)
                                    Y[i][j][k][t] = cplex.getValue(Yx[i][j][k][t]);
                                else
                                    Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t]));
                                double aux;
                                aux = cplex.getValue(Yx[i][j][k][t]);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    //printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][k][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t]));
                        }
                }
            }
            Z[i][k][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][k][T+1]));
        }

    return re;
}

void RF::getFOptValue()
{
    int pool;
    pool = cplex.getSolnPoolNsolns();
    //printf("Qtd de sol encontradas %d!\n",pool);
    for(int p=0;p<pool;p++)
    {
        for(int i=1; i<=N; i++)
            for(int k=1; k<=M; k++)
            {
                for(int t=1; t<=T; t++)
                {
                    for(int u=t; u<=T; u++)
                        X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u],p);
                    Z[i][k][t] = IloRound(cplex.getValue(Zx[i][k][t],p));

                    for(int j=1; j<=N; j++)
                        Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t],p));
                }
                Z[i][k][T+1] = IloRound(cplex.getValue(Zx[i][k][T+1],p));
            }
        solutionRF = cplex.getObjValue(p);
        store_solution();
        //printf("Solution stored in FixAndOpt %.4lf!\n",solutionRF);
        if(BestSolution > solutionRF)
            BestSolution = solutionRF;
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                Z[i][k][t] = IloRound(cplex.getValue(Zx[i][k][t]));
                G[i][k][t] = IloRound(cplex.getValue(Gx[i][k][t]));
                F[0][i][k][t] = cplex.getValue(Fx[0][i][k][t]);
                for(int j=1; j<=N; j++)
                {
                    F[i][j][k][t] = cplex.getValue(Fx[i][j][k][t]);
                    Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t]));
                }
            }
            Z[i][k][T+1] = IloRound(cplex.getValue(Zx[i][k][T+1]));
        }
}

bool RF::verify_duplication()
{
    for(int k = 0; k < nm; k++)
    {
        bool duplication = true;
        if(vet[k].sol != solutionRF) continue;
        for(int i = 1; i <= N; i++)
            for(int l = 1; l <= M; l++)
                for(int t = 1; t <= T; t++)
                    for(int j = 1; j <= N; j++)
                    {
                        if(vet[k].ystore[i][j][l][t] != Y[i][j][l][t])
                        {
                            duplication = false;
                            i = N+1;
                            j = N+1;
                            t = T+1;
                            l = M+1;
                        }
                    }

        for(int i = 1; i <= N; i++)
            for(int l = 1; l <= M; l++)
                for(int t = 1; t <= T; t++)
                    for(int u = t; u <= T; u++)
                    {
                        if(vet[k].xstore[i][l][t][u] != X[i][l][t][u])
                        {
                            duplication = false;
                            i = N+1;
                            u = T+1;
                            t = T+1;
                            l = M+1;
                        }
                    }

        for(int i = 1; i <= N; i++)
            for(int l = 1; l <= M; l++)
                for(int t = 1; t <= T; t++)
                    if(vet[k].zstore[i][l][t] != Z[i][l][t])
                    {
                        duplication = false;
                        i = N+1;
                        l = M+1;
                        t = T+1;
                    }

        if(duplication) return true;
    }
    return false;
}

void RF::store_solution()
{
    //printf("Posicao %d\n",pos);
    bool duplication = false;
    duplication = verify_duplication();
    if(duplication) printf("Solucao repetida\n");

    if(!duplication)
    {
        if((vet[0].sol < solutionRF) || (vet[0].sol == INF))
        {
            vet[0].sol = solutionRF;
            for(int i = 1; i <= N; i++)
                for(int k = 1; k <= M; k++)
                {
                    for(int t = 1; t <= T; t++)
                    {
                        vet[0].zstore[i][k][t] = Z[i][k][t];
                        for(int u = t; u <= T; u++)
                            vet[0].xstore[i][k][t][u] = X[i][k][t][u];
                        for(int j = 1; j <= N; j++)
                            vet[0].ystore[i][j][k][t] = Y[i][j][k][t];
                    }
                    vet[0].zstore[i][k][T+1] = Z[i][k][T+1];
                }
        }

        if(vet[1].sol > solutionRF)
        {
            vet[1].sol = solutionRF;
            for(int i = 1; i <= N; i++)
                for(int k = 1; k <= M; k++)
                {
                    for(int t = 1; t <= T; t++)
                    {
                        vet[1].zstore[i][k][t] = Z[i][k][t];
                        for(int u = t; u <= T; u++)
                            vet[1].xstore[i][k][t][u] = X[i][k][t][u];
                        for(int j = 1; j <= N; j++)
                            vet[1].ystore[i][j][k][t] = Y[i][j][k][t];
                    }
                    vet[1].zstore[i][k][T+1] = Z[i][k][T+1];
                }
        }

        //printf("Deslocou solucoes\n");


        //if(nm < vetKS) nm++;
    }
    //printf("Guardou!\n");
}

void RF::AllocateVariables()
{
    X.resize(N+1);
    Z.resize(N+1);
    F.resize(N+1);
    G.resize(N+1);
    Y.resize(N+1);
    firstItem.resize(M+1);
    reducedCostY.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        F[i].resize(N+1);
        G[i].resize(M+1);
        X[i].resize(M+1);
        Z[i].resize(M+1);
        Y[i].resize(N+1);
        reducedCostY[i].resize(N+1);
        for(int k=0;k<=M;k++)
        {
            X[i][k].resize(T+1);
            Z[i][k].resize(T+2,0.0);
            G[i][k].resize(T+1,0.0);
            for(int t=0;t<=T;t++)
                X[i][k][t].resize(T+1,0.0);
        }
        for(int j=0;j<=N;j++)
        {
            reducedCostY[i][j].resize(M+1);
            Y[i][j].resize(M+1);
            F[i][j].resize(M+1);
            for(int k=0;k<=M;k++)
            {
                reducedCostY[i][j][k].resize(T+1,0.0);
                Y[i][j][k].resize(T+1,0);
                F[i][j][k].resize(T+1,0.0);
            }
        }
    }

    for(int k=0; k<vetKS; k++)
    {
        vet[k].sol = INF;
        vet[k].xstore.resize(N+1);
        vet[k].ystore.resize(N+1);
        vet[k].zstore.resize(N+1);
        for(int i = 0; i <= N; i++)
        {
            vet[k].xstore[i].resize(M+1);
            vet[k].ystore[i].resize(N+1);
            vet[k].zstore[i].resize(M+1);
            for(int l = 0; l <= M; l++)
            {
                vet[k].xstore[i][l].resize(T+1);
                vet[k].zstore[i][l].resize(T+2,0.0);
                for(int t = 0; t <= T; t++)
                    vet[k].xstore[i][l][t].resize(T+1,0.0);
            }
            for(int j = 0; j <= N; j++)
            {
                vet[k].ystore[i][j].resize(M+1);
                for(int l = 0; l <= M; l++)
                    vet[k].ystore[i][j][l].resize(T+1,0);
            }
        }

    }
}

void RF::FreeVariables()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            for(int k=0;k<=M;k++)
            {
                reducedCostY[i][j][k].clear();
                Y[i][j][k].clear();
                F[i][j][k].clear();
            }
            reducedCostY[i][j].clear();
            Y[i][j].clear();
            F[i][j].clear();
        }
        for(int k=0;k<=M;k++)
        {
            for(int t=0;t<=T;t++)
                X[i][k][t].clear();
            X[i][k].clear();
            Z[i][k].clear();
            G[i][k].clear();
        }
        X[i].clear();
        Z[i].clear();
        F[i].clear();
        G[i].clear();
        Y[i].clear();
        reducedCostY[i].clear();
    }
    firstItem.clear();
    X.clear();
    reducedCostY.clear();
    Y.clear();
    F.clear();
    G.clear();
    Z.clear();

    for(int k=vetKS-1;k>=0;k--)
    {
        for(int i = 0;i <= N; i++)
        {
            for(int j = 0;j <= N; j++)
            {
                for(int l = 0;l <= M; l++)
                    vet[k].ystore[i][j][l].clear();
                vet[k].ystore[i][j].clear();
            }
            for(int l = 0;l <= M; l++)
            {
                for(int t = 0; t <= T; t++)
                    vet[k].xstore[i][l][t].clear();
                vet[k].zstore[i][l].clear();
                vet[k].xstore[i][l].clear();
            }

            vet[k].xstore[i].clear();
            vet[k].ystore[i].clear();
            vet[k].zstore[i].clear();
        }
        vet[k].xstore.clear();
        vet[k].ystore.clear();
        vet[k].zstore.clear();
    }
}

void RF::SetParameters()
{
    cplex.setParam(IloCplex::Threads,4);
    cplex.setParam(IloCplex::TiLim, 150);
    cplex.setParam(IloCplex::EpGap, 0.000000001);
    cplex.setParam(IloCplex::TreLim, 14000.0);
    cplex.setOut(env.getNullStream());
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
