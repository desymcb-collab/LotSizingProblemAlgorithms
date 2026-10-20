#include "scf_fl.h"

void Models::RodaGuimaraes()
{
    bool re = false;
    best_bound = INF;
    BestSolution = INF;

    //ins->show_instance();
    //getchar();
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    aloca_Guimaraes();

    env = IloEnv();
    re = SolveGuimaraes();
    env.end();

    funcao_objetivo();
    printSolution();
    //printSeq();
    libera_Guimaraes();

    cout << "Best Solution: " << BestSolution << " --- Best Bound: " << best_bound;
        //cout << "\nBest Solution: " << BestSolution << "\tBest Bound: " << best_bound << "\tTempo: " << time << "\n";

}

bool Models::SolveGuimaraes()
{
    bool re = false;
    IloModel model(env);
    Guimaraes_Model(model);                                         //Constroi modelo classico

    set_parameters();

    re = cplex.solve();
    BestSolution = re?cplex.getObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    best_bound = re?cplex.getBestObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF

    ///Get solution of the current problem
    if(re)
    {
        getValue();
        /*while(!Algoritmo_FordFulkerson(model) && re)
        {

            re = cplex.solve();
            if(re)
            {
                getValue();
                BestSolution = cplex.getObjValue();
                printf("BestSolution: %.2lf\n",BestSolution);

            }else re = false;
        }
        BestSolution = re?cplex.getObjValue():min(BestSolution,INF);*/

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
        for(int k=1; k<=M; k++)
            if(ins->ProdTime[i][k] < menor)
                menor = ins->ProdTime[i][k];

    Xx = Num4D(env,N+1);
    X_a = Num3D(env,N+1);
    X_b = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = Num3D(env,M+1);
        X_a[i] = Num2D(env,M+1);
        X_b[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Xx[i][k] = Num2D(env,T+1);
            X_a[i][k] = IloNumVarArray(env,T+1);
            X_b[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Xx[i][k][t] = IloNumVarArray(env,T+1);

                X_a[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_a%d_%d_%d",i,k,t);
                X_a[i][k][t].setName(name);

                X_b[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_b%d_%d_%d",i,k,t);
                X_b[i][k][t].setName(name);

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

    SALEx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        SALEx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            SALEx[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"S%d_%d",i,t);
            SALEx[i][t].setName(name);
        }
    }

    Yx = Num5D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num4D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = Num3D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                Yx[i][j][k] = Num2D(env,T+1);
                for(int t=0; t<=T; t++)
                {
                    Yx[i][j][k][t] = IloNumVarArray(env,T+1);
                    for(int u=0; u<=T; u++)
                    {
                        Yx[i][j][k][t][u] = IloNumVar(env,0,IloInfinity,ILOINT);
                        sprintf(name,"Yx%d_%d_%d_%d_%d",i,j,k,t,u);
                        Yx[i][j][k][t][u].setName(name);
                    }
                }
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Zx%d_%d_%d",i,k,t);
                Zx[i][k][t].setName(name);
            }
        }
    }

    Lx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Lx[i] = Num2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            Lx[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Lx[i][k][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Lx%d_%d_%d",i,k,t);
                Lx[i][k][t].setName(name);
            }
        }
    }

    Ix = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Ix[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            Ix[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            if(t==0)
                Ix[i][t].setBounds(ins->istock[i][t+1],ins->istock[i][t+1]);
            sprintf(name,"I%d_%d",i,t);
            Ix[i][t].setName(name);
        }
    }

    Fx = Num4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Fx[i] = Num3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Fx[i][j] = Num2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                Fx[i][j][k] = IloNumVarArray(env,T+1);
                for(int t=0; t<=T; t++)
                {
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
        for(int k=1; k<=M; k++)
                obj1 += ins->iscost[i][k]*Zx[i][k][1];

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    for(int u=t; u<=T; u++)
                    {
                        //if(i!=j)
                            obj1 += ins->SetupCost[i][j][k]*Yx[i][j][k][t][u];
                    }
    /*cost of last produced item per period*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                obj1 += ins->lscost[i][k]*Lx[i][k][t];
            }

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    /*market price*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            obj1 -= ins->MP[i][t]*SALEx[i][t];

    /*initial holding cost*/
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
            expr1[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int u=1; u<=t; u++)
                expr1[i][t] += Ix[i][u];
            constraints.add(expr1[i][t] == ins->istock[i][0]);
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
            }
            constraints.add(expr2[k][t] <= ins->Capacity[k][t]);
        }

    expr3 = Exp2D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        expr3[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr3[k][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr3[i][t] += SALEx[i][t];
            constraints.add(expr3[i][t] >= ins->Demand[i][t]);
        }

    expr4 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr4[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr4[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int k=1; k<=M; k++)
            {
                for(int u=1; u<=t; u++)
                    expr4[i][t] += Xx[i][k][u][t];
            }
            expr4[i][t] += Ix[i][t-1];
            expr4[i][t] -= Ix[i][t];
            expr4[i][t] -= SALEx[i][t];
            constraints.add(expr4[i][t] == 0);
        }

    expr5 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr5[k] = IloExprArray(env,T+2);
        for(int t=0; t<=T+1; t++)
            expr5[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T+1; t++)
        {
            for(int i=1; i<=N; i++)
                expr5[k][t] += Zx[i][k][t];
            constraints.add(expr5[k][t] == 1);
        }

    expr25 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr25[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr25[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
                expr25[k][t] += Lx[i][k][t];
            constraints.add(expr25[k][t] == 1);
        }

    expr6 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr6[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr6[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr6[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr6[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr6[i][k][t] += Yx[j][i][k][t][t];
                    expr6[i][k][t] -= Yx[i][j][k][t][t];
//                    for(int u=t+1; u<=T; u++)
//                    {
//                        expr6[i][k][t] -= Yx[i][j][k][t][u];
//                        expr6[i][k][t] += Zx[j][k][u];
//                    }
                }
                expr6[i][k][t] -= Lx[i][k][t];

                constraints.add(expr6[i][k][t] == 0);
            }

    expr26 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr26[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr26[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr26[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<T; t++)
            {
                for(int j=1; j<=N; j++)
                {
                        expr26[i][k][t] += Yx[j][i][k][t][t+1];

                }
                expr26[i][k][t] -= Zx[i][k][t+1];
                constraints.add(expr26[i][k][t] == 0);
            }

    expr27 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr27[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr27[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr27[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<T; t++)
            {
                for(int j=1; j<=N; j++)
                {
                        expr27[i][k][t] += Yx[i][j][k][t][t+1];

                }
                expr27[i][k][t] -= Lx[i][k][t];
                constraints.add(expr27[i][k][t] == 0);
            }



    /*expr21 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr21[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr21[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr21[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr21[i][k][t] += Yx[j][i][k][t];
                constraints.add(expr21[i][k][t] <= 1);
            }

    expr22 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr22[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr22[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr22[i][k][t] = IloExpr(env);
        }
    }

    for(int t=1; t<=T; t++)
        for(int k=1; k<=M; k++)
            for(int i=1; i<=N; i++)
            {
                for(int j=1; j<=N; j++)
                    expr22[i][k][t] += Yx[i][j][k][t];
                constraints.add(expr22[i][k][t] <= 1);
            }*/

    expr8 = Exp2D(env,M+1);
    for(int k=0; k<=M; k++)
    {
        expr8[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr8[k][t] = IloExpr(env);
    }

    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                expr8[k][t] += Fx[0][j][k][t];
                expr8[k][t] -= Gx[j][k][t];
            }
            constraints.add(expr8[k][t] == 0);
        }

    expr9 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr9[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr9[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr9[i][k][t] -= Gx[i][k][t];
                expr9[i][k][t] += Fx[0][i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr9[i][k][t] += Fx[j][i][k][t];
                    expr9[i][k][t] -= Fx[i][j][k][t];
                }
                constraints.add(expr9[i][k][t] == 0);
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

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
            {
                expr10[i][k][t] += Fx[0][i][k][t];
                expr10[i][k][t] -= N*Zx[i][k][t];
                constraints.add(expr10[i][k][t] <= 0);
            }

    expr11 = Exp4D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = Exp3D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr11[i][j] = Exp2D(env,M+1);
            for(int k=0; k<=M; k++)
            {
                expr11[i][j][k] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr11[i][j][k][t] = IloExpr(env);
            }
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                {
                    expr11[i][j][k][t] += Fx[i][j][k][t];
                    expr11[i][j][k][t] -= N*Yx[i][j][k][t][t];
                    constraints.add(expr11[i][j][k][t] <= 0);
                }


    expr12 = Exp3D(env,N+1);
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
                expr12[i][k][t] += Zx[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr12[i][k][t] += Yx[j][i][k][t][t];
                }
                expr12[i][k][t] -= Gx[i][k][t];
                constraints.add(expr12[i][k][t] >= 0);
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
                {
                    expr13[i][k][t] += Yx[j][i][k][t][t];
                }
                expr13[i][k][t] += Zx[i][k][t];
                expr13[i][k][t] -= BigM*Gx[i][k][t];
                constraints.add(expr13[i][k][t] <= 0);
            }

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
                expr14[i][k][t] -= X_b[i][k][t];
                expr14[i][k][t] -= X_a[i][k][t];
                constraints.add(expr14[i][k][t] == 0);
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
            for(int t=1; t<T; t++)
            {
                expr16[i][k][t] += X_a[i][k][t];
                for(int j=1; j<=N; j++)
                    expr16[i][k][t] -= ins->minl[i][k]*Yx[j][i][k][t][t];


                expr16[i][k][t] += ins->minl[i][k]*Zx[i][k][t+1];
                constraints.add(expr16[i][k][t] >= 0);
            }

    expr17 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr17[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr17[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<T; t++)
                {
                expr17[i][k][t] += X_a[i][k][t];
                expr17[i][k][t] += X_b[i][k][t+1];
                for(int j=1; j<=N; j++)
                    expr17[i][k][t] -= ins->minl[i][k]*Yx[j][i][k][t][t];
                constraints.add(expr17[i][k][t] >= 0);
            }

    expr23 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr23[i] = IloExprArray(env,M+1);
        for(int k=0; k<=M; k++)
            expr23[i][k] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            expr23[i][k] += X_a[i][k][T];
            for(int j=1; j<=N; j++)
            {
                expr23[i][k] -= ins->minl[i][k]*Yx[j][i][k][T][T];
            }
            //expr23[i][k] -= ins->maxl[i][k]*Zx[i][k][T];
            constraints.add(expr23[i][k] >= 0);
        }

    expr18 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr18[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr18[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr18[i][k][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            //for(int t=1; t<=T; t++)
            {
//                for(int u=t; u<=T; u++)
//                    expr18[i][k][t] += Xx[i][k][t][u];
//                for(int j=1; j<=N; j++)
//                    expr18[i][k][t] -= ins->maxl[i][k]*Yx[j][i][k][t];
//                expr18[i][k][t] -= ins->maxl[i][k]*Zx[i][k][t];
//                expr18[i][k][t] -= ins->maxl[i][k]*Zx[i][k][t+1];
                expr18[i][k][1] += X_b[i][k][1];
                expr18[i][k][1] -= ins->minl[i][k]*Zx[i][k][1];
                constraints.add(expr18[i][k][1] >= 0);
            }

    expr19 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr19[i] = Exp2D(env,M+1);
        for(int k=0; k<=M; k++)
        {
            expr19[i][k] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr19[i][k][t] = IloExpr(env);
        }
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
                expr15[i][k][t] += X_b[i][k][t];
                expr15[i][k][t] -= ins->maxl[i][k]*Zx[i][k][t];
                constraints.add(expr15[i][k][t] <= 0);
            }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<T; t++)
            {
                expr19[i][k][t] += X_a[i][k][t];
                expr19[i][k][t] += X_b[i][k][t+1];
                for(int j=1; j<=N; j++)
                    expr19[i][k][t] -= ins->maxl[i][k]*Yx[j][i][k][t][t];

                constraints.add(expr19[i][k][t] <= 0);
            }

    expr20 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr20[i] = IloExprArray(env,M+1);
        for(int k=0; k<=M; k++)
            expr20[i][k] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
        {
            expr20[i][k] += X_a[i][k][T];
            for(int j=1; j<=N; j++)
            {
                expr20[i][k] -= ins->maxl[i][k]*Yx[j][i][k][T][T];
            }
            //expr20[i][k] -= ins->maxl[i][k]*Zx[i][k][T];
            constraints.add(expr20[i][k] <= 0);
        }


    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        for(int k=1; k<=M; k++)
        {
            Zx[i][k][0].setBounds(0,0);
            Ix[i][T].setBounds(0,0);

            if(ins->MachineProd[i][k]==0)
            {
                for(int t=1; t<=T; t++)
                {
                    Zx[i][k][t].setBounds(0,0);
                    Gx[i][k][t].setBounds(0,0);
                    for(int j=1; j<=N; j++)
                        Yx[j][i][k][t][t].setBounds(0,0);
                }
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t][t].setBounds(0,0);
        }
    }
    for(int k=1; k<=M; k++)
        if(ins->istate[k] != 0)
            Zx[ins->istate[k]][k][1].setBounds(ins->istate[k],ins->istate[k]);

    cplex = IloCplex(model);
    cplex.exportModel("lp/fl.lp");
}

void Models::CleanGM()
{
    obj.end();
    constraints.end();

    Fx.end();
    Gx.end();
    Xx.end();
    X_b.end();
    X_a.end();
    Yx.end();
    Zx.end();
    I.end();
    SALE.end();
    R.end();

    expr1.end();
    expr2.end();
    expr3.end();
    expr4.end();
    expr5.end();
    expr6.end();
    //expr7.end();
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
    expr18.end();
    expr19.end();
    expr20.end();
//    expr21.end();
//    expr22.end();
    expr23.end();
    expr24.end();
    expr25.end();
    expr26.end();

    cplex.end();
}

void Models::printSeq()
{
    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            bool visited=false;
            printf("M - %d --- T - %d\n",k,t);
            for(int i=1; i<=N; i++)
            {
                int a=0;
                if(Z[i][k][t] > 0)
                {
                    printf("%d ->\n",i);
                    a++;
                }
                for(int j=1; j<=N; j++)
                {

                    if(Y[i][j][k][t][t] > 0)
                    {
                        printf("%d -> %d\n",i,j);
                        a++;
                    }
                }
                if(G[i][k][t] > 0 && a == 0)
                    printf("-> %d\n",i);


            }

        }
}

void Models::printSolution()
{
    for(int t=1; t<=T; t++)
        for(int i=1; i<=N; i++)
        {
            for(int k=1; k<=M; k++)
            {
                for(int u=t; u<=T; u++)
                    printf("X_%d_%d_%d_%d - %.2lf\n",i,k,t,u,X[i][k][t][u]);

                printf("Z_%d_%d_%d - %d -- ",i,k,t,Z[i][k][t]);

                printf("L_%d_%d_%d - %d\n",i,k,t,L[i][k][t]);

                //printf("G_%d_%d_%d - %d -- ",i,k,t,G[i][k][t]);

                //printf("F_0_%d_%d_%d - %.2lf\n",i,k,t,F[0][i][k][t]);
                for(int j=1; j<=N; j++)
                {

                    //printf("F_%d_%d_%d_%d - %.2lf -- ",i,j,k,t,F[i][j][k][t]);
                    for(int u=t; u<=T; u++)
                        if(Y[i][j][k][t][u] > 0)
                            printf("Y_%d_%d_%d_%d_%d - %d\n",i,j,k,t,u,Y[i][j][k][t][u]);
                }
            }
            printf("%.4lf\n",I[i][t]);
            getchar();
        }

//    for(int k=1; k<=M; k++)
//        for(int i=1; i<=N; i++)
//            printf("%d\n",Z[i][k][T+1]);

}


void Models::getValue()
{
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int k=1; k<=M; k++)
            {
                for(int u=t; u<=T; u++)
                    X[i][k][t][u] = cplex.getValue(Xx[i][k][t][u]);
                Z[i][k][t] = IloRound(cplex.getValue(Zx[i][k][t]));
                L[i][k][t] = IloRound(cplex.getValue(Lx[i][k][t]));

                G[i][k][t] = IloRound(cplex.getValue(Gx[i][k][t]));
                F[0][i][k][t] = cplex.getValue(Fx[0][i][k][t]);
                for(int j=1; j<=N; j++)
                {
                    F[i][j][k][t] = cplex.getValue(Fx[i][j][k][t]);
                    for(int u=t; u<=T; u++)
                            Y[i][j][k][t][u] = IloRound(cplex.getValue(Yx[i][j][k][t][u]));
                }
            }
            I[i][t] = IloRound(cplex.getValue(Ix[i][t]));
            SALE[i][t] = IloRound(cplex.getValue(SALEx[i][t]));
        }

//    for(int k=1; k<=M; k++)
//        for(int j=1; j<=N; j++)
//            Z[j][k][T+1] = IloRound(cplex.getValue(Zx[j][k][T+1]));
}

void Models::funcao_objetivo()
{
    double fo,aux;
    aux = 0;

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    for(int u=t; u<=T; u++)
                    {
                        aux += ins->SetupCost[i][j][k]*double(Y[i][j][k][t][u]);
                        //printf("%.6lf %.2lf %.2lf \n",fo,ins->SetupCost[i][j],Y[i][j][t]);
                    }

    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
                aux += ins->iscost[i][k]*double(Z[i][k][1]);

    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
            for(int t=1; t<=T; t++)
            {
                aux += ins->lscost[i][k]*double(L[i][k][t]);
            }

    fo = aux;
    printf("Setup cost: %.6lf \n",aux);
    aux = 0;
    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                {
                    aux += (u-t)*ins->HoldingCost[i]*X[i][k][t][u];
                    //printf("%.6lf %.2lf %.2lf \n",fo,ins->HoldingCost[i],I[i][t]);
                }

    /*initial holding cost*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            aux += ins->HoldingCost[i]*I[i][t];

    fo += aux;
    printf("Holding cost: %.6lf \n",aux);
    aux = 0;

    /*market price*/
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            aux += ins->MP[i][t]*SALE[i][t];

    fo -= aux;
    printf("Market price cost: %.6lf \n",aux);

    printf("Funcao Objetivo: %.3lf\n",fo);

}

void Models::aloca_Guimaraes()
{
    X.resize(N+1);
    Z.resize(N+1);
    L.resize(N+1);
    I.resize(N+1);
    SALE.resize(N+1);
    F.resize(N+1);
    G.resize(N+1);
    Y.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        F[i].resize(N+1);
        G[i].resize(M+1);
        X[i].resize(M+1);
        Z[i].resize(M+1);
        L[i].resize(M+1);
        I[i].resize(T+1,0.0);
        SALE[i].resize(T+1,0.0);
        Y[i].resize(N+1);
        for(int k=0;k<=M;k++)
        {
            X[i][k].resize(T+1);
            Z[i][k].resize(T+2,0);
            L[i][k].resize(T+1,0);
            G[i][k].resize(T+1,0);
            for(int t=0;t<=T;t++)
                X[i][k][t].resize(T+1,0.0);
        }
        for(int j=0;j<=N;j++)
        {
            Y[i][j].resize(M+1);
            F[i][j].resize(M+1);
            for(int k=0;k<=M;k++)
            {
                Y[i][j][k].resize(T+1);
                for(int t=0; t<=T; t++)
                    Y[i][j][k][t].resize(T+1,0);
                F[i][j][k].resize(T+1,0);
            }
        }
    }
}

void Models::libera_Guimaraes()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            for(int k=0;k<=M;k++)
            {
                for(int t=0;t<=T;t++)
                    Y[i][j][k][t].clear();
                Y[i][j][k].clear();
                F[i][j][k].clear();
            }
            Y[i][j].clear();
            F[i][j].clear();
        }
        for(int k=0;k<=M;k++)
        {
            for(int t=0;t<=T;t++)
                X[i][k][t].clear();
            X[i][k].clear();
            Z[i][k].clear();
            L[i][k].clear();
            G[i][k].clear();
        }
        X[i].clear();
        Z[i].clear();
        L[i].clear();
        F[i].clear();
        G[i].clear();
        Y[i].clear();
        I[i].clear();
        SALE[i].clear();
    }
    X.clear();
    Y.clear();
    F.clear();
    G.clear();
    Z.clear();
    L.clear();
    I.clear();
    SALE.clear();
}


void Models::set_parameters()
{
    cplex.setParam(IloCplex::Threads,4);
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.setParam(IloCplex::EpGap, 0.0000000000001);
    cplex.setParam(IloCplex::TreLim, 16000.0);
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
