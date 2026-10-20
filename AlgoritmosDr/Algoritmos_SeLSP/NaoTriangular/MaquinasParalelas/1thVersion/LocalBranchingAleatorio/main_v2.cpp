#include <iostream>
#include "raiz.h"
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


ILOSTLBEGIN

typedef IloArray<IloNumVarArray> DoisD;
typedef IloArray<IloArray<IloNumVarArray> > TresD;
typedef IloArray<IloArray<IloArray <IloNumVarArray> > > QuatroD;

typedef IloArray<IloExprArray> Exp2D;
typedef IloArray<Exp2D> Exp3D;
typedef IloArray<Exp3D> Exp4D;
typedef IloArray<Exp4D> Exp5D;

long seed = 0;


int geraAleatorio(int min,int max)
{
    int r;
    // srand((unsigned)time(NULL));
    if(!seed)
    {
        time(&seed);
        srand((unsigned)seed);
    }

    r = min+rand()%(max-min+1);

    return r;

}

double resolve_subproblema(Instance *ins)
{
    int N,M,T;
    N=ins->N;M=ins->M;T=ins->T;
    char name[20];

    IloEnv env;
    IloExprArray expr1;
    Exp2D expr2;
    IloExprArray expr3;
    Exp2D expr4;
    Exp3D expr5;
    Exp2D expr6;
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
    Exp4D expr17;
//********************************************* Variáveis *********************************************//

    QuatroD x(env, N+1); // quantidade do item i produzido na máquina j no período t
    for(int i = 0; i <= N; i++)
    {
        x[i] = TresD(env,M+1);
        for(int k = 0; k <= M; k++)
        {
            x[i][k] = DoisD(env,T+1);
            for(int t = 0; t <= T; t++)
            {
                x[i][k][t] = IloNumVarArray(env,T+1);
                for(int u = 0; u <= T; u++)
                {
                    x[i][k][t][u] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"x%d_%d_%d_%d",i,k,t,u);
                    x[i][k][t][u].setName(name);
                }
            }
        }
    }

    //FILE *teste;
    //teste = fopen("relax.txt", "w");

    QuatroD   y(env, N+1); // binária: 1 se o item i é produzido na máquina j no período t e 0 caso contrário
    for(int i = 0; i <= N; i++)
    {
        y[i] = TresD(env, N+1);
        for(int j = 0; j <= N; j++)
        {
            y[i][j] = DoisD(env, M+1);
            for(int k = 0; k <= M; k++)
            {
                y[i][j][k] = IloNumVarArray(env, T+1);
                for(int t = 0; t <= T; t++)
                {
                    if(ins->relax[i][j][k][t] == -1)
                    {
                        y[i][j][k][t] = IloNumVar(env,0,1,ILOFLOAT);
                        sprintf(name,"y%d_%d_%d_%d",i,j,k,t);
                        y[i][j][k][t].setName(name);
                        //fprintf(teste,"y[%d][%d][%d][%d] RELAXADO\n", i, j, k, t);
                    }
                    else
                    {
                        y[i][j][k][t] = IloNumVar(env,0,1,ILOINT);
                        sprintf(name,"y%d_%d_%d_%d",i,j,k,t);
                        y[i][j][k][t].setName(name);
                        //fprintf(teste,"y[%d][%d][%d][%d] INTEIRO\n", i, j, k, t);
                    }
                }
            }
        }
    }

    TresD   z(env, N+1);   // binaria: 1 se o setup carry over eh realizado, o cc
    for(int i = 0; i <= N; i++)
    {
        z[i] = DoisD(env, M+1);
        for(int k = 0; k <= M; k++)
        {
            z[i][k] = IloNumVarArray(env, T+2);
            for(int t = 0; t <= T+1; t++)
            {
                if(ins->relaxz[i][k][t] == -1)
                {
                    z[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"z%d_%d_%d",i,k,t);
                    z[i][k][t].setName(name);
                    //fprintf(teste,"z[%d][%d][%d] RELAXADO\n", i, k, t);
                }
                else
                {
                    z[i][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"z%d_%d_%d",i,k,t);
                    z[i][k][t].setName(name);
                    //fprintf(teste,"z[%d][%d][%d] INTEIRO\n", i, k, t);
                }
            }
        }
    }

    //getchar();
    //fclose(teste);

    QuatroD Fx(env, N+1);     // quantidade de fluxo de setup de maquina
    for(int i = 0; i <= N; i++)
    {
        Fx[i] = TresD(env, N+1);
        for(int j = 0; j <= N; j++)
        {
            Fx[i][j] = DoisD(env, M+1);
            for(int k = 0; k <= M; k++)
            {
                Fx[i][j][k] = IloNumVarArray(env, T+1);
                for(int t = 0; t <= T; t++)
                {
                    Fx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                    sprintf(name,"Fx%d_%d_%d_%d",i,j,k,t);
                    Fx[i][j][k][t].setName(name);
                }
            }
        }
    }

    TresD   Gx(env, N+1);     // binaria: 1 se a maquina foi preparada para produzir o item i no periodo t, o cc
    for(int i=0; i<=N; i++)
    {
        Gx[i] =DoisD(env,M+1);
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

    TresD X_ini(env,N+1);
    for(int i=0; i<=N; i++)
    {
        X_ini[i] = DoisD(env,M+1);
        for(int k=0; k<=M; k++)
        {
            X_ini[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                X_ini[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_ini%d_%d_%d",i,k,t);
                X_ini[i][k][t].setName(name);
            }
        }
    }

    TresD X_co(env,N+1);
    for(int i=0; i<=N; i++)
    {
        X_co[i] = DoisD(env,M+1);
        for(int k=0; k<=M; k++)
        {
            X_co[i][k] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                X_co[i][k][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"X_co%d_%d_%d",i,k,t);
                X_co[i][k][t].setName(name);
            }
        }
    }

    DoisD R(env,M+1);
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

    IloExpr SomZero(env);
    IloExpr SomUm(env);
    IloExpr SomZeroz(env);
    IloExpr SomUmz(env);

    double soma, resul;

    if(T == 10)
    {
      ins->CLB = (10*N)/M;//6*(N/10)*ins->it;//(T)*(N/10);
      ins->CLBCo = (10*N)/M;//6*(N/10)*ins->it;//(T)*(N/10);

    }
    else
    {
      ins->CLB = (5*N)/M;//(T)*(N/10);
      ins->CLBCo = (5*N)/M;//(T)*(N/10);
    }


    double sd[N+1][T+1][T+2];

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int u=1; u<=T+1; u++)
            {
                sd[i][t][u] = 0;
            }
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
            {
                for(int l=t; l<=u; l++)
                    sd[i][t][u] += ins->Demand[i][l];
                //printf("%.2lf\t",sd[i][t][u]);
            }

        }
        //sd[i][T][T] = ins->MaxProd[i][T];
        //printf("\n");
    }

    IloExpr cost(env);
    IloModel model(env);
    IloCplex cplex;
    IloRangeArray constraints(env);
    int conta = 1;


    //********************************************* Funcao Objetivo *********************************************//

    IloExpr obj(env);//initialize expression for objective function

    /*setup cost*/
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    obj += ins->SetupCost[i][j][k]*y[i][j][k][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<T; t++)
                for(int u=t+1; u<=T; u++)
                {
                    double aux = 0.0;
                    for(int s=t+1; s<=u; s++)
                        for(int l=t; l<=s-1; l++)
                            aux += ins->Demand[i][s]*ins->HoldingCost[i];

                    obj += aux*x[i][k][t][u];
                }

    model.add(IloMinimize(env, obj));

    //********************************************* Expressões *********************************************//

    if(ins->beta < ins->it)
    {
        conta = ins->it - ins->beta;

    }

    for (int i = 1; i <= N; i++)
    {
        for(int k = 1; k <= M; k++)
        {
            for(int t = (ins->it - 1); t >= conta; t--)
            {
                for(int j = 1; j <= M; j++)
                {
                    if(ins->ZO[i][j][k][ins->vetq[t]] == 2)
                    {
                        SomZero += y[i][j][k][ins->vetq[t]];
                    }

                    if(ins->Z1[i][j][k][ins->vetq[t]] == 2)
                    {
                        SomUm += (1 - y[i][j][k][ins->vetq[t]]);
                    }
                }

                if(ins->YO[i][k][ins->vetq[t]] == 2)
                {
                    SomZeroz += z[i][k][ins->vetq[t]];
                }

                if(ins->Y1[i][k][ins->vetq[t]] == 2)
                {
                    SomUmz += (1 - z[i][k][ins->vetq[t]]);
                }
            }

        }

    }
    //********************************************* Restrições *********************************************//

    //////////////////////////////////////////////////////////////
    expr1 = IloExprArray(env,N+1);
    for(int i=0; i<=N; i++)
        expr1[i] = IloExpr(env);

    for(int i=1; i<=N; i++)
    {
        for(int u=1; u<=T; u++)
            for(int k=1; k<=M; k++)
                expr1[i] += x[i][k][1][u];
        constraints.add(expr1[i] == 1);
    }

    expr2 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr2[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr2[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=2; t<=T; t++)
        {
            for(int k=1; k<=M; k++)
            {
                for(int u=1; u<=t-1; u++)
                    expr2[i][t] += x[i][k][u][t-1];
                for(int u=t; u<=T; u++)
                    expr2[i][t] -= x[i][k][t][u];
            }
                constraints.add(expr2[i][t] == 0);
        }

    expr3 = IloExprArray(env,N+1);
    for(int i=0; i<=N; i++)
        expr3[i] = IloExpr(env);

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
            for(int k=1; k<=M; k++)
                expr3[i] += x[i][k][t][T];
            constraints.add(expr3[i] == 1);
    }

    /////////////////////////////////////////////////////////////

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
            {
                for(int u=t; u<=T; u++)
                    expr4[k][t] += ins->ProdTime[i][k]*x[i][k][t][u]*sd[i][t][u];
                for(int l=1; l<=N; l++)
                    expr4[k][t] += ins->SetupTime[l][i][k]*y[l][i][k][t];
            }
            constraints.add(expr4[k][t] <= ins->Capacity[k][t]);
        }

    ///////////////////////////////////////////////////////////////

    if(ins->it > ins->beta)
    {
        for(int i = 1; i <= N; i++) // Custo de produção + Custo de preparação
        {
            for(int k = 1; k <= M; k++)
            {
                for(int t = 1; t < (ins->it-ins->beta); t++)
                {
                    for(int j = 1; j <= N; j++)
                        model.add(y[i][j][k][ins->vetq[t]] == ins->relax[i][j][k][ins->vetq[t]]);
                    model.add(z[i][k][ins->vetq[t]] == ins->relaxz[i][k][ins->vetq[t]]);

                }
            }
        }

    }

    ////////////////////////////////////////////////////////////////

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
                for(int u=t; u<=T; u++)
                    expr5[i][k][t] += x[i][k][t][u];

                for(int j=1; j<=N; j++)
                        expr5[i][k][t] -= y[j][i][k][t];
                    expr5[i][k][t] -= z[i][k][t];
                constraints.add(expr5[i][k][t] <= 0);
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
            for(int i=1; i<=N; i++)
                expr6[k][t] += z[i][k][t];
            constraints.add(expr6[k][t] == 1);
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
                expr7[i][k][t] += z[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr7[i][k][t] += y[j][i][k][t];
                    expr7[i][k][t] -= y[i][j][k][t];
                }
                expr7[i][k][t] -= z[i][k][t+1];
                constraints.add(expr7[i][k][t] == 0);
            }

    /////////////////////////////////////////////////////////////

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
                expr10[i][k][t] -= N*z[i][k][t];
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
                    expr11[i][j][k][t] -= N*y[i][j][k][t];
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
                expr12[i][k][t] += z[i][k][t];
                for(int j=1; j<=N; j++)
                {
                    expr12[i][k][t] += y[j][i][k][t];
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
                    expr13[i][k][t] += y[j][i][k][t];
                }
                expr13[i][k][t] += z[i][k][t];
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
                    expr14[i][k][t] += sd[i][t][u]*x[i][k][t][u];
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
                expr15[i][k][t] -= BigM*z[i][k][t];
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
                    expr16[i][k][t] -= 25*y[j][i][k][t];

                expr16[i][k][t] += 25*z[i][k][t+1];
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
                        expr17[i][k][t][s] -= 25*y[j][i][k][t];
                    for(int u=t+1; u<=s-1; u++)
                        expr17[i][k][t][s] += BigM*R[k][u];
                    expr17[i][k][t][s] -= BigM*R[k][s];
                    constraints.add(expr17[i][k][t][s] >= -BigM);
                }


    /////////////////////////////////////////////////////////////////

    if (ins->it != 1)
    {
        model.add(SomZero + SomUm <= ins->CLB);
        model.add(SomZeroz + SomUmz <= ins->CLBCo);
    }

    //********************************************* Função Objetivo *********************************************/

    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        for(int k=1; k<=M; k++)
        {
            z[i][k][0].setBounds(0,0);
            if(ins->MachineProd[i][k]==0)
            {
                for(int t=1; t<=T; t++)
                {
                    z[i][k][t].setBounds(0,0);
                    Gx[i][k][t].setBounds(0,0);
                    for(int j=1; j<=N; j++)
                        y[j][i][k][t].setBounds(0,0);
                }
            }
            for(int t=1; t<=T; t++)
                y[i][i][k][t].setBounds(0,0);
        }
    }

    cplex = IloCplex(model);

    if(T == 5)
            cplex.setParam(IloCplex::TiLim, 200);
    if(T == 10)
            cplex.setParam(IloCplex::TiLim, 100);
    cplex.setParam(IloCplex::Threads,4);
    cplex.setOut(env.getNullStream());

    //********************************************* Resultados *********************************************/
    bool re = false;
    re = cplex.solve();
    //printf("Foi\n");
    if (re)
    {
        //printf("Foi2\n");
        resul = re?cplex.getObjValue():INF;
        //printf("Foi3\n");
        if(ins->bol == 1)
            ins->Solu = re?cplex.getObjValue():INF;
        //printf("Foi4\n");

        for(int i = 1; i <= N; i++)
        {
            for(int k = 1; k <= M; k++)
            {
                for(int t = 1; t <= ins->it; t++)
                {
                    for(int j = 1; j <= N; j++)
                        ins->relax[i][j][k][ins->vetq[t]] = cplex.getValue(y[i][j][k][ins->vetq[t]]);
                    ins->relaxz[i][k][ins->vetq[t]] = cplex.getValue(z[i][k][ins->vetq[t]]);

                }
            }
        }
        //printf("Foi5\n");
        cout << " Funcao Objetivo = " << cplex.getObjValue() << endl;
        //getchar();

        /*if (ins->bol == 1)
        {
            FILE *saida2;

            saida2 = fopen("saida2.txt","w");

            for(int i = 1; i <= N; i++)
            {
                //fprintf(saida2, " Produto %d ", i);
                for(int j = 1; j <= N; j++)
                {
                    fprintf(saida2, "Produto %d -> %d \n\n", i, j);

                    for(int k = 1; k <= M; k++)
                    {
                        fprintf(saida2, " Máquina %d \n", k);
                        for(int t = 1; t <= T; t++)
                        {
                            fprintf (saida2, "y[%d][%d][%d][%d] = %.2lf\n",i,j,k,t,cplex.getValue(y[i][j][k][t]));

                        }
                        fprintf (saida2, "\n");
                    }
                }

                fprintf (saida2, "\n\n");


            }

            for(int i = 1; i <= N; i++)
            {
                fprintf(saida2, " Produto %d \n\n", i);
                for(int j = 1; j <= M; j++)
                {
                    fprintf(saida2, " Máquina %d \n", j);

                    for(int t = 1; t <= T; t++)
                    {

                        fprintf (saida2, "z[%d][%d][%d] = %.2lf\n",i,j,t,cplex.getValue(z[i][j][t]));

                    }
                    fprintf (saida2, "\n");
                }

                fprintf (saida2, "\n\n");


            }
            fclose(saida2);
        }*/

        env.end();

        return resul;

    }
    else
    {
        cout << " No solution found" << endl;
        env.end();
        return -1;
    }

}

void particao(int p, Instance *ins)
{
    int N,M,T;
    N=ins->N;M=ins->M;T=ins->T;

    for(int i = 1; i <= N; i++)
    {
        for(int k = 1; k <= M; k++)
        {
            for(int j = 1; j <= N; j++)
                ins->relax[i][j][k][ins->vetq[p]] = -2; //Inteira
            ins->relaxz[i][k][ins->vetq[p]] = -2; //Inteira
        }
    }
}

void AtualizaConjuntos(Instance *ins)
{
    int N,M,T;
    N=ins->N;M=ins->M;T=ins->T;

    for(int i = 1; i <= N; i++ )
    {
        for(int k = 1; k <= M; k++)
        {
            for(int t = 1; t <= ins->it; t++)
            {
                for(int j = 1; j <= N; j++)
                {

                    if(ins->relax[i][j][k][ins->vetq[t]] == 0)
                    {
                        ins->ZO[i][j][k][ins->vetq[t]] = 2;

                    }
                    if (ins->relax[i][j][k][ins->vetq[t]] >= 1 && (ins->relax[i][j][k][ins->vetq[t]] - abs(ins->relax[i][j][k][ins->vetq[t]])==0))
                    {
                        ins->Z1[i][j][k][ins->vetq[t]] = 2;
                    }
                }

                if(ins->relaxz[i][k][ins->vetq[t]] == 0)
                {
                    ins->YO[i][k][ins->vetq[t]] = 2;

                }
                if (ins->relaxz[i][k][ins->vetq[t]] >= 1 && (ins->relaxz[i][k][ins->vetq[t]] - abs(ins->relaxz[i][k][ins->vetq[t]])==0))
                {
                    ins->Y1[i][k][ins->vetq[t]] = 2;
                }

            }
        }
    }
}

int verifica(Instance * ins, int n)
{
    int c, r = 0;

    for(int j = 1; j < n; j++)
    {
        if(ins->vetq[n] == ins->vetq[j])
        {
            r = 1;
            j = n + 1;
        }
    }

    return r;
}

void RelaxAndFix(Instance * ins)
{
    int N,M,T;
    N=ins->N;M=ins->M;T=ins->T;
    int qtdade, qtot, ini, veri, ca;
    double resp, alfa;
    double cond;
    double lim;

    lim = 0.05;

    ins->bol = 0;
    resp = 0;
    ins->totalite = 0;
    ins->Solu = 99999999;
    ins->MelhorSolu = 99999999;
    ins->tempomeio = 0;

    cond = double((abs(ins->solz - 1)));

    while(cond > lim && ins->tempomeio < 3600)
    {
        ins->bol = 0;
        resp = 0;
        ins->totalite = ins->totalite + 1;
//         printf("cond = %lf");
//         getchar();

        // srand((unsigned)time(NULL));
        //   alfa =((double)(rand())/((double)(RAND_MAX)*1)); // gera número aleatório entre 0 e 1
        if(T == 10)
            qtdade = geraAleatorio(2, T/2);
        if(T == 5)
            qtdade = geraAleatorio(1, T-2);

        /* if(T < 18)
             qtdade = ceil(alfa*(ins->c)/1.5); // qtdade = maior inteiro alfa*nbTemp (qtade vai de 1 até nbTemp)
         else
             qtdade = ceil(alfa*(T)/1.5); */


        ini = 2;
        ins->vetq[1] = 1;

        for(int t = 2; t <= T; t++)
        {
            ins->vetq[t] = -1;
        }

        if(qtdade < T)
            qtot = qtdade + 1;
        else
            qtot = qtdade;

        veri = 1;
        srand((unsigned)time(NULL));

        do
        {
            for(int i = ini; i <= qtot; i ++)
            {
                do
                {
                    ins->vetq[i] = geraAleatorio(ini,qtot); //gera aleatório entre ini e qtot
                    //printf("%d - %d\n",ins->vetq[i],i);
                    //getchar();
                    veri = verifica(ins, i); //verifica se vetq[i] já foi gerado, se sim gera outro numero
                }
                while (veri == 1);
//                printf("%d - %d\n",ins->vetq[i],i);
            }

            ini = ini + qtdade;
            qtot = qtot + qtdade;

            if(qtot > T)
            {
                qtot = T;
            }

        }
        while(ins->vetq[T] == -1);


        ins->auxp = qtdade;
        ins->auxfi = qtdade;
        ins->pp = 0;
        srand((unsigned)time(NULL));
        if(qtdade  <= T)
            ins->beta = geraAleatorio(qtdade, T-1);
        else
            ins->beta = T;//geraAleatorio(qtdade, T-1);


//         printf("qtdade = %d ", qtdade);
//         printf("beta = %d\n", ins->beta);
//        for(int t = 1; t <= T; t++ )
//             {
//               printf("v[%d] = %d   ", t, ins->vetq[t]);
//             }
//        printf("\n");
         //getchar();


        for(int i = 0; i <= N; i++)
        {
            for(int k = 0; k <= M; k++)
            {
                for(int t = 0; t <= T; t++ )
                {
                    for(int j = 0; j <= N; j++)
                    {
                        ins->ZO[i][j][k][t] = -1;
                        ins->Z1[i][j][k][t] = -1;
                        ins->relax[i][j][k][t] = -1;
                    }
                    ins->YO[i][k][t] = -1;
                    ins->Y1[i][k][t] = -1;
                    ins->relaxz[i][k][t] = -1;
                }
            }
        }

        resp = 0;
        /*gettimeofday(&(ins->tempo_meio),NULL);
        ins->tm = (double)(ins->tempo_meio).tv_usec + ((double)(ins->tempo_meio).tv_sec * (1000000.0));
        ins->tempomeio = (ins->tm - ins->ti) / 1000;
        ins->tempomeio =  ins->tempomeio/1000; */
        int t = 1;

        while(t <= T)
        {
            ins->it = t;

            if(resp != -1 && ins->tempomeio < 3600)
            {

                ins->tempomeio = 0;
                if(t == T)
                {
                    ins->bol = 1;
                }

                particao(t, ins);
                resp = resolve_subproblema(ins);
                //printf("Aqui %d!!\n",t);
                AtualizaConjuntos(ins);
                //printf("Aqui %d!!\n",t);
                t++;
                //printf("ENTREIIIII");
                //getchar();

            }
            else
            {
                t = T+1;
                printf("ite =  %d\n", ins->it);
                printf("SOLUÇÃO NÃO ENCONTRADA!");
                //getchar();

            }

            gettimeofday(&(ins->tempo_meio),NULL);
            ins->tm = (double)(ins->tempo_meio).tv_usec + ((double)(ins->tempo_meio).tv_sec * (1000000.0));
            ins->tempomeio = (ins->tm - ins->ti) / 1000;
            ins->tempomeio =  ins->tempomeio/1000;
            printf("tempo = %lf\n", ins->tempomeio);
            //getchar();

        }

        if(ins->Solu < ins->MelhorSolu && resp != -1)
            ins->MelhorSolu = ins->Solu;

        if(resp != -1)
            cond = double((abs(ins->solz - ins->Solu)))/ abs(ins->Solu);


    }

    printf("iteracoes = %d\n", ins->totalite);
    gettimeofday(&ins->tempo_fim,NULL);
    ins->tf = (double)ins->tempo_fim.tv_usec + ((double)ins->tempo_fim.tv_sec * (1000000.0));
    ins->tempo = (ins->tf - ins->ti) / 1000;
    ins->tempo = ins->tempo/1000;
    // printf("Tempo gasto em segundos %.3f\n",ins->tempo); */


}

int main()
{
    try
    {
        //********************************************* Fixos *********************************************//
        Instance ins;
        double elapsedTime;
        struct timespec start,endIt;
        int qt_period,index_items,qt,id,id2,id3,id4,id5;
        int qt_items[4];
        int p[2];
        int cd[2],prob1[2],prob2[2];

        double cl1[2];
        double cl2 = 0.5;

        qt_items[0] = 15; qt_items[1] = 20; qt_items[2] = 15; qt_items[3] = 400;
        cl1[0] = 0.6; cl1[1] = 0.8;
        cd[0] = 50; cd[1] = 100;
        prob1[0] = 60; prob1[1] = 80;
        prob2[0] = 10; prob2[1] = 20;
        p[0] = 2; p[1] = 3;

        for(id5 = 1; id5 <= 1; id5++) // 0 a 1 p
            for(id = 0; id <= 0; id++) // 0 a 1 cl1
                for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
                    for(id3 = 0; id3 <= 0; id3++) // 0 a 1 prob1
                        for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                            for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                                for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                                    for(qt = 1; qt <= 10; qt++)  // 0 a 9
                                        //for(int i = 0; i < 1; i++) // 5 execucoes por instancia
                                        {
                                            printf("\nData%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\n",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);

                                            sprintf(ins.FileInfo,"../../Instancias/MultiMachineCapVar/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);
                                            sprintf(ins.FileSaida,"resultados/Results-4.txt");
                                            //sprintf(ins.FileInfo2,"resultados/Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d_saida.dat",itensmaq[i][0],itensmaq[i][1],contador,periods[t]);

                                            ins.car = 1;
                                            //printf("%s\n",ins.FileInfo);
                                            //  getchar();

                                            ins.le_dados_SM(ins.FileInfo);
                                            //ins.le_dados(ins.FileInfo);
                                            //ins.le_dados2(ins.FileInfo2);
                                            // printf("Capacidade = %lf", ins.d[49][17]);
                                            // getchar();

                                            if(ins.car == 1)
                                            {
                                                ins.ti = ins.tf = ins.tempo = ins.tm = ins.tempomeio  = 0;

                                                gettimeofday(&ins.tempo_inicio,NULL);


                                                ins.ti = (double)ins.tempo_inicio.tv_usec + ((double)ins.tempo_inicio.tv_sec * (1000000.0));

                                                RelaxAndFix(&ins);
                                                FILE *saida;
                                                saida = fopen(ins.FileSaida, "a");

                                                fprintf(saida, "Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod\t",p[id5],qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt); //fo
                                                fprintf(saida, "%lf\t", ins.MelhorSolu); //fo
                                                fprintf(saida, "%lf\n", ins.tempo);//tempo
                                                //fprintf(saida, "%d\n", ins.totalite);//iterações

                                                fclose(saida);
                                            }

                                            ins.desaloca_variaveis();
                                        }
    }
    catch (IloException & ex)
    {
        cerr << "Error: " << ex << endl;
    }
    catch (...)
    {
        cerr << "Error" << endl;
    }


    return 0;
}
