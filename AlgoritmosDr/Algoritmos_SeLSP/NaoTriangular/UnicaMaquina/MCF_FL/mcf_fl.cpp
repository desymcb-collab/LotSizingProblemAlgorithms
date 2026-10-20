#include "mcf_fl.h"

void Models::RodaMCF()
{
    double inicio,fim;
    bool re = false;
    best_bound = INF;
    BestSolution = INF;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    aloca_MCF();
    inicio = time(NULL);
    env = IloEnv();
    re = SolveMCF();
    env.end();
    fim = time(NULL);
    tempo = difftime(fim,inicio);
    //funcao_objetivo();
    //printSolution();
    //printSeq();
    libera_MCF();

    cout << "Best Solution: " << BestSolution << " --- Best Bound: " << best_bound << "\tTempo: " << tempo << "\n";
        //cout << "\nBest Solution: " << BestSolution << "\tBest Bound: " << best_bound << "\tTempo: " << time << "\n";

    ///Exit file with cplex results

    sprintf(ins->pathname,"resultados/MCF_FL.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf %.4lf\n",BestSolution,best_bound,tempo);
    fclose(saida);
}

bool Models::SolveMCF()
{
    bool re = false;
    IloModel model(env);
    MCF_Model(model);                                         //Constroi modelo classico

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

void Models::MCF_Model(IloModel model)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    char name[20];

    double menor = INF;
    for(int i=1; i<=N; i++)
        if(ins->ProdTime[i] < menor)
            menor = ins->ProdTime[i];

    Xx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Xx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Xx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Xx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
                sprintf(name,"Xx%d_%d_%d",i,j,t);
                Xx[i][j][t].setName(name);
            }
        }
    }

    X_ini = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        X_ini[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            X_ini[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"X_ini%d_%d",i,t);
            X_ini[i][t].setName(name);
        }
    }

    X_co = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        X_co[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            X_co[i][t] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
            sprintf(name,"X_co%d_%d",i,t);
            X_co[i][t].setName(name);
        }
    }

    R = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        R[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"R%d",t);
        R[t].setName(name);
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
        Zx[i] = IloNumVarArray(env,T+2);
        for(int t=0; t<=T+1; t++)
        {
            Zx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Zx%d_%d",i,t);
            Zx[i][t].setName(name);
        }
    }

    Fx = Num4D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        Fx[k] = Num3D(env,N+1);
        for(int i=0; i<=N; i++)
        {
            Fx[k][i] = Num2D(env,N+1);
            for(int j=0; j<=N; j++)
            {
                Fx[k][i][j] = IloNumVarArray(env,T+1);
                for(int t=0; t<=T; t++)
                {
                    Fx[k][i][j][t] = IloNumVar(env,0,1,ILOINT);
                    sprintf(name,"Fx%d_%d_%d_%d",k,i,j,t);
                    Fx[k][i][j][t].setName(name);
                }
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
            for(int u=t; u<=T; u++)
                obj1 += (u-t)*ins->HoldingCost[i]*Xx[i][t][u];

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
        for(int u=1; u<=T; u++)
        {
            for(int t=1; t<=u; t++)
                expr1[i][u] += Xx[i][t][u];
            constraints.add(expr1[i][u] == ins->Demand[i][u]);
        }

    expr2 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr2[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int k=1; k<=N; k++)
        {
            for(int u=t; u<=T; u++)
                expr2[t] += ins->ProdTime[k]*Xx[k][t][u];
            for(int l=1; l<=N; l++)
                expr2[t] += ins->SetupTime[l][k]*Yx[l][k][t];
        }
        constraints.add(expr2[t] <= ins->Capacity[t]);
    }

    expr3 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr3[i] = Exp2D(env,T+1);
        for(int t=0; t<=T; t++)
        {
            expr3[i][t] = IloExprArray(env,T+1);
            for(int u=0; u<=T; u++)
                expr3[i][t][u] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int u=t; u<=T; u++)
            {
                int cap;
                cap = int((ins->Capacity[t])/menor);
                if(cap > ins->Demand[i][u])
                    cap = ins->Demand[i][u];
                expr3[i][t][u] += Xx[i][t][u];
                for(int j=1; j<=N; j++)
                    expr3[i][t][u] -= ins->MaxProd[i][t]*Yx[j][i][t];
                expr3[i][t][u] -= ins->MaxProd[i][t]*Zx[i][t];
                constraints.add(expr3[i][t][u] <= 0);
            }

    expr13 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr13[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                expr13[i][t] += Xx[i][t][u];
            for(int j=1; j<=N; j++)
                expr13[i][t] -= ins->minl[i]*Yx[j][i][t];
            expr13[i][t] -= ins->minl[i]*Zx[i][t];
            constraints.add(expr13[i][t] >= 0);
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
    }

    expr6 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr6[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[i][t] = IloExpr(env);
    }

    for(int k=1; k<=N; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
                expr6[k][t] += Fx[k][0][j][t];

            expr6[k][t] -= Gx[k][t];
            constraints.add(expr6[k][t] == 0);
        }

    expr7 = Exp3D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        expr7[k] = Exp2D(env,N+1);
        for(int i=0; i<=N; i++)
        {
            expr7[k][i] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr7[k][i][t] = IloExpr(env);
        }
    }

    for(int k=1; k<=N; k++)
        for(int i=1; i<=N; i++)
            if(i!=k)
                for(int t=1; t<=T; t++)
                {
                    expr7[k][i][t] += Fx[k][0][i][t];
                    for(int j=1; j<=N; j++)
                    {
                        if(j!=i)
                        {
                            if(j!=k)
                                expr7[k][i][t] += Fx[k][j][i][t];
                            expr7[k][i][t] -= Fx[k][i][j][t];
                        }
                    }
                    constraints.add(expr7[k][i][t] == 0);
                }

    expr12 = Exp2D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        expr12[k] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr12[k][t] = IloExpr(env);
    }

    for(int k=1; k<=N; k++)
        for(int t=1; t<=T; t++)
        {
            for(int j=0; j<=N; j++)
                expr12[k][t] += Fx[k][j][k][t];
            expr12[k][t] -= Gx[k][t];

            constraints.add(expr12[k][t] == 0);
        }


    expr8 = Exp3D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        expr8[k] = Exp2D(env,N+1);
        for(int i=0; i<=N; i++)
        {
            expr8[k][i] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[k][i][t] = IloExpr(env);
        }
    }

    for(int k=1; k<=N; k++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                expr8[k][j][t] += Fx[k][0][j][t];
                expr8[k][j][t] -= Zx[j][t];
                constraints.add(expr8[k][j][t] <= 0);
            }

    expr9 = Exp4D(env,N+1);
    for(int k=0; k<=N; k++)
    {
        expr9[k] = Exp3D(env,N+1);
        for(int i=0; i<=N; i++)
        {
            expr9[k][i] = Exp2D(env,N+1);
            for(int j=0; j<=N; j++)
            {
                expr9[k][i][j] = IloExprArray(env,T+1);
                for(int t=0; t<=T; t++)
                    expr9[k][i][j][t] = IloExpr(env);
            }
        }
    }

    for(int k=1; k<=N; k++)
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
                for(int t=1; t<=T; t++)
                {
                    expr9[k][i][j][t] += Fx[k][i][j][t];
                    expr9[k][i][j][t] -= Yx[i][j][t];
                    constraints.add(expr9[k][i][j][t] <= 0);
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
            expr11[i][t] -= BigM*Gx[i][t];
            constraints.add(expr11[i][t] <= 0);
        }

    }

    /*expr14 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr14[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                expr14[i][t] += Xx[i][t][u];
            expr14[i][t] -= X_co[i][t];
            expr14[i][t] -= X_ini[i][t];
            constraints.add(expr14[i][t] == 0);
        }

    expr15 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr15[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr15[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr15[i][t] += X_co[i][t];
            expr15[i][t] -= BigM*Z[i][t];
            constraints.add(expr15[i][t] <= 0);
        }

    expr16 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr16[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr16[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr16[i][t] += X_ini[i][t];
            for(int j=1; j<=N; j++)
                expr16[i][t] -= ins->minl[i]*Yx[j][i][t];

            expr16[i][t] += ins->minl[i]*Zx[i][t+1];
            constraints.add(expr16[i][t] >= 0);
        }

    expr17 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr17[i] = Exp2D(env,T+1);
        for(int j=0; j<=T; j++)
        {
            expr17[i][j] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr17[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int s=t+1; s<=T; s++)
            {
                expr17[i][t][s] += X_ini[i][t];
                for(int k=t+1; k<=s; k++)
                    expr17[i][t][s] += X_co[i][k];
                for(int j=1; j<=N; j++)
                    expr17[i][t][s] -= ins->minl[i]*Yx[j][i][t];
                for(int k=t+1; k<=s-1; k++)
                    expr17[i][t][s] += BigM*R[k];
                expr17[i][t][s] -= BigM*R[s];
                constraints.add(expr17[i][t][s] >= -BigM);
            }*/

    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
        {
            for(int k=1; k<=N; k++)
            {
                Fx[k][i][i][t].setBounds(0,0);
                Fx[k][k][i][t].setBounds(0,0);
            }
            Yx[i][i][t].setBounds(0,0);
        }
    }

    cplex = IloCplex(model);
    //cplex.exportModel("lp/SM_NTS.lp");
}

void Models::CleanGM()
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
    /*expr14.end();
    expr15.end();
    expr16.end();
    expr17.end();*/
    cplex.end();
}

///Busca em largura dos elementos///

bool Models::bfs(int st, int lt, int parent[])
{
    bool R[N+2];
    memset(R,0,sizeof(R));

    queue <int> NR;
    NR.push(st);
    R[st] = true;
    parent[st] = -1;

    while(!NR.empty())
    {
        int u = NR.front();
        NR.pop();

        for(int i=1;i<=N+1;i++)
            if(R[i]==false && NetG[u][i] > zero)
            {
                NR.push(i);
                parent[i] = u;
                R[i] = true;
            }
    }
    return (R[lt] == true);
}

///Busca em profundidade dos elementos visitados a partir de um noh raiz///

void Models::dfs(int st, bool visitedS[], bool visited[])
{
    visitedS[st] = true;
    visited[st] = true;
    for(int i=1;i<=N+1;i++)
        if(NetG[st][i] > zero && !visitedS[i])
            dfs(i,visitedS,visited);
}


bool Models::Algoritmo_FordFulkerson(IloModel model)
{
    double flow = 0.0, mincut;
    bool ro = true;
    int st,z;

    NetG.resize(N+2);
    for(int i=0;i<=N+1;i++)
        NetG[i].resize(N+2,0.0);
    z = N+1;

    for(int t=1;t<=T;t++)
    {
        bool visited[N+2], visit[N+2];
        memset(visit, false, sizeof(visit));
        for(int s=1;s<=N;s++)
        {
                NetworkG(t);
                int parent[N+1];
                memset(parent,0,sizeof(parent));
                mincut=0.0;
                //printf("s-z:%d-%d\n",s,z);
                ///Print residual graph
//                for(int i=1;i<=N;i++)
//                    for(int j=1;j<=N;j++)
//                    {
//                        if(NetG[i][j]>zero)
//                            cout << i << " - " << j << " : " << NetG[i][j] << endl;
//                    }

                while(bfs(s,z,parent))
                {
                    int u;
                    double path = INF*1.0;
                    for(int l=z;l!=s;l=parent[l])
                    {
                        u = parent[l];
                        path = min(path,NetG[u][l]);
                    }
                    for(int l=z;l!=s;l=parent[l])
                    {
                        u = parent[l];
                        NetG[u][l] -= path;
                        NetG[l][u] += path;
                    }

                    mincut += path;
                }
                //printf("\n");
                ///Print residual graph
//                for(int i=1;i<=N;i++)
//                    for(int j=1;j<=N;j++)
//                    {
//                        if(NetG[i][j]>zero)
//                            cout << i << " - " << j << " : " << NetG[i][j] << endl;
//                    }

                int conj_S[N+1], sizeS = 1;
                bool visitedS[N+1],visitedZ[N+1];
                memset(visitedS, false, sizeof(visitedS));

                dfs(s,visitedS,visit);

                for (int i = 1; i <= N; i++)
                    if(visitedS[i])
                    {
                        conj_S[i] = 1;
                        sizeS++;
                    }
                    else
                        conj_S[i] = 2;

                /*for(int i=1;i<=N;i++)
                {
                    if(conj_S[i]==1)
                    {
                        soma += double(Z[i][t+1]);
                        for(int j=1;j<=N;j++)
                            if(conj_S[j]==2)
                                soma += double(Y[i][j][t]);
                    }
                }

                if(soma != mincut)
                {
                    printf("Corte diferente (%d,%d) %.6lf x %.6lf! ERRO\n",s,t,soma,mincut);
                    printErro(t);
                    for (int i = 1; i <= N; i++)
                        cout << i << " - " << conj_S[i] << endl;
                    getchar();
                }*/

                if(sizeS > 0)
                {
                    mincut = 0.0;
                    for(int i=1;i<=N;i++)
                    {
                        if(conj_S[i]==1)
                        {
                            mincut += double(G[i][t]);
                            for(int j=1;j<=N;j++)
                                if(conj_S[j]==1)
                                    mincut -= double(Y[i][j][t]);
                        }
                    }
                    //printf("Soma de %.6lf\n",mincut);
                    for(int k=1;k<=N;k++)
                    {
                        double soma = 0.0;
                        if(conj_S[k]==1)
                        {
                            soma -= (1.0/double(BigM))*double(G[k][t]);
                            soma += mincut;
                            if(soma < -zero)
                            {
                                ro = false;
                                ///Insert new constraint to the model
                                addMinCutConstraints(model,conj_S,k,t,mincut);
                                //printf("Inequacao invalida: %d,%d,%d - %.7lf. (Almada-Lobo)\n",s,k,t,soma);
                                //getchar();
                                return ro;
                            }
                        }
                    }
                    //getchar();
                    //printf("Inequacao valida: %d,%d - %.2lf. (Almada-Lobo)\n",s,t,soma);
                }


        }
    }

    for(int i=0;i<=N;i++)
        NetG[i].clear();
    NetG.clear();

    return ro;
}

void Models::NetworkG(int t)
{
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            NetG[i][j] = 0.0;
            NetG[i][j] = double(Y[i][j][t]);
        }
        NetG[i][N+1] = 0.0;
        NetG[i][N+1] += double(Z[i][t+1]);
    }
}

void Models::addMinCutConstraints(IloModel model, int conj_S[], int k, int t, double mincut)
{
    expr18 = IloExpr(env);
    int aux = 0;
    for(int i=1;i<=N;i++)
    {
        if(conj_S[i]==1)
        {
            expr18 += Gx[i][t];
            for(int j=1;j<=N;j++)
                if(conj_S[j]==1)
                    expr18 -= Yx[i][j][t];
        }
    }

    expr18 -= (1.0/double(BigM))*double(G[k][t]);

    model.add(expr18 >= 0);

    expr18.end();

}


void Models::printSeq()
{
    for(int t=1; t<=T; t++)
    {
        bool visited=false;
        printf("\n");
        for(int i=1; i<=N; i++)
        {
            int a=0;
            if(Z[i][t] > 0)
            {
                printf("%d ->\n",i);
                a++;
            }
            for(int j=1; j<=N; j++)
            {

                if(Y[i][j][t] > 0)
                {
                    printf("%d -> %d\n",i,j);
                    a++;
                }
            }
            if(G[i][t] > 0 && a == 0)
                printf("-> %d\n",i);


        }

    }
}

void Models::printSolution()
{
    for(int t=1; t<=T; t++)
        for(int i=1; i<=N; i++)
        {

            for(int u=t; u<=T; u++)
            {
                printf("X_%d_%d_%d - %.2lf\n",i,t,u,X[i][t][u]);
            }

            printf("Z_%d_%d - %d -- ",i,t,Z[i][t]);

            printf("G_%d_%d - %d -- ",i,t,G[i][t]);

            for(int k=1; k<=N; k++)
                printf("F_%d_0_%d_%d - %d\n",k,i,t,F[k][0][i][t]);
            for(int j=1; j<=N; j++)
            {
                for(int k=1; k<=N; k++)
                    printf("F_%d_%d_%d_%d - %d \n ",k,i,j,t,F[k][i][j][t]);

                //printf("Y_%d_%d_%d - %d\n",i,j,t,Y[i][j][t]);
            }
            //getchar();
        }
        for(int i=1; i<=N; i++)
            printf("%d\n",Z[i][T+1]);

}


void Models::getValue()
{
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
            {
                X[i][t][u] = cplex.getValue(Xx[i][t][u]);
            }
            Z[i][t] = IloRound(cplex.getValue(Zx[i][t]));
            G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
            for(int k=1; k<=N; k++)
                F[k][0][i][t] = IloRound(cplex.getValue(Fx[k][0][i][t]));
            for(int j=1; j<=N; j++)
            {
                for(int k=1; k<=N; k++)
                    if(i!=j && k!=i)
                        F[k][i][j][t] = IloRound(cplex.getValue(Fx[k][i][j][t]));
                Y[i][j][t] = IloRound(cplex.getValue(Yx[i][j][t]));
            }
        }

    for(int j=1; j<=N; j++)
        Z[j][T+1] = IloRound(cplex.getValue(Zx[j][T+1]));
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
            for(int u=t; u<=T; u++)
                {
                    fo += (u-1)*ins->HoldingCost[i]*X[i][t][u];
                    //printf("%.6lf %.2lf %.2lf \n",fo,ins->HoldingCost[i],I[i][t]);
                }

    printf("Funcao Objetivo: %.2lf\n",fo);
}

void Models::aloca_MCF()
{
    X.resize(N+1);
    Z.resize(N+1);
    F.resize(N+1);
    G.resize(N+1);
    Y.resize(N+1);
    for(int k=0; k<=N; k++)
    {
        F[k].resize(N+1);
        for(int i=0;i<=N;i++)
        {
            F[k][i].resize(N+1);
            for(int j=0;j<=N;j++)
                F[k][i][j].resize(T+1,0);
        }
    }
    for(int i=0;i<=N;i++)
    {
        X[i].resize(N+1);
        G[i].resize(T+1,0);
        Z[i].resize(T+2,0);
        Y[i].resize(N+1);
        for(int j=0;j<=N;j++)
            X[i][j].resize(T+1,0);
        for(int j=0;j<=N;j++)
            Y[i][j].resize(T+1,0);

    }
}

void Models::libera_MCF()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
            Y[i][j].clear();
        for(int j=0;j<=N;j++)
            X[i][j].clear();
        X[i].clear();
        Z[i].clear();
        G[i].clear();
        Y[i].clear();
    }
    for(int k=0; k<=N; k++)
    {
        for(int i=0;i<=N;i++)
        {
            for(int j=0;j<=N;j++)
                F[k][i][j].clear();
            F[k][i].clear();
        }
        F[k].clear();
    }
    X.clear();
    Y.clear();
    F.clear();
    G.clear();
    Z.clear();
}


void Models::set_parameters()
{
    cplex.setParam(IloCplex::Threads,4);
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.setParam(IloCplex::EpGap, 0.0000001);
    cplex.setParam(IloCplex::TreLim, 4000.0);
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
