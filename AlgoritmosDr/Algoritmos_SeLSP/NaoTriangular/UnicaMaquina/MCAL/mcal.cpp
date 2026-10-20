#include "mcal.h"
ILOSTLBEGIN

void MCAL::RunMCAL()
{
/****Variables initialization****/

    bool re = false;
    best_bound = -INF;
    BestSolution = INF;
    alocate_MCAL();

/****Solution method (CPLEX) and processing time****/

    SolveMCAL();
    env.end();

/********************Final results and memory dealocation************************/

    cout << "\nBest Solution: " << BestSolution << " --- Best Bound: " << best_bound << "\n";
    //obj_function();
    //print_solution();
    free_MCAL();

/***********************Exit file with cplex results*****************************/

    sprintf(ins->pathname,"resultados/MCAL2.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf",BestSolution,best_bound);
    fclose(saida);

/********************************************************************************/
}

void MCAL::SolveMCAL()
{
    bool ro = true, re = false;
    double tempo;
    struct timespec curT;
    IloModel model(env);
    MCAL_Model(model);                                         //Constroi modelo classico
    set_parameters();

    re = cplex.solve();
    BestSolution = re?cplex.getObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    clock_gettime(CLOCK_MONOTONIC, &curT);
    tempo = curT.tv_sec - start;
    ///Get solution of the current problem
    if(re && tempo <3300)
    {
        getValue();
        while(!Algoritmo_FordFulkerson(model) && re && tempo < 3300)
        {
            re = cplex.solve();
            if(re)
                getValue();
                //BestSolution = cplex.getObjValue();
                //printf("BestSolution: %.2lf\n",BestSolution);
            clock_gettime(CLOCK_MONOTONIC, &curT);
            tempo = curT.tv_sec - start;
        }
        if(Algoritmo_FordFulkerson(model) && re)
        {
            BestSolution = re?cplex.getObjValue():min(BestSolution,INF);
            best_bound = re?cplex.getBestObjValue():-INF;
        }
        ///Relaxed problem
        /*IloModel relax(env);
        relax.add(model);
        for(int i=0; i<=N; i++)
        {
            relax.add(IloConversion(env, Zx[i], ILOFLOAT));
            relax.add(IloConversion(env, Gx[i], ILOFLOAT));
            for(int j=0; j<=N; j++)
            {
                relax.add(IloConversion(env, Yx[i][j], ILOFLOAT));
                relax.add(IloConversion(env, Qx[i][j], ILOFLOAT));
            }
        }

        IloCplex cplexRelaxed(relax);
        re = cplexRelaxed.solve();
        best_bound = re?cplexRelaxed.getObjValue():-INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
        relax.end();*/
        model.end();
    }
}

void MCAL::MCAL_Model(IloModel model)
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

    X_ini = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        X_ini[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
            int cap = 0;
            cap = int(ins->Capacity[t]/menor);
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
            int cap = 0;
            cap = int(ins->Capacity[t]/menor);
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

    Vx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Vx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T; t++)
        {
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
                Yx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Yx%d_%d_%d",i,j,t);
                Yx[i][j][t].setName(name);
            }
        }
    }

    Qx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Qx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Qx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Qx[i][j][t] = IloNumVar(env,0,1,ILOINT);
                sprintf(name,"Qx%d_%d_%d",i,j,t);
                Qx[i][j][t].setName(name);
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

    Gx = Num2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Gx[i] = IloNumVarArray(env,T+1);
        for(int t=0; t<=T+1; t++)
        {
            Gx[i][t] = IloNumVar(env,0,1,ILOINT);
            sprintf(name,"Gx%d_%d",i,t);
            Gx[i][t].setName(name);
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
    model.add(obj);
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
            int cap;
            //cap = int((ins->Capacity[t])/menor);
            //if(cap > ins->MaxProd[i][t])
            cap = ins->MaxProd[i][t];
            expr3[i][t] += Xx[i][t];
            for(int j=1; j<=N; j++)
                expr3[i][t] -= cap*Yx[j][i][t];
            expr3[i][t] -= cap*Zx[i][t];
            constraints.add(expr3[i][t] <= 0);
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
                    for(int k=1; k<=N; k++)
                        expr6[i][j][t] -= BigM*Yx[k][i][t];
                    expr6[i][j][t] += 2*BigM*Qx[i][j][t];
                    expr6[i][j][t] -= BigM*Zx[i][t];
                    expr6[i][j][t] -= Vx[j][t];
                    constraints.add(expr6[i][j][t] <= (BigM-1));
                }

    expr7 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = Exp2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr7[i][j] = IloExprArray(env,N+1);
            for(int t=0; t<=T; t++)
                expr7[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                expr7[i][j][t] += Yx[i][j][t];
                expr7[i][j][t] -= Qx[i][j][t];
                constraints.add(expr7[i][j][t] >= 0);
            }

    expr8 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr8[i] = Exp2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr8[i][j] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr8[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
            {
                expr8[i][j][t] += Yx[i][j][t];
                expr8[i][j][t] -= BigM*Qx[i][j][t];
                constraints.add(expr8[i][j][t] <= 0);
            }

    expr9 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr9[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr9[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
                expr9[i][t] += Yx[j][i][t];
            expr9[i][t] += Zx[i][t];
            expr9[i][t] -= Gx[i][t];
            constraints.add(expr9[i][t] >= 0);
        }

    expr10 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr10[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr10[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
                expr10[i][t] += Yx[j][i][t];
            expr10[i][t] += Zx[i][t];
            expr10[i][t] -= BigM*Gx[i][t];
            constraints.add(expr10[i][t] <= 0);
        }

    expr11 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr11[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr11[i][t] += Xx[i][t];
            expr11[i][t] -= X_co[i][t];
            expr11[i][t] -= X_ini[i][t];
            constraints.add(expr11[i][t] == 0);
        }

    expr12 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr12[i][t] = IloExpr(env);
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            expr12[i][t] += X_co[i][t];
            expr12[i][t] -= BigM*Z[i][t];
            constraints.add(expr12[i][t] <= 0);
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
            expr13[i][t] += X_ini[i][t];
            for(int j=1; j<=N; j++)
                expr13[i][t] -= ins->minl[i]*Yx[j][i][t];

            expr13[i][t] += ins->minl[i]*Zx[i][t+1];
            constraints.add(expr13[i][t] >= 0);
        }

    expr14 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = Exp2D(env,T+1);
        for(int j=0; j<=T; j++)
        {
            expr14[i][j] = IloExprArray(env,T+1);
            for(int t=0; t<=T; t++)
                expr14[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int s=t+1; s<=T; s++)
            {
                expr14[i][t][s] += X_ini[i][t];
                for(int k=t+1; k<=s; k++)
                    expr14[i][t][s] += X_co[i][k];
                for(int j=1; j<=N; j++)
                    expr14[i][t][s] -= ins->minl[i]*Yx[j][i][t];
                for(int k=t+1; k<=s-1; k++)
                    expr14[i][t][s] += BigM*R[k];
                expr14[i][t][s] -= BigM*R[s];
                constraints.add(expr14[i][t][s] >= -BigM);
            }

    model.add(constraints);

    //Zx[1][1].setBounds(1,1);
    for(int i=1; i<=N; i++)
    {
        Ix[i][0].setBounds(0,0);
        Ix[i][T].setBounds(0,0);
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
    }

    cplex = IloCplex(model);
//    cplex.exportModel("lp/SM_NTS_MCAL.lp");
}


///Busca em largura dos elementos///

bool MCAL::bfs(int st, int lt, int parent[])
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

void MCAL::dfs(int st, bool visitedS[], bool visited[])
{
    visitedS[st] = true;
    visited[st] = true;
    for(int i=1;i<=N+1;i++)
        if(NetG[st][i] > zero && !visitedS[i])
            dfs(i,visitedS,visited);
}


bool MCAL::Algoritmo_FordFulkerson(IloModel model)
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

                int conj_S[N+1], sizeS = 0;
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

                if(sizeS > 2)
                {
                    mincut = 0.0;
                    for(int i=1;i<=N;i++)
                    {
                        if(conj_S[i]==1)
                        {
                            mincut += double(Z[i][t]);
                            for(int j=1;j<=N;j++)
                                if(conj_S[j]==2)
                                    mincut += double(Y[j][i][t]);
                        }
                    }
                    double soma = 0.0;
                    for(int k=1;k<=N;k++)
                        if(conj_S[k]==1)
                        {

                            soma += double(BigM*G[k][t]);
                            soma -= double(BigM);
                        }
                    soma = 1 - soma;
                    soma -= mincut;
//                    printf("Soma de %.6lf\n",soma);
                    if(soma > zero)
                    {
                        ro = false;
                        ///Insert new constraint to the model
                        addMinCutConstraints(model,conj_S,t,mincut);
                        //printf("Inequacao invalida: %d,%d - %.2lf. (Almada-Lobo)\n",s,t,soma);
                        //getchar();
                        return ro;
                    }
                    //printf("Inequacao valida: %d,%d - %.2lf. (Almada-Lobo)\n",s,t,soma);
                }
        }
    }

    for(int i=0;i<=N;i++)
        NetG[i].clear();
    NetG.clear();

    return ro;
}

void MCAL::NetworkG(int t)
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

void MCAL::addMinCutConstraints(IloModel model, int conj_S[], int t, double mincut)
{
    int aux = 0;
    expr15 = IloExpr(env);

    for(int i=1; i<=N; i++)
    {
        if(conj_S[i]==1)
        {
            expr15 += Zx[i][t];
            for(int j=1; j<=N; j++)
                if(conj_S[j]==2)
                    expr15 += Yx[j][i][t];
            aux += BigM;
            expr15 -= BigM*Gx[i][t];
        }
    }

    aux = 1 - aux;

    model.add(expr15 >= aux);

    expr15.end();

}


void MCAL::getValue()
{
    for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                X[i][t] = cplex.getValue(Xx[i][t]);
                Z[i][t] = cplex.getValue(Zx[i][t]);
                G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
                I[i][t] = cplex.getValue(Ix[i][t]);
                for(int j=1; j<=N; j++)
                {
                    Q[i][j][t] = IloRound(cplex.getValue(Qx[i][j][t]));
                    Y[i][j][t] = IloRound(cplex.getValue(Yx[i][j][t]));
                }
            }
        }
        for(int i=1; i<=N; i++)
            Z[i][T+1] = cplex.getValue(Zx[i][T+1]);
}

void MCAL::obj_function()
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

void MCAL::printErro(int t)
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

void MCAL::print_solution()
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

void MCAL::alocate_MCAL()
{
    X.resize(N+1);
    V.resize(N+1);
    Z.resize(N+1);
    G.resize(N+1);
    I.resize(N+1);
    Y.resize(N+1);
    Q.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        X[i].resize(T+1,0);
        V[i].resize(T+1,0);
        G[i].resize(T+1,0);
        Z[i].resize(T+2,0);
        I[i].resize(T+1,0);
        Y[i].resize(N+1);
        Q[i].resize(N+1);
        for(int j=0;j<=N;j++)
        {
            Y[i][j].resize(T+1,0);
            Q[i][j].resize(T+1,0);
        }
    }
}

void MCAL::free_MCAL()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            Y[i][j].clear();
            Q[i][j].clear();
        }
        V[i].clear();
        X[i].clear();
        Z[i].clear();
        G[i].clear();
        I[i].clear();
        Y[i].clear();
        Q[i].clear();
    }
    V.clear();
    X.clear();
    Y.clear();
    Q.clear();
    G.clear();
    Z.clear();
    I.clear();
}

void MCAL::set_parameters()
{
    cplex.setParam(IloCplex::Threads,4);
    cplex.setParam(IloCplex::TiLim, 300);
    //cplex.setParam(IloCplex::EpGap, 0.000001);
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
