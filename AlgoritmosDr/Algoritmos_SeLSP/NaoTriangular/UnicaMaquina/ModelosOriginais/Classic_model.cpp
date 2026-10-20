#include "Classic_model.h"
ILOSTLBEGIN

void CM::RodaClassic()
{
    double inicio,fim;
    bool re = false;
    best_bound = INF;
    BestSolution = INF;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    aloca_Classic();
    inicio = time(NULL);
    //env = IloEnv();
    SolveClassic();
    env.end();

    fim = time(NULL);
    tempo = difftime(fim,inicio);
    //funcao_objetivo();

//    for(int t=1; t<=T; t++)
//    {
//        //printf("\nG %d\n",G[t]);
//        for(int i=0; i<=N; i++)
//        {
//            printf("\n(i,t)%d %d\n",i,t);
//            printf("(X,Z,I)%.2lf %d %.2lf\n",X[i][t],Z[i][t],I[i][t]);
//            for(int j=1; j<=N; j++)
//                printf("\n(i,j,t) %d %d %d - Y: %d",i,j,t,Y[i][j][t]);
//            printf("\n");
//        }
//    }

    libera_Classic();

    cout << "\nBest Solution: " << BestSolution << "\tTempo: " << tempo << "\n";

    ///Exit file with cplex results

    /*sprintf(ins->pathname,"resultados/Classic_SingleMachine.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf\n",BestSolution,tempo);
    fclose(saida);*/
}


///Busca em largura dos elementos///

bool CM::bfs(int st, int lt, int parent[])
{
    bool R[N+1];
    memset(R,0,sizeof(R));

    queue <int> NR;
    NR.push(st);
    R[st] = true;
    parent[st] = -1;

    while(!NR.empty())
    {
        int u = NR.front();
        NR.pop();

        for(int i=1;i<=N;i++)
            if(R[i]==false && WG[u][i] > 0)
            {
                NR.push(i);
                parent[i] = u;
                R[i] = true;
            }
    }
    return (R[lt] == true);
}

///Busca em profundidade dos elementos visitados a partir de um noh raiz///

void CM::dfs(int st, bool visited[])
{
    visited[st] = true;
    for(int i=1;i<=N;i++)
        if(WG[st][i] > 0 && !visited[i])
            dfs(i,visited);
}

int CM::stZ(int t)
{
    int first=0;
    int st=-1;
    for(int i=1;i<=N;i++)
        if(first < Z[i][t])
        {
            first = Z[i][t];
            st = i;
        }
    return st;
}

int CM::dfsY(int lt)
{
    for(int i=1;i<=N;i++)
        if(cY[lt][i] > 0)
            lt = dfsY(i);
    return lt;
}

bool CM::Algoritmo_FordFulkerson(IloModel model)
{
    double flow = 0.0;
    bool ro = true;

    cY.resize(N+1);
    WG.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        cY[i].resize(N+1,0);
        WG[i].resize(N+1,0);
    }

    for(int t=1;t<=T;t++)
    {
        int parent[N+1],st=0,lt=0;
        memset(parent,0,sizeof(parent));
        copiaY(t);

        ///Find the first item produced in period t
        st = stZ(t);

        ///Find the last item produced in period t
        lt = dfsY(st);

        printf("\nPrimeiro: %d --- Ultimo: %d\n",st,lt);

        WeightedGraph(t); ///Grafo com pesos
        ///Print residual graph

        if((st!=lt) && Y[st][lt][t]!=1)
        {
            while(bfs(st,lt,parent))
            {
                int u;
                int path = INF;
                for(int i=lt;i!=st;i=parent[i])
                {
                    u = parent[i];
                    path = min(path,WG[u][i]);
                }
                for(int i=lt;i!=st;i=parent[i])
                {
                    u = parent[i];
                    WG[u][i] -= path;
                    WG[i][u] += path;
                }
            }
            ///Print residual graph
//            for(int i=1;i<=N;i++)
//            {
//                for(int j=1;j<=N;j++)
//                {
//                    printf("%d\t",WG[i][j]);
//                    //getchar();
//                }
//                printf("\n");
//            }
        }

        ///Find sets S and N\S

        int conj_S[N+1], sizeS = 0;
        bool visited[N+1];
        memset(visited, false, sizeof(visited));
        dfs(st,visited);

        for(int i=1;i<=N;i++)
            if(visited[i])
            {
               conj_S[i] = 1;
               sizeS += 1;
//               printf("existe caminho entre %d e %d\n",st,i);
            }else
                conj_S[i] = 0;

        printf("Size %d\n",sizeS);

        if(sizeS > 0)
        {
            ///Verify if constraint is satisfied

            for(int k=1;k<=N;k++)
            {
                if(conj_S[k]==1)
                {
                    int soma = 0;

                    for(int i=1; i<=N; i++)
                    {
                        if(conj_S[i]==1)
                        {
                            soma += Z[i][t+1];
                            for(int j=1; j<=N; j++)
                                if(conj_S[j]==0)
                                    soma += Y[i][j][t];
                        }
                    }
                    for(int j=1; j<=N; j++)
                        soma -= Y[j][k][t];

                    printf("Soma de %d\n",soma);

                    if(soma < 0)
                    {
                        ro = false;
                        getchar();
                        ///Insert new constraint to the model
                        insereRestricoes1(model,conj_S,k,t);
                        printf("Nao satisfaz restricoes. (Almada-Lobo)\n");
                        getchar();
                    }
                }
            }
        }
        //getchar();
    }

    for(int i=0;i<=N;i++)
    {
        cY[i].clear();
        WG[i].clear();
    }

    cY.clear();
    WG.clear();

    return ro;
}

bool CM::Algoritmo_Barany(IloModel model)
{
    bool ro = true;
    int conj_S[T+1][T+1];

    for(int l=1; l<=T; l++)
    {
        for(int t=1; t<=l; t++)
        {
            for(int i=0; i<=N; i++)
            {
                int yz=0;

                double soma=0.0;
                for(int u=t; u<=l; u++)
                    soma += ins->Demand[i][u];

                for(int j=0; j<=N; j++)
                    yz += Y[j][i][t];
                yz += Z[i][t];
                soma = double(soma*yz);

                if(X[i][t] <= soma)
                    conj_S[l][t] = 1;
                else
                    conj_S[l][t] = 0;
            }
        }
    }

    ///Verify if constraint is satisfied
    for(int l=1; l<=T; l++)
        for(int i=0; i<=N; i++)
        {
            double somaT=0.0,Tdemand=0.0;
            int somaX=0;
            for(int t=1; t<=l; t++)
            {
                Tdemand += ins->Demand[i][t];
                if(conj_S[l][t] = 0)
                {
                    double soma=0.0;
                    int yz = 0;
                    for(int u=t; u<=l; u++)
                        soma += ins->Demand[i][u];
                    for(int j=0; j<=N; j++)
                        yz += Y[j][i][t];
                    yz += Z[i][t];
                    somaT = double(soma*yz);
                }else
                    somaX += X[i][t];
            }
            somaT = double(somaT + somaX - Tdemand);
            if(somaT < 0.0)
            {
                ro = false;
                printf("Nao satisfaz restricoes. (metodo do Barani)\n");
                ///Insert new constraint to the model
                insereRestricoes2(model,conj_S[l],l,i);
            }
        }

    return ro;
}

void CM::SolveClassic()
{
    bool ro = false, re = false;

    IloModel model(env);
    Classic_Model(model);                                         //Constroi modelo classico
    set_parameters();

    ///Get solution of the current problem
//    while(!ro && re)
//    {
        re = cplex.solve();
        BestSolution = re?cplex.getObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF

        if(re)
        {
            for(int t=1; t<=T; t++)
            {
//                G[t] = cplex.getValue(Gx[t]);
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
//
////            if(Algoritmo_FordFulkerson(model))// && Algoritmo_Barany(model))
////                ro = true;

            IloModel relax(env);
            relax.add(model);
            for(int i=0; i<=N; i++)
            {
                relax.add(IloConversion(env, Zx[i], ILOFLOAT));
                for(int j=0; j<=N; j++)
                    relax.add(IloConversion(env, Yx[i][j], ILOFLOAT));
            }
            relax.add(IloConversion(env, Gx, ILOFLOAT));

            IloCplex cplexRelaxed(relax);
            cplexRelaxed.exportModel("lp/relaxedProblem.lp");

        if(cplexRelaxed.solve())
            cplexRelaxed.end();

            CleanCM();
            model.end();
        }
//    }

//    if(re)
//    {
//
//    }
}

void CM::Classic_Model(IloModel model)
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

    Gx = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        Gx[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"Gx%d",t);
        Gx[t].setName(name);
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
    {

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



    expr13 = Exp3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr13[i] = Exp2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            expr13[i][j] = IloExprArray(env,N+1);
            for(int t=0; t<=T; t++)
                expr13[i][j][t] = IloExpr(env);
        }
    }

    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=N; j++)
        {
            if(i!=j)
            {
                for(int t=1; t<=T; t++)
                {
                    expr13[i][j][t] += Vx[i][t];
                    expr13[i][j][t] += N*Yx[i][j][t];
                    expr13[i][j][t] -= N*Zx[i][t];
                    expr13[i][j][t] -= Vx[j][t];
                    constraints.add(expr13[i][j][t] <= (N-1));
                }
            }
        }
    }


    /*expr14 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr14[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr14[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
                expr14[i][t] += Yx[j][i][t];
            constraints.add(expr14[i][t] <= 1);
        }
    }*/

    expr6 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr6[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr6[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            expr6[i][t] += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr6[i][t] -= Yx[i][j][t];
            }
            expr6[i][t] -= Gx[t];
            constraints.add(expr6[i][t] <= 0);
        }
    }

    expr7 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr7[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr7[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                expr7[i][t] -= Yx[j][i][t-1];
            }
            expr7[i][t] += Zx[i][t];
            expr7[i][t] -= Gx[t-1];
            constraints.add(expr7[i][t] <= 0);
        }
    }

    expr8 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr8[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
                expr8[t] += Yx[i][j][t];
        expr8[t] += Gx[t];
        constraints.add(expr8[t] >= 1);
    }

    expr9 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr9[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
                expr9[t] += Yx[i][j][t]/double(N);
        expr9[t] += Gx[t];
        constraints.add(expr9[t] <= 1);
    }


    model.add(constraints);

    for(int t=1; t<=T; t++)
        Gx[t].setBounds(0,0);
    Zx[1][1].setBounds(1,1);
    for(int i=1; i<=N; i++)
    {
        Ix[i][0].setBounds(0,0);
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
    }

    cplex = IloCplex(model);
    cplex.exportModel("lp/SM_NTS_CM.lp");
}

void CM::insereRestricoes1(IloModel model, int conj_S[], int k, int t)
{
    expr10 = IloExpr(env);

    for(int i=1; i<=N && conj_S[i]==1; i++)
    {
        expr10 += Zx[i][t+1];
        for(int j=1; j<=N && conj_S[j]==0; j++)
            expr10 += Yx[i][j][t];
    }

    for(int j=1; j<=N && conj_S[j]==0; j++)
        expr10 -= Yx[j][k][t];

    constraints.add(expr10 >= 0);
    model.add(constraints);

    cplex.exportModel("lp/SM_NTS_CM.lp");

    expr10.end();

}

void CM::insereRestricoes2(IloModel model, int conj_S[], int l, int i)
{
    expr10 = IloExpr(env);
    expr11 = IloExprArray(env,l+1);
    for(int t=0;t<=l;t++)
        expr11[t] = IloExpr(env);

    for(int t=1; t<=l; t++)
    {
        expr10 -= ins->Demand[i][t];
        if(conj_S[t] = 0)
        {
            double soma = 0.0;
            for(int u=t; u<=l; u++)
                soma += ins->Demand[i][u];
            for(int j=0; j<=N; j++)
                expr11[t] += Y[j][i][t];
            expr11[t] += Z[i][t];
            expr10 += soma*expr11[t];
        }else
            expr10 += X[i][t];
    }

    constraints.add(expr10 >= 0);

    model.add(constraints);
    expr10.end();
    expr11.end();
}

void CM::CleanCM()
{
    obj.end();
    constraints.end();
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
//    expr6.end();
//    expr7.end();
//    expr8.end();
//    expr9.end();
    expr13.end();
    //expr14.end();

    cplex.end();
}


void CM::funcao_objetivo()
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

void CM::aloca_Classic()
{
    X.resize(N+1);
    V.resize(N+1);
    Z.resize(N+1);
    I.resize(N+1);
    G.resize(T+1,0);
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

void CM::libera_Classic()
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
    G.clear();
}


void CM::set_parameters()
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

void CM::copiaY(int t)
{
    for(int i=1;i<=N;i++)
        for(int j=1;j<=N;j++)
            cY[i][j] = Y[i][j][t];
}

void CM::WeightedGraph(int t)
{
    for(int i=1;i<=N;i++)
        for(int j=1;j<=N;j++)
            WG[i][j] = int(ins->SetupCost[i][j])*Y[i][j][t];
}

SetList CM::inicializaList(SetList conj)
{
    conj = (SetList)malloc(sizeof(s_list));
    conj->n_sets = 0;
    conj->S = NULL;
    conj->NS = NULL;
    return conj;
}

void CM::insereSet(SetList cj)
{
    Set_S n = (Set_S)malloc(sizeof(ISet));
    if(cj->n_sets == 0)
    {
        n->n_items = 0;
        n->IS = NULL;
        cj->S = n;
        cj->n_sets = 1;
    }else
    {
        n->n_items = 0;
        n->IS = NULL;
        cj->NS = n;
        cj->n_sets += 1;
    }
}

void CM::insereItem(Set_S st,int a)
{
    tipoItem it = (tipoItem)malloc(sizeof(Item));
    if(st->n_items == 0)
    {
        it->itm = a;
        it->I_prox = NULL;
        st->IS = it;
        st->n_items += 1;
    }else
    {
        it->itm = a;
        it->I_prox = NULL;
        it->I_prox = st->IS;
        st->IS = it;
        st->n_items += 1;
    }
}

void CM::libera_conjunto(Set_S st)
{
    for(int i=1;i<=st->n_items;i++)
    {
        tipoItem it = (tipoItem)malloc(sizeof(Item));
        it = st->IS;
        st->IS = st->IS->I_prox;
        free(it);
    }
    free(st);
}

SetList CM::libera_conjuntos(SetList cj)
{

    free(cj->S->IS);
    free(cj->S);
    free(cj->NS->IS);
    free(cj->NS);
    return cj;
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
