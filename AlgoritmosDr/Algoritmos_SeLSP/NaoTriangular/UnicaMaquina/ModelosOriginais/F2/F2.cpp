#include "F2.h"
ILOSTLBEGIN

void F2::RunF2()
{
/****Variables initialization****/

    double inicio,fim;
    bool re = false;
    best_bound = -INF;
    BestSolution = INF;
    alocate_F2();

/****Solution method (CPLEX) and processing time****/

    inicio = time(NULL);
    SolveF2();
    env.end();
    fim = time(NULL);
    tempo = difftime(fim,inicio);

/********************Final results and memory dealocation************************/

    //cout << "\nBest Solution: " << BestSolution << "\tTempo: " << tempo << "\n";
    //obj_function();
    //print_solution();
    free_F2();

/***********************Exit file with cplex results*****************************/

    /*sprintf(ins->pathname,"resultados/F2_SingleMachine.txt");
    saida = fopen(ins->pathname,"a");
    if(saida == NULL) printf("Error in saida file!!!\n");
    fprintf(saida,"%.4lf %.4lf\n",BestSolution,tempo);
    fclose(saida);*/

/********************************************************************************/
}

void F2::SolveF2()
{
    bool ro = false, re = false;

    IloModel model(env);
    F2_Model(model);                                      //Constroi modelo classico
    set_parameters();

    re = cplex.solve();
    BestSolution = re?cplex.getObjValue():INF;

    ///Get solution of the current problem
    if(re)
    {
        getValue();
        while(!Algoritmo_FordFulkerson(model) && re)
        {
            //cplex.exportModel("lp/F2-FF.lp");
            if(cplex.solve())
            {
                getValue();
                BestSolution = re?max(cplex.getObjValue(),BestSolution):min(BestSolution,INF);
            }else
            {
                re = false;
                printf("Solucao infactivel!!\n");
                getchar();
            }
        }
        printf("BestSolution (Fulkerson): %.2lf\n",BestSolution);

        GInequalities(model);
        re = cplex.solve();
        BestSolution = re?max(cplex.getObjValue(),BestSolution):min(BestSolution,INF);
        printf("BestSolution (GInequalities): %.2lf\n",BestSolution);

        model.remove(CapacityConstraints);
        //model.remove(FulkersonsConstraints);
        //model.remove(constraints);
        //BaranysInequalities(model);
        //cplex.exportModel("lp/F2_Barany.lp");

        re = cplex.solve();
        getValue();
        while(!Algoritmo_Barany(model) && re)
        {
            //cplex.exportModel("lp/F2_Barany.lp");
            re = cplex.solve();
            if(re)
            {
                getValue();
            }else
            {
                re = false;
                printf("Solucao infactivel!!\n");
                getchar();
            }
        }

        //model.add(constraints);
        model.add(CapacityConstraints);
        //model.remove(BaranysConstraints);

        re = cplex.solve();
        BestSolution = re?max(cplex.getObjValue(),BestSolution):min(BestSolution,INF);
        printf("BestSolution (Barany): %.2lf\n",BestSolution);

        model.end();
        //print_solution();
    }
}

void F2::F2_Model(IloModel model)
{
    ///CPLEX variables inicialization
    constraints = IloRangeArray(env);
    CapacityConstraints = IloRangeArray(env);
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

    Yx = Num3D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        Yx[i] = Num2D(env,N+1);
        for(int j=0; j<=N; j++)
        {
            Yx[i][j] = IloNumVarArray(env,T+1);
            for(int t=0; t<=T; t++)
            {
                Yx[i][j][t] = IloNumVar(env,0,1,ILOFLOAT);
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
            Zx[i][t] = IloNumVar(env,0,1,ILOFLOAT);
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
        Gx[t] = IloNumVar(env,0,1,ILOFLOAT);
        sprintf(name,"Gx%d",t);
        Gx[t].setName(name);
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
    expr1.end();

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
        CapacityConstraints.add(expr2[t] <= ins->Capacity[t]);
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
            double cap;
            cap = ins->Capacity[t]/ins->ProdTime[i];
            if(cap > ins->MaxProd[i][t])
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
            for(int j=1; j<=N; j++)
                expr10[i][t] += Yx[j][i][t];
            constraints.add(expr10[i][t] <= 1);
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
                expr11[i][t] += Yx[i][j][t];
            constraints.add(expr11[i][t] <= 1);
        }
    }

    model.add(CapacityConstraints);
    model.add(constraints);

    Zx[1][1].setBounds(1,1);
    for(int i=1; i<=N; i++)
    {
        Ix[i][0].setBounds(0,0);
        Zx[i][0].setBounds(0,0);
        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
    }

    cplex = IloCplex(model);
//    UPcplex.exportModel("lp/SM_NTS_F2.lp");
}


void F2::GInequalities(IloModel model)
{
    GConstraints = IloRangeArray(env);

    expr = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                expr[i][t] -= Yx[j][i][t-1];
            }
            expr[i][t] += Zx[i][t];
            expr[i][t] -= Gx[t-1];
            GConstraints.add(expr[i][t] <= 0);
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
        GConstraints.add(expr8[t] >= 1);
    }

    expr9 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr9[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
                expr9[t] += (Yx[i][j][t]/N);
        expr9[t] += Gx[t];
        GConstraints.add(expr9[t] <= 1);
    }

    expr12 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr12[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr12[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
        {
            expr12[i][t] += Zx[i][t];
            for(int j=1; j<=N; j++)
            {
                expr12[i][t] -= Yx[i][j][t];
            }
            expr12[i][t] -= Gx[t];
            GConstraints.add(expr12[i][t] <= 0);
        }
    }

    model.add(GConstraints);

    for(int t=1; t<=T; t++)
        Gx[t].setBounds(0,0);
}

void F2::BaranysInequalities(IloModel model)
{
    BaranysConstraints = IloRangeArray(env);

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
            double aux = 0.0;
            for(int u=1; u<=t; u++)
            {
                expr1[i][t] += Xx[i][u];
                aux += ins->Demand[i][u];
            }
            if(t<T)
                BaranysConstraints.add(expr1[i][t] >= aux);
            else
                BaranysConstraints.add(expr1[i][t] == aux);
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
            BaranysConstraints.add(expr3[i][t] <= 0);
        }
    expr3.end();

    expr4 = IloExprArray(env,T+1);
    for(int t=0; t<=T; t++)
        expr4[t] = IloExpr(env);

    for(int t=1; t<=T; t++)
    {
        for(int i=1; i<=N; i++)
            expr4[t] += Zx[i][t];
        BaranysConstraints.add(expr4[t] == 1);
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
            BaranysConstraints.add(expr5[i][t] == 0);
        }
    expr5.end();

    model.add(BaranysConstraints);

}

bool F2::Algoritmo_Barany(IloModel model)
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
                    somaT += soma*yz;
                    //cout << soma << " * " << yz << " = " << somaT << endl;
                }else
                    somaT += X[i][t];

            }
            soma = somaT - Tdemand;
            //cout << somaT << " - " << Tdemand << " = " << soma << endl;
            if(soma < -zero)
            {
                ro = false;
                printf("Inequacao invalida. (metodo do Barani)\n");
                ///Insert new constraint to the model
                addBaranysConstraints(model,conj_S[i][l],i,l);
                return ro;
            }
        }

    return ro;
}

void F2::addBaranysConstraints(IloModel model, int conj_S[], int i, int l)
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

///Busca em largura dos elementos///

bool F2::bfs(int st, int lt, int parent[])
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

void F2::dfs(int st, bool visitedS[], bool visited[])
{
    visitedS[st] = true;
    visited[st] = true;
    for(int i=1;i<=N+1;i++)
        if(NetG[st][i] > zero && !visitedS[i])
            dfs(i,visitedS,visited);
}

void F2::dfsRev(int st, bool visitedS[], bool visited[])
{
    visitedS[st] = true;
    visited[st] = true;
    for(int i=1;i<=N+1;i++)
        if(NetG[i][st] > zero && !visited[i])
            dfs(i,visitedS,visited);
}

int F2::stZ(int t)
{
    double first=0.0;
    int st=-1;
    for(int i=1;i<=N;i++)
        if(first < Z[i][t])
        {
            first = Z[i][t];
            st = i;
        }
    return st;
}

int F2::dfsY(int lt, bool visited[], int t)
{
    visited[lt] = true;
    for(int i=1;i<=N+1;i++)
        if(Y[lt][i][t] > zero && !visited[i])
            lt = dfsY(i,visited,t);

    return lt;
}

int F2::dfsG(int lt, bool visited[], int z)
{
    visited[lt] = true;
    for(int i=1;i<=N+1;i++)
        if(NetG[lt][i] > zero && !visited[i])
        {
            if(i == z)
            {
                lt = z;
                i = N+2;
            }
            else
                lt = dfsG(i,visited,z);

        }
    return lt;
}

bool F2::Algoritmo_FordFulkerson(IloModel model)
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
                int parent[N+2];
                memset(parent,0,sizeof(parent));
                mincut=0.0;
//                printf("N: %d\n",N);
                ///Print residual graph
//                for(int i=1;i<=N+1;i++)
//                {
//                    for(int j=1;j<=N+1;j++)
//                    {
//                        printf("%.2lf\t",NetG[i][j]);
//                        //if(NetG[i][j]>zero)
//                            //cout << NetG[i][j] << " " << endl;
//                    }
//                    printf("\n");
//                }
//                printf("\nCorte\n");

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
//                    for(int i=1;i<=N+1;i++)
//                    {
//                        for(int j=1;j<=N+1;j++)
//                        {
//                            printf("%.2lf\t",NetG[i][j]);
//                            //if(NetG[i][j]>zero)
//                                //cout << NetG[i][j] << " " << endl;
//                        }
//                        printf("\n");
//                    }
//                    printf("\n");
//                    getchar();
                }
//                printf("\n");
                ///Print residual graph
//                for(int i=1;i<=N+1;i++)
//                    for(int j=1;j<=N+1;j++)
//                    {
//                        if(NetG[i][j]>zero)
//                            cout << i << " - " << j << " : " << NetG[i][j] << endl;
//                    }

                int conj_S[N+1], sizeS = 0;
                bool visitedS[N+2],visitedZ[N+2];
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

//                getchar();
//                printf("\n");
//                for (int i = 1; i <= N; i++)
//                    if(conj_S[i] == 1)
//                        cout << i << " " << endl;
//                printf("\n");
//                for (int i = 1; i <= N; i++)
//                    if(conj_S[i] == 2)
//                        cout << i << " " << endl;
//                printf("\n");

                double soma = 0.0;

                for(int i=1;i<=N;i++)
                    if(conj_S[i]==1)
                    {
                        soma += Z[i][t+1];
                        for(int j=1;j<=N;j++)
                            if(conj_S[j]==2)
                                soma += Y[i][j][t];

                    }

//                if(soma != mincut)
//                {
//                    printf("Corte diferente %.6lf x %.6lf! ERRO\n",soma,mincut);
//                    //getchar();
//                }else
//                    printf("Corte correto %.6lf x %.6lf! OK\n",soma,mincut);


                if(sizeS > 0)
                    for(int k=1;k<=N;k++)
                        if(conj_S[k]==1)
                        {
                            soma = 0.0;
                            for(int j=1; j<=N; j++)
                                soma += Y[j][k][t];

//                            printf("Period %d : %.6lf - %.6lf = %.6lf\n",t,soma,mincut,soma-mincut);
                            soma -= mincut;

                            if(soma > zero)
                            {
                                ro = false;
                                ///Insert new constraint to the model
                                addMinCutConstraints(model,conj_S,k,t,mincut);
                                printf("Inequacao invalida: %d,%d - %.2lf. (Almada-Lobo)\n",k,t,soma);
                                return ro;
//                                getchar();
                            }
                        }
//                    getchar();


        }
    }

    for(int i=0;i<=N;i++)
        NetG[i].clear();
    NetG.clear();

    return ro;
}

void F2::NetworkG(int t)
{
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            NetG[i][j] = 0.0;
            NetG[i][j] = Y[i][j][t];
        }
        NetG[N+1][i] = 0.0;
        NetG[i][N+1] = Z[i][t+1];
    }
}

void F2::addMinCutConstraints(IloModel model, int conj_S[], int k, int t, double mincut)
{
    expr13 = IloExpr(env);

    for(int i=1; i<=N; i++)
    {
        if(conj_S[i]==1)
        {
            expr13 += Zx[i][t+1];
            for(int j=1; j<=N; j++)
                if(conj_S[j]==2)
                    expr13 += Yx[i][j][t];
        }
    }

    for(int j=1; j<=N ; j++)
        expr13 -= Yx[j][k][t];

    FulkersonsConstraints = IloRangeArray(env);
    FulkersonsConstraints.add(expr13 >= 0);
    model.add(FulkersonsConstraints);

    expr13.end();

}

void F2::getValue()
{
    for(int t=1; t<=T; t++)
        {
            for(int i=1; i<=N; i++)
            {
                X[i][t] = cplex.getValue(Xx[i][t]);
                Z[i][t] = cplex.getValue(Zx[i][t]);
                I[i][t] = cplex.getValue(Ix[i][t]);
                for(int j=1; j<=N; j++)
                    Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
            }
        }
        for(int i=1; i<=N; i++)
            Z[i][T+1] = cplex.getValue(Zx[i][T+1]);
}

void F2::obj_function()
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

void F2::print_solution()
{
    for(int t=1; t<=T; t++)
    {
        for(int i=0; i<=N; i++)
        {
            printf("\n(i,t)%d %d\n",i,t);
            printf("(X,Z,I)%.2lf %.2lf %.2lf\n",X[i][t],Z[i][t],I[i][t]);
            for(int j=1; j<=N; j++)
                printf("\n(i,j,t) %d %d %d - Y: %.6lf",i,j,t,Y[i][j][t]);
            printf("\n");
        }
    }
}

void F2::alocate_F2()
{
    X.resize(N+1);
    //V.resize(N+1);
    Z.resize(N+1);
    I.resize(N+1);
    Y.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        X[i].resize(T+1,0);
        //V[i].resize(T+1,0);
        Z[i].resize(T+2,0);
        I[i].resize(T+1,0);
        Y[i].resize(N+1);
        for(int j=0;j<=N;j++)
            Y[i][j].resize(T+1,0);
    }
}

void F2::free_F2()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
            Y[i][j].clear();
        //V[i].clear();
        X[i].clear();
        Z[i].clear();
        I[i].clear();
        Y[i].clear();
    }
    //V.clear();
    X.clear();
    Y.clear();
    Z.clear();
    I.clear();
}

void F2::set_parameters()
{
    cplex.setParam(IloCplex::TiLim, 3600);
    cplex.setParam(IloCplex::EpGap, 0.000001);
    cplex.setParam(IloCplex::TreLim, 4000.0);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::Threads,4);
    //cplex.setParam(IloCplex::WorkMem, 128.0);
    //cplex.setParam(IloCplex::NumericalEmphasis, 1);
    cplex.setParam(IloCplex::EpRHS,1e-9);
    cplex.setParam(IloCplex::EpOpt,1e-9);
    cplex.setParam(IloCplex::EpInt,1e-9);
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
