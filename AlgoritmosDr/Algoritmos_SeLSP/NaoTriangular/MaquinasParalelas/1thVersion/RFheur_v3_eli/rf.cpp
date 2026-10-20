#include "rf.h"

void RF::RelaxAndFix()
{
    struct timespec fim;
    double elapsedTime;
    bool fo = false, re = false, lt = true;
    best_bound = -INF;
    BestSolution = INF;

    AllocateVariables();
    env = IloEnv();
    nm = 0;

    int jump = 2;
    int jump_fo = 1;
    int t=2;
    binary = 2;
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
        }
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
            if(elapsedTime < 3000)
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
                if(elapsedTime < 3000)
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
    if(lt)
        PathRelinking();

    env.end();
    FreeVariables();
    cout << "Best solution: " << BestSolution;
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
        for(int t=1; t<=T; t++)
        {
            for(int j=1; j<=N; j++)
            {
                start.add(Yx[i][j][t]);
                val.add(vet[0].ystore[i][j][t]);
            }

            for(int u=t; u<=T; u++)
            {
                start.add(Xx[i][t][u]);
                val.add(vet[0].xstore[i][t][u]);
            }
            start.add(Zx[i][t]);
            val.add(vet[0].zstore[i][t]);
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
        getFOptValue();                                       //Constroi modelo classico
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

            X_ini[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_ini%d_%d",i,t);
            X_ini[i][t].setName(name);

            X_co[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_co%d_%d",i,t);
            X_co[i][t].setName(name);

            for(int u=0; u<=T; u++)
            {
                Xx[i][t][u] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Xx%d_%d_%d",i,t,u);
                Xx[i][t][u].setName(name);
            }
        }

    }

    R = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        R[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"R%d",t);
        R[t].setName(name);
    }

    he = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        he[t] = IloNumVar(env,0,16,ILOINT);
        sprintf(name,"R%d",t);
        he[t].setName(name);
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
                sprintf(name,"Yx%d_%d_%d",i,j,t);
                Yx[i][j][t].setName(name);
                if(t>0 && (t>l || t<l-1))
                    Yx[i][j][t].setBounds(Y[i][j][t],Y[i][j][t]);

                Fx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Fx%d_%d_%d",i,j,t);
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
            sprintf(name,"Gx%d_%d",i,t);
            Gx[i][t].setName(name);
        }
    }


    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int t=1; t<=T; t++)
        obj1 += ins->ExtraCost[t]*he[t];

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
            expr1[i][t] = IloExpr(env);
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
        for(int i=1; i<=N; i++)
        {
            for(int u=t; u<=T; u++)
                expr2[t] += ins->ProdTime[i]*Xx[i][t][u];
            for(int l=1; l<=N; l++)
                expr2[t] += ins->SetupTime[l][i]*Yx[l][i][t];
        }
        expr2[t] -= he[t];
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
                for(int l=t; l<=T; l++)
                {
                    for(int j=1; j<=N; j++)
                        expr3[i][t][u] -= ins->Demand[i][l]*Yx[j][i][l];
                    expr3[i][t][u] -= ins->Demand[i][l]*Zx[i][l];
                }
                constraints.add(expr3[i][t][u] <= 0);
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

    expr11 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr11[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
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
                expr16[i][t] -= 25*Yx[j][i][t];

            expr16[i][t] += 25*Zx[i][t+1];
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
                for(int u=t+1; u<=s; u++)
                    expr17[i][t][s] += X_co[i][u];
                for(int j=1; j<=N; j++)
                    expr17[i][t][s] -= 25*Yx[j][i][t];
                for(int u=t+1; u<=s-1; u++)
                    expr17[i][t][s] += BigM*R[u];
                expr17[i][t][s] -= BigM*R[s];
                constraints.add(expr17[i][t][s] >= -BigM);
            }*/

    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Zx[i][0].setBounds(0,0);

        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
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
        //cout << "RF: " << BestSolution << "\n";
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

            X_ini[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_ini%d_%d",i,t);
            X_ini[i][t].setName(name);

            X_co[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_co%d_%d",i,t);
            X_co[i][t].setName(name);

            for(int u=0; u<=T; u++)
            {
                Xx[i][t][u] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Xx%d_%d_%d",i,t,u);
                Xx[i][t][u].setName(name);
            }
        }

    }

    R = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        R[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"R%d",t);
        R[t].setName(name);
    }

    he = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        he[t] = IloNumVar(env,0,16,ILOINT);
        sprintf(name,"he%d",t);
        he[t].setName(name);
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
                    if(t<=l)
                    {
                        Yx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                        sprintf(name,"Yx%d_%d_%d",i,j,t);
                        Yx[i][j][t].setName(name);
                    }
                    else
                    {
                        Yx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                        sprintf(name,"Yx%d_%d_%d",i,j,t);
                        Yx[i][j][t].setName(name);
                    }

                    Fx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"Fx%d_%d_%d",i,j,t);
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
            sprintf(name,"Gx%d_%d",i,t);
            Gx[i][t].setName(name);
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int t=1; t<=T; t++)
        obj1 += ins->ExtraCost[t]*he[t];

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
            expr1[i][t] = IloExpr(env);
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
        for(int i=1; i<=N; i++)
        {
            for(int u=t; u<=T; u++)
                expr2[t] += ins->ProdTime[i]*Xx[i][t][u];
            for(int l=1; l<=N; l++)
                expr2[t] += ins->SetupTime[l][i]*Yx[l][i][t];
        }
        expr2[t] -= he[t];
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
                for(int l=t; l<=T; l++)
                {
                    for(int j=1; j<=N; j++)
                        expr3[i][t][u] -= ins->Demand[i][l]*Yx[j][i][l];
                    expr3[i][t][u] -= ins->Demand[i][l]*Zx[i][l];
                }
                constraints.add(expr3[i][t][u] <= 0);
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

    expr11 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr11[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
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
                expr16[i][t] -= 25*Yx[j][i][t];

            expr16[i][t] += 25*Zx[i][t+1];
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
                for(int u=t+1; u<=s; u++)
                    expr17[i][t][s] += X_co[i][u];
                for(int j=1; j<=N; j++)
                    expr17[i][t][s] -= 25*Yx[j][i][t];
                for(int u=t+1; u<=s-1; u++)
                    expr17[i][t][s] += BigM*R[u];
                expr17[i][t][s] -= BigM*R[s];
                constraints.add(expr17[i][t][s] >= -BigM);
            }*/

    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Zx[i][0].setBounds(0,0);

        for(int t=1; t<=T; t++)
            Yx[i][i][t].setBounds(0,0);
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF_%d.lp",l);
    //cplex.exportModel(name);
}

void RF::PathRelinking()
{
    bool re = false;
    IloModel model(env);
    PRModel(model);                                         //Constroi modelo classico

    SetParameters();

    re = cplex.solve();
    solutionRF = re?cplex.getObjValue():INF;                        //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //best_bound = re?cplex.getBestObjValue():INF;                   //If re = true then s->lp_value = cplex.getObjValue(), otherwise, s->lp_value = INF
    //cout << "RF solution period " << t << ": " << solutionRF << " -- Quantity of fea sol: " << cplex.getSolnPoolNsolns() << "\n";
    ///Get solution of the current problem
    if(re)
    {
        getFOptValue();
        CleanModel();
        model.end();
        cout << "PR1: " << BestSolution << "\n";
    }

}

void RF::PRModel(IloModel model)
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

            X_ini[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_ini%d_%d",i,t);
            X_ini[i][t].setName(name);

            X_co[i][t] = IloNumVar(env,0,IloInfinity,ILOINT);
            sprintf(name,"X_co%d_%d",i,t);
            X_co[i][t].setName(name);

            for(int u=0; u<=T; u++)
            {
                Xx[i][t][u] = IloNumVar(env,0,IloInfinity,ILOINT);
                sprintf(name,"Xx%d_%d_%d",i,t,u);
                Xx[i][t][u].setName(name);
            }
        }

    }

    R = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        R[t] = IloNumVar(env,0,1,ILOINT);
        sprintf(name,"R%d",t);
        R[t].setName(name);
    }

    he = IloNumVarArray(env,T+1);
    for(int t=0; t<=T; t++)
    {
        he[t] = IloNumVar(env,0,16,ILOINT);
        sprintf(name,"he%d",t);
        he[t].setName(name);
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
                    sprintf(name,"Yx%d_%d_%d",i,j,t);
                    Yx[i][j][t].setName(name);

                    Fx[i][j][t] = IloNumVar(env,0,IloInfinity,ILOINT);
                    sprintf(name,"Fx%d_%d_%d",i,j,t);
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
            sprintf(name,"Gx%d_%d",i,t);
            Gx[i][t].setName(name);
        }
    }

    IloExpr obj1(env);//initialize expression for objective function

    /*setup cost*/
    for(int t=1; t<=T; t++)
        obj1 += ins->ExtraCost[t]*he[t];

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
            expr1[i][t] = IloExpr(env);
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
        for(int i=1; i<=N; i++)
        {
            for(int u=t; u<=T; u++)
                expr2[t] += ins->ProdTime[i]*Xx[i][t][u];
            for(int l=1; l<=N; l++)
                expr2[t] += ins->SetupTime[l][i]*Yx[l][i][t];
        }
        expr2[t] -= he[t];
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
                for(int l=t; l<=T; l++)
                {
                    for(int j=1; j<=N; j++)
                        expr3[i][t][u] -= ins->Demand[i][l]*Yx[j][i][l];
                    expr3[i][t][u] -= ins->Demand[i][l]*Zx[i][l];
                }
                constraints.add(expr3[i][t][u] <= 0);
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

    expr11 = Exp2D(env,N+1);
    for(int i=0; i<=N; i++)
    {
        expr11[i] = IloExprArray(env,T+1);
        for(int t=0; t<=T; t++)
            expr11[i][t] = IloExpr(env);
    }

    for(int t=1; t<=T; t++)
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
                expr16[i][t] -= 25*Yx[j][i][t];

            expr16[i][t] += 25*Zx[i][t+1];
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
                for(int u=t+1; u<=s; u++)
                    expr17[i][t][s] += X_co[i][u];
                for(int j=1; j<=N; j++)
                    expr17[i][t][s] -= 25*Yx[j][i][t];
                for(int u=t+1; u<=s-1; u++)
                    expr17[i][t][s] += BigM*R[u];
                expr17[i][t][s] -= BigM*R[s];
                constraints.add(expr17[i][t][s] >= -BigM);
            }*/

    model.add(constraints);

    for(int i=1; i<=N; i++)
    {
        Zx[i][0].setBounds(0,0);

        for(int t=1; t<=T; t++)
        {
            Yx[i][i][t].setBounds(0,0);
            for(int j=1; j<=N; j++)
                if(vet[0].ystore[i][j][t]==0 && vet[1].ystore[i][j][t]==0 && vet[0].zstore[i][t]==0 && vet[1].zstore[i][t]==0)
                    Yx[i][j][t].setBounds(0,0);
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
        {
            if(reducedCostY[i][j][t] <= 0)
                Y[i][j][t] = -1;
        }
}

void RF::getFirstItems()
{
    double initItem;
    //printf("Qtd de sol encontradas %d!\n",pool);


        initItem = 0.0;
        for(int i=1; i<=N; i++)
        {
            Z[i][1] = cplex.getValue(Zx[i][1]);
            if(initItem <= Z[i][1])
            {
                initItem = Z[i][1];
                firstItem = i;
            }
        }
        printf("Primeiro item: %d\n",firstItem);

}

void RF::getRelaxedValue()
{
    bool re = true;
    double initItem;
    //printf("Qtd de sol encontradas %d!\n",pool);

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                X[i][t][u] = cplex.getValue(Xx[i][t][u]);
            Z[i][t] = cplex.getValue(Zx[i][t]);

            for(int j=1; j<=N; j++)
            {

                Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
                //printf("\nVariable %d_%d_%d_%d : %.4lf\n",i,j,k,t,Y[i][j][k][t]);
                //getchar();
                double aux;
                aux = cplex.getValue(Yx[i][j][t]);
                aux = aux - floor(aux+0.1);
                if(fabs(aux) > 0.01 && re==true)
                {
                    re = false;
                    printf("\nContinuous value %.4lf %.4lf %.4lf\n",Y[i][j][t],aux,fabs(aux));
                }
            }
        }
        Z[i][T+1] = cplex.getValue(Zx[i][T+1]);
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
        for(int j=1; j<=N; j++)
            for(int t=1; t<=T; t++)
                reducedCostY[i][j][t] = cplex.getReducedCost(Yx[i][j][t]);
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
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    X[i][t][u] = cplex.getValue(Xx[i][t][u],p);

                Z[i][t] = 1.0*IloRound(cplex.getValue(Zx[i][t],p));

                for(int j=1; j<=N; j++)
                {
                    if(id==1)
                    {
                        if(t > l)
                        {
                            Y[i][j][t] = cplex.getValue(Yx[i][j][t],p);
                            double aux;
                            aux = cplex.getValue(Yx[i][j][t],p);
                            aux = aux - floor(aux+0.1);
                            if(fabs(aux) > 0.01 && re==true)
                            {
                                re = false;
                                //printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][k][t],aux,fabs(aux));
                            }
                        }else
                            Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t],p));
                    }else
                    {
                        if(id==2)
                        {
                            if(t < l)
                            {
                                Y[i][j][t] = cplex.getValue(Yx[i][j][t],p);
                                double aux;
                                aux = cplex.getValue(Yx[i][j][t],p);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t],p));
                        }else{
                            if(t > l && t < (T-l+1))
                            {
                                Y[i][j][t] = cplex.getValue(Yx[i][j][t],p);
                                double aux;
                                aux = cplex.getValue(Yx[i][j][t],p);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t],p));
                        }
                    }
                }
            }
            Z[i][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][T+1],p));
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
    {
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                X[i][t][u] = cplex.getValue(Xx[i][t][u]);
            Z[i][t] = 1.0*IloRound(cplex.getValue(Zx[i][t]));
            G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
            F[0][i][t] = cplex.getValue(Fx[0][i][t]);
            for(int j=1; j<=N; j++)
            {
                F[i][j][t] = cplex.getValue(Fx[i][j][t]);
                if(id==1)
                    {
                        if(t > l)
                        {
                            Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
                            double aux;
                            aux = cplex.getValue(Yx[i][j][t]);
                            aux = aux - floor(aux+0.1);
                            if(fabs(aux) > 0.01 && re==true)
                            {
                                re = false;
                                //printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][k][t],aux,fabs(aux));
                            }
                        }else
                            Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t]));
                    }else
                    {
                        if(id==2)
                        {
                            if(t < l)
                            {
                                Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
                                double aux;
                                aux = cplex.getValue(Yx[i][j][t]);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t]));
                        }else{
                            if(t > l && t < (T-l+1))
                            {
                                Y[i][j][t] = cplex.getValue(Yx[i][j][t]);
                                double aux;
                                aux = cplex.getValue(Yx[i][j][t]);
                                aux = aux - floor(aux+0.1);
                                if(fabs(aux) > 0.01 && re==true)
                                {
                                    re = false;
                                    printf("Continuous value %.4lf %.4lf %.4lf\n",Y[i][j][t],aux,fabs(aux));
                                }
                            }else
                                Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t]));
                        }
                    }
            }
        }
        Z[i][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][T+1]));
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
        {
            for(int t=1; t<=T; t++)
            {
                for(int u=t; u<=T; u++)
                    X[i][t][u] = cplex.getValue(Xx[i][t][u],p);
                Z[i][t] = IloRound(cplex.getValue(Zx[i][t],p));

                for(int j=1; j<=N; j++)
                    Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t],p));
            }
            Z[i][T+1] = IloRound(cplex.getValue(Zx[i][T+1],p));
        }

        solutionRF = cplex.getObjValue(p);
        store_solution();
        //printf("Solution stored in FixAndOpt %.4lf!\n",solutionRF);
        if(BestSolution > solutionRF)
            BestSolution = solutionRF;

    }

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            for(int u=t; u<=T; u++)
                X[i][t][u] = cplex.getValue(Xx[i][t][u]);
            Z[i][t] = IloRound(cplex.getValue(Zx[i][t]));
            G[i][t] = IloRound(cplex.getValue(Gx[i][t]));
            F[0][i][t] = cplex.getValue(Fx[0][i][t]);
            for(int j=1; j<=N; j++)
            {
                F[i][j][t] = cplex.getValue(Fx[i][j][t]);
                Y[i][j][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][t]));
            }
        }
        Z[i][T+1] = IloRound(cplex.getValue(Zx[i][T+1]));
    }
}

bool RF::verify_duplication()
{
    for(int k = 0; k < nm; k++)
    {
        bool duplication = true;
        if(vet[k].sol != solutionRF) continue;
        for(int i = 1; i <= N; i++)
            for(int t = 1; t <= T; t++)
                for(int j = 1; j <= N; j++)
                {
                    if(vet[k].ystore[i][j][t] != Y[i][j][t])
                    {
                        duplication = false;
                        i = N+1;
                        j = N+1;
                        t = T+1;
                    }
                }

        for(int i = 1; i <= N; i++)
            for(int t = 1; t <= T; t++)
                for(int u = t; u <= T; u++)
                {
                    if(vet[k].xstore[i][t][u] != X[i][t][u])
                    {
                        duplication = false;
                        i = N+1;
                        u = T+1;
                        t = T+1;
                    }
                }

        for(int i = 1; i <= N; i++)
            for(int t = 1; t <= T; t++)
                if(vet[k].zstore[i][t] != Z[i][t])
                {
                    duplication = false;
                    i = N+1;
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

            {
                for(int t = 1; t <= T; t++)
                {
                    vet[0].zstore[i][t] = Z[i][t];
                    for(int u = t; u <= T; u++)
                        vet[0].xstore[i][t][u] = X[i][t][u];
                    for(int j = 1; j <= N; j++)
                        vet[0].ystore[i][j][t] = Y[i][j][t];
                }
                vet[0].zstore[i][T+1] = Z[i][T+1];
            }
        }

        if(vet[1].sol > solutionRF)
        {
            vet[1].sol = solutionRF;
            for(int i = 1; i <= N; i++)
            {
                for(int t = 1; t <= T; t++)
                {
                    vet[1].zstore[i][t] = Z[i][t];
                    for(int u = t; u <= T; u++)
                        vet[1].xstore[i][t][u] = X[i][t][u];
                    for(int j = 1; j <= N; j++)
                        vet[1].ystore[i][j][t] = Y[i][j][t];
                }
                vet[1].zstore[i][T+1] = Z[i][T+1];
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
    reducedCostY.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        F[i].resize(N+1);
        G[i].resize(T+1,0);
        X[i].resize(T+1);
        Z[i].resize(T+2,0.0);
        Y[i].resize(N+1);
        reducedCostY[i].resize(N+1);

        for(int t=0;t<=T;t++)
            X[i][t].resize(T+1,0.0);

        for(int j=0;j<=N;j++)
        {
            reducedCostY[i][j].resize(T+1,0.0);
            Y[i][j].resize(T+1,0.0);
            F[i][j].resize(T+1,0.0);
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
            vet[k].xstore[i].resize(T+1);
            vet[k].ystore[i].resize(N+1);
            vet[k].zstore[i].resize(T+2,0.0);

            for(int t = 0; t <= T; t++)
                vet[k].xstore[i][t].resize(T+1,0.0);

            for(int j = 0; j <= N; j++)
                vet[k].ystore[i][j].resize(T+1,0.0);

        }

    }
}

void RF::FreeVariables()
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            reducedCostY[i][j].clear();
            Y[i][j].clear();
            F[i][j].clear();
        }

        for(int t=0;t<=T;t++)
            X[i][t].clear();

        X[i].clear();
        Z[i].clear();
        F[i].clear();
        G[i].clear();
        Y[i].clear();
        reducedCostY[i].clear();
    }
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
                vet[k].ystore[i][j].clear();

            for(int t = 0; t <= T; t++)
                vet[k].xstore[i][t].clear();

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
