#include "rfks.h"

void RFKS::KernelSearch()
{
    env = IloEnv();

    BestSolution = ins->vet[1].sol;
    solutionKS = INF;
    
    LRC_ks kernel,bucket;

    Phase1(&kernel,&bucket);

    Phase2(&kernel,&bucket);
    
    env.end();
    //cout << "Number of feasible solutions: " << nm << "\n";
    cout << "RFKS: " << BestSolution << "\n";
}

/******************************************************************************
                        Kernel Search - Auxiliary functions
******************************************************************************/

int RFKS::busca_binaria_ks(double LS, LRC_ks vetor)
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

void RFKS::Insere_ks(LRC_ks *lrc, double campo1, int campo2, int campo3, int campo4, int campo5, int kernel)
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

int RFKS::busca_elemento(LRC_ks vetor,int t,int j,int i,int k)
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

void RFKS::retira_elemento(LRC_ks *lrc, int campo2, int campo3,int campo4, int campo5)
{ //campos da estrutura
    int ii, aux;

    aux = busca_elemento(*lrc,campo2,campo3,campo4,campo5);
    lrc->campo[aux].kernel = 1;
//    for(ii=aux; ii<=lrc->nelementos-1; ii++)
//        lrc->campo[ii] = lrc->campo[ii+1];
//    lrc->nelementos--;
}

/******************************************************************************
                        Kernel Search - Inicialization Phase
******************************************************************************/

void RFKS::Define_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket)
{
    bool rc=false;

    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    if(ins->reducedCostY[i][j][k][t] > 0.001)
                    {
                        //printf("Custo reduzido: %.2lf!!!\n",ins->reducedCostY[i][j][k][t]);
                        rc = true;
                        i = N+1;
                        j = N+1;
                        k = M+1;
                        t = T+1;                        
                    }//else printf("Custos reduzidos iguais a 0!!!\n");

    for(int i = 1; i<=N; i++)
        for(int j = 1; j<=N; j++)
            for(int k = 1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    if(i!=j)
                        if(rc)
                        {
                            if((ins->reducedCostY[i][j][k][t] <= 0) || (ins->vet[0].ystore[i][j][k][t] > 0.1) || (ins->vet[1].ystore[i][j][k][t] >= 0.5) || (ins->vet[2].ystore[i][j][k][t] == 1)  || (i==ins->nontri[k-1][0]) || (i== ins->nontri[k-1][1]) || (j==ins->nontri[k-1][0]) || (j== ins->nontri[k-1][1]))    //criterio para criacao do kernel (a variavel y tem valor zero em mais de taxa% das solucoes)
                            {
                                Insere_ks(kernel,ins->vet[0].ystore[i][j][k][t],t,j,i,k,1);  //as variaveis que obtiverem solucao menor que o criterio durante mais de 5 solucao serao excluidas do kernel
                            }else{
                                Insere_ks(bucket,ins->reducedCostY[i][j][k][t],t,j,i,k,0);   // 0 indica se a variavel do bucket foi inserida no kernel durante alguma iteracao
                            }
                        }else{
                            if((ins->vet[0].ystore[i][j][k][t] > 0.1) || (ins->vet[1].ystore[i][j][k][t] >= 0.5) || (ins->vet[2].ystore[i][j][k][t] == 1)  || (i==ins->nontri[k-1][0]) || (i== ins->nontri[k-1][1]) || (j==ins->nontri[k-1][0]) || (j== ins->nontri[k-1][1]))    //criterio para criacao do kernel (a variavel y tem valor zero em mais de taxa% das solucoes)
                            {
                                Insere_ks(kernel,ins->vet[0].ystore[i][j][k][t],t,j,i,k,1);  //as variaveis que obtiverem solucao menor que o criterio durante mais de 5 solucao serao excluidas do kernel
                            }else{
                                Insere_ks(bucket,ins->vet[0].ystore[i][j][k][t],t,j,i,k,0);   // 0 indica se a variavel do bucket foi inserida no kernel durante alguma iteracao
                            }
                        }
    
}

void RFKS::Update_Kernel_Bucket(LRC_ks *kernel, LRC_ks *bucket)
{

    int l;
    for(int i = 1; i<=N; i++) {
        for(int j = 1; j<=N; j++) {
            for(int k = 1; k<=M; k++) {
                for(int t=1; t<=T; t++) {
                    if(ins->Y[i][j][k][t] == 1.0){
                        l=busca_elemento(*kernel,t,j,i,k);
                        if(l<0)
                        {
                            Insere_ks(kernel,ins->Y[i][j][k][t],t,j,i,k,1);
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

void RFKS::Phase1(LRC_ks *kernel, LRC_ks *bucket)
{
    struct timespec fim;
    
    total_time = 1800 + additionalTimeKS;
        
    kernel->nelementos = 0;
    bucket->nelementos = 0;

    Define_Kernel_Bucket(kernel,bucket);

    qtd_var_bucket = kernel->nelementos;
    qtd_buckets = ceil(double(bucket->nelementos)/double(qtd_var_bucket));            //number of bucket
    printf("Quantidade de elementos no kernel: %d\n",kernel->nelementos);
    printf("Quantidade de elementos no bucket: %d\n",bucket->nelementos);
    printf("Quantidade de bucket: %d\n",qtd_buckets);

    clock_gettime(CLOCK_MONOTONIC, &fim);
    elapsedTime = (fim.tv_sec - start);
    subKS_time = floor((total_time - elapsedTime)/(qtd_buckets+1));
    //printf("Tempo subP: %.2lf\n",subKS_time);

    if(kernel->nelementos == 0)
        printf("Deu ruim: kernel vazio!!!\n");

    Prob_Reduzido_kernel(*bucket);
    printf("Kernel --- Sol: %.2lf\n",solutionKS);
}

void RFKS::Prob_Reduzido_kernel(LRC_ks bucket)
{
    bool re = false;

    env = IloEnv();
    re = Prob_Reduzido(bucket);                       //Variante do modelo de Sambasivan para a utilizacao de carry-over
    env.end();

    //cout << "Solucao Prob_Reduzido: " << LS->BestSolution << "\tLimitante Inferior: " << LS->LowerBound << "\n";
}

bool RFKS::Prob_Reduzido(LRC_ks bucket)
{
    bool re = false;

    IloModel model(env);

    FLP_Reduzido(model,bucket);   //Constroi modelo localizacao de facilidades

    SetParameters();

    IloNumArray val(env);
    IloNumVarArray start(env);
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
            for(int k=1; k<=M; k++)
            {
                for(int j=1; j<=N; j++)
                {
                    start.add(Yx[i][j][k][t]);
                    val.add(ins->Y[i][j][k][t]);
                }
                for(int u=t; u<=T; u++)
                {
                    start.add(Xx[i][k][t][u]);
                    val.add(ins->X[i][k][t][u]);
                }
                start.add(Zx[i][k][t]);
                    val.add(ins->Z[i][k][t]);
            }

    cplex.addMIPStart(start,val);
    start.end();
    val.end();
    re = cplex.solve();
    ///Get solution of the current problem
    solutionKS = re?cplex.getObjValue():INF;

    if(re)
    {   
        getFOptValue();
        CleanModel();
        model.end();
    }

    return re;
}

void RFKS::FLP_Reduzido(IloModel model, LRC_ks bucket)
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
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
                Zx[i][k][T+1].setBounds(0,0);
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);

        }
    }

    cplex = IloCplex(model);
    //sprintf(name,"lp/RF_%hu.lp",l);
    //cplex.exportModel(name);
}

/******************************************************************************
                        Kernel Search - Extension Phase
******************************************************************************/
void RFKS::Phase2(LRC_ks *kernel,LRC_ks *bucket)
{
    struct timespec fim;

    for(int ii=1;ii<=qtd_buckets;ii++)
    {
        clock_gettime(CLOCK_MONOTONIC, &fim);
        elapsedTime = (fim.tv_sec - start);
        //tempo = difftime(tempo2,start);
        if(elapsedTime < total_time)
        {
            subKS_time = floor((total_time - elapsedTime)/(qtd_buckets-ii+1));
            Prob_Reduzido_bucket(*kernel,*bucket,ii);
            printf("Bucket: %d --- Sol: %.2lf\n",ii,solutionKS);
            if(solutionKS < INF)
            {
                Update_Kernel_Bucket(kernel,bucket);
            }
        }else ii = qtd_buckets+1;
    }
}

void RFKS::Prob_Reduzido_bucket(LRC_ks kernel, LRC_ks bucket, int id)
{
    bool re = false;

    env = IloEnv();
    re = Prob_bucket(kernel,bucket,id);                                    //Variante do modelo de Sambasivan para a utilizacao de carry-over
    env.end();
    //cout << "SolucaoFLP: " << LS->BestSolution << "\tLimitante Inferior: " << LS->LowerBound << "\n";
}

bool RFKS::Prob_bucket(LRC_ks kernel, LRC_ks bucket, int id)
{
    bool re = false;

    IloModel model(env);

    FLP_Reduzido_bucket(model,kernel,bucket,id);   //Constroi modelo localizacao de facilidades

    SetParameters();

    if(ins->vet[6].sol < INF)
    {
        IloNumArray val(env);
        IloNumVarArray start(env);
        for(int i=1; i<=N; i++)
            for(int t=1; t<=T; t++)
                for(int k=1; k<=M; k++)
                {
                    for(int j=1; j<=N; j++)
                    {
                        start.add(Yx[i][j][k][t]);
                        val.add(ins->vet[6].ystore[i][j][k][t]);
                    }
                    for(int u=t; u<=T; u++)
                    {
                        start.add(Xx[i][k][t][u]);
                        val.add(ins->vet[6].xstore[i][k][t][u]);
                    }
                    start.add(Zx[i][k][t]);
                        val.add(ins->vet[6].zstore[i][k][t]);
                }

        cplex.addMIPStart(start,val);
        start.end();
        val.end();
    }

    re = cplex.solve();
    ///Get solution of the current problem
    solutionKS = re?cplex.getObjValue():INF;

    if(re)
    {
        getFOptValue();
        CleanModel();
        model.end();
    }

    return re;
}

void RFKS::FLP_Reduzido_bucket(IloModel model, LRC_ks kernel, LRC_ks bucket, int id)
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
            R[i][t] = IloNumVar(env,0,1,ILOINT);
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
                Zx[i][k][t] = IloNumVar(env,0,1,ILOINT);
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
                    Yx[i][j][k][t] = IloNumVar(env,0,IloInfinity,ILOINT);
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
                Gx[i][k][t] = IloNumVar(env,0,1,ILOINT);
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
        for(int j=1; j<=N; j++)
            for(int k=1; k<=M; k++)
                for(int t=1; t<=T; t++)
                    expr += ins->SetupCost[i][j][k]*Yx[i][j][k][t];

    /*holding cost*/
    for(int i=1; i<=N; i++)
        for(int k=1; k<=M; k++)
            for(int t=1; t<=T; t++)
                for(int u=t; u<=T; u++)
                    expr += (u-t)*ins->HoldingCost[i]*Xx[i][k][t][u];

    constraints.add(expr <= (BestSolution -0.0001));
    expr.clear();

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

    int aux;
    aux = id*qtd_var_bucket;
    if(aux > bucket.nelementos)
        aux = bucket.nelementos;
    for(int i=(id-1)*qtd_var_bucket; i < aux; i++)
            expr += Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1];

    //for(int i=0; i < kernel.nelementos; i++)
      //  expr_ks_8 += z[kernel.campo[i].item][kernel.campo[i].j1][kernel.campo[i].t1];

    constraints.add(expr >= 1);
    expr.clear();

    model.add(constraints);

    for(int i=0; i < (id-1)*qtd_var_bucket; i++)
        if(bucket.campo[i].kernel == 0)
            Yx[bucket.campo[i].i1][bucket.campo[i].j1][bucket.campo[i].k1][bucket.campo[i].t1].setBounds(0,0);
    for(int i=id*qtd_var_bucket; i < bucket.nelementos; i++)
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
                Zx[i][k][T+1].setBounds(0,0);
            }
            for(int t=1; t<=T; t++)
                Yx[i][i][k][t].setBounds(0,0);
        }
    }

    cplex = IloCplex(model);
    
    //cplex.exportModel("resultados/KS.lp");
}

/*****************************************************************************/
/*****************************************************************************/

void RFKS::CleanModel()
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

void RFKS::getFOptValue()
{
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
                    ins->Y[i][j][k][t] = 1.0*IloRound(cplex.getValue(Yx[i][j][k][t]));
                }
            }
            ins->Z[i][k][T+1] = 1.0*IloRound(cplex.getValue(Zx[i][k][T+1]));
        }

    solutionKS = cplex.getObjValue();

    if(BestSolution > solutionKS)
    {
        BestSolution = solutionKS;
        store_solution();
    }
}

void RFKS::store_solution()
{
    if(ins->vet[6].sol > BestSolution)
    {
        ins->vet[6].sol = BestSolution;
        for(int i = 1; i <= N; i++)
            for(int k = 1; k <= M; k++)
            {
                for(int t = 1; t <= T; t++)
                {
                    ins->vet[6].zstore[i][k][t] = ins->Z[i][k][t];
                    for(int u = t; u <= T; u++)
                        ins->vet[6].xstore[i][k][t][u] = ins->X[i][k][t][u];
                    for(int j = 1; j <= N; j++)
                        ins->vet[6].ystore[i][j][k][t] = ins->Y[i][j][k][t];
                }
                ins->vet[6].zstore[i][k][T+1] = ins->Z[i][k][T+1];
            }
    }
}

void RFKS::SetParameters()
{
    cplex.setParam(IloCplex::Param::Threads,1);
    cplex.setParam(IloCplex::Param::TimeLimit,subKS_time);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.000001);
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
