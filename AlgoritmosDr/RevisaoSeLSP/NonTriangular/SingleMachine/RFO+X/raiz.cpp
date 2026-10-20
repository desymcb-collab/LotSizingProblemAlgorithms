#include "raiz.h"

/************Inicializa variaveis do problema************/
void Instance::inicializa_variaveis(void)
{
    //Parameter
    Capacity.resize(T+1,0.0);
    HoldingCost.resize(N+1,0);
    SetupTime.resize(N+1);
    SetupCost.resize(N+1);
    ProdCost.resize(N+1,0.0);
    ProdTime.resize(N+1,0.0);
    Demand.resize(N+1);
    MaxProd.resize(N+1);
    MachineProd.resize(N+1,0.0);
    minl.resize(N+1,0.0);
    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0.0);
        MaxProd[i].resize(T+2,0.0);        
        SetupTime[i].resize(N+1,0.0);
        SetupCost[i].resize(N+1,0.0);
    }

    //Variables
    X.resize(N+1);
    Z.resize(N+1);
    F.resize(N+1);
    G.resize(N+1);
    Y.resize(N+1);
    reducedCostY.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        F[i].resize(N+1);
        Y[i].resize(N+1);
        reducedCostY[i].resize(N+1);
        
        X[i].resize(T+1);
        Z[i].resize(T+2,0.0);
        G[i].resize(T+1,0.0);
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
            {
                vet[k].ystore[i][j].resize(T+1,0.0);
            }
        }

    }

}

/************Inicializa variaveis do problema************/
void Instance::desaloca_variaveis(void)
{
    //Parameter
    for(int i=0;i<=N;i++)
    {
        Demand[i].clear();
        MaxProd[i].clear();
        SetupTime[i].clear();
        SetupCost[i].clear();
    }

    Capacity.clear();
    HoldingCost.clear();
    SetupTime.clear();
    SetupCost.clear();
    ProdCost.clear();
    ProdTime.clear();
    Demand.clear();
    MaxProd.clear();
    MachineProd.clear();
    minl.clear();

    //Variables
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
            {
                vet[k].ystore[i][j].clear();
            }
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

/********************Read instances**********************/
/********************Read instances**********************/
void Instance::le_dados_SM(const char * input_data)
{
    FILE   *arq_ent;

    double  faux1, faux2, faux3, faux4;
    double  aux1;
    int Iaux,IIaux;

    arq_ent = fopen(input_data,"r");

    int r3;
    char aux[20];
    /* read general data of the problem */
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&N);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&T);

    inicializa_variaveis();
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        HoldingCost[i] = aux1;
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        ProdTime[i] = aux1;
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupTime[i][j] = aux1;
        }

    nontri[0] = 0;
    nontri[1] = 0;
    nontri[2] = 0;
    for(int i=1; i<=N; i++)
    {
        Iaux = 0;
        for(int j=1; j<=N; j++)
        {
            if(SetupTime[i][j] > 4)
            {
                Iaux += 1;
                if(Iaux == 1)
                {
                    nontri[0] = i;
                    nontri[1] = j;
                }else
                {
                    if(Iaux == 2 && N > 20)
                    {
                        nontri[2] = j;
                    }else
                    {
                        j = N+1;
                        nontri[0] = 0;
                        nontri[1] = 0;
                        nontri[2] = 0;
                    }
                }
            }
        }
        if(nontri[0] > 0)
            i = N+1;
    }
    printf("Custos n-triangular para %d, %d e %d \n",nontri[0],nontri[1],nontri[2]);
    
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupCost[i][j] = aux1;
        }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int t=1; t<=T; t++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        Capacity[t] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            Demand[i][t] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=T; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            MaxProd[i][j] = aux1;
        }
        MaxProd[i][T+1] = 0;
    }
    
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        MachineProd[i] = Iaux;
        //printf("%d ",MachineProd[i]);
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        MachineProd[i] = Iaux;
        //printf("%d ",MachineProd[i]);
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        minl[i] = Iaux;
    }

    fclose(arq_ent);
}


/**************Mostra dados das instancias***************/
void Instance::show_instance(void)
{
    printf("%d\n",N);
    printf("%d\n",T);

    for(int i=1; i<=N; i++)
        printf("%.1lf ",HoldingCost[i]);

    printf("\n");
    for(int i=1; i<=N; i++)
        printf("%.1lf ",ProdTime[i]);

    printf("\n");

    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=N; j++)
            printf("%.1lf ",SetupTime[i][j]);
        printf("\n");
    }

    for(int t=1; t<=T; t++)
        printf("%.1lf ",Capacity[t]);
    

    printf("\n");
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
            printf("%.1lf ",Demand[i][t]);
        printf("\n");
    }

    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=T; j++)
            printf("%.1lf ",MaxProd[i][j]);
        printf("\n");
    }

    
    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=N; j++)
            printf("%.1lf ",SetupCost[i][j]);
        printf("\n");
    }

    for(int i=1; i<=N; i++)
    {
        printf("%d ",MachineProd[i]);
    }
}

void Instance::show_solution()
{
    int id;

    if(vet[1].sol >= vet[2].sol)
    {
        id = 2;
    }else id = 1;

    printf("\nsol: %.2lf\n",vet[id].sol);
    for(int i = 1; i <= N; i++)
        for(int t = 1; t <= T; t++)
        {
            if(vet[id].zstore[i][t] > 0)
                printf("\nz_%d_%d: %.2lf\n",i,t,vet[id].zstore[i][t]);
            for(int u = t; u <= T; u++)
            {
                if(vet[id].xstore[i][t][u] > 0)
                    printf("x_%d_%d_%d: %.2lf\n",i,t,u,vet[id].xstore[i][t][u]);
            }
            
            for(int j = 1; j <= N; j++)
            {
                if(vet[id].ystore[i][j][t] > 0)
                printf("y_%d_%d_%d:%.2lf\n",i,j,t,vet[id].ystore[i][j][t]);
            }                
        }
        
}