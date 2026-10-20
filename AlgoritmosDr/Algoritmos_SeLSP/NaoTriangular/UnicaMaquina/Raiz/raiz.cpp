#include "raiz.h"

/************Inicializa variaveis do problema************/
void Instance::inicializa_variaveis(void)
{
    Capacity.resize(T+1,0);
    //for(int j=0;j<=M;j++)
        //Capacity[j].resize(T+1,0.0);

    HoldingCost.resize(N+1,0);
    SetupTime.resize(N+1);
    SetupCost.resize(N+1);
    ProdCost.resize(N+1);
    ProdTime.resize(N+1,0);
    Demand.resize(N+1);
    MaxProd.resize(N+1);
    isetup.resize(N+1,0);
    istock.resize(N+1,0);
    maxl.resize(N+1,0);
    minl.resize(N+1,0);
    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0);
        MaxProd[i].resize(T+1,0);
        SetupTime[i].resize(N+1,0);
        SetupCost[i].resize(N+1,0);
        ProdCost[i].resize(N+1,0);
    }

}

/************Inicializa variaveis do problema************/
void Instance::desaloca_variaveis(void)
{
    for(int i=0;i<=N;i++)
    {
        Demand[i].clear();
        MaxProd[i].clear();
        SetupTime[i].clear();
        SetupCost[i].clear();
        ProdCost[i].clear();
    }

    Capacity.clear();
    HoldingCost.clear();
    SetupTime.clear();
    SetupCost.clear();
    ProdCost.clear();
    ProdTime.clear();
    Demand.clear();
    MaxProd.clear();
    isetup.clear();
    istock.clear();
    maxl.clear();
    minl.clear();
}

/********************Read instances**********************/
void Instance::le_dados_SM(char * input_data)
{
    FILE   *arq_ent;
    char aux[10];

    double  faux1, faux2, faux3, faux4;
    double  aux1;

    arq_ent = fopen(input_data,"r");

    int r3;
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
        for(int j=1; j<=T; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            MaxProd[i][j] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        isetup[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        istock[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        minl[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        maxl[i] = aux1;
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
            printf("%.1lf ",Demand[i][t] );
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
}
