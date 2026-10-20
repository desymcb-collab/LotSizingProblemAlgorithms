#include "raiz.h"

/************Inicializa variaveis do problema************/
void Instance::inicializa_variaveis(void)
{
    Capacity.resize(M+1);
    for(int j=0;j<=M;j++)
        Capacity[j].resize(T+1,0.0);

    HoldingCost.resize(N+1,0);
    SetupTime.resize(N+1);
    SetupCost.resize(N+1);
    ProdCost.resize(N+1);
    ProdTime.resize(N+1);
    Demand.resize(N+1);
    MaxProd.resize(N+1);
    MachineProd.resize(N+1);
    minl.resize(N+1);
    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0.0);
        MaxProd[i].resize(T+2,0.0);
        MachineProd[i].resize(M+1,0);
        ProdTime[i].resize(M+1,0.0);
        ProdCost[i].resize(M+1,0.0);
        minl[i].resize(M+1,0);
        SetupTime[i].resize(N+1);
        SetupCost[i].resize(N+1);
        for(int j=0;j<=N;j++)
        {
            SetupTime[i][j].resize(M+1,0.0);
            SetupCost[i][j].resize(M+1,0.0);
        }
    }

}

/************Inicializa variaveis do problema************/
void Instance::desaloca_variaveis(void)
{
    for(int i=0;i<=N;i++)
    {
        for(int j=0;j<=N;j++)
        {
            SetupTime[i][j].clear();
            SetupCost[i][j].clear();
        }
        Demand[i].clear();
        MaxProd[i].clear();
        MachineProd[i].clear();
        SetupTime[i].clear();
        SetupCost[i].clear();
        ProdCost[i].clear();
        minl[i].clear();
    }

    for(int j=0;j<=M;j++)
        Capacity[j].clear();

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
}

/********************Read instances**********************/
void Instance::le_dados_SM(char * input_data)
{
    FILE   *arq_ent;

    double  faux1, faux2, faux3, faux4;
    double  aux1;
    int Iaux;

    arq_ent = fopen(input_data,"r");

    int r3;
    char aux[20];
    /* read general data of the problem */
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&N);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&T);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&M);

    inicializa_variaveis();
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        HoldingCost[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            ProdTime[i][k] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%lf",&aux1);
                SetupTime[i][j][k] = aux1;
            }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            Capacity[k][t] = aux1;
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
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%lf",&aux1);
                SetupCost[i][j][k] = aux1;
            }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%d",&Iaux);
            MachineProd[i][k] = Iaux;
            //printf("%d ",MachineProd[i][k]);
        }
        //printf("\n");
    }

    /*r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%d",&Iaux);
            minl[i][k] = Iaux;
        }*/

    fclose(arq_ent);
}

/**************Mostra dados das instancias***************/
void Instance::show_instance(void)
{
    printf("%d\n",N);
    printf("%d\n",T);
    printf("%d\n",M);

    for(int i=1; i<=N; i++)
        printf("%.1lf ",HoldingCost[i]);

    printf("\n");
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
            printf("%.1lf ",ProdTime[i][k]);
        printf("\n");
    }

    printf("\n");

    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
                printf("%.1lf ",SetupTime[i][j][k]);
            printf("\n");
        }

    for(int k=1; k<=M; k++)
    {
        for(int t=1; t<=T; t++)
            printf("%.1lf ",Capacity[k][t]);
        printf("\n");
    }

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

    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
                printf("%.1lf ",SetupCost[i][j][k]);
            printf("\n");
        }
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            printf("%d ",MachineProd[i][k]);
        }
        printf("\n");
    }
}
