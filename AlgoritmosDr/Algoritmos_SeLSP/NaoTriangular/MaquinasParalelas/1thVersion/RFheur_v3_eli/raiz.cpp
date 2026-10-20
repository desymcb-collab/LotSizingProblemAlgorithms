#include "raiz.h"

/************Inicializa variaveis do problema************/
void Instance::inicializa_variaveis(void)
{
    Capacity.resize(T+1,0.0);
    HoldingCost.resize(N+1,0.0);
    ProdTime.resize(N+1,0.0);
    HExtra.resize(N+1,0.0);
    ExtraCost.resize(N+1,0.0);
    maxl.resize(N+1,0.0);

    SetupTime.resize(N+1);
    Demand.resize(N+1);

    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0.0);
        SetupTime[i].resize(N+1,0.0);
    }

}

/************Inicializa variaveis do problema************/
void Instance::desaloca_variaveis(void)
{
    for(int i=0;i<=N;i++)
    {
        Demand[i].clear();
        SetupTime[i].clear();
    }

    Capacity.clear();
    HoldingCost.clear();
    SetupTime.clear();
    ProdTime.clear();
    Demand.clear();
    HExtra.clear();
    ExtraCost.clear();
    maxl.clear();
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
    for(int i=1; i<=T; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        HExtra[i] = aux1;
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int t=1; t<=T; t++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        ExtraCost[t] = aux1;
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
