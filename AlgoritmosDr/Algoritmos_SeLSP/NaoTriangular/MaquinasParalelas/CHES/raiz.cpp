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
    MP.resize(N+1);
    MaxProd.resize(N+1);
    MachineProd.resize(N+1);
    minl.resize(N+1);
    maxl.resize(N+1);
    iscost.resize(N+1);
    lscost.resize(N+1);
    istock.resize(N+1);
    istate.resize(M+1,0);
    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0.0);
        MP[i].resize(T+1,0.0);
        MaxProd[i].resize(T+2,0.0);
        MachineProd[i].resize(M+1,0);
        ProdTime[i].resize(M+1,0.0);
        ProdCost[i].resize(M+1,0.0);
        minl[i].resize(M+1,0.0);
        maxl[i].resize(M+1,0.0);
        SetupTime[i].resize(N+1);
        SetupCost[i].resize(N+1);
        iscost[i].resize(M+1,0.0);
        lscost[i].resize(M+1,0.0);
        istock[i].resize(T+1,0.0);
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
        MP[i].clear();
        MaxProd[i].clear();
        MachineProd[i].clear();
        SetupTime[i].clear();
        SetupCost[i].clear();
        ProdCost[i].clear();
        minl[i].clear();
        maxl[i].clear();
        iscost[i].clear();
        lscost[i].clear();
        istock[i].clear();
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
    MP.clear();
    MaxProd.clear();
    MachineProd.clear();
    minl.clear();
    maxl.clear();
    iscost.clear();
    lscost.clear();
    istock.clear();
    istate.clear();
}

/********************Read instances**********************/
void Instance::le_dados_SM(char * input_data)
{
    FILE   *arq_ent;

    double  faux1, faux2, faux3, faux4;
    double  aux1;
    int Iaux;

    arq_ent = fopen(input_data,"r");
    if(arq_ent == NULL)
        printf("Erro ao abrir documento!!\n");

    int r3;
    char aux[20];
    /* read general data of the problem */
    //printf("fopen\n");
    r3 = fscanf(arq_ent,"%s",aux);
    //printf("fopen %s\n",aux);
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
//    r3 = fscanf(arq_ent,"%s",aux);
//    for(int k=1; k<=M; k++)
//        for(int i=1; i<=N; i++)
//            for(int j=1; j<=N; j++)
//            {
//                r3 = fscanf(arq_ent,"%lf",&aux1);
//                SetupTime[i][j][k] = aux1;
//            }
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
    //r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            MaxProd[i][t] = 0;
            for(int u=t; u<=T; u++)
                MaxProd[i][t] += Demand[i][u];
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
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            iscost[i][k] = aux1;
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

    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            minl[i][k] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            maxl[i][k] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        istate[k] = Iaux;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        istock[i][1] = aux1;
//        for(int t=2; t<=T; t++)
//        {
//            double dem;
//            dem = istock[i][t-1] - Demand[i][t-1];
//            if(dem > 0.0)
//                istock[i][t] = dem;
//            else
//                istock[i][t] = 0.0;
//        }
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            MP[i][t] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            lscost[i][k] = aux1;
        }


    fclose(arq_ent);
}

/**************Mostra dados das instancias***************/
void Instance::show_instance(void)
{
    printf("%d\n",N);
    printf("%d\n",T);
    printf("%d\n",M);

    for(int i=1; i<=N; i++)
        printf("%.2lf ",HoldingCost[i]);

    printf("\n");
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
            printf("%.5lf ",ProdTime[i][k]);
        printf("\n");
    }

    printf("\n");

//    for(int k=1; k<=M; k++)
//        for(int i=1; i<=N; i++)
//        {
//            for(int j=1; j<=N; j++)
//                printf("%.1lf ",SetupTime[i][j][k]);
//            printf("\n");
//        }

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
                printf("%.2lf ",SetupCost[i][j][k]);
            printf("\n");
        }
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
                printf("%.1lf ",iscost[i][k]);

        }
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
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
                printf("%.1lf ",minl[i][k]);

        }
         printf("\n");
    }

    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
                printf("%.1lf ",maxl[i][k]);

        }
         printf("\n");
    }
    for(int k=1; k<=M; k++)
        printf("%d ",istate[k]);
    printf("\n");

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
            printf("%.1lf ",istock[i][t] );
        printf("\n");
    }

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
            printf("%.1lf ",MP[i][t] );
        printf("\n");
    }

    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
                printf("%.1lf ",lscost[i][k]);

        }
         printf("\n");
    }
}
