#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

char input_arq[200],output_arq[200];

void convert_arq(int factor);

int main()
{
    int qt_period,index_items,qt,id,id2,id3,id4;
    int qt_items[4];
    int p = 3;
    int cd[2],prob1[2],prob2[2];

    double cl1[2];
    double cl2 = 0.5;

    qt_items[0] = 15; qt_items[1] = 20; qt_items[2] = 15; qt_items[3] = 400;
    cl1[0] = 0.6; cl1[1] = 0.8;
    cd[0] = 50; cd[1] = 100;
    prob1[0] = 60; prob1[1] = 80;
    prob2[0] = 10; prob2[1] = 20;

    for(id = 1; id <= 1; id++) // 0 a 1 cl1
        for(id2 = 1; id2 <= 1; id2++) // 0 a 1 cd
            for(id3 = 0; id3 <= 0; id3++) // 0 a 1 prob1
                for(id4 = 1; id4 <= 1; id4++) // 0 a 1 prob2
                    for(qt_period = 10; qt_period <= 10; qt_period = qt_period +5) // 5 a 15
                        for(index_items = 0; index_items <= 0; index_items ++)  // 0 a 1
                            for(qt = 1; qt <= 10; qt++)  // 0 a 9
                            {
                                sprintf(input_arq,"../Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d.txt",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);
                                sprintf(output_arq,"../Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);
                                printf("Data%d-%d-%d-%.1lf-%.1lf-%d-%d-%d-%d-mod.dat\n",p,qt_items[index_items],qt_period,cl1[id],cl2,cd[id2],prob1[id3],prob2[id4],qt);
                                convert_arq(cd[id2]);
                            }
    return 0;
}

/********************Convert arq instances**********************/
void convert_arq(int factor)
{
    FILE   *arq_ent,*new_arq;
    srand(time(NULL));
    int  faux1, faux2, faux3, faux4;
    int  aux1,N,M,T;

    arq_ent = fopen(input_arq,"r");
    new_arq = fopen(output_arq,"w+");

    int r3;
    char aux[20];
    /* read general data of the problem */

    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"N\n");
    printf("Entrou\n");
    r3 = fscanf(arq_ent,"%d",&N);
    fprintf(new_arq,"%d\n",N);
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"T\n");
    r3 = fscanf(arq_ent,"%d",&T);
    fprintf(new_arq,"%d\n",T);
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"M\n");
    r3 = fscanf(arq_ent,"%d",&M);
    fprintf(new_arq,"%d\n",M);
    printf("Saiua\n");

    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&aux1);
        fprintf(new_arq,"%d ",aux1);
    }
    fprintf(new_arq,"\n");
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%d",&aux1);
            fprintf(new_arq,"%d ",aux1);
        }
        fprintf(new_arq,"\n");
    }

    ///Setup Time
    int SetupTime[N+1][N+1][M+1],SetupCost[N+1][N+1][M+1];
    int rand1[M+1][2];
    //printf("Saiu1c\n");

    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int k=1; k<=M; k++)
    {
        rand1[k][0] = 1 + rand()%N;
        rand1[k][1] = 1 + rand()%N;
        while(rand1[k][0] == rand1[k][1])
            rand1[k][1] = 1 + rand()%N;
        printf("%d %d\n",rand1[k][0],rand1[k][1]);
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%d",&aux1);

                SetupTime[i][j][k] = aux1;

                if(i== rand1[k][0] && i!=j && j!=rand1[k][1])
                    SetupTime[i][j][k] = 2 + rand()%3;
                if(j== rand1[k][0] && i!=j && i!=rand1[k][1])
                    SetupTime[i][j][k] = 2 + rand()%3;
                if(i== rand1[k][1] && i!=j && j!=rand1[k][0])
                    SetupTime[i][j][k] = 2 + rand()%3;
                if(j== rand1[k][1] && i!=j && i!=rand1[k][0])
                    SetupTime[i][j][k] = 2 + rand()%3;
                //printf("Leu %d %d: %d\n",i,j,SetupTime[i][j][k]);
                fprintf(new_arq,"%d",SetupTime[i][j][k]);
                //printf("Leu %d %d: %d\n",i,j,SetupTime[i][j][k]);
                if(j!=N)
                    fprintf(new_arq," ");

            }
            //printf("Pulou linha\n");
            fprintf(new_arq,"\n");
        }
        //fprintf(new_arq,"\n");
    }
    //printf("Saiu1\n");
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int k=1; k<=M; k++)
    {
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%d",&aux1);
            fprintf(new_arq,"%d ",aux1);
        }
        fprintf(new_arq,"\n");
    }

    //printf("Saiu2\n");
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%d",&aux1);
            fprintf(new_arq,"%d ",aux1);
        }
        fprintf(new_arq,"\n");
    }
    //printf("Saiu3\n");

    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=T; j++)
        {
            r3 = fscanf(arq_ent,"%d",&aux1);
            fprintf(new_arq,"%d ",aux1);
        }
        fprintf(new_arq,"\n");
    }
    //printf("Saiu4\n");
    ///Setup Cost
    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%d",&aux1);
                SetupCost[i][j][k] = aux1;
                if(i== rand1[k][0] && i!=j && j!=rand1[k][1])
                    SetupCost[i][j][k] = (SetupTime[i][j][k]*factor) + (rand()%25 - rand()%25);
                if(j== rand1[k][0] && i!=j && i!=rand1[k][1])
                    SetupCost[i][j][k] = (SetupTime[i][j][k]*factor) + (rand()%25 - rand()%25);
                if(i== rand1[k][1] && i!=j && j!=rand1[k][0])
                    SetupCost[i][j][k] = (SetupTime[i][j][k]*factor) + (rand()%25 - rand()%25);
                if(j== rand1[k][1] && i!=j && i!=rand1[k][0])
                    SetupCost[i][j][k] = (SetupTime[i][j][k]*factor) + (rand()%25 - rand()%25);
                fprintf(new_arq,"%d ",SetupCost[i][j][k]);
            }
            fprintf(new_arq,"\n");
        }
        //fprintf(new_arq,"\n");
    }
    //printf("Saiu5\n");


    r3 = fscanf(arq_ent,"%s",aux);
    fprintf(new_arq,"%s\n",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%d",&aux1);
            fprintf(new_arq,"%d ",aux1);
        }
        fprintf(new_arq,"\n");
    }
    //printf("Saiu6\n");

    fprintf(new_arq,"minl\n",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            fprintf(new_arq,"%d ",25);
        }
        fprintf(new_arq,"\n");
    }

    fclose(new_arq);
    fclose(arq_ent);
    printf("Saiu7\n");
}
