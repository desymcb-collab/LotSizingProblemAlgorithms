#include <iostream>
#include "raiz.h"
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


ILOSTLBEGIN

typedef IloArray<IloNumVarArray> DoisD;
typedef IloArray<IloArray<IloNumVarArray> > TresD;
typedef IloArray<IloArray<IloArray <IloNumVarArray> > > QuatroD;

typedef IloArray<IloExprArray> Exp2D;
typedef IloArray<Exp2D> Exp3D;
typedef IloArray<Exp3D> Exp4D;
typedef IloArray<Exp4D> Exp5D;

long seed = 0;


int geraAleatorio(int min,int max)
{
    int r;
    // srand((unsigned)time(NULL));
    if(!seed)
    {
        time(&seed);
        srand((unsigned)seed);
    }

    r = min+rand()%(max-min+1);

    return r;

}

double resolve_subproblema(Instance *ins)
{
    IloEnv env;
//********************************************* Variáveis *********************************************//

    QuatroD x(env, ins->nbProdutos); // quantidade do item i produzido na máquina j no período t
    TresD   y(env, ins->nbProdutos); // binária: 1 se o item i é produzido na máquina j no período t e 0 caso contrário
    TresD   z(env, ins->nbProdutos);


    Exp2D sumProdTotal(env, ins->nbProdutos); // produção total por produto por período
    Exp2D usoMaqTotal(env,ins->nbMaquinas); // uso total de cada máquina por período
    Exp2D SomZ(env,ins->nbMaquinas); // uso total de cada máquina por período
    Exp3D sumProdMaq(env, ins->nbProdutos); // produção total por produto por período

    IloExpr SomZero(env);
    IloExpr SomUm(env);
    IloExpr SomZeroz(env);
    IloExpr SomUmz(env);

    int i, j, t, u;
    double soma, resul;

    if(ins->nbTemp == 18)
    {
      ins->CLB = 18*50;//6*(ins->nbProdutos/10)*ins->it;//(ins->nbTemp)*(ins->nbProdutos/10);
      ins->CLBCo = 18*50;//6*(ins->nbProdutos/10)*ins->it;//(ins->nbTemp)*(ins->nbProdutos/10);

    }
    else
    {
      ins->CLB = 6*(ins->nbProdutos/10);//(ins->nbTemp)*(ins->nbProdutos/10);
      ins->CLBCo = 6*(ins->nbProdutos/10);//(ins->nbTemp)*(ins->nbProdutos/10);
    }

     FILE *teste;
     teste = fopen("relax.txt", "w");

    for(i = 0; i < ins->nbProdutos; i++)
    {
        x[i]  = TresD(env,ins->nbMaquinas);
        y[i] = DoisD(env, ins->nbMaquinas);
        z[i] = DoisD(env, ins->nbMaquinas);
        sumProdTotal[i] =  IloExprArray(env, ins->nbTemp);
        sumProdMaq[i] = Exp2D(env, ins->nbMaquinas);


        for(j = 0; j < ins->nbMaquinas; j++)
        {
            x[i][j] = DoisD(env, ins->nbTemp);
            y[i][j] = IloNumVarArray(env, ins->nbTemp);
            usoMaqTotal[j] = IloExprArray(env, ins->nbTemp);
            sumProdMaq[i][j] = IloExprArray(env, ins->nbTemp);
            z[i][j] = IloNumVarArray(env, ins->nbTemp);
            SomZ[j] = IloExprArray(env, ins->nbTemp);


            for(t = 0; t < ins->nbTemp; t++)
            {
                x[i][j][t] = IloNumVarArray(env, ins->nbTemp);
                SomZ[j][t] =  IloExpr(env);


                if(ins->relax[i][j][t] == -1)
                {
                    y[i][j][t] = IloNumVar(env);
                     fprintf(teste,"y[%d][%d][%d] RELAXADO\n", i, j, t);
                }
                else
                {
                    y[i][j][t] = IloNumVar(env,0,1, ILOINT);
                    fprintf(teste,"y[%d][%d][%d] INTEIRO\n", i, j, t);
                }

                if(ins->relaxz[i][j][t] == -1)
                {
                    z[i][j][t] = IloNumVar(env);
                    fprintf(teste,"z[%d][%d][%d] RELAXADO\n", i, j, t);
                }
                else
                {
                    z[i][j][t] = IloNumVar(env,0,1, ILOINT);
                     fprintf(teste,"z[%d][%d][%d] INTEIRO\n", i, j, t);
                }



                sumProdTotal[i][t] = IloExpr(env);
                usoMaqTotal[j][t] = IloExpr(env);
                sumProdMaq[i][j][t] = IloExpr(env);

                for(u = 0; u < ins->nbTemp; u++)
                {
                    x[i][j][t][u] = IloNumVar(env);

                }

            }
        }
    }

      getchar();
      fclose(teste);

    IloExpr cost(env);
    IloModel model(env);
    IloCplex cplex(model);
    int conta = 0;

    //********************************************* Expressões *********************************************//

    for (i = 0; i < ins->nbProdutos; i++)
    {
        for (t = 0; t < ins->nbTemp; t++)
        {

            for (j = 0; j < ins->nbMaquinas; j++)
            {
                for (u = 0; u <= t; u++)
                {
                    sumProdTotal[i][t] += x[i][j][u][t];
                }
            }
        }
    }

    for (j = 0; j < ins->nbMaquinas; j++)
    {
        for (t = 0; t < ins->nbTemp; t++)
        {

            for (i = 0; i < ins->nbProdutos; i++)
            {
                usoMaqTotal[j][t] +=  ins->f[i][j]*y[i][j][t];

                SomZ[j][t] += z[i][j][t];

                for (u = t; u < ins->nbTemp; u++)
                {
                    usoMaqTotal[j][t] += ins->b[i][j]*x[i][j][t][u];

                }
            }

        }
    }

    for (i = 0; i < ins->nbProdutos; i++)
    {
        for (j = 0; j < ins->nbMaquinas; j++)
        {
            for (t = 0; t < ins->nbTemp; t++)
            {

                for (u = t; u < ins->nbTemp; u++)
                {
                    sumProdMaq[i][j][t] += x[i][j][t][u];
                }
            }
        }
    }

    if(ins->beta < ins->it)
    {
        conta = ins->it - ins->beta;

    }

    for (i = 0; i < ins->nbProdutos; i++)
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {

            for(t = (ins->it - 1); t >= conta; t--)
            {

                if(ins->ZO[i][j][ins->vetq[t]] == 2)
                {
                    SomZero += y[i][j][ins->vetq[t]];
                }

                if(ins->Z1[i][j][ins->vetq[t]] == 2)
                {
                    SomUm += (1 - y[i][j][ins->vetq[t]]);
                }
                if(ins->YO[i][j][ins->vetq[t]] == 2)
                {
                    SomZeroz += z[i][j][ins->vetq[t]];
                }

                if(ins->Y1[i][j][ins->vetq[t]] == 2)
                {
                    SomUmz += (1 - z[i][j][ins->vetq[t]]);
                }

            }

        }

    }
    //********************************************* Restrições *********************************************//
    for (i = 0; i < ins->nbProdutos; i++) // equações de balanço de estoque
    {

        for (t = 0; t < ins->nbTemp; t++)
        {
            model.add(sumProdTotal[i][t]  == ins->d[i][t]);
        }

    }

    for (j = 0; j < ins->nbMaquinas; j++) // limitação da capacidade de cada máquina em cada período
    {

        for (t = 0; t < ins->nbTemp; t++)
        {
            model.add(usoMaqTotal[j][t] <= ins->Cap[j][t]);
            model.add(SomZ[j][t] == 1);
        }

    }

    if(ins->it > ins->beta)
    {
        for(i = 0; i < ins->nbProdutos; i++) // Custo de produção + Custo de preparação
        {
            for(j = 0; j < ins->nbMaquinas; j++)
            {
                for(t = 0; t < (ins->it-ins->beta); t++)
                {
                    model.add(y[i][j][ins->vetq[t]] == ins->relax[i][j][ins->vetq[t]]);
                    model.add(z[i][j][ins->vetq[t]] == ins->relaxz[i][j][ins->vetq[t]]);

                }
            }
        }

    }

    for(i = 0; i < ins->nbProdutos; i++)
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {
            for(t = 0; t < ins->nbTemp; t++ )
            {

                if(t==0)
                {
                    model.add(z[i][j][t] - y[i][j][t] <= 0);
                }
                else
                {
                    model.add(z[i][j][t] - y[i][j][t] - z[i][j][t-1] <= 0);
                }

                for(u = t; u < ins->nbTemp; u++ )
                {
                    if(t != 0)
                    {
                        model.add(x[i][j][t][u] <= ins->mini[i][j][t][u]*(y[i][j][t] + z[i][j][t-1]));
                    }
                    else
                    {
                        model.add(x[i][j][t][u] <= 999999999*y[i][j][t]);//ins->mini[i][j][t][u]*y[i][j][t]);
                    }
                }

            }
        }

    }

    int q;

    for(i = 0; i < ins->nbProdutos; i++) // Configuraço setup carry over
    {

        for(j = 0; j < ins->nbMaquinas; j++)
        {
            for(t = 0; t < ins->nbTemp; t++ )
            {

                for(q = 0; q < ins->nbProdutos; q++)
                {

                    if(q != i)
                    {
                        if(t!=0)
                            model.add(z[i][j][t] + z[i][j][t-1] - y[i][j][t] + y[q][j][t] <= 2);
                        else
                            model.add(z[i][j][t]  - y[i][j][t] + y[q][j][t] <= 2);

                    }

                }


            }
        }
    }

    if (ins->it != 0)
    {
        model.add(SomZero + SomUm <= ins->CLB);
        model.add(SomZeroz + SomUmz <= ins->CLBCo);
    }

    //********************************************* Custos da função objetivo *********************************************/
    for(i = 0; i < ins->nbProdutos; i++) // Custo de produção + Custo de preparação
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {
            for(t = 0; t < ins->nbTemp; t++)
            {

                cost+=  ins->s[i][j][t]*y[i][j][t];

            }
        }
    }


    for(i = 0; i < ins->nbProdutos; i++) // Custos de estoque
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {

            for(t = 0; t < ins->nbTemp; t++)
            {
                for(u = 0; u < ins->nbTemp; u++)
                {

                    cost+= ins->CustoAux[i][j][t][u]*x[i][j][t][u];
                }
            }
        }

    }

    //********************************************* Função Objetivo *********************************************/
    model.add(IloMinimize(env, cost));


    if(ins->nbTemp == 18)
            cplex.setParam(IloCplex::TiLim, 100);
    if(ins->nbTemp == 12)
            cplex.setParam(IloCplex::TiLim, 200);
    if(ins->nbTemp == 6)
            cplex.setParam(IloCplex::TiLim, 200);


    //********************************************* Resultados *********************************************/
    if (cplex.solve())
    {
        resul = cplex.getObjValue();

        if(ins->bol == 1)
            ins->Solu = cplex.getObjValue();


        for(i = 0; i < ins->nbProdutos; i++)
        {
            for(j = 0; j < ins->nbMaquinas; j++)
            {
                for(t = 0; t <= ins->it; t++)
                {
                    ins->relax[i][j][ins->vetq[t]] = cplex.getValue(y[i][j][ins->vetq[t]]);
                    ins->relaxz[i][j][ins->vetq[t]] = cplex.getValue(z[i][j][ins->vetq[t]]);

                }
            }
        }

        cout << " Funcao Objetivo = " << cplex.getObjValue() << endl;
        getchar();

        if (ins->bol == 1)
        {
           FILE *saida2;

          saida2 = fopen("saida2.txt","w");

            for(i = 0; i < ins->nbProdutos; i++)
            {
                fprintf(saida2, " Produto %d \n\n", i);
                for(j = 0; j < ins->nbMaquinas; j++)
                {
                    fprintf(saida2, " Máquina %d \n", j);

                    for(t = 0; t < ins->nbTemp; t++)
                    {
                        fprintf (saida2, "y[%d][%d][%d] = %.2lf\n",i,j,t,cplex.getValue(y[i][j][t]));

                    }
                    fprintf (saida2, "\n");
                }

                fprintf (saida2, "\n\n");


            }

            for(i = 0; i < ins->nbProdutos; i++)
            {
                fprintf(saida2, " Produto %d \n\n", i);
                for(j = 0; j < ins->nbMaquinas; j++)
                {
                    fprintf(saida2, " Máquina %d \n", j);

                    for(t = 0; t < ins->nbTemp; t++)
                    {

                        fprintf (saida2, "z[%d][%d][%d] = %.2lf\n",i,j,t,cplex.getValue(z[i][j][t]));

                    }
                    fprintf (saida2, "\n");
                }

                fprintf (saida2, "\n\n");


            }
        }

        env.end();

        return resul;

    }

    else
    {
        cout << " No solution found" << endl;
        env.end();
        return -1;
    }

}

void particao(int p, Instance *ins)
{
    int i, j;
    for(i = 0; i < ins->nbProdutos; i++)
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {

            ins->relax[i][j][ins->vetq[p]] = -2; //Inteira
            ins->relaxz[i][j][ins->vetq[p]] = -2; //Inteira

        }

    }
}

void AtualizaConjuntos(Instance *ins)
{
    int i, j, t;
    for(i = 0; i < ins->nbProdutos; i++ )
    {
        for(j = 0; j < ins->nbMaquinas; j++)
        {
            for(t = 0; t <= ins->it; t++)
            {

                if(ins->relax[i][j][ins->vetq[t]] == 0)
                {
                    ins->ZO[i][j][ins->vetq[t]] = 2;

                }
                if (ins->relax[i][j][ins->vetq[t]] == 1)
                {
                    ins->Z1[i][j][ins->vetq[t]] = 2;
                }
                if(ins->relaxz[i][j][ins->vetq[t]] == 0)
                {
                    ins->YO[i][j][ins->vetq[t]] = 2;

                }
                if (ins->relaxz[i][j][ins->vetq[t]] == 1)
                {
                    ins->Y1[i][j][ins->vetq[t]] = 2;
                }


            }
        }
    }


}

int verifica(int vet[], int n)
{
    int j, c, r = 0;

    for(j = 0; j < n; j++)
    {
        if(vet[n] == vet[j])
        {
            r = 1;
            j = n + 1;
        }

    }

    return r;
}

void RelaxAndFix(Instance * ins)
{
    int i,j,t,k, qtdade, qtot, ini, veri, ca;
    double resp, alfa;
    double cond;
    double lim;

    lim = 0.05;

    ins->bol = 0;
    resp = 0;
    ins->totalite = 0;
    ins->Solu = 99999999;
    ins->MelhorSolu = 99999999;
    ins->tempomeio = 0;

    cond = double((abs(ins->solz - 1)));

    while(cond > lim && ins->tempomeio < 1800)
    {
        ins->bol = 0;
        resp = 0;
        ins->totalite = ins->totalite + 1;
        // printf("cond = %lf");
        // getchar();

        // srand((unsigned)time(NULL));
        //   alfa =((double)(rand())/((double)(RAND_MAX)*1)); // gera número aleatório entre 0 e 1
        if(ins->nbTemp == 18)
            qtdade = geraAleatorio(2, ins->nbTemp/3);
        if(ins->nbTemp == 12)
            qtdade = geraAleatorio(4, ins->nbTemp/2);
        if(ins->nbTemp == 6)
            qtdade = geraAleatorio(1, ins->nbTemp/2);

        /* if(ins->nbTemp < 18)
             qtdade = ceil(alfa*(ins->c)/1.5); // qtdade = maior inteiro alfa*nbTemp (qtade vai de 1 até nbTemp)
         else
             qtdade = ceil(alfa*(ins->nbTemp)/1.5); */


        ini = 1;
        ins->vetq[0] = 0;

        for(t = 1; t < ins->nbTemp; t++)
        {
            ins->vetq[t] = -1;
        }

        if(qtdade < ins->nbTemp)
            qtot = qtdade + 1;
        else
            qtot = qtdade;



        veri = 1;
        srand((unsigned)time(NULL));

        do
        {
            for(i = ini; i < qtot; i ++)
            {
                do
                {
                    ins->vetq[i] = geraAleatorio(ini,qtot-1); //gera aleatório entre ini e qtot
                    veri = verifica(ins->vetq, i); //verifica se vetq[i] já foi gerado, se sim gera outro numero
                }
                while (veri == 1);

            }

            ini = ini + qtdade;
            qtot = qtot + qtdade;

            if(qtot > ins->nbTemp)
            {
                qtot = ins->nbTemp;
            }

        }
        while(ins->vetq[ins->nbTemp-1] == -1);


        ins->auxp = qtdade;
        ins->auxfi = qtdade;
        ins->pp = 0;
        srand((unsigned)time(NULL));
        if(qtdade  < ins->nbTemp)
            ins->beta = qtdade ;//;= geraAleatorio(qtdade, 2*qtdade);
        else
            ins->beta = ins->nbTemp-1;//geraAleatorio(qtdade, ins->nbTemp-1);


         printf("qtdade = %d", qtdade);
         printf("beta = %d", ins->beta);
            for(t = 0; t < ins->nbTemp; t++ )
                 {
                   printf("v[%d] = %d   ", t, ins->vetq[t]);
                 }
         getchar();


        for(i = 0; i < ins->nbProdutos; i++)
        {
            for(j = 0; j < ins->nbMaquinas; j++)
            {
                for(t = 0; t < ins->nbTemp; t++ )
                {
                    ins->ZO[i][j][t] = -1;
                    ins->Z1[i][j][t] = -1;
                    ins->relax[i][j][t] = -1;
                    ins->YO[i][j][t] = -1;
                    ins->Y1[i][j][t] = -1;
                    ins->relaxz[i][j][t] = -1;
                }
            }
        }

        resp = 0;
        /*gettimeofday(&(ins->tempo_meio),NULL);
        ins->tm = (double)(ins->tempo_meio).tv_usec + ((double)(ins->tempo_meio).tv_sec * (1000000.0));
        ins->tempomeio = (ins->tm - ins->ti) / 1000;
        ins->tempomeio =  ins->tempomeio/1000; */
        t = 0;

        while(t < ins->nbTemp)
        {
            ins->it = t;


            if(resp != -1 && ins->tempomeio < 1800)
            {

                ins->tempomeio = 0;
                if(t == ins->nbTemp-1)
                {
                    ins->bol = 1;
                }

                particao(t, ins);
                resp = resolve_subproblema(ins);
                AtualizaConjuntos(ins);
                t++;
                printf("ENTREIIIII");
                getchar();

            }

            else
            {
                t = ins->nbTemp;
                printf("ite =  %d\n", ins->it);
                printf("SOLUÇÃO NÃO ENCONTRADA!");
                getchar();

            }

            gettimeofday(&(ins->tempo_meio),NULL);
            ins->tm = (double)(ins->tempo_meio).tv_usec + ((double)(ins->tempo_meio).tv_sec * (1000000.0));
            ins->tempomeio = (ins->tm - ins->ti) / 1000;
            ins->tempomeio =  ins->tempomeio/1000;
            printf("tempo = %lf", ins->tempomeio);
            //getchar();

        }

        if(ins->Solu < ins->MelhorSolu && resp != -1)
            ins->MelhorSolu = ins->Solu;

        if(resp != -1)
            cond = double((abs(ins->solz - ins->Solu)))/ abs(ins->Solu);


    }

    printf("iteracoes = %d", ins->totalite);
    gettimeofday(&ins->tempo_fim,NULL);
    ins->tf = (double)ins->tempo_fim.tv_usec + ((double)ins->tempo_fim.tv_sec * (1000000.0));
    ins->tempo = (ins->tf - ins->ti) / 1000;
    ins->tempo = ins->tempo/1000;
    // printf("Tempo gasto em segundos %.3f\n",ins->tempo); */


}

int main()
{
    try
    {
        //********************************************* Fixos *********************************************//
        Instance ins;
        int itensmaq[4][2],periods[3];

        itensmaq[4][2];

        itensmaq[0][0] = 2;
        itensmaq[0][1] = 6;
        itensmaq[1][0] = 4;
        itensmaq[1][1] = 12;
        itensmaq[2][0] = 6;
        itensmaq[2][1] = 50;
        //itensmaq[3][0] = 6;
        //itensmaq[3][1] = 50;

        /*items[0] = 6;
        items[1] = 12;
        items[2] = 25;
        items[3] = 50;

        plants[0] = 2;
        plants[1] = 4;
        plants[2] = 6; */

        periods[0] = 6;
        periods[1] = 12;
        periods[2] = 18;


        for(int i = 2; i < 3; i++)
            // for(int j = 0; j < 2; j++)
            for(int t = 2; t < 3; t++)
                for(int contador = 0; contador <= 9; contador++)
                {
                    /* sprintf(ins.FileInfo,"exemplos/%d-%d/NBB0/NBB0%d%d.dat",plants[j],items[i],contador,periods[t]);
                     sprintf(ins.FileSaida,"exemplosSaida/%d-%d/NBB0/NBB0%d%d_saida.dat",plants[j],items[i],contador,periods[t]);
                     sprintf(ins.FileInfo2,"exemplosSaidaRela/%d-%d/NBB0/NBB0%d%d_saida.dat",plants[j],items[i],contador,periods[t]); */
                    sprintf(ins.FileInfo,"exemplos/%d-%d/AAA0/AAA0%d%d.dat",itensmaq[i][0],itensmaq[i][1],contador,periods[t]);
                    sprintf(ins.FileSaida,"exemplosSaida/%d-%d/AAA0/AAA0%d%d_saida.dat",itensmaq[i][0],itensmaq[i][1],contador,periods[t]);
                    sprintf(ins.FileInfo2,"exemplosSaidaRela/%d-%d/AAA0/AAA0%d%d_saida.dat",itensmaq[i][0],itensmaq[i][1],contador,periods[t]);

                    ins.car = 1;
                    //printf("%s\n",ins.FileInfo);
                    //  getchar();


                    ins.le_dados(ins.FileInfo);
                    ins.le_dados2(ins.FileInfo2);
                   // printf("Capacidade = %lf", ins.d[49][17]);
                   // getchar();

                    if(ins.car == 1)
                    {
                        ins.ti = ins.tf = ins.tempo = ins.tm = ins.tempomeio  = 0;

                        gettimeofday(&ins.tempo_inicio,NULL);


                        ins.ti = (double)ins.tempo_inicio.tv_usec + ((double)ins.tempo_inicio.tv_sec * (1000000.0));

                        RelaxAndFix(&ins);
                         FILE *saida;
                         saida = fopen(ins.FileSaida, "w");


                        fprintf(saida, "%lf\n", ins.MelhorSolu); //fo
                        fprintf(saida, "%lf\n", ins.tempo);//tempo
                        fprintf(saida, "%d\n", ins.totalite);//iterações

                        fclose(saida);
                    }

                    ins.Liberar_matriz4D(ins.nbProdutos, ins.nbMaquinas, ins.nbTemp, ins.CustoAux);
                    ins.Liberar_matriz4D(ins.nbProdutos, ins.nbMaquinas, ins.nbTemp, ins.mini);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.s);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.c);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.relax);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.relaxz);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.ZO);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.Z1);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.YO);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbMaquinas, ins.Y1);
                    ins.Liberar_matriz3D(ins.nbProdutos, ins.nbTemp, ins.hs);
                    ins.Liberar_matriz2D(ins.nbProdutos, ins.nbTemp, ins.h);
                    ins.Liberar_matriz2D(ins.nbProdutos, ins.nbTemp, ins.d);
                    ins.Liberar_matriz2D(ins.nbProdutos, ins.nbMaquinas, ins.b);
                    ins.Liberar_matriz2D(ins.nbProdutos, ins.nbMaquinas, ins.f);
                    ins.Liberar_matriz2D(ins.nbMaquinas, ins.nbTemp, ins.Cap);
                    ins.Liberar_matriz2D(ins.nbProdutos, ins.nbMaquinas, ins.aux);
                    free(ins.vetq);
                }
    }

    catch (IloException & ex)
    {
        cerr << "Error: " << ex << endl;
    }
    catch (...)
    {
        cerr << "Error" << endl;
    }


    return 0;
}
