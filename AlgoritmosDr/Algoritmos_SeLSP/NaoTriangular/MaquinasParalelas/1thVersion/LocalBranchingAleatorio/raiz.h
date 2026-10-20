#ifndef RAIZ_H_INCLUDED
#define RAIZ_H_INCLUDED

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define	INF    999999999.0
#define zero        0.00001
#define Max_T   16
#define Max_M   3
#define BigM    25

using namespace std;

class Instance
{
public:
    Instance() {}
    ~Instance() {}

    char FileInfo[150];
    char FileInfo2[150];
    char FileSaida[150];
    int T  ; // unidades de tempo
    int N ; // variedade de produtos
    int M ;// quantidade de máquinas
    int C; //capacidade
    int it;
    int contador;
    int bol;
    int iteracao ;
    int auxp;
    int pp;
    int auxfi;
    double solz;
    int totalite;
    double Solu;
    double MelhorSolu;

    std::vector< std::vector< std::vector< std::vector< double > > > > relax;
    std::vector< std::vector< std::vector< double > > > relaxz;
    std::vector< std::vector< std::vector< std::vector< double > > > > ZO; //recebe 2 se o valor yijkt é zero na solução, e -1 caso contrario
    std::vector< std::vector< std::vector< std::vector< double > > > > Z1; //recebe 2 se o valor yijkt é 1 na solução, e -1 caso contrario
    std::vector< std::vector< std::vector< double > > > YO; //recebe 2 se o valor zijt é zero na solução, e -1 caso contrario
    std::vector< std::vector< std::vector< double > > > Y1; //recebe 2 se o valor zijt é 1 na solução, e -1 caso contrario
    std::vector< int > vetq;
    int beta;

    double CLB;
    double CLBCo;
    double part;

    int car;

    /************************************************************************
                                   Parameters
    ************************************************************************/
    ///Paralell machine
    std::vector< std::vector< double > > Capacity;                      //plant capacity MxT
    std::vector< double > HoldingCost;                                  //holding cost N
    std::vector< std::vector< std::vector< double > > > SetupTime;      //setup time NxNxM
    std::vector< std::vector< std::vector< double > > > SetupCost;      //setup cost NxNxM
    std::vector< std::vector< double > > ProdTime;                      //production time NxM
    std::vector< std::vector< double > > ProdCost;                      //production time NxM
    std::vector< std::vector< double > > Demand;                        //production demand NxT
    std::vector< std::vector< double > > MaxProd;                       //production cost NxT
    std::vector< std::vector< int > > MachineProd;                   //binary variable - 1 indicates that machine k produced item i, 0 otherwise MxN
    std::vector< double> isetup;    //estado inicial da maquina para cada produto
    std::vector< double> istock;    //estoque inicial de cada produto
    std::vector< std::vector< int > > minl;      //tamanho minimo do lote
    std::vector< double> maxl;      //tamanho maximo do lote

    FILE *saida;
    //FILE *teste;
    //FILE *teste2;

    double ti,tf,tempo, tm, tempomeio;  // ti = tempo inicial // tf = tempo final


    timeval tempo_inicio,tempo_fim, tempo_meio;


    /********************************Functions********************************/
    void le_dados_SM(char * input_data);

    void inicializa_variaveis(void);

    void desaloca_variaveis(void);
};

#endif // RAIZ_H_INCLUDED
