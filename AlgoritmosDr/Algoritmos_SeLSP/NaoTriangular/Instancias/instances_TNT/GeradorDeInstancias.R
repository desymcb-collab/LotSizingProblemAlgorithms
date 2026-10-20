geraTodas <- function(tipo1,tipo2,tipo3,tipo4,tipo5)
{
  for(p in seq(12,18,by=6))
    for(r in (0:4)){
      for(j in seq(2,2,by=1)){
        M = j
        for(i in seq(30,50,by=20)){
          N = i
          n = c(i,T)
          fname = sprintf("Data%d-%d-%d-%d-%d-%d-%d-%d-%d.dat",j,i,p,tipo1,tipo2,tipo3,tipo4,tipo5,r)
          write(n,fname,ncolumns = (M*N),append = TRUE,sep = " ")
          write(j,fname,append = TRUE,sep = " ")
          custos = matrix(0,(M*N),4)
          hpdc = c(1:(M*N))
          rpdc = matrix(0,M,N)
          dpdc = matrix(0,T,(M*N))
          cap = c(1:(M))
          
            for(jj in (1:j)){
              for(ii in (1:i)){
                custos[(((jj-1)*i)+ii),1] = round(runif(1,1.0,5.0),digits=1)
                custos[(((jj-1)*i)+ii),2] = round(runif(1,10.0,50.0) * tipo1,digits=1)
                custos[(((jj-1)*i)+ii),3] = round(runif(1,20.0,100.0) * tipo2,digits=1)
                custos[(((jj-1)*i)+ii),4] = round(runif(1,1.5,2.5),digits=1) 
              }
            }
          
          for(jj in (1:j)){
            for(ii in (1:i)){
             hpdc[(((jj-1)*i)+ii)] = round(runif(1,0.2,0.4),digits=1)
            }
          }
          
          for(tt in (1:T)){  
            for(jj in (1:j)){
              for(ii in (1:i)){
                dpdc[tt,(((jj-1)*i)+ii)] = ceiling(runif(1,0,100))
              }
            }
          }
          
          for(jj in (1:j))
            for(ii in (1:jj))
              if(ii != jj)
                rpdc[jj,ii] = round(runif(1,0.2,0.4),digits=1)
              
          soma =0
          for(jj in (1:j)){
            for(tt in (1:T)){  
              for(ii in (1:i)){
                soma = soma + ((custos[(((jj-1)*i)+ii),1]*(dpdc[tt,(((jj-1)*i)+ii)])) + custos[(((jj-1)*i)+ii),2])
              }
            }
            cap[jj] = ceiling((soma/T)*tipo3)
          }
          
          for(jj in (1:j))
            write(cap[jj],fname,append = TRUE,sep = " ")
          
          for(jj in (1:j))
            for(ii in (1:i) )
              write(custos[(((jj-1)*i)+ii),],fname,append = TRUE,sep = " ")
          
          write(hpdc,fname,ncolumns = (M*N),append = TRUE,sep = " ")
          
          for(tt in (1:T)) 
            write(dpdc[tt,],ncolumns = (M*N),fname,append = TRUE,sep = " ")
            
            for(jj in (1:j))
              for(ii in (1:jj))
                if(ii != jj)
                  write(rpdc[jj,ii],fname,append = TRUE,sep = " ")
        }
        
      }
    }
}

tipo1 = 0.8      #capacidade de utilizacao
tipo2 = 0.5      #
tipo3 = 0.75
tipo4 = 0
tipo5 = 0
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "AAB"
tipo1 = 1
tipo2 = 10
tipo3 = 0.75
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "ABA"
tipo1 = 10
tipo2 = 1
tipo3 = 0.75
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "ABB"
tipo1 = 1
tipo2 = 1
tipo3 = 0.75
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "NAA"
tipo1 = 10
tipo2 = 10
tipo3 = 0.85
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "NAB"
tipo1 = 1
tipo2 = 10
tipo3 = 0.85
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "NBA"
tipo1 = 10
tipo2 = 1
tipo3 = 0.85
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

classe = "NBB"
tipo1 = 1
tipo2 = 1
tipo3 = 0.85
geraTodas(tipo1,tipo2,tipo3,tipo4,tipo5)

