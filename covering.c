#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void printint64(long long int n)  {
    long long int M=n;
    int digit[20],i,len;
    
    if(M<0)  printf("-"),M=-M;
    for(i=0;i<20;i++)  digit[i]=M%10,M/=10;
    len=19;
    while((len>0)&&(digit[len]==0))  len--;
    while(len>=0) printf("%d",digit[len]),len--;
    printf("\n");
 
    return;
}

int gcd(int a,int b)
{
  int c;

  while(b>0)  {
     if(a>=b)  {
        a-=b;
        if(a>=b)  {
           a-=b;
           if(a>=b)  {
              a-=b;
              if(a>=b)  {
                 a-=b;
                 if(a>=b)  {
                    a-=b;
                    if(a>=b)  {
                       a-=b;
                       if(a>=b)  {
                          a-=b;
                          if(a>=b)  {
                             a-=b;
                             if(a>=b)  a%=b;
              }}}}}}}}
     c=a,a=b,b=c;
  }
  return a;
}

int single_modinv (int a,int modulus)

{ /* start of single_modinv */

 a%=modulus;

 int ps1, ps2, dividend, divisor, rem, q, t;

 int parity;

 q = 1;
 rem = a;
 dividend = modulus;
 divisor = a;
 ps1 = 1;
 ps2 = 0;
 parity = 0;

 while (divisor > 1)
 {
   rem = dividend - divisor;
   t = rem - divisor;
   if (t >= 0) {
     q += ps1;
     rem = t;
     t -= divisor;
     if (t >= 0) {
       q += ps1;
       rem = t;
       t -= divisor;
       if (t >= 0) {
         q += ps1;
         rem = t;
         t -= divisor;
         if (t >= 0) {
           q += ps1;
           rem = t;
           t -= divisor;
           if (t >= 0) {
             q += ps1;
             rem = t;
             t -= divisor;
             if (t >= 0) {
               q += ps1;
               rem = t;
               t -= divisor;
               if (t >= 0) {
                 q += ps1;
                 rem = t;
                 t -= divisor;
                 if (t >= 0) {
                   q += ps1;
                   rem = t;
                   if (rem >= divisor) {
                     q = dividend/divisor;
                     rem = dividend - q * divisor;
                     q *= ps1;
                   }}}}}}}}}
   q += ps2;
   parity = ~parity;
   dividend = divisor;
   divisor = rem;
   ps2 = ps1;
   ps1 = q;
 }

 if(parity==0)
   return (ps1);
 else
   return (modulus - ps1);
} /* end of single_modinv from Mersenneforum.org*/

int main()  {

    int N,b,w[32],A[32],L,LC,limit,temp,i,j,p,I,h,n,up,add,found;
    int pos,test,ord,all,p2,possible,all2,*S,*R,*pr,*primes,*isprime,np,ct,stored_LC[32];
    long long int best,res,C,P,**RES,**RES2;
    long long int alpha,alpha2,beta,v[64],M,K,E;
    
    pr=(int*)(malloc)(6542*sizeof(int));
    primes=(int*)(malloc)(6542*sizeof(int));
    isprime=(int*)(malloc)(65536*sizeof(int));

    for(i=0;i<65536;i++)  isprime[i]=1;
    isprime[0]=0;
    isprime[1]=0;
    for(n=2;n<256;n++)  {
        if(isprime[n])  {
           for(j=n*n;j<65536;j+=n)  isprime[j]=0;
        }
    }
    np=0;
    for(n=0;n<65536;n++)
        if(isprime[n])  primes[np]=n,np++;
    
    scanf("%d%d%I64d%d%I64d",&N,&b,&C,&limit,&best);
    
    R=(int*) (malloc) (N*sizeof(int));
    S=(int*) (malloc) (N*sizeof(int));
    
    L=0;
    for(I=0;I<limit;I+=65536)  {
        if(I==0)  {
           for(i=0;i<np;i++)  pr[i]=primes[i];
           ct=np;
        }
        else {
           up=I+65536;
           for(i=0;i<65536;i++)  isprime[i]=1;
           for(i=0;primes[i]*primes[i]<=up;i++)  {
               p=primes[i];
               for(j=((I+p-1)/p)*p-I;j<65536;j+=p)  isprime[j]=0;
           }
           ct=0;
           for(i=0;i<65536;i++)
               if(isprime[i])  pr[ct]=i+I,ct++;
        }
        for(h=0;(h<ct)&&(pr[h]<limit);h++)  {
        p=pr[h];
        res=1;
        ord=0;
        P=p;
        for(i=1;i<=N;i++)  {
            res=(long long int) res*b;
            res%=P;
            if((ord==0)&&(res==1))  ord=i;
        }
        if((res==1)&&(b%p!=1))  {
           test=1;
           for(i=2;i*i<=p;i++)  if(p%i==0)  test=0;
           if(test)  v[L]=p,w[L]=ord,L++;
        }
        }
    }
    for(i=0;i<L;i++)  {
        for(j=0;i+j+1<L;j++)  {
            if(w[j]>w[j+1])  {
               temp=v[j],v[j]=v[j+1],v[j+1]=temp;
               temp=w[j],w[j]=w[j+1],w[j+1]=temp;
    }}}

    RES=(long long int**) (malloc) (L*sizeof(long long int*));
    RES2=(long long int**) (malloc) (L*sizeof(long long int*));    
    
    printf("Checking k%c%d%cn",'*',b,'^');
    if(C>0) printf("%c%I64d",'+',C);
    else    printf("%I64d",C);
    printf(" sequence for exponent=%d, bound for primes in the covering set=%d, bound for k is ",N,limit),printint64(best);
    printf("Examining primes in the covering set: ");
    for(i=0;i<L;i++)  {
         printf("%I64d",v[i]);
         if(i!=L-1) printf(",");
         else printf("\n");
    }
    printf("And their orders: ");
    for(i=0;i<L;i++)  {
         printf("%d",w[i]);
         if(i!=L-1) printf(",");
         else printf("\n");
    }
   
    // i=0..L-1, j=0..N-1
    // RES[i][j]=Mod(b,v[i])^(-j)
    // RES2[i][j]=Mod(b,v[i])^j
    
    for(i=0;i<L;i++)  {
        RES[i]=(long long int*) (malloc) ((N+1)*sizeof(long long int));
        RES2[i]=(long long int*)(malloc) ((N+1)*sizeof(long long int));
        res=1;
        for(j=0;j<=N;j++)  {
            RES2[i][j]=res;
            RES[i][j]=single_modinv(RES2[i][j],v[i]);
            res=((long long int) res*b)%v[i];
        }
    }
    
    stored_LC[L-1]=w[L-1];
    for(i=L-2;i>=0;i--)  stored_LC[i]=gcd(w[i],stored_LC[i+1]);
    
    pos=0;
    for(i=0;i<L;i++)  A[i]=-1;
    while(pos>=0)  {
          for(i=0;i<N;i++)  R[i]=0;
          for(i=0;i<=pos;i++)  {
              if(A[i]>=0)  {
                 for(j=A[i];j<N;j+=w[i])  R[j]=1;
              }
          }
          all=0;
          for(i=0;i<N;i++)  all+=R[i];
          possible=0;
          for(i=pos+1;i<L;i++)  possible+=N/w[i];
          if(pos!=L-1)  {
             LC=stored_LC[pos+1];
             for(i=0;i<LC;i++)  S[i]=0;
             for(i=0;i<N;i++)
                 if(R[i]==0)  S[i%LC]=1;
             p2=0;
             for(i=0;i<LC;i++)  p2+=S[i];
          }
          else p2=0;
          
          if((possible+all<N)||(p2>L-pos))  {
              while((pos>=0)&&(A[pos]==w[pos]-1))  pos--;
              if(pos>=0)  A[pos]++;
          }
          else {
              alpha=0;
              beta=1;
              for(i=0;i<=pos;i++)  {
                  if(A[i]>=0)  {
                     alpha2=(RES[i][A[i]]*(-C))%v[i];
                     E=(((alpha2-alpha)%v[i])*single_modinv(beta%v[i],v[i]))%v[i];
                     if(E<0)  E+=v[i];
                     alpha=E*beta+alpha;
                     beta*=v[i];
                  }
              }
              K=alpha;
              if(all==N)  {
                 if(K<best)  {
                    while((K<best)&&(gcd((K+C)%(b-1),b-1)>1))  K+=beta;
                    if(K<best)  {
                        best=K;
                        found=1;
                        printf("**************** Solution found ****************\n");
                        printint64(K);
                    }
                 }
                 while((pos>=0)&&(A[pos]==w[pos]-1))  pos--;
                 if(pos>=0)  A[pos]++;
              }   
              else  {
                 M=1;
                 for(i=0;i<=pos;i++)
                      if(A[i]>=0)  M*=v[i];
                 if(best<M)  {
                     all2=all;
                     while(best>K)  {
                     all=all2;
                     for(i=0;i<N;i++)  {
                         if(R[i]==0)  {
                            add=0;
                            for(j=pos+1;j<L;j++)  {
                                E=((K%v[j])*RES2[j][i]+C)%v[j];
                                if(E<0) E+=v[j];
                                if(E==0)  {
                                   all++;
                                   add=1;
                                   break;
                                }
                            }
                            if(add==0)  break;
                         }
                     }
                     if((all==N)&&(K<best)&&(gcd((K+C)%(b-1),b-1)==1))  {
                         best=K;
                         found=1;
                         printf("**************** Solution found ****************\n");
                         printint64(K);
                     }
                     K+=M;
                     }
                     while((pos>=0)&&(A[pos]==w[pos]-1))  pos--;
                     if(pos>=0)  A[pos]++;
                 }
                 else {
                   if(pos!=L-1)  {
                      pos++;
                      A[pos]=-1;
                   }
                   else {
                      while((pos>=0)&&(A[pos]==w[pos]-1))  pos--;
                      if(pos>=0)  A[pos]++;
                   }
                 }           
             }
          }
    }
    for(i=0;i<L;i++)  free(RES[i]);
    free(RES);
    for(i=0;i<L;i++)  free(RES2[i]);
    free(RES2);
    free(R);
    free(S);
    
    return 0;
}
