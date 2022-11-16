#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"


int gcd(int a,int b)
{
  if(a<0)  a=-a;
  if(b<0)  b=-b;
  if(a==0)  return b;
  if(b==0)  return a;
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
 
 if(modulus<0)  modulus=-modulus;
 a%=modulus;
 if(a<0)  a+=modulus;

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

    int e,n,v[256],w[256],a[256],l,lc,temp,temp2,add,g,h,i,j,p,up,found,inv,limit,**res,**res2;
    int pos,ord,all,p2,possible,all2,*s,*r,*pr,*primes,*isprime,np,ct,stored_lc[64];
    
    mpz_t B;
    mpz_t C;
    mpz_t E;
    mpz_t G;
    mpz_t K;
    mpz_t M;
    mpz_t P;
    mpz_t RES;
    mpz_t BEST;
    mpz_t ALPHA;
    mpz_t ALPHA2;
    mpz_t BETA;
    mpz_t TEMP;
    mpz_t MINUS_C;
    mpz_t BMINUSONE;
    mpz_init(B);
    mpz_init(C);
    mpz_init(E);
    mpz_init(G);
    mpz_init(K);
    mpz_init(M);
    mpz_init(P);
    mpz_init(RES);
    mpz_init(BEST);
    mpz_init(ALPHA);
    mpz_init(ALPHA2);
    mpz_init(BETA);
    mpz_init(TEMP);
    mpz_init(MINUS_C);
    mpz_init(BMINUSONE);
    
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
    
    gmp_scanf("%d%Zd%Zd%d%Zd",&n,&B,&C,&limit,&BEST);
    mpz_neg(MINUS_C,C);
    mpz_sub_ui(BMINUSONE,B,1);
    
    r=(int*) (malloc) (n*sizeof(int));
    s=(int*) (malloc) (n*sizeof(int));
    
    l=0;
    for(g=0;g<limit;g+=65536)  {
        if(g==0)  {
           for(i=0;i<np;i++)  pr[i]=primes[i];
           ct=np;
        }
        else {
           up=g+65536;
           for(i=0;i<65536;i++)  isprime[i]=1;
           for(i=0;primes[i]*primes[i]<=up;i++)  {
               p=primes[i];
               for(j=((g+p-1)/p)*p-g;j<65536;j+=p)  isprime[j]=0;
           }
           ct=0;
           for(i=0;i<65536;i++)
               if(isprime[i])  pr[ct]=i+g,ct++;
        }
        for(h=0;(h<ct)&&(pr[h]<limit);h++)  {
        p=pr[h];
        mpz_set_ui(P,p);
        mpz_powm_ui(RES,B,n,P);
        if(mpz_cmp_ui(RES,1)==0)  {
           temp=mpz_mod_ui(RES,B,p);
           if(temp!=1)  {
              mpz_mod(RES,B,P);
              ord=1;
              while(mpz_cmp_ui(RES,1)!=0)  {
                    ord++;
                    mpz_mul(RES,RES,B);
                    mpz_mod(RES,RES,P);
             }
             v[l]=p;
             w[l]=ord;
             l++;             
          }
        }
        }
    }
    for(i=0;i<l;i++)  {
        for(j=0;i+j+1<l;j++)  {
            if(w[j]>w[j+1])  {
               temp=v[j],v[j]=v[j+1],v[j+1]=temp;
               temp=w[j],w[j]=w[j+1],w[j+1]=temp;
    }}}

    res=(int**) (malloc) (l*sizeof(int*));
    res2=(int**) (malloc) (l*sizeof(int*));    
    
    gmp_printf("Checking k%c%Zd%cn",'*',B,'^');
    if(mpz_sgn(C)>0) gmp_printf("%c%Zd",'+',C);
    else             gmp_printf("%Zd",C);
    gmp_printf(" sequence for exponent=%d, bound for primes in the covering set=%d, bound for k is %Zd\n",n,limit,BEST);
    printf("Examining primes in the covering set: ");
    for(i=0;i<l;i++)  {
         printf("%d",v[i]);
         if(i!=l-1) printf(",");
         else printf("\n");
    }
    printf("And their orders: ");
    for(i=0;i<l;i++)  {
         printf("%d",w[i]);
         if(i!=l-1) printf(",");
         else printf("\n");
    }
   
    // i=0..L-1, j=0..N-1
    // RES[i][j]=Mod(b,v[i])^(-j)
    // RES2[i][j]=Mod(b,v[i])^j
    
    for(i=0;i<l;i++)  {
        res[i]=(int*) (malloc) ((n+1)*sizeof(int));
        res2[i]=(int*)(malloc) ((n+1)*sizeof(int));
        mpz_set_ui(RES,1);
        mpz_set_ui(P,v[i]);
        for(j=0;j<=n;j++)  {
            res2[i][j]=mpz_get_ui(RES);
            res[i][j]=single_modinv(res2[i][j],v[i]);
            mpz_mul(RES,RES,B);
            mpz_mod(RES,RES,P);
        }
    }
    
    stored_lc[l-1]=w[l-1];
    for(i=l-2;i>=0;i--)  stored_lc[i]=gcd(w[i],stored_lc[i+1]);
    
    pos=0;
    for(i=0;i<l;i++)  a[i]=-1;
    while(pos>=0)  {
          for(i=0;i<n;i++)  r[i]=0;
          for(i=0;i<=pos;i++)  {
              if(a[i]>=0)  {
                 for(j=a[i];j<n;j+=w[i])  r[j]=1;
              }
          }
          all=0;
          for(i=0;i<n;i++)  all+=r[i];
          possible=0;
          for(i=pos+1;i<l;i++)  possible+=n/w[i];
          if(pos!=l-1)  {
             lc=stored_lc[pos+1];
             for(i=0;i<lc;i++)  s[i]=0;
             for(i=0;i<n;i++)
                 if(r[i]==0)  s[i%lc]=1;
             p2=0;
             for(i=0;i<lc;i++)  p2+=s[i];
          }
          else p2=0;
          
          if((possible+all<n)||(p2>l-pos))  {
              while((pos>=0)&&(a[pos]==w[pos]-1))  pos--;
              if(pos>=0)  a[pos]++;
          }
          else {
              mpz_set_ui(ALPHA,0);
              mpz_set_ui(BETA,1);
              for(i=0;i<=pos;i++)  {
                  if(a[i]>=0)  {
                     mpz_mul_ui(ALPHA2,MINUS_C,res[i][a[i]]);
                     temp=mpz_mod_ui(ALPHA2,ALPHA2,v[i]);
                     mpz_sub(E,ALPHA2,ALPHA);
                     temp=mpz_mod_ui(E,E,v[i]);
                     
                     temp2=mpz_mod_ui(TEMP,BETA,v[i]);
                     inv=single_modinv(temp2,v[i]);
                     
                     mpz_mul_ui(E,E,inv);
                     e=mpz_mod_ui(E,E,v[i]);
                     
                     mpz_addmul_ui(ALPHA,BETA,e);
                     mpz_mul_ui(BETA,BETA,v[i]);
                  }
              }
              mpz_set(K,ALPHA);
              if(all==n)  {
                 if(mpz_cmp(K,BEST)<0)  {
                    while(1)  {
                          if(mpz_cmp(K,BEST)>=0)  break;
                          mpz_add(TEMP,K,C);
                          mpz_gcd(G,TEMP,BMINUSONE);
                          if(mpz_cmp_ui(G,1)==0)  break;
                          mpz_add(K,K,BETA);
                    }
                    if(mpz_cmp(K,BEST)<0)  {
                        mpz_add(TEMP,K,C);
                        mpz_gcd(G,TEMP,BMINUSONE);
                        if(mpz_cmp_ui(G,1)==0)  {
                            mpz_set(BEST,K);
                            found=1;
                            printf("**************** Solution found ****************\n");
                            gmp_printf("%Zd\n",K);
                        }
                    }
                 }
                 while((pos>=0)&&(a[pos]==w[pos]-1))  pos--;
                 if(pos>=0)  a[pos]++;
              }   
              else  {
                 mpz_set_ui(M,1);
                 for(i=0;i<=pos;i++)
                      if(a[i]>=0)  mpz_mul_ui(M,M,v[i]);
                 if(mpz_cmp(BEST,M)<0)  {
                     all2=all;
                     while(mpz_cmp(BEST,K)>0)  {
                     all=all2;
                     for(i=0;i<n;i++)  {
                         if(r[i]==0)  {
                            add=0;
                            for(j=pos+1;j<l;j++)  {
                                temp=mpz_mod_ui(RES,K,v[j]);
                                mpz_mul_ui(RES,RES,res2[j][i]);
                                mpz_add(RES,RES,C);
                                if(mpz_divisible_ui_p(RES,v[j])>0)  {
                                   all++;
                                   add=1;
                                   break;
                                }
                            }
                            if(add==0)  break;
                         }
                     }
                   if(all==n)  {
                      mpz_add(TEMP,K,C);
                      mpz_gcd(G,TEMP,BMINUSONE);
                      if(mpz_cmp_ui(G,1)==0)  {
                         mpz_set(BEST,K);
                         found=1;
                         printf("**************** Solution found ****************\n");
                         gmp_printf("%Zd\n",K);
                      }
                    }
                    mpz_add(K,K,M);
                    }
                    while((pos>=0)&&(a[pos]==w[pos]-1))  pos--;
                    if(pos>=0)  a[pos]++;
                 }
                 else {
                   if(pos!=l-1)  {
                      pos++;
                      a[pos]=-1;
                   }
                   else {
                      while((pos>=0)&&(a[pos]==w[pos]-1))  pos--;
                      if(pos>=0)  a[pos]++;
                   }
                 }           
             }
          }
    }
    for(i=0;i<l;i++)  free(res[i]);
    free(res);
    for(i=0;i<l;i++)  free(res2[i]);
    free(res2);
    free(r);
    free(s);
    mpz_clear(B);
    mpz_clear(C);
    mpz_clear(E);
    mpz_clear(G);
    mpz_clear(K);
    mpz_clear(M);
    mpz_clear(P);
    mpz_clear(RES);
    mpz_clear(BEST);
    mpz_clear(ALPHA);
    mpz_clear(ALPHA2);
    mpz_clear(BETA);
    mpz_clear(TEMP);
    mpz_clear(MINUS_C);
    mpz_clear(BMINUSONE);
    
    return 0;
}
