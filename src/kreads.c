/* kreads patterns.txt < fasta : prints the single-line FASTA records whose sequence contains any pattern (all length K<=31, ACGT). Patterns should include both strands. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define H 65536
static uint64_t key[H]; static int idx[H]; static int K;
static int code(int c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;}return -1;}
static unsigned hh(uint64_t x){x^=x>>33;x*=0xff51afd7ed558ccdULL;x^=x>>33;return x&(H-1);}
int main(int c,char**v){FILE*f=fopen(v[1],"r");char b[256];int n=0;memset(idx,-1,sizeof idx);
 while(fgets(b,256,f)){b[strcspn(b,"\r\n")]=0;if(!*b)continue;if(!K)K=strlen(b);uint64_t x=0;for(int i=0;i<K;i++)x=(x<<2)|code(b[i]);unsigned h=hh(x);while(idx[h]>=0&&key[h]!=x)h=(h+1)&(H-1);if(idx[h]<0){key[h]=x;idx[h]=n;}n++;}
 uint64_t mask=((1ULL<<(2*K))-1);size_t cap=1<<20;char*line=malloc(cap),*hd=malloc(1<<16);hd[0]=0;ssize_t L;
 while((L=getline(&line,&cap,stdin))>0){if(line[0]=='>'){strncpy(hd,line,(1<<16)-1);continue;}
  uint64_t x=0;int run=0,ok=0;for(ssize_t i=0;i<L&&!ok;i++){int d=code(line[i]);if(d<0){run=0;x=0;continue;}x=((x<<2)|d)&mask;if(++run>=K){unsigned h=hh(x);while(idx[h]>=0){if(key[h]==x){ok=1;break;}h=(h+1)&(H-1);}}}
  if(ok){fputs(hd,stdout);fputs(line,stdout);}}
 return 0;}
