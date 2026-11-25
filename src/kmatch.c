/* kmatch patterns.txt < fasta : counts exact occurrences of each pattern (all length K<=31, ACGT only) in FASTA sequence lines; each occurrence counted (reads are single lines). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define H 65536
static uint64_t key[H]; static int idx[H]; static char *pat[4096]; static uint64_t cnt[4096]; static int K;
static int code(int c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;}return -1;}
static unsigned hh(uint64_t x){x^=x>>33;x*=0xff51afd7ed558ccdULL;x^=x>>33;return x&(H-1);}
static int enc(const char*s,uint64_t*o){uint64_t v=0;for(int i=0;i<K;i++){int c=code(s[i]);if(c<0)return 0;v=(v<<2)|c;}*o=v;return 1;}
int main(int c,char**v){FILE*f=fopen(v[1],"r");char b[256];int n=0;memset(idx,-1,sizeof idx);
 while(fgets(b,256,f)){b[strcspn(b,"\r\n")]=0;if(!*b)continue;if(!K)K=strlen(b);pat[n]=strdup(b);uint64_t x;if(!enc(b,&x)){fprintf(stderr,"bad pattern\n");return 2;}unsigned h=hh(x);while(idx[h]>=0&&key[h]!=x)h=(h+1)&(H-1);if(idx[h]<0){key[h]=x;idx[h]=n;}n++;}
 uint64_t mask=(K==32)?~0ULL:((1ULL<<(2*K))-1);size_t cap=1<<20;char*line=malloc(cap);ssize_t L;while((L=getline(&line,&cap,stdin))>0){if(line[0]=='>')continue;uint64_t x=0;int run=0;for(ssize_t i=0;i<L;i++){int d=code(line[i]);if(d<0){run=0;x=0;continue;}x=((x<<2)|d)&mask;if(++run>=K){unsigned h=hh(x);while(idx[h]>=0){if(key[h]==x){cnt[idx[h]]++;break;}h=(h+1)&(H-1);}}}}
 for(int i=0;i<n;i++)if(cnt[i])printf("%llu %s\n",(unsigned long long)cnt[i],pat[i]);return 0;}
