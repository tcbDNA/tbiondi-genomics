/* tyscan tRNA.fa [minscore] < single-line-fasta : POST-HOC (cold review 2, check 4). Finds places in reads that resemble the tRNA body (no flanks, relaxed identity).
   Seeds: 10-mers of the tRNA (both strands), two seeds on a consistent diagonal (+-6) within 90 nt (an 11-mer, or two separate 10-mers). Each seeded window (+-110 nt) is aligned to the tRNA by local
   Smith-Waterman (+2 match, -2 mismatch, -3 gap); windows with score >= minscore (default 46) are printed in tRNA-sense orientation: >read:pos:strand:score */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define K 10
#define N (1<<(2*K))
static int16_t pos[2][N]; static char T[2][512]; static int L;
static int code(int c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;}return -1;}
static char comp(char c){switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';}return 'N';}
static int sw(const char*a,int n,const char*b,int m){static int p[4096],q[4096];int best=0;memset(p,0,sizeof(int)*(m+1));for(int i=1;i<=n;i++){q[0]=0;for(int j=1;j<=m;j++){int s=p[j-1]+(a[i-1]==b[j-1]&&a[i-1]!='N'?2:-2);int d=p[j]-3,e=q[j-1]-3;if(d>s)s=d;if(e>s)s=e;if(s<0)s=0;q[j]=s;if(s>best)best=s;}memcpy(p,q,sizeof(int)*(m+1));}return best;}
int main(int c,char**v){FILE*f=fopen(v[1],"r");int min=c>2?atoi(v[2]):46;char b[4096];while(fgets(b,4096,f)){b[strcspn(b,"\r\n")]=0;if(*b!='>')strcat(T[0],b);}L=strlen(T[0]);for(int i=0;i<L;i++)T[1][i]=comp(T[0][L-1-i]);
 memset(pos,-1,sizeof pos);for(int s=0;s<2;s++)for(int i=0;i+K<=L;i++){uint32_t x=0;int ok=1;for(int j=0;j<K;j++){int d=code(T[s][i+j]);if(d<0)ok=0;x=(x<<2)|(d&3);}if(ok&&pos[s][x]<0)pos[s][x]=i;}
 size_t cap=1<<20;char*line=malloc(cap),hd[1<<12];hd[0]=0;ssize_t len;uint32_t mask=N-1;long nread=0,nw=0;
 while((len=getline(&line,&cap,stdin))>0){if(line[0]=='>'){strncpy(hd,line+1,4000);hd[strcspn(hd," \r\n")]=0;continue;}
  while(len>0&&(line[len-1]=='\n'||line[len-1]=='\r'))line[--len]=0; nread++;
  for(ssize_t i=0;i<len;i++)if(line[i]>='a')line[i]-=32;
  uint32_t x=0;int run=0;long lp[2]={-1000,-1000},ld[2]={0,0},emit[2]={-1000,-1000};
  for(ssize_t i=0;i<len;i++){int d=code(line[i]);if(d<0){run=0;continue;}x=((x<<2)|d)&mask;if(++run<K)continue;
   for(int s=0;s<2;s++){int t=pos[s][x];if(t<0)continue;long st=i-K+1,dg=st-t;
    if(st-lp[s]<=90&&labs(dg-ld[s])<=6&&st!=lp[s]&&dg-emit[s]>40){long a=dg-110,e=dg+L+110;if(a<0)a=0;if(e>len)e=len;int sc=sw(line+a,e-a,T[s],L);if(sc>=min){emit[s]=dg;nw++;printf(">%s:%ld:%c:%d\n",hd,a+1,s?'-':'+',sc);if(!s)printf("%.*s\n",(int)(e-a),line+a);else{for(long k=e-1;k>=a;k--)putchar(comp(line[k]));putchar('\n');}}else if(sc>=min-10){emit[s]=dg-39;}}
    lp[s]=st;ld[s]=dg;}}}
 fprintf(stderr,"reads %ld windows %ld\n",nread,nw);return 0;}
