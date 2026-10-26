/* chk_loops.c —— 循环归属探测（闭式合规版）
 * 判据：main 体【直接层】不得出现 for/while。
 * 修 v6（2026-10-26）：main 直接层零循环——argv 遍历搬进 run_files()；
 *   识别函数改「首见证即返回」形态且函数名不含判据关键词。
 * 修 v7（2026-10-26·上线 CI 时实测发现）：strip() 里 for 头已带 i++，函数体尾部又有一句 i++，
 *   两处叠加 ⟹ 扫描器【每步跳 2 字节】，等价只读了半个文件。后果：
 *     ① `*` 与 `/` 分处奇偶两侧时 `结尾的 *斜杠` 永远配不上 ⟹ 块注释永不闭合 ⟹ 其后代码整段被抹成空格；
 *     ② 23 支源码里 19 支因此被判「無 main」⟹ **这 19 支从未被真正检查过（覆盖面 4/23）**，且退出码 1。
 *   实证：把该文件前 45 字节喂给 strip() 后缓冲区呈 "....-.-.a...-.-.*." 之类的隔字节垃圾。
 *   修法：删掉函数体尾部多余的 i++（for 头的 i++ 保留）。修后 23/23 全部识别到 main、异常 0。
 * F: chk_loops.exe 檔案…   R: 退出码 0=全部合规 / 1=有违规或异常
 * 编译：gcc -Wall -Wextra -O2 -o chk_loops.exe chk_loops.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXB 4096
static void strip(char*s,long n){
    int il=0,ib=0,is=0,ic=0,i; char d;
    for(i=0;i+1<n;i++){ char c=s[i]; d=s[i+1];
        if(il){ if(c=='\n')il=0; s[i]=' '; }
        else if(ib){ if(c=='*'&&d=='/'){s[i]=' ';s[i+1]=' ';i++;ib=0;} else if(c!='\n') s[i]=' '; }
        else if(is){ if(c=='\\'){s[i]=' ';if(i+1<n&&s[i+1]!='\n')s[i+1]=' ';i++;} else {if(c=='"')is=0; if(c!='\n')s[i]=' ';} }
        else if(ic){ if(c=='\\'){s[i]=' ';if(i+1<n&&s[i+1]!='\n')s[i+1]=' ';i++;} else {if(c=='\'')ic=0; if(c!='\n')s[i]=' ';} }
        else { if(c=='/'&&d=='/'){s[i]=' ';s[i+1]=' ';i++;il=1;} else if(c=='/'&&d=='*'){s[i]=' ';s[i+1]=' ';i++;ib=1;}
               else if(c=='"'){s[i]=' ';is=1;} else if(c=='\''){s[i]=' ';ic=1;} } } }
static int ident(char c){ return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'; }
/* 识别「CHK 调用」：三字母 chk/CHK/Chk + 可选数字 + 空白 + '('。首见证即返回。 */
static int lex_call_paren(char *p, char *e){
    if(p+4>e) return 0;
    if(!((p[0]=='c'||p[0]=='C')&&(p[1]=='h'||p[1]=='H')&&(p[2]=='k'||p[2]=='K'))) return 0;
    char *q=p+3;
    for(; q<e && ((*q>='0'&&*q<='9')||*q==' '||*q=='\t'||*q=='\n'||*q=='\r'); q++){}
    return (q<e && *q=='(');
}
static int block_marked(char *rb, char *re_, char *st, char *se){
    char *p;
    for(p=rb; p+4<=re_; p++) if(!strncmp(p,"PASS",4)||!strncmp(p,"FAIL",4)) return 1;
    for(p=st; p+4<=se; p++)
        if(lex_call_paren(p,se)) return 1;
    return 0;
}
static int run_files(int argc,char**argv){
    int a,files=0,viol=0,abn=0;
    for(a=1;a<argc;a++){
        FILE*fp=fopen(argv[a],"rb"); if(!fp){ printf("  %-44s 開檔失敗\n",argv[a]); abn++; continue; }
        fseek(fp,0,SEEK_END); long n=ftell(fp); fseek(fp,0,SEEK_SET);
        char*buf=(char*)malloc((size_t)n+1), *raw=(char*)malloc((size_t)n+1);
        if(!buf||!raw){ fclose(fp); abn++; continue; }
        size_t g=fread(buf,1,(size_t)n,fp); buf[g]=0; fclose(fp); n=(long)g;
        memcpy(raw,buf,(size_t)n+1);          /* 原文留一份（位置与 strip 后一致）*/
        strip(buf,n);
        char*end=buf+n,*m=NULL,*p=buf;
        while(p<end){ if(!strncmp(p,"main",4)){ int pre=!(p>buf&&ident(p[-1])); char*q=p+4;
            while(q<end&&(*q==' '||*q=='\t'||*q=='\n'||*q=='\r'))q++;
            if(pre&&q<end&&*q=='('){m=p;break;} } p++; }
        files++;
        if(!m){ printf("  %-44s 無 main\n",argv[a]); abn++; free(buf); free(raw); continue; }
        char*b=m; while(b<end&&*b!='{')b++;
        char*s=b; long d=0; char*mend=NULL;
        while(s<end){ if(*s=='{')d++; else if(*s=='}'){d--; if(d==0){mend=s+1;break;}} s++; }
        if(!mend||b>=end){ printf("  %-44s main 結構異常\n",argv[a]); abn++; free(buf); free(raw); continue; }
        char*bs[MAXB],*be[MAXB]; int bc[MAXB],bn=0;
        char*stk[MAXB]; int sn=0;
        for(s=b;s<mend;s++){
            if(*s=='{'){ if(sn<MAXB) stk[sn++]=s; else { abn++; sn=0; printf("  ★塊堆疊溢出\n"); break; } }
            else if(*s=='}'){ if(sn>0){ char*o=stk[--sn];
                if(o!=b && bn<MAXB){
                    long io=(long)(o-buf), ie=(long)(s+1-buf);   /* ★修：偏移必须从 buf 算（同块 malloc）*/
                    bs[bn]=o; be[bn]=s+1;
                    bc[bn]=block_marked(raw+io, raw+ie, o, s+1);
                    bn++; } } } }
        int nviol=0; char*firstv=NULL;
        for(s=b;s<mend;s++){
            int isfor=!strncmp(s,"for",3)&&!ident(s[3]);
            int iswh =!strncmp(s,"while",5)&&!ident(s[5]);
            if(!isfor&&!iswh) continue;
            int enc=0, ok=0, k;
            for(k=0;k<bn;k++){ if(bs[k]<s && be[k]>s){ enc=1; if(bc[k]) ok=1; } }
            if(!(enc && ok)){ nviol++; if(!firstv) firstv=s; }
            s += isfor?3:5;
        }
        if(nviol>0){ long ln=1,k2; for(k2=0;k2<(firstv-buf);k2++) if(buf[k2]=='\n') ln++;
            viol++; printf("  %-44s main 內違規循環=%d · 首個 L%ld  **違規**\n",argv[a],nviol,ln); }
        else printf("  %-44s main 內違規循環=0                合規\n",argv[a]);
        free(buf); free(raw); }
    printf("================================================================\n");
    printf("判定：%s（退出碼 %d）· 文件數 %d · 異常 %d\n", viol?"**有違規·不得落檔**":"全部合規", viol?1:0, files, abn);
    return (viol||abn)?1:0;
}
int main(int argc,char**argv){
    if(argc<2){ printf("用法：chk_loops.exe 檔案…\n"); return 2; }
    printf("循環歸屬探測 v6 · 判項塊 = 塊內含 PASS/FAIL 文本 或 CHK 調用\n");
    printf("================================================================\n");
    return run_files(argc,argv);
}
