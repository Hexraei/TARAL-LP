// TARAL-LP C++ engine, M0: from-scratch MPS reader, sparse LU (left-looking, threshold pivoting) + eta file,
// bounded-variable primal revised simplex (composite phase 1). C++ standard library only.
// CLI: taral MODEL.mps --time-limit S --sol OUT.sol --json OUT.json
#include <bits/stdc++.h>
using namespace std;
typedef double D;
static const D INF=numeric_limits<D>::infinity();
static D now(){return chrono::duration<D>(chrono::steady_clock::now().time_since_epoch()).count();}

struct Model{
  int m=0,n=0; vector<string> cname; vector<int> cs,ci; vector<D> cv; vector<D> c,lo,up; // lo/up size n+m (structurals then row activities)
  D off=0; bool maxi=false;
};
static bool readMPS(const string&path,Model&M,string&err){
  ifstream f(path); if(!f){err="cannot open";return false;}
  bool fixedMode=false; { ifstream g(path); string l,sc; while(getline(g,l)){ if(l.empty()||l[0]=='*')continue; if(l[0]!=' '&&l[0]!='\t'){istringstream is(l);is>>sc;continue;} istringstream is(l);vector<string> q;string w;while(is>>w)q.push_back(w);
      if(sc=="ROWS"&&q.size()!=2)fixedMode=true; if(sc=="COLUMNS"&&q.size()!=3&&q.size()!=5)fixedMode=true; if((sc=="RHS"||sc=="RANGES")&&q.size()!=3&&q.size()!=5)fixedMode=true; if(fixedMode)break; } }
  string line,sec; unordered_map<string,int> rid,cid; vector<char> rtype; vector<D> rhs,rng; string objrow; bool haveObj=false;
  vector<vector<pair<int,D>>> cols; vector<D> cost; vector<D> blo,bup; vector<string> cn;
  D objrhs=0; bool objsense=false;
  while(getline(f,line)){
    if(line.empty()||line[0]=='*')continue;
    if(line[0]!=' '&&line[0]!='\t'){ istringstream is(line);string w;is>>w; sec=w; if(w=="OBJSENSE"){ string r; if(is>>r){ if(r=="MAX"||r=="MAXIMIZE")M.maxi=true; sec="";} } continue; }
    vector<string> t; if(fixedMode&&(sec=="ROWS"||sec=="COLUMNS"||sec=="RHS"||sec=="RANGES"||sec=="BOUNDS")){
      auto P=[&](size_t a,size_t b){ string z=line.size()>a?line.substr(a,b-a):string(); size_t i=z.find_first_not_of(" \t\r"),j=z.find_last_not_of(" \t\r"); return i==string::npos?string():z.substr(i,j-i+1); };
      if(sec=="ROWS"){t={P(1,3),P(4,12)};}
      else if(sec=="BOUNDS"){t={P(1,3),P(4,12),P(14,22)}; string v=P(24,36); if(!v.empty())t.push_back(v);}
      else {t={P(4,12),P(14,22),P(24,36)}; string a=P(39,47); if(!a.empty()){t.push_back(a);t.push_back(P(49,61));}}
    } else { istringstream is(line);string w;while(is>>w)t.push_back(w); } if(t.empty())continue;
    if(sec=="OBJSENSE"){ if(t[0]=="MAX"||t[0]=="MAXIMIZE")M.maxi=true; continue;}
    if(sec=="ROWS"){ char ty=t[0][0]; if(ty=='N'){ if(!haveObj){haveObj=true;objrow=t[1];} rid[t[1]]=-1; } else { rid[t[1]]=rtype.size(); rtype.push_back(ty);} }
    else if(sec=="COLUMNS"){
      if(t.size()>=3&&t[1]=="'MARKER'")continue;
      auto it=cid.find(t[0]); int j; if(it==cid.end()){j=cn.size();cid[t[0]]=j;cn.push_back(t[0]);cols.emplace_back();cost.push_back(0);} else j=it->second;
      for(size_t k=1;k+1<t.size();k+=2){ auto r=rid.find(t[k]); if(r==rid.end())continue; D v=atof(t[k+1].c_str()); if(r->second==-1){ if(t[k]==objrow)cost[j]+=v; } else if(v!=0)cols[j].push_back({r->second,v}); }
    } else if(sec=="RHS"||sec=="RANGES"){
      if(rhs.size()<rtype.size()){rhs.assign(rtype.size(),0);rng.assign(rtype.size(),0);} 
      size_t k=(t.size()%2==1)?1:0;
      for(;k+1<t.size();k+=2){ auto r=rid.find(t[k]); if(r==rid.end())continue; D v=atof(t[k+1].c_str());
        if(r->second==-1){ if(sec=="RHS"&&t[k]==objrow)objrhs=v; } else if(sec=="RHS")rhs[r->second]=v; else rng[r->second]=v; }
    } else if(sec=="BOUNDS"){
      if(blo.size()<cn.size()){blo.resize(cn.size(),0);bup.resize(cn.size(),INF);} string ty=t[0]; string cname_; D v=0;
      bool noval=(ty=="FR"||ty=="MI"||ty=="PL"||ty=="BV");
      if(noval){cname_=t.size()>=3?t[2]:t[1];} else { if(t.size()>=4){cname_=t[2];v=atof(t[3].c_str());} else {cname_=t[1];v=atof(t[2].c_str());} }
      auto it=cid.find(cname_); if(it==cid.end())continue; int j=it->second;
      if(ty=="UP"){bup[j]=v; if(v<0&&blo[j]==0)blo[j]=-INF;} else if(ty=="LO")blo[j]=v; else if(ty=="FX"){blo[j]=bup[j]=v;} else if(ty=="FR"){blo[j]=-INF;bup[j]=INF;}
      else if(ty=="MI")blo[j]=-INF; else if(ty=="PL")bup[j]=INF; else if(ty=="BV"){blo[j]=0;bup[j]=1;} else if(ty=="LI")blo[j]=v; else if(ty=="UI")bup[j]=v;
    }
  }
  int m=rtype.size(),n=cn.size(); if(rhs.size()<(size_t)m){rhs.assign(m,0);rng.assign(m,0);} blo.resize(n,0);bup.resize(n,INF);
  M.m=m;M.n=n;M.cname=cn; M.cs.assign(n+1,0); for(int j=0;j<n;j++){ sort(cols[j].begin(),cols[j].end()); for(auto&p:cols[j]){M.ci.push_back(p.first);M.cv.push_back(p.second);} M.cs[j+1]=M.ci.size(); }
  M.c=cost; if(M.maxi)for(auto&x:M.c)x=-x; M.lo.assign(n+m,0);M.up.assign(n+m,INF);
  for(int j=0;j<n;j++){M.lo[j]=blo[j];M.up[j]=bup[j];}
  for(int i=0;i<m;i++){ D l,u,r=rhs[i],R=rng[i]; char ty=rtype[i];
    if(ty=='E'){ if(R>0){l=r;u=r+R;} else if(R<0){l=r+R;u=r;} else l=u=r; } else if(ty=='L'){l=-INF;u=r; if(R!=0)l=r-fabs(R);} else {l=r;u=INF; if(R!=0)u=r+fabs(R);} M.lo[n+i]=l;M.up[n+i]=u; }
  M.off=-objrhs; return true;
}

struct LU{
  int m=0; vector<int> prow,colpos; vector<int> Ls,Li; vector<D> Lv; vector<int> Us,Ui; vector<D> Uv; vector<D> piv;
  // eta file (PFI)
  vector<int> er; vector<int> es,ei; vector<D> ev; vector<D> epiv;
  vector<int> failed;
  template<class F> void factor(int m_,const vector<int>&nnzc,F getcol){ // getcol(pos, idx, val)
    m=m_; prow.assign(m,-1);colpos.assign(m,-1);Ls.assign(1,0);Li.clear();Lv.clear();Us.assign(1,0);Ui.clear();Uv.clear();piv.assign(m,0);failed.clear();
    er.clear();es.assign(1,0);ei.clear();ev.clear();epiv.clear();
    vector<int> order(m);iota(order.begin(),order.end(),0); stable_sort(order.begin(),order.end(),[&](int a,int b){return nnzc[a]<nnzc[b];});
    vector<int> rowcnt(m,0); vector<int> pi; vector<D> pv; for(int p=0;p<m;p++){pi.clear();pv.clear();getcol(p,pi,pv);for(int r:pi)rowcnt[r]++;}
    vector<char> used(m,0); vector<D> w(m,0); int k=0; vector<int> pidx(m,-1); // row->pivot idx
    vector<int> nz; vector<char> mark(m,0);
    for(int oi=0;oi<m;oi++){ int p=order[oi]; pi.clear();pv.clear();getcol(p,pi,pv); for(size_t t=0;t<pi.size();t++)w[pi[t]]=pv[t];
      for(int j=0;j<k;j++){ D x=w[prow[j]]; if(x!=0){ for(int e=Ls[j];e<Ls[j+1];e++)w[Li[e]]-=Lv[e]*x; } }
      D mx=0; for(int i=0;i<m;i++) if(!used[i]) mx=max(mx,fabs(w[i]));
      if(mx<1e-9){ failed.push_back(p); for(int i=0;i<m;i++)w[i]=0; continue; }
      int best=-1;int bc=INT_MAX; for(int i=0;i<m;i++) if(!used[i]&&fabs(w[i])>=0.1*mx){ if(rowcnt[i]<bc||(rowcnt[i]==bc&&fabs(w[i])>fabs(w[best]))){bc=rowcnt[i];best=i;} }
      D pv0=w[best]; prow[k]=best;colpos[k]=p;piv[k]=pv0;used[best]=1;pidx[best]=k;
      for(int i=0;i<m;i++){ D x=w[i]; if(x==0)continue; if(used[i]&&i!=best){ Ui.push_back(pidx[i]);Uv.push_back(x);} else if(!used[i]){ Li.push_back(i);Lv.push_back(x/pv0);} w[i]=0; }
      w[best]=0; Ls.push_back(Li.size());Us.push_back(Ui.size()); k++; }
    // note: Ls/Us have k+1 entries
  }
  void ftranLU(vector<D>&a,vector<D>&x){ // a indexed by row (destroyed), x by basis position
    int m_=m; for(int j=0;j<m_;j++){ if(prow[j]<0)continue; D v=a[prow[j]]; if(v!=0)for(int e=Ls[j];e<Ls[j+1];e++)a[Li[e]]-=Lv[e]*v; }
    for(int k=m_-1;k>=0;k--){ if(prow[k]<0)continue; D z=a[prow[k]]/piv[k]; x[colpos[k]]=z; if(z!=0)for(int e=Us[k];e<Us[k+1];e++)a[prow[Ui[e]]]-=Uv[e]*z; }
  }
  void ftran(vector<D>&a,vector<D>&x){ ftranLU(a,x); for(size_t t=0;t<er.size();t++){ int r=er[t]; D xr=x[r]/epiv[t]; if(xr!=0){ for(int e=es[t];e<es[t+1];e++)x[ei[e]]-=ev[e]*xr; } x[r]=xr; } }
  void btran(vector<D>&c,vector<D>&y){ // c indexed by basis position (destroyed); y by row
    for(int t=(int)er.size()-1;t>=0;t--){ int r=er[t]; D s=c[r]; for(int e=es[t];e<es[t+1];e++)s-=ev[e]*c[ei[e]]; c[r]=s/epiv[t]; }
    vector<D> v(m); for(int k=0;k<m;k++){ D s=c[colpos[k]]; for(int e=Us[k];e<Us[k+1];e++)s-=Uv[e]*v[Ui[e]]; v[k]=s/piv[k]; }
    for(int k=0;k<m;k++)y[prow[k]]=v[k];
    for(int j=m-1;j>=0;j--){ D s=0; for(int e=Ls[j];e<Ls[j+1];e++)s+=Lv[e]*y[Li[e]]; y[prow[j]]-=s; }
  }
  void addEta(int r,const vector<D>&alpha){ er.push_back(r);epiv.push_back(alpha[r]); for(int i=0;i<m;i++) if(i!=r&&alpha[i]!=0){ei.push_back(i);ev.push_back(alpha[i]);} es.push_back(ei.size()); }
};

struct Res{string status="failed";D obj=0;long it=0;string msg;vector<D> x;};
static Res solve(const Model&M,D tlim,D t0){
  Res R; int m=M.m,n=M.n,N=n+m; vector<D> lo=M.lo,up=M.up; bool perturbed=false,pertDone=false,wantPert=false; int pertTh=getenv("TARAL_PERT")?atoi(getenv("TARAL_PERT")):3000; mt19937_64 rng(12345); const D ftol=1e-8,dtol=1e-9;
  vector<int> head(m),pos(N,-1); vector<D> xv(N,0),xB(m); vector<char> atUp(N,0);
  for(int i=0;i<m;i++){head[i]=n+i;pos[n+i]=i;}
  for(int j=0;j<n;j++){ if(lo[j]>-INF){xv[j]=lo[j];} else if(up[j]<INF){xv[j]=up[j];atUp[j]=1;} else xv[j]=0; }
  LU lu; auto colget=[&](int var,vector<int>&idx,vector<D>&val){ if(var<n){ for(int e=M.cs[var];e<M.cs[var+1];e++){idx.push_back(M.ci[e]);val.push_back(M.cv[e]);} } else {idx.push_back(var-n);val.push_back(-1);} };
  auto refactor=[&]()->bool{
    for(int tries=0;tries<5;tries++){
      vector<int> nnzc(m); for(int p=0;p<m;p++){int v=head[p];nnzc[p]=v<n?M.cs[v+1]-M.cs[v]:1;}
      lu.factor(m,nnzc,[&](int p,vector<int>&i,vector<D>&v){colget(head[p],i,v);});
      if(lu.failed.empty())return true;
      // replace failed basic columns with slacks of unpivoted rows
      vector<char> used(m,0); for(int k=0;k<m;k++)if(lu.prow[k]>=0)used[lu.prow[k]]=1;
      vector<int> freeRows; for(int i=0;i<m;i++)if(!used[i])freeRows.push_back(i);
      for(size_t t=0;t<lu.failed.size()&&t<freeRows.size();t++){ int p=lu.failed[t]; int v=head[p]; pos[v]=-1; xv[v]=(lo[v]>-INF)?lo[v]:(up[v]<INF?up[v]:0); atUp[v]=!(lo[v]>-INF)&&up[v]<INF; int s=n+freeRows[t]; head[p]=s;pos[s]=p; }
    } return false; };
  auto computeXB=[&](){ vector<D> rhs(m,0); for(int j=0;j<N;j++) if(pos[j]<0&&xv[j]!=0){ if(j<n){ for(int e=M.cs[j];e<M.cs[j+1];e++)rhs[M.ci[e]]-=M.cv[e]*xv[j]; } else rhs[j-n]+=xv[j]; } lu.ftran(rhs,xB); };
  if(!refactor()){R.msg="singular start";return R;} computeXB();
  auto perturb=[&](){ for(int j=0;j<N;j++){ D r1=(rng()%1000000)/1e6; if(lo[j]>-INF&&lo[j]!=up[j]){D d=(1e-7+1e-6*r1)*(1+fabs(lo[j]));lo[j]-=d;} if(up[j]<INF&&lo[j]!=up[j]){D d=(1e-7+1e-6*r1)*(1+fabs(up[j]));up[j]+=d;} }
    for(int j=0;j<N;j++) if(pos[j]<0){ if(atUp[j]&&up[j]<INF)xv[j]=up[j]; else if(lo[j]>-INF)xv[j]=lo[j]; else if(up[j]<INF)xv[j]=up[j]; } perturbed=true;pertDone=true; computeXB(); };
  auto unperturb=[&](){ lo=M.lo;up=M.up; for(int j=0;j<N;j++) if(pos[j]<0){ if(atUp[j]&&up[j]<INF)xv[j]=up[j]; else if(lo[j]>-INF)xv[j]=lo[j]; else if(up[j]<INF)xv[j]=up[j]; else xv[j]=0; } perturbed=false; computeXB(); };
  vector<D> cB(m),y(m),alpha(m),col(m); long it=0; int degen=0; bool bland=false; int sinceRef=0; bool verified=false;
  auto cost=[&](int v)->D{return v<n?M.c[v]:0;};
  vector<D> dwork(N);
  while(true){
    if((it&15)==0&&now()-t0>tlim){R.status="timeout";R.it=it;return R;}
    if(wantPert){ wantPert=false; perturb(); }
    bool ph1=false; for(int i=0;i<m;i++){int v=head[i];D x=xB[i]; if(x<lo[v]-ftol){cB[i]=-1;ph1=true;} else if(x>up[v]+ftol){cB[i]=1;ph1=true;} else cB[i]=0;}
    if(!ph1)for(int i=0;i<m;i++)cB[i]=cost(head[i]);
    { vector<D> c2=cB; lu.btran(c2,y); }
    int q=-1;D bd=0;int dir=0;
    for(int j=0;j<N;j++){ if(pos[j]>=0)continue; if(lo[j]==up[j])continue; D cj=ph1?0:cost(j); D dj; if(j<n){D s=0;for(int e=M.cs[j];e<M.cs[j+1];e++)s+=y[M.ci[e]]*M.cv[e];dj=cj-s;} else dj=cj+y[j-n];
      int dd=0; bool canUp=xv[j]<up[j]-0&&!(atUp[j]&&up[j]<INF&&lo[j]>-INF&&false); bool free_=lo[j]==-INF&&up[j]==INF;
      if(free_){ if(dj<-dtol)dd=1; else if(dj>dtol)dd=-1; }
      else if(!atUp[j]&&dj<-dtol&&xv[j]<up[j])dd=1; else if(atUp[j]&&dj>dtol&&xv[j]>lo[j])dd=-1; (void)canUp;
      if(!dd)continue; D sc=fabs(dj); if(bland){q=j;dir=dd;bd=sc;break;} if(sc>bd){bd=sc;q=j;dir=dd;} }
    if(q<0){
      if(perturbed){ unperturb(); sinceRef=0; if(!refactor()){R.msg="singular";return R;} computeXB(); bland=false; degen=0; continue; }
      if(sinceRef>0){ if(!refactor()){R.msg="singular";return R;} computeXB();sinceRef=0;continue; }
      R.it=it; if(ph1){R.status="infeasible";return R;} R.status="optimal"; break; }
    fill(col.begin(),col.end(),0); { vector<int> ii;vector<D> vv;colget(q,ii,vv);for(size_t t=0;t<ii.size();t++)col[ii[t]]=vv[t]; }
    lu.ftran(col,alpha);
    // ratio test (Harris two-pass)
    D tmax=INF; const D hz=1e-9;
    for(int i=0;i<m;i++){ D a=alpha[i]; if(fabs(a)<1e-9)continue; D rate=-dir*a; int v=head[i];D x=xB[i];D lim;
      if(rate>0){ if(x<lo[v]-ftol)lim=lo[v]; else if(x>up[v]+ftol)continue; else lim=up[v]; } else { if(x>up[v]+ftol)lim=up[v]; else if(x<lo[v]-ftol)continue; else lim=lo[v]; }
      if(lim==INF||lim==-INF)continue; D r=(lim+(rate>0?hz:-hz)-x)/rate; if(r<tmax)tmax=r; }
    int r=-1;D t=INF;D besta=0;
    if(tmax<INF){ for(int i=0;i<m;i++){ D a=alpha[i]; if(fabs(a)<1e-9)continue; D rate=-dir*a; int v=head[i];D x=xB[i];D lim;
      if(rate>0){ if(x<lo[v]-ftol)lim=lo[v]; else if(x>up[v]+ftol)continue; else lim=up[v]; } else { if(x>up[v]+ftol)lim=up[v]; else if(x<lo[v]-ftol)continue; else lim=lo[v]; }
      if(lim==INF||lim==-INF)continue; D rr=(lim-x)/rate; if(rr<=tmax&&fabs(a)>besta){besta=fabs(a);r=i;t=max(rr,0.0);} } }
    D flip=(lo[q]>-INF&&up[q]<INF)?up[q]-lo[q]:INF;
    if(r<0&&flip==INF){ if(ph1){R.msg="phase1 unbounded ray";R.status="failed";R.it=it;return R;} R.status="unbounded";R.it=it;return R; }
    it++;
    if(flip<=t){ // bound flip
      for(int i=0;i<m;i++)xB[i]-=dir*alpha[i]*flip; xv[q]=dir>0?up[q]:lo[q]; atUp[q]=dir>0; degen=0; continue; }
    if(t<1e-12){ ++degen; if(!pertDone&&(degen>pertTh||(degen>150&&now()-t0>0.4*tlim))){ wantPert=true; degen=0; } else if(degen>5000)bland=true; } else {degen=0;bland=false;}
    int lv=head[r]; D rate=-dir*alpha[r]; for(int i=0;i<m;i++)xB[i]-=dir*alpha[i]*t; 
    D xq=xv[q]+dir*t; // leaving var goes to the bound it hit
    bool toUp= rate>0 ? !(xB[r]<lo[lv]-ftol&&false) : false; D lx;
    { D x0=xB[r]; // after update it should sit at a bound
      if(rate>0){ lx=(x0>up[lv]-1e-6*(1+fabs(up[lv]))&&up[lv]<INF)?up[lv]:lo[lv]; if(lo[lv]>-INF&&fabs(x0-lo[lv])<fabs(x0-lx))lx=lo[lv]; }
      else { lx=(lo[lv]>-INF)?lo[lv]:up[lv]; if(up[lv]<INF&&fabs(x0-up[lv])<fabs(x0-lx))lx=up[lv]; }
      toUp=(lx==up[lv]&&up[lv]<INF&&!(lx==lo[lv])); }
    pos[lv]=-1;xv[lv]=lx;atUp[lv]=toUp; head[r]=q;pos[q]=r;xB[r]=xq;atUp[q]=0;
    lu.addEta(r,alpha); sinceRef++;
    if(sinceRef>=100){ if(!refactor()){R.msg="singular";R.it=it;return R;} computeXB();sinceRef=0; }
  }
  // final accurate recompute
  if(!refactor()){R.msg="singular";return R;} computeXB();
  R.x.assign(n,0); for(int j=0;j<n;j++) R.x[j]=pos[j]>=0?xB[pos[j]]:xv[j];
  D o=M.off; for(int j=0;j<n;j++)o+=M.c[j]*R.x[j]; R.obj=M.maxi?-o:o; R.it=it; return R;
}
int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage\n");return 2;} string mps=argv[1],sol,js; D tl=60; for(int i=2;i<argc;i++){string a=argv[i]; if(a=="--time-limit"&&i+1<argc)tl=atof(argv[++i]); else if(a=="--sol"&&i+1<argc)sol=argv[++i]; else if(a=="--json"&&i+1<argc)js=argv[++i];}
  D t0=now(); Model M; string err; Res R; if(!readMPS(mps,M,err)){R.msg=err;} else R=solve(M,tl,t0);
  if(!sol.empty()&&R.status=="optimal"){FILE*f=fopen(sol.c_str(),"w");for(int j=0;j<M.n;j++)fprintf(f,"%s %.17g\n",M.cname[j].c_str(),R.x[j]);fclose(f);}
  if(!js.empty()){FILE*f=fopen(js.c_str(),"w");fprintf(f,"{\"status\":\"%s\",\"objective\":%.17g,\"iterations\":%ld,\"wall_s\":%.4f,\"message\":\"%s\"}\n",R.status.c_str(),R.obj,R.it,now()-t0,R.msg.c_str());fclose(f);}
  return 0; }
