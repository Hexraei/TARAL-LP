#include "../src/structural_integer.hpp"
#include <cassert>
int main(){
 Model m;m.cols={{{0,31013}},{{0,-41014}},{{0,-51015}}};m.row_lo={0};m.row_up={0};m.col_lo={1,0,0};m.col_up={kInf,kInf,kInf};m.cost={1,0,0};m.is_int={1,1,1};
 auto d=std::chrono::steady_clock::time_point::max();auto p=structural_integer::three_integer(m,d);assert(p.kind==structural_integer::Proof::Optimal && p.x==std::vector<double>({25508,1,15506}));
 for(int k=0;k<8;++k){auto a=m;
  if(k==0)a.cols[0][0].value=31013.000001;
  if(k==1)a.is_int[1]=0;
  if(k==2)a.cols[2][0].value=-41014;
  if(k==3)a.row_up[0]=1;
  if(k==4)a.col_up[1]=1;
  if(k==5)a.cost[1]=1;
  if(k==6)a.maximize=true;
  if(k==7)a.col_lo[0]=0;
  assert(structural_integer::three_integer(a,d).kind==structural_integer::Proof::None);
 }
 assert(structural_integer::three_integer(m,std::chrono::steady_clock::time_point::min()).kind==structural_integer::Proof::Timeout);
 // Exact known optimum tiny case.
 m.cols={{{0,3}},{{0,-5}},{{0,-7}}};auto q=structural_integer::three_integer(m,d);assert(q.kind==structural_integer::Proof::Optimal && q.x[0]==4);
 assert(structural_integer::three_integer(m,std::chrono::steady_clock::time_point::min()).kind==structural_integer::Proof::Timeout);
}
