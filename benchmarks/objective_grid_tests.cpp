#include "../src/structural_integer.hpp"
#include <cassert>
int main() {
 Model m;m.cols={{{0,1}},{{0,-1},{1,2}},{{1,-1}}};m.row_lo={-kInf,0};m.row_up={0,0};m.col_lo={0,0,0};m.col_up={kInf,kInf,1};m.cost={-1,0,0};m.is_int={0,0,1};
 auto g=structural_integer::objective_grid(m);assert(g.step==0.5);assert(structural_integer::lattice_bound(-1.8,g)==-1.5);
 assert(structural_integer::lattice_bound(-1.50000000001,g)==-1.5);
 for(int kind=0;kind<9;++kind){auto a=m;
  if(kind==0)a.cols[2][0].value=-1.00000001;
  if(kind==1)a.is_int[2]=0;
  if(kind==2)a.cost[2]=1;
  if(kind==3)a.cols[0].push_back({1,1});
  if(kind==4)a.col_lo[1]=-kInf;
  if(kind==5)a.col_up[1]=4;
  if(kind==6)a.row_up[1]=1;
  if(kind==7)a.col_up[0]=1.25;
  if(kind==8)a.cols[1][1].value=2.1;
  assert(structural_integer::objective_grid(a).step==0);
 }
 g.offset=0.25;assert(structural_integer::lattice_bound(-1.8,g)==-1.75);
}
