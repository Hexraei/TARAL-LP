// Includes the real implementation to test canonicalization, not a mirror.
#include "../src/milp.cpp"
#include <cassert>
int main(){
 std::vector<Change> a{{2,0,8},{1,-4,9},{2,3,8},{1,-4,4},{2,3,5}};auto b=a;
 canonical_changes(b);assert(b.size()==2 && b[0].j==1 && b[1].j==2);
 assert(b[0].lo==-4 && b[0].up==4 && b[1].lo==3 && b[1].up==5);
 for(unsigned seed=0;seed<100;++seed){std::mt19937 r(seed);
  std::vector<Change> hist;std::vector<double>lo(20,0),up(20,1000);
  for(int k=0;k<2000;++k){int j=r()%20;if(r()%2)lo[j]=std::min(up[j],lo[j]+double(r()%3));else up[j]=std::max(lo[j],up[j]-double(r()%3));hist.push_back({j,lo[j],up[j]});}
  auto map=hist;canonical_changes(map);assert(map.size()<=20);
  std::vector<double>clo(20,0),cup(20,1000);for(auto c:map)clo[c.j]=c.lo,cup[c.j]=c.up;assert(clo==lo && cup==up);
  auto down=map,upnode=map;down.push_back({0,lo[0],lo[0]});upnode.push_back({0,up[0],up[0]});canonical_changes(down);canonical_changes(upnode);
  assert(map[0].lo==lo[0] && map[0].up==up[0]);assert(down[0].up==lo[0] && upnode[0].lo==up[0]);
 }
}
