// Exact, restricted structural integer proofs. No solver library dependencies.
// All arithmetic decisions use exactly represented integral doubles with bounded
// magnitude. Unsupported structures fall through without changing the model.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <vector>
#include "taral.hpp"

namespace structural_integer {
using Bits = std::vector<uint64_t>;
inline bool bit(const Bits& a, int k) { return (a[k / 64] >> (k % 64)) & 1; }
inline void toggle(Bits& a, int k) { a[k / 64] ^= uint64_t(1) << (k % 64); }
inline void add(Bits& a, const Bits& b) { for (size_t i = 0; i < a.size(); ++i) a[i] ^= b[i]; }
inline bool exact_int(double x) { return std::isfinite(x) && std::abs(x) <= 1099511627776.0 && x == std::trunc(x); }
struct Proof {
    enum Kind { None, Infeasible, Optimal, Timeout } kind = None;
    std::vector<double> x;
    std::string json;
};
inline std::string rows_json(const Bits& b, const std::vector<int>& eq) {
    std::ostringstream s; s << '['; bool comma = false;
    for (size_t i=0;i<eq.size();++i) if (bit(b,int(i))) { if(comma)s<<',';s<<eq[i];comma=true; }
    s << ']'; return s.str();
}
// At most 2048 equalities/variables. A parity contradiction is valid even when
// the continuous model is feasible. A unique binary assignment is NOT by itself
// a feasible point: remaining integer auxiliaries and every original row/bound
// must be checked exactly before declaring optimality.
inline Proof parity(const Model& m, std::chrono::steady_clock::time_point deadline) {
    Proof out; const int n=int(m.cols.size()), nr=int(m.row_lo.size());
    if(n==0 || n>2048 || nr>4096 || !m.qobj.empty() || int(m.is_int.size())!=n) return out;
    for(int j=0;j<n;++j) if(!m.is_int[j]) return out;
    std::vector<std::vector<Entry>> rows(nr);
    for(int j=0;j<n;++j) for(const auto& e:m.cols[j]) rows[e.index].push_back({j,e.value});
    std::vector<int> eq;
    for(int i=0;i<nr;++i) if(m.row_lo[i]==m.row_up[i] && exact_int(m.row_lo[i])) {
        bool safe=true;for(const auto& e:rows[i]) safe &= exact_int(e.value);
        if(safe)eq.push_back(i);
    }
    if(eq.empty() || eq.size()>2048) return out;
    const int ne=int(eq.size());
    std::vector<Bits> a(ne,Bits((n+1+63)/64)), why(ne,Bits((ne+63)/64));
    for(int i=0;i<ne;++i) {
        for(const auto& e:rows[eq[i]]) if(int64_t(e.value)%2)toggle(a[i],e.index);
        if(int64_t(m.row_lo[eq[i]])%2)toggle(a[i],n);
        toggle(why[i],i);
    }
    int rank=0;std::vector<int> piv;
    for(int j=0;j<n;++j) {
        if(std::chrono::steady_clock::now()>=deadline) {out.kind=Proof::Timeout;return out;}
        int k=rank;while(k<ne && !bit(a[k],j))++k;
        if(k==ne)continue;
        std::swap(a[k],a[rank]);std::swap(why[k],why[rank]);
        for(int i=0;i<ne;++i)if(i!=rank && bit(a[i],j)){add(a[i],a[rank]);add(why[i],why[rank]);}
        piv.push_back(j);++rank;
    }
    for(int i=rank;i<ne;++i) if(bit(a[i],n)) {
        // Independently recompute XOR from original equations. Never trust only
        // the eliminated matrix when producing an infeasibility status.
        Bits check((n+1+63)/64);
        for(int k=0;k<ne;++k)if(bit(why[i],k)) {
            for(const auto& e:rows[eq[k]]) if(int64_t(e.value)%2)toggle(check,e.index);
            if(int64_t(m.row_lo[eq[k]])%2)toggle(check,n);
        }
        bool ok=bit(check,n);for(int j=0;j<n;++j)ok &= !bit(check,j);
        if(!ok)return out;
        out.kind=Proof::Infeasible;
        out.json="{\"type\":\"integer_parity_contradiction\",\"rows\":"+rows_json(why[i],eq)+"}";
        return out;
    }
    // Solve only the fully determined binary + one-row auxiliary pattern.
    // General integers have a residue, not a value, so never fix them from GF(2).
    std::vector<char> binary(n,0),known(n,0);std::vector<double>x(n,0);
    for(int j=0;j<n;++j)binary[j]=m.col_lo[j]==0.0 && m.col_up[j]==1.0;
    std::ostringstream cert;cert<<"{\"type\":\"integer_parity_unique\",\"fixed_binary\":[";bool comma=false;
    for(int i=0;i<rank;++i) {
        int j=piv[i];if(!binary[j])continue;
        bool unique=true;for(int k=0;k<n;++k)if(k!=j && bit(a[i],k))unique=false;
        if(!unique)continue;
        x[j]=bit(a[i],n)?1.0:0.0;known[j]=1;
        if(comma) cert<<',';
        comma=true;
        cert<<"{\"col\":"<<j<<",\"value\":"<<int(x[j])<<",\"rows\":"<<rows_json(why[i],eq)<<'}';
    }
    cert<<"]}";
    for(int j=0;j<n;++j)if(binary[j] && !known[j])return out;
    for(int j=0;j<n;++j)if(!known[j]) {
        if(m.cols[j].size()!=1)return out;
        auto e=m.cols[j][0];int i=e.index;
        if(m.row_lo[i]!=m.row_up[i] || !exact_int(m.row_lo[i]) || !exact_int(e.value) || e.value==0)return out;
        int64_t total=0;
        for(const auto& t:rows[i]) if(t.index!=j) {
            if(!known[t.index] || !binary[t.index] || !exact_int(t.value))return out;
            // Each value is binary here. 2048 * 2^40 stays safely in int64_t.
            total+=int64_t(t.value)*int64_t(x[t.index]);
        }
        int64_t rhs=int64_t(m.row_lo[i])-total,c=int64_t(e.value);
        if(rhs%c)return out; // valid necessary congruence but unsupported full solution
        int64_t v=rhs/c;
        if(std::abs(double(v))>1099511627776.0)return out;
        x[j]=double(v);known[j]=1;
    }
    // EXACT row/bound validation for this restricted small-integer input.
    for(int j=0;j<n;++j) {
        if(!exact_int(x[j]) || x[j]<m.col_lo[j] || x[j]>m.col_up[j])return out;
        // Restrict auxiliaries further so sum of products fits int64_t.
        if(std::abs(x[j])>1048576.0)return out;
    }
    for(int i=0;i<nr;++i) {
        int64_t sum=0;
        for(const auto& e:rows[i]) {
            if(!exact_int(e.value) || std::abs(e.value)>1048576.0)return out;
            sum+=int64_t(e.value)*int64_t(x[e.index]);
        }
        if(std::isfinite(m.row_lo[i]) && (!exact_int(m.row_lo[i]) || double(sum)<m.row_lo[i]))return out;
        if(std::isfinite(m.row_up[i]) && (!exact_int(m.row_up[i]) || double(sum)>m.row_up[i]))return out;
    }
    // All variables are unique: objective coefficients need not be integral.
    out.kind=Proof::Optimal;out.x=std::move(x);out.json=cert.str();return out;
}
} // namespace structural_integer
