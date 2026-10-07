// Standalone experimental runner. Not linked into the production taral target.
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include "miqp.hpp"
int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: r2-miqp model.mps [--time-limit seconds] [--json path]\n"; return 2; }
    double limit=60; const char* json=nullptr;
    for (int i=2;i<argc;++i) {
        std::string flag=argv[i];
        if (i+1>=argc) return 2;
        if (flag=="--time-limit") limit=std::stod(argv[++i]);
        else if (flag=="--json") json=argv[++i]; else return 2;
    }
    if (!std::isfinite(limit) || limit <= 0) return 2;
    try {
        Model md=read_mps(argv[1]);
        auto r=miqp_proto::solve_miqp(md,limit,100000);
        std::ostream* out=&std::cout; std::ofstream file;
        if(json) { file.open(json); if(!file) return 2; out=&file; }
        auto number=[&](double v) { if(std::isfinite(v)) *out << std::setprecision(17) << v; else *out << "null"; };
        *out << "{\"status\":\"" << r.status << "\",\"has_solution\":" << (r.has_solution?"true":"false") << ",\"objective\":";
        number(r.objective); *out << ",\"best_bound\":"; number(r.best_bound);
        *out << ",\"nodes\":" << r.nodes << ",\"unresolved\":" << r.unresolved << "}\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
