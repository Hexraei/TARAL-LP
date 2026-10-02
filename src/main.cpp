// CLI per the engine contract:  taral MODEL.mps --time-limit S --sol OUT.sol --json OUT.json
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>
#include <string>

#include "ipm.hpp"
#include "taral.hpp"

namespace {
std::string json_escape(const std::string& s) {
    std::string o;
    for (char ch : s) {
        if (ch == '"' || ch == '\\') o += '\\', o += ch;
        else if (static_cast<unsigned char>(ch) < 0x20) o += ' ';
        else o += ch;
    }
    return o;
}

void write_json(const char* path, const char* status, const Result* r, double wall, const std::string& msg) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    std::fprintf(f, "{\"status\": \"%s\", \"objective\": ", status);
    if (r && r->status == Status::Optimal) std::fprintf(f, "%.17g", r->objective);
    else std::fprintf(f, "null");
    std::fprintf(f, ", \"iterations\": %ld, \"wall_s\": %.6f, \"message\": \"%s\"}\n", r ? r->iterations : 0L, wall,
                 json_escape(msg).c_str());
    std::fclose(f);
}

void num_or_null(std::FILE* f, double v) {
    if (std::isfinite(v)) std::fprintf(f, "%.17g", v);
    else std::fprintf(f, "null");
}

void write_milp_json(const char* path, const MilpResult& r, double wall) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    std::fprintf(f, "{\"status\": \"%s\", \"objective\": ", r.status.c_str());
    num_or_null(f, r.has_solution && r.status != "unbounded" ? r.objective : kInf);
    std::fprintf(f, ", \"best_bound\": ");
    num_or_null(f, r.best_bound);
    std::fprintf(f, ", \"gap\": ");
    num_or_null(f, r.has_solution && r.status != "unbounded" ? r.gap : kInf);
    std::fprintf(f, ", \"nodes\": %ld, \"unresolved_nodes\": %ld, \"has_solution\": %s", r.nodes, r.unresolved_nodes,
                 r.has_solution ? "true" : "false");
    std::fprintf(f, ", \"iterations\": %ld, \"wall_s\": %.6f, \"message\": \"%s\"}\n", r.lp_iterations, wall,
                 json_escape(r.message).c_str());
    std::fclose(f);
}

void write_ipm_json(const char* path, const IpmResult& r, double wall) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    bool ok = r.status == IpmStatus::Optimal;
    std::fprintf(f, "{\"status\": \"%s\", \"objective\": ", ipm_status_name(r.status));
    num_or_null(f, ok ? r.objective : kInf);
    std::fprintf(f, ", \"dual_objective\": ");
    num_or_null(f, ok ? r.dual_objective : kInf);
    std::fprintf(f, ", \"primal_res\": %.3e, \"dual_res\": %.3e, \"gap\": %.3e, \"max_row_viol\": %.3e, "
                 "\"max_bound_viol\": %.3e, \"iterations\": %ld, \"wall_s\": %.6f, \"message\": \"%s\"}\n",
                 r.primal_res, r.dual_res, r.gap, r.max_row_viol, r.max_bound_viol, r.iterations, wall,
                 json_escape(r.message).c_str());
    std::fclose(f);
}

void write_sol(const char* path, const Model& md, const std::vector<double>& x) {
    if (std::FILE* f = std::fopen(path, "w")) {
        for (size_t j = 0; j < md.col_names.size(); ++j) std::fprintf(f, "%s %.17g\n", md.col_names[j].c_str(), x[j]);
        std::fclose(f);
    }
}
// Exit codes: 0 definitive answer (optimal, infeasible, unbounded, ...), 2 usage, 3 parse error,
// 4 stopped at a time/iteration/node limit, 5 any other failure (numerical, unsupported, nonconvex).
int exit_code(const std::string& status) {
    if (status == "optimal" || status == "infeasible" || status == "unbounded" || status == "dual_infeasible" ||
        status == "unbounded_relaxation")
        return 0;
    if (status == "parse_error") return 3;
    if (status == "time_limit" || status == "iteration_limit" || status == "node_limit") return 4;
    return 5;
}
}  // namespace

int main(int argc, char** argv) {
    const char *model = nullptr, *sol = nullptr, *json = nullptr;
    double limit = 60;
    long node_limit = std::numeric_limits<long>::max();
    std::string method = "simplex";  // "ipm": interior point; "dual": dual simplex (LPs)
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--time-limit") && i + 1 < argc) limit = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--sol") && i + 1 < argc) sol = argv[++i];
        else if (!std::strcmp(argv[i], "--json") && i + 1 < argc) json = argv[++i];
        else if (!std::strcmp(argv[i], "--node-limit") && i + 1 < argc) node_limit = std::atol(argv[++i]);
        else if (!std::strcmp(argv[i], "--method") && i + 1 < argc) method = argv[++i];
        else model = argv[i];
    }
    if (!model) {
        std::fprintf(stderr, "usage: taral MODEL.mps [--time-limit S] [--node-limit N] [--sol OUT.sol] [--json OUT.json]\n");
        return 2;
    }
    auto t0 = std::chrono::steady_clock::now();
    auto wall = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    Model md;
    try {
        md = read_mps(model);
    } catch (const ParseError& e) {
        if (json) write_json(json, "parse_error", nullptr, wall(), e.what());
        std::printf("status parse_error: %s\n", e.what());
        return exit_code("parse_error");
    }
    if ((!md.qobj.empty() || method == "ipm") && !md.has_integers()) {  // convex QP, or LP by interior point
        IpmOptions opt;
        opt.time_limit = limit - wall();
        IpmResult r = ipm_solve(md, opt);
        double w = wall();
        if (json) write_ipm_json(json, r, w);
        if (sol && r.status == IpmStatus::Optimal) write_sol(sol, md, r.x);
        std::printf("status %s objective %.12g iterations %ld wall %.3fs %s\n", ipm_status_name(r.status), r.objective,
                    r.iterations, w, r.message.c_str());
        return exit_code(ipm_status_name(r.status));
    }
    if (!md.qobj.empty()) {  // quadratic objective with integer variables: never a silent relaxation
        std::string why = "mixed-integer quadratic models are not supported";
        if (json) write_json(json, "unsupported", nullptr, wall(), why);
        std::printf("status unsupported: %s\n", why.c_str());
        return exit_code("unsupported");
    }
    if (md.has_integers()) {
        MilpResult r = solve_milp(md, limit - wall(), node_limit);
        double w = wall();
        if (json) write_milp_json(json, r, w);
        if (sol && r.has_solution) write_sol(sol, md, r.x);
        std::printf("status %s objective %.12g best_bound %.12g gap %.3g nodes %ld iterations %ld wall %.3fs %s\n",
                    r.status.c_str(), r.has_solution ? r.objective : NAN, r.best_bound, r.gap, r.nodes, r.lp_iterations,
                    w, r.message.c_str());
        return exit_code(r.status);
    }
    Result r = method == "dual" ? solve_lp_dual(md, md.col_lo, md.col_up, nullptr, limit - wall())
                                  : solve(md, limit - wall());
    double w = wall();
    if (json) write_json(json, status_name(r.status), &r, w, r.message);
    if (sol && r.status == Status::Optimal) write_sol(sol, md, r.x);
    std::printf("status %s objective %.12g iterations %ld wall %.3fs %s\n", status_name(r.status), r.objective,
                r.iterations, w, r.message.c_str());
    return exit_code(status_name(r.status));
}
