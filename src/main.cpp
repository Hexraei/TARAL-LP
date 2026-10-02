// CLI per the engine contract:  taral MODEL.mps --time-limit S --sol OUT.sol --json OUT.json
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include <string>

#include "ipm.hpp"
#include "taral.hpp"

namespace {
// Escapes for a JSON string. Bytes that are not part of well-formed UTF-8 (model names are arbitrary bytes)
// become '?', so the file is always valid JSON.
std::string json_escape(const std::string& s) {
    std::string o;
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char ch = static_cast<unsigned char>(s[i]);
        if (ch == '"' || ch == '\\') o += '\\', o += static_cast<char>(ch);
        else if (ch < 0x20) o += ' ';
        else if (ch < 0x80) o += static_cast<char>(ch);
        else {
            size_t len = ch >= 0xF0 ? 4 : ch >= 0xE0 ? 3 : ch >= 0xC2 ? 2 : 0;
            bool ok = len > 0 && ch <= 0xF4 && i + len <= s.size();
            for (size_t k = 1; ok && k < len; ++k) ok = (static_cast<unsigned char>(s[i + k]) & 0xC0) == 0x80;
            if (ok && len > 2) {  // no overlong forms, surrogates or values above U+10FFFF
                unsigned char c1 = static_cast<unsigned char>(s[i + 1]);
                ok = !(ch == 0xE0 && c1 < 0xA0) && !(ch == 0xED && c1 > 0x9F) && !(ch == 0xF0 && c1 < 0x90) &&
                     !(ch == 0xF4 && c1 > 0x8F);
            }
            if (ok) o.append(s, i, len), i += len - 1;
            else o += '?';
        }
    }
    return o;
}

void num_or_null(std::FILE* f, double v);

void write_json(const char* path, const char* status, const Result* r, double wall, const std::string& msg) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    std::fprintf(f, "{\"status\": \"%s\", \"objective\": ", status);
    if (r && r->status == Status::Optimal) std::fprintf(f, "%.17g", r->objective);
    else std::fprintf(f, "null");
    std::fprintf(f, ", \"certificate_quality\": \"%s\", \"certificate_tolerance\": 1e-8",
                 r ? r->certificate_quality.c_str() : "unknown");
    if (r && (r->status == Status::Optimal || !r->x.empty())) {
        auto scalar = [&](const char* name, double value) {
            std::fprintf(f, ", \"%s\": ", name);
            num_or_null(f, value);
        };
        auto array = [&](const char* name, const std::vector<double>& values) {
            std::fprintf(f, ", \"%s\": [", name);
            for (size_t k = 0; k < values.size(); ++k) {
                if (k) std::fprintf(f, ", ");
                num_or_null(f, values[k]);
            }
            std::fprintf(f, "]");
        };
        scalar("dual_objective", r->dual_objective);
        scalar("primal_res", r->primal_res);
        scalar("dual_res", r->dual_res);
        scalar("gap", r->gap);
        scalar("complementarity", r->complementarity);
        scalar("max_row_viol", r->max_row_viol);
        scalar("max_bound_viol", r->max_bound_viol);
        scalar("max_row_violation_magnitude_scaled", r->max_row_violation_magnitude_scaled);
        array("row_violation_abs", r->row_violation_abs);
        array("row_violation_magnitude_scaled", r->row_violation_magnitude_scaled);
        array("row_term_magnitude", r->row_term_magnitude);
        array("x", r->x);
        array("row_activity", r->row_activity);
        array("row_dual", r->row_dual);
        array("reduced_cost", r->reduced_cost);
    }
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

void sci_or_null(std::FILE* f, double v) {
    if (std::isfinite(v)) std::fprintf(f, "%.3e", v);
    else std::fprintf(f, "null");
}

void write_ipm_json(const char* path, const IpmResult& r, double wall) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    bool ok = r.status == IpmStatus::Optimal;
    std::fprintf(f, "{\"status\": \"%s\", \"objective\": ", ipm_status_name(r.status));
    num_or_null(f, ok ? r.objective : kInf);
    std::fprintf(f, ", \"dual_objective\": ");
    num_or_null(f, ok ? r.dual_objective : kInf);
    std::fprintf(f, ", \"primal_res\": ");
    sci_or_null(f, r.primal_res);
    std::fprintf(f, ", \"dual_res\": ");
    sci_or_null(f, r.dual_res);
    std::fprintf(f, ", \"gap\": ");
    sci_or_null(f, r.gap);
    std::fprintf(f, ", \"max_row_viol\": ");
    sci_or_null(f, r.max_row_viol);
    std::fprintf(f, ", \"max_bound_viol\": ");
    sci_or_null(f, r.max_bound_viol);
    std::fprintf(f, ", \"iterations\": %ld, \"wall_s\": %.6f, \"message\": \"%s\"}\n", r.iterations, wall,
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
    const char *warm_sol = nullptr, *warm_dual = nullptr;
    double cross_tol = 1e-3;
    std::string method = "simplex";  // "ipm": interior point; "dual": dual simplex (LPs)
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--time-limit") && i + 1 < argc) limit = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--sol") && i + 1 < argc) sol = argv[++i];
        else if (!std::strcmp(argv[i], "--json") && i + 1 < argc) json = argv[++i];
        else if (!std::strcmp(argv[i], "--node-limit") && i + 1 < argc) node_limit = std::atol(argv[++i]);
        else if (!std::strcmp(argv[i], "--method") && i + 1 < argc) method = argv[++i];
        else if (!std::strcmp(argv[i], "--warm-sol") && i + 1 < argc) warm_sol = argv[++i];
        else if (!std::strcmp(argv[i], "--warm-dual") && i + 1 < argc) warm_dual = argv[++i];
        else if (!std::strcmp(argv[i], "--cross-tol") && i + 1 < argc) cross_tol = std::atof(argv[++i]);
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
        if (r.status == IpmStatus::Optimal && !std::isfinite(r.objective))
            r.status = IpmStatus::NumericalFailure, r.message = "non-finite objective";
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
        if (r.status == "optimal" && !std::isfinite(r.objective))
            r.status = "numerical_failure", r.has_solution = false, r.message = "non-finite objective";
        double w = wall();
        if (json) write_milp_json(json, r, w);
        if (sol && r.has_solution) write_sol(sol, md, r.x);
        std::printf("status %s objective %.12g best_bound %.12g gap %.3g nodes %ld iterations %ld wall %.3fs %s\n",
                    r.status.c_str(), r.has_solution ? r.objective : NAN, r.best_bound, r.gap, r.nodes, r.lp_iterations,
                    w, r.message.c_str());
        return exit_code(r.status);
    }
    std::vector<char> warm;  // crossover: start from the basis of an approximate (PDHG) primal point
    if (warm_sol) {
        std::unordered_map<std::string, double> val;
        std::ifstream in(warm_sol);
        std::string name;
        double v;
        while (in >> name >> v) val[name] = v;
        std::vector<double> x(md.col_names.size(), 0.0);
        for (size_t j = 0; j < x.size(); ++j) x[j] = val.count(md.col_names[j]) ? val[md.col_names[j]] : 0.0;
        std::vector<double> y;  // optional row multipliers, in row order
        if (warm_dual) {
            std::unordered_map<std::string, double> yv;
            std::ifstream din(warm_dual);
            while (din >> name >> v) yv[name] = v;
            for (const std::string& r : md.row_names) y.push_back(yv.count(r) ? yv[r] : 0.0);
        }
        int interior = 0;
        warm = basis_from_point(md, x, y, cross_tol, &interior);
        std::fprintf(stderr, "crossover: %d variables farther than %g from a bound, %zu rows (%.3fs)\n", interior, cross_tol,
                     md.row_names.size(), wall());
    }
    const std::vector<char>* wp = warm.empty() ? nullptr : &warm;
    Result r = method == "dual" ? solve_lp_dual(md, md.col_lo, md.col_up, wp, limit - wall())
                                  : solve_lp(md, md.col_lo, md.col_up, wp, limit - wall());
    if (r.status == Status::Optimal && !std::isfinite(r.objective))
        r.status = Status::NumericalFailure, r.message = "non-finite objective";
    double w = wall();
    if (json) write_json(json, status_name(r.status), &r, w, r.message);
    if (sol && r.status == Status::Optimal) write_sol(sol, md, r.x);
    std::printf("status %s objective %.12g iterations %ld wall %.3fs %s\n", status_name(r.status), r.objective,
                r.iterations, w, r.message.c_str());
    return exit_code(status_name(r.status));
}
