// CLI per the engine contract:  taral MODEL.mps --time-limit S --sol OUT.sol --json OUT.json
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
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

std::string g_json_extra;  // optional ", \"key\": value" fragments appended to the LP JSON (--work-limit)
// --presolve reports the ORIGINAL-space primal point only. Dual postsolve is out of scope, so dual/KKT fields of the
// reduced model are written as null, never as placeholders next to an original-space x, and reduced-model
// infeasible/unbounded certificates are not mapped back (no original-space proof is claimed).
static bool g_primal_only = false;
static void presolve_sanitize(Result& r, const Model& orig) {
    g_primal_only = true;
    if (r.status == Status::Infeasible || r.status == Status::Unbounded) {
        r.certificate_quality = "presolve_reduced_model_only";
        r.certificate_verified = false;
        r.farkas_row_lower.clear(); r.farkas_row_upper.clear(); r.farkas_col_lower.clear(); r.farkas_col_upper.clear();
        r.ray.clear(); r.x.clear();
        r.message = "reduced-model result after presolve; the certificate is not mapped back and is not an original-space proof. " + r.message;
        return;
    }
    r.certificate_quality = "presolve_primal_only";
    r.row_activity.assign(orig.row_lo.size(), 0.0);
    for (size_t j = 0; j < orig.cols.size() && j < r.x.size(); ++j)
        for (const Entry& e : orig.cols[j]) r.row_activity[e.index] += e.value * r.x[j];
}

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
        if (g_primal_only) {
            std::fprintf(f, ", \"original_space_primal_only\": true");
            for (const char* nm : {"dual_objective", "primal_res", "dual_res", "gap", "complementarity", "max_row_viol", "max_bound_viol",
                                   "max_row_violation_magnitude_scaled", "row_violation_abs", "row_violation_magnitude_scaled",
                                   "row_term_magnitude", "row_dual", "reduced_cost"})
                std::fprintf(f, ", \"%s\": null", nm);
            array("x", r->x);
            array("row_activity", r->row_activity);
        } else {
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
    }
    if (r && (r->status == Status::Infeasible || r->status == Status::Unbounded ||
              !r->farkas_row_lower.empty() || !r->farkas_col_lower.empty() || !r->ray.empty())) {
        std::fprintf(f, ", \"certificate_verified\": %s, \"certificate_residual\": %.17g, \"certificate_margin\": %.17g",
                     r->certificate_verified ? "true" : "false", r->certificate_residual, r->certificate_margin);
        auto array = [&](const char* name, const std::vector<double>& values) {
            std::fprintf(f, ", \"%s\": [", name);
            for (size_t k = 0; k < values.size(); ++k) {
                if (k) std::fprintf(f, ", ");
                if (std::isfinite(values[k])) std::fprintf(f, "%.17g", values[k]);
                else std::fprintf(f, "null");
            }
            std::fprintf(f, "]");
        };
        array("farkas_row_lower", r->farkas_row_lower); array("farkas_row_upper", r->farkas_row_upper);
        array("farkas_col_lower", r->farkas_col_lower); array("farkas_col_upper", r->farkas_col_upper);
        if (!(r->status == Status::Optimal || !r->x.empty())) array("x", r->x);  // the certificate block above already wrote x
        array("ray", r->ray);
    }
    std::fprintf(f, ", \"iterations\": %ld, \"wall_s\": %.6f, \"message\": \"%s\"%s}\n", r ? r->iterations : 0L, wall,
                 json_escape(msg).c_str(), g_json_extra.c_str());
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
    std::fprintf(f, ", \"prop_tightened\": %ld, \"prop_crossed\": %ld, \"prop_crossed_lp_infeasible\": %ld, \"prop_pruned\": %ld",
                 r.prop_tightened, r.prop_crossed, r.prop_crossed_lp_infeasible, r.prop_pruned);
    std::fprintf(f, ", \"rc_fixed\": %ld, \"rc_skipped\": %ld", r.rc_fixed, r.rc_skipped);
    std::fprintf(f, ", \"audit\": [");
    for (size_t k = 0; k < r.audit.size(); ++k) std::fprintf(f, "%s%ld", k ? ", " : "", r.audit[k]);
    std::fprintf(f, "]");
    if (r.persistent_nodes_used) std::fprintf(f, ", \"peak_trail_records\": %ld, \"peak_trail_nodes\": %ld",r.peak_trail_records,r.peak_trail_nodes);
    if (r.compact_nodes_used) std::fprintf(f, ", \"peak_open_nodes\": %ld, \"peak_open_changes\": %ld, \"max_depth\": %ld", r.peak_open_nodes,r.peak_open_changes,r.max_depth);
    if (!r.structural_certificate.empty()) {
        std::fprintf(f, ", \"structural_certificate\": %s", r.structural_certificate.c_str());
        if (r.has_solution) {
            std::fprintf(f, ", \"x\": [");
            for (size_t j=0;j<r.x.size();++j) { if(j)std::fprintf(f, ",");num_or_null(f,r.x[j]); }
            std::fprintf(f, "]");
        }
    }
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


bool write_explanation_json(const char* path, const Model& md, const InfeasibilityExplanation& e) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return false;
    auto arr = [&](const std::vector<double>& v) {
        std::fprintf(f, "[");
        for (size_t k = 0; k < v.size(); ++k) { if (k) std::fprintf(f, ", "); num_or_null(f, v[k]); }
        std::fprintf(f, "]");
    };
    std::fprintf(f, "{\"scope\": \"lp_relaxation_row_irreducible_retained_column_bounds\", \"status\": \"%s\", \"message\": \"%s\"",
                 e.status.c_str(), json_escape(e.message).c_str());
    std::fprintf(f, ", \"witness_tolerance\": 1e-7, \"certificate_verified\": %s, \"certificate_margin\": ",
                 e.certificate_verified ? "true" : "false");
    num_or_null(f, e.certificate_margin);
    std::fprintf(f, ", \"certificate_residual\": ");
    num_or_null(f, e.certificate_residual);
    std::fprintf(f, ", \"rows\": [");
    for (size_t k = 0; k < e.rows.size(); ++k) {
        int i = e.rows[k];
        std::fprintf(f, "%s{\"index\": %d, \"name\": \"%s\", \"lower\": ", k ? ", " : "", i, json_escape(md.row_names[i]).c_str());
        num_or_null(f, md.row_lo[i]);
        std::fprintf(f, ", \"upper\": ");
        num_or_null(f, md.row_up[i]);
        std::fprintf(f, ", \"multiplier_lower\": ");
        num_or_null(f, e.farkas_row_lower[i]);
        std::fprintf(f, ", \"multiplier_upper\": ");
        num_or_null(f, e.farkas_row_upper[i]);
        std::fprintf(f, ", \"removal_witness_violation\": ");
        num_or_null(f, e.witness[k].empty() ? kInf : e.witness_violation[k]);
        std::fprintf(f, ", \"removal_witness\": ");
        if (e.witness[k].empty()) std::fprintf(f, "null"); else arr(e.witness[k]);
        std::fprintf(f, "}");
    }
    std::fprintf(f, "], \"unproven_rows\": [");
    for (size_t k = 0; k < e.unproven_rows.size(); ++k) std::fprintf(f, "%s%d", k ? ", " : "", e.unproven_rows[k]);
    std::fprintf(f, "], \"bound_columns\": [");
    for (size_t k = 0; k < e.bound_columns.size(); ++k) {
        int j = e.bound_columns[k];
        std::fprintf(f, "%s{\"index\": %d, \"name\": \"%s\", \"lower\": ", k ? ", " : "", j, json_escape(md.col_names[j]).c_str());
        num_or_null(f, md.col_lo[j]);
        std::fprintf(f, ", \"upper\": ");
        num_or_null(f, md.col_up[j]);
        std::fprintf(f, ", \"multiplier_lower\": ");
        num_or_null(f, e.farkas_col_lower[j]);
        std::fprintf(f, ", \"multiplier_upper\": ");
        num_or_null(f, e.farkas_col_upper[j]);
        std::fprintf(f, "}");
    }
    std::fprintf(f, "], \"inconsistent_bound_columns\": [");
    for (size_t k = 0; k < e.inconsistent_bound_columns.size(); ++k)
        std::fprintf(f, "%s%d", k ? ", " : "", e.inconsistent_bound_columns[k]);
    std::fprintf(f, "]");
    std::fprintf(f, ", \"farkas_row_lower\": "); arr(e.farkas_row_lower);
    std::fprintf(f, ", \"farkas_row_upper\": "); arr(e.farkas_row_upper);
    std::fprintf(f, ", \"farkas_col_lower\": "); arr(e.farkas_col_lower);
    std::fprintf(f, ", \"farkas_col_upper\": "); arr(e.farkas_col_upper);
    std::fprintf(f, ", \"relaxation\": {\"status\": \"%s\", \"objective\": ", e.relaxation_status.c_str());
    num_or_null(f, e.relaxation_status == "optimal" ? e.relaxation_objective : kInf);
    std::fprintf(f, ", \"objective_definition\": \"sum_i w_i*(lower_relax_i+upper_relax_i), w_i=1/(1+max finite |row bound|), all rows, column bounds retained\"");
    std::fprintf(f, ", \"scale_factor\": "); num_or_null(f, e.relaxation_scale > 0 ? e.relaxation_scale : kInf);  // null when the scaled LP was never attempted
    std::fprintf(f, ", \"objective_note\": \"reporting-only, not claimed minimal; objective is the solver value on the scaled LP divided by scale_factor\"");
    std::fprintf(f, ", \"certificate_quality\": \"%s\", \"certificate_quality_model\": \"scaled_lp\"", relaxation_quality_label(e));
    std::fprintf(f, ", \"kkt_gap_scaled_model\": "); num_or_null(f, e.relaxation_gap);
    std::fprintf(f, ", \"kkt_gap_abs_unscaled\": "); num_or_null(f, e.relaxation_gap_abs_unscaled);
    std::fprintf(f, ", \"kkt_gap\": "); num_or_null(f, e.relaxation_gap);
    std::fprintf(f, ", \"max_violation\": "); num_or_null(f, e.relaxation_violation);
    std::fprintf(f, ", \"rows\": [");
    bool first = true;
    for (size_t i = 0; i < e.relax_lower.size(); ++i) {
        double lo = e.relax_lower[i], up = e.relax_upper[i];
        if (!(lo > 1e-12 * (1 + std::abs(md.row_lo[i])) || up > 1e-12 * (1 + std::abs(md.row_up[i])))) continue;
        std::fprintf(f, "%s{\"index\": %zu, \"name\": \"%s\", \"lower_relaxed_by\": ", first ? "" : ", ", i,
                     json_escape(md.row_names[i]).c_str());
        first = false;
        num_or_null(f, lo);
        std::fprintf(f, ", \"upper_relaxed_by\": ");
        num_or_null(f, up);
        std::fprintf(f, ", \"weight\": ");
        num_or_null(f, e.relax_weight[i]);
        std::fprintf(f, "}");
    }
    std::fprintf(f, "], \"x\": "); arr(e.relaxation_x);
    std::fprintf(f, "}, \"lp_solves\": %ld, \"wall_s\": %.6f}\n", e.lp_solves, e.wall_s);
    bool ok = !std::ferror(f);
    if (std::fclose(f) != 0) ok = false;
    return ok;
}

// FNV-1a 64 over a canonical text of (status, iterations, objective, x) with %.17g. Excludes wall time and the
// message (which carries timings). Bit-exact: equal only when the arithmetic result is identical.
unsigned long long result_hash(const char* status, const Result& r) {
    unsigned long long h = 14695981039346656037ULL;
    auto feed = [&](const std::string& t) {
        for (unsigned char c : t) h = (h ^ c) * 1099511628211ULL;
        h = (h ^ 0xff) * 1099511628211ULL;
    };
    char b[48];
    feed(status);
    feed(std::to_string(r.iterations));
    std::snprintf(b, sizeof b, "%.17g", r.status == Status::Optimal ? r.objective : 0.0);
    feed(b);
    for (double v : r.x) { std::snprintf(b, sizeof b, "%.17g", v); feed(b); }
    return h;
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
const char* kUsage =
    "usage: taral MODEL.mps [--time-limit S] [--node-limit N] [--method simplex|dual|ipm]\n"
    "             [--sol OUT.sol] [--json OUT.json] [--warm-sol F [--warm-dual F] [--cross-tol T]]\n"
    "             [--no-fallback] [--audit-prop] [--no-prop-prune] [--explain-infeasible OUT.json]\n"
    "             [--work-limit N]\n"
    "             [--presolve] [--presolve-log OUT.json]\n";

void print_help() {
    std::printf("%s"
                "\n"
                "options:\n"
                "  --time-limit S    wall-clock limit in seconds; finite and > 0 (default 60)\n"
                "  --node-limit N    branch-and-bound node limit; integer >= 0 (default unlimited)\n"
                "  --method M        simplex (default), dual (dual simplex, LPs), or ipm (interior point)\n"
                "  --warm-sol FILE   crossover: start from the basis of an approximate primal point (name value per line)\n"
                "  --warm-dual FILE  optional row multipliers for --warm-sol (name value per line)\n"
                "  --cross-tol T     crossover: distance from a bound that counts as interior (default 1e-3)\n"
                "  --no-fallback     LP simplex: do not hand a stalled primal run to the dual simplex\n"
                "  --work-limit N    LP (simplex/dual): deterministic iteration budget; routing never uses the wall clock (an explicit --time-limit that fires still stops the run, non-deterministically)\n"
                "  --presolve        verified presolve (linear LP/MILP): replayable log, original-space audit of the result\n"
                "  --presolve-log F  write the presolve log and audit to F (implies --presolve)\n"
                "  --persistent-nodes MILP: shared checkpointed domains (opt-in)\n"
                "  --compact-nodes   MILP: canonical interval storage and bounded dives (opt-in)\n"
                "  --integer-structure  MILP: opt-in exact structural parity proof\n"
                "  --audit-prop      MILP: re-check propagation prunes (diagnostic)\n"
                "  --no-prop-prune   MILP: disable propagation pruning\n"
                "  --explain-infeasible FILE  write a verified LP-relaxation infeasibility explanation (irreducible rows + minimum relaxation) and exit\n"
                "  --sol FILE        write the solution (name value per line) when optimal\n"
                "  --json FILE       write a JSON result summary\n"
                "  -h, --help        show this help and exit\n",
                kUsage);
}

int usage_error(const std::string& msg) {
    std::fprintf(stderr, "taral: %s\n%sTry 'taral --help'.\n", msg.c_str(), kUsage);
    return 2;
}

bool parse_double(const char* s, double& out) {
    if (!*s) return false;
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s, &end);
    if (*end || errno == ERANGE || !std::isfinite(v)) return false;
    out = v;
    return true;
}

bool parse_long(const char* s, long& out) {
    if (!*s) return false;
    char* end = nullptr;
    errno = 0;
    long v = std::strtol(s, &end, 10);
    if (*end || errno == ERANGE) return false;
    out = v;
    return true;
}
}  // namespace

int main(int argc, char** argv) {
    const char *model = nullptr, *sol = nullptr, *json = nullptr;
    double limit = 60;
    long node_limit = std::numeric_limits<long>::max();
    const char *warm_sol = nullptr, *warm_dual = nullptr;
    double cross_tol = 1e-3;
    MilpOptions mopt;
    std::string method = "simplex";  // "ipm": interior point; "dual": dual simplex (LPs)
    bool fallback = true;            // primal simplex that stalls hands the rest of the time to the dual (--no-fallback: off)
    bool only_files = false;
    const char* explain_path = nullptr;  // --explain-infeasible FILE
    long work_limit = 0;       // --work-limit N: deterministic iteration budget (LP simplex/dual)
    bool time_limit_given = false;
    bool do_presolve = false;
    const char* presolve_log = nullptr;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!only_files && !std::strcmp(a, "--")) { only_files = true; continue; }
        if (!only_files && (!std::strcmp(a, "--help") || !std::strcmp(a, "-h"))) { print_help(); return 0; }
        if (!only_files && a[0] == '-' && a[1]) {
            auto value = [&](const char*& v) {
                if (i + 1 >= argc) return false;
                v = argv[++i];
                return true;
            };
            const char* v = nullptr;
            if (!std::strcmp(a, "--time-limit")) {
                if (!value(v)) return usage_error("--time-limit requires a value");
                time_limit_given = true;
                if (!parse_double(v, limit) || !(limit > 0))
                    return usage_error(std::string("--time-limit must be a finite number > 0, got '") + v + "'");
            } else if (!std::strcmp(a, "--work-limit")) {
                if (!value(v)) return usage_error("--work-limit requires a value");
                if (!parse_long(v, work_limit) || work_limit < 1)
                    return usage_error(std::string("--work-limit must be an integer >= 1, got '") + v + "'");
            } else if (!std::strcmp(a, "--node-limit")) {
                if (!value(v)) return usage_error("--node-limit requires a value");
                if (!parse_long(v, node_limit) || node_limit < 0)
                    return usage_error(std::string("--node-limit must be an integer >= 0, got '") + v + "'");
            } else if (!std::strcmp(a, "--sol")) {
                if (!value(sol)) return usage_error("--sol requires a value");
            } else if (!std::strcmp(a, "--json")) {
                if (!value(json)) return usage_error("--json requires a value");
            } else if (!std::strcmp(a, "--method")) {
                if (!value(v)) return usage_error("--method requires a value");
                method = v;
                if (method != "simplex" && method != "dual" && method != "ipm")
                    return usage_error("unknown --method '" + method + "' (expected simplex, dual or ipm)");
            } else if (!std::strcmp(a, "--warm-sol")) {
                if (!value(warm_sol)) return usage_error("--warm-sol requires a value");
            } else if (!std::strcmp(a, "--warm-dual")) {
                if (!value(warm_dual)) return usage_error("--warm-dual requires a value");
            } else if (!std::strcmp(a, "--explain-infeasible")) {
                if (!value(explain_path)) return usage_error("--explain-infeasible requires a value");
            } else if (!std::strcmp(a, "--presolve")) {
                do_presolve = true;
            } else if (!std::strcmp(a, "--presolve-log")) {
                if (!value(presolve_log)) return usage_error("--presolve-log requires a value");
                do_presolve = true;
            } else if (!std::strcmp(a, "--no-fallback")) {
                fallback = false;
            } else if (!std::strcmp(a, "--persistent-nodes")) {
                mopt.persistent_nodes = true;
            } else if (!std::strcmp(a, "--compact-nodes")) {
                mopt.compact_nodes = true;
            } else if (!std::strcmp(a, "--integer-structure")) {
                mopt.integer_structure = true;
            } else if (!std::strcmp(a, "--audit-prop")) {
                mopt.audit_prop = true;
            } else if (!std::strcmp(a, "--no-prop-prune")) {
                mopt.no_prop_prune = true;
            } else if (!std::strcmp(a, "--cross-tol")) {
                if (!value(v)) return usage_error("--cross-tol requires a value");
                if (!parse_double(v, cross_tol) || !(cross_tol > 0))
                    return usage_error(std::string("--cross-tol must be a finite number > 0, got '") + v + "'");
            } else {
                return usage_error(std::string("unknown option '") + a + "'");
            }
        } else {
            if (model) return usage_error(std::string("unexpected extra argument '") + a + "'");
            model = a;
        }
    }
    if (!model) {
        std::fprintf(stderr, "usage: taral MODEL.mps [--time-limit S] [--node-limit N] [--sol OUT.sol] [--json OUT.json]\n");
        return 2;
    }
    if (mopt.integer_structure && do_presolve)
        return usage_error("--integer-structure and --presolve cannot be combined until certificate postsolve is supported");
    if (work_limit) {
        if (do_presolve) return usage_error("--work-limit and --presolve cannot be combined (presolve work accounting is not supported)");
        if (method == "ipm") return usage_error("--work-limit supports --method simplex or dual only");
        if (explain_path) return usage_error("--work-limit does not apply to --explain-infeasible");
        work_budget().cap = work_limit;
        if (!time_limit_given) limit = 1e9;  // the work limit, not the wall clock, ends the solve
    }
    auto t0 = std::chrono::steady_clock::now();
    auto wall = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    Model md;
    try {
        md = read_mps(model, limit < 1e8 ? t0 + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(limit))
                              : std::chrono::steady_clock::time_point::max());
    } catch (const ParseTimeLimit& e) {
        if (json) write_json(json, "time_limit", nullptr, wall(), e.what());
        std::printf("status time_limit: %s\n", e.what());
        return exit_code("time_limit");
    } catch (const ParseError& e) {
        if (json) write_json(json, "parse_error", nullptr, wall(), e.what());
        std::printf("status parse_error: %s\n", e.what());
        return exit_code("parse_error");
    }
    if (explain_path) {  // LP-relaxation infeasibility explanation; ignores integrality and the objective
        InfeasibilityExplanation ex = explain_infeasibility(md, limit - wall());
        if (!write_explanation_json(explain_path, md, ex)) {
            std::fprintf(stderr, "error: cannot write explanation to '%s'\n", explain_path);
            std::printf("status output_error: cannot write '%s'\n", explain_path);
            return 5;
        }
        std::printf("explain %s rows %zu unproven %zu solves %ld wall %.3fs %s\n", ex.status.c_str(), ex.rows.size(),
                    ex.unproven_rows.size(), ex.lp_solves, wall(), ex.message.c_str());
        if (ex.status == "no_verified_proof") return 5;
        return ex.status == "reduced_unproven" && ex.wall_s >= limit - 1e-3 ? 4 : 0;
    }
    if (work_limit && !md.qobj.empty()) return usage_error("--work-limit supports linear LPs only");
    Model orig_md;  // --presolve: solve a reduced model, report and audit in original space
    PresolveResult pre;
    if (do_presolve) {
        if (!md.qobj.empty() || method == "ipm")
            return usage_error("--presolve supports linear models with --method simplex or dual only");
        orig_md = md;
        auto log_fail = [&]() {
            std::fprintf(stderr, "error: cannot write presolve log to '%s'\n", presolve_log);
            std::printf("status output_error: cannot write '%s'\n", presolve_log);
            return 5;
        };
        const auto dl = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                                                              std::chrono::duration<double>(std::max(0.0, std::min(limit - wall(), 1e8))));
        long hook = -1;  // test hook: TARAL_PRESOLVE_TEST_TIMEOUT_AFTER_OPS=K behaves as if the deadline passed after K reductions
        if (const char* h = std::getenv("TARAL_PRESOLVE_TEST_TIMEOUT_AFTER_OPS")) hook = std::atol(h);
        pre = presolve_model(orig_md, dl, hook);
        if (pre.timed_out) {  // honest incomplete state: the partial reduction is not used
            if (presolve_log && !write_presolve_log(presolve_log, orig_md, pre, nullptr)) return log_fail();
            if (json) write_json(json, "time_limit", nullptr, wall(), "presolve did not finish within the time limit");
            std::printf("status time_limit (presolve incomplete after %zu reductions) wall %.3fs\n", pre.log.size(), wall());
            return exit_code("time_limit");
        }
        if (pre.infeasible) {
            if (presolve_log && !write_presolve_log(presolve_log, orig_md, pre, nullptr)) return log_fail();
            if (json) write_json(json, "infeasible", nullptr, wall(), "presolve: " + pre.infeasible_reason);
            std::printf("status infeasible (presolve: %s) wall %.3fs\n", pre.infeasible_reason.c_str(), wall());
            return exit_code("infeasible");
        }
        std::printf("presolve: %zu -> %zu rows, %zu -> %zu cols, %zu reductions, %d passes\n", orig_md.row_lo.size(),
                    pre.kept_rows.size(), orig_md.cols.size(), pre.kept_cols.size(), pre.log.size(), pre.passes);
        if (presolve_log && !write_presolve_log(presolve_log, orig_md, pre, nullptr)) return log_fail();  // rewritten with the audit when a point is reported
        md = pre.reduced;
        if (md.cols.empty()) {  // everything fixed: the answer is the substitution itself
            std::vector<double> x = presolve_expand(orig_md, pre, {});
            PresolveAudit au = presolve_audit(orig_md, x, md.obj_const);
            if (presolve_log && !write_presolve_log(presolve_log, orig_md, pre, &au)) return log_fail();
            const char* st = au.ok ? "optimal" : "numerical_failure";
            Result full;  // original-space result, so the JSON carries the objective and point
            full.status = au.ok ? Status::Optimal : Status::NumericalFailure;
            full.x = x;
            full.objective = au.objective;
            full.certificate_quality = "presolve_substitution";
            full.message = "presolve solved the model";
            presolve_sanitize(full, orig_md);
            if (json) write_json(json, st, &full, wall(), "presolve solved the model");
            if (sol && au.ok) write_sol(sol, orig_md, x);
            std::printf("status %s objective %.12g (presolve) wall %.3fs\n", st, au.objective, wall());
            return exit_code(st);
        }
    }
    auto presolve_report = [&](std::vector<double>& x, double obj, bool& ok) {  // expand + audit; x becomes original-space
        PresolveAudit au;
        if (do_presolve) {
            x = presolve_expand(orig_md, pre, x);
            au = presolve_audit(orig_md, x, obj);
            ok = au.ok;
            if (presolve_log && !write_presolve_log(presolve_log, orig_md, pre, &au)) {
                std::fprintf(stderr, "error: cannot write presolve log to '%s'\n", presolve_log);
                ok = false;
            }
            std::printf("presolve audit %s: row %.3g bound %.3g int %.3g objective diff %.3g\n", au.ok ? "ok" : "FAILED",
                        au.max_row_violation, au.max_bound_violation, au.max_int_violation, au.objective_diff);
        }
    };
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
    if (work_limit && md.has_integers()) {  // B&B node LPs share one counter; no determinism claim for MILP
        std::string why = "--work-limit supports continuous LPs only (use --node-limit for MILP)";
        if (json) write_json(json, "unsupported", nullptr, wall(), why);
        std::printf("status unsupported: %s\n", why.c_str());
        return exit_code("unsupported");
    }
    if (md.has_integers()) {
        MilpResult r = solve_milp(md, limit - wall(), node_limit, mopt);
        if (r.status == "optimal" && !std::isfinite(r.objective))
            r.status = "numerical_failure", r.has_solution = false, r.message = "non-finite objective";
        if (r.has_solution && do_presolve) {
            bool ok = true;
            presolve_report(r.x, r.objective, ok);
            if (!ok) r.status = "numerical_failure", r.has_solution = false, r.message = "presolve audit or log write failed";
        }
        double w = wall();
        if (json) write_milp_json(json, r, w);
        if (sol && r.has_solution) write_sol(sol, do_presolve ? orig_md : md, r.x);
        std::printf("status %s objective %.12g best_bound %.12g gap %.3g nodes %ld iterations %ld prop_crossed %ld wall %.3fs %s\n",
                    r.status.c_str(), r.has_solution ? r.objective : NAN, r.best_bound, r.gap, r.nodes, r.lp_iterations,
                    r.prop_crossed, w, r.message.c_str());
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
    Result r = solve_lp_gated(md, md.col_lo, md.col_up, wp, limit - wall(), method == "dual", fallback);
    if (r.status == Status::Optimal && !std::isfinite(r.objective))
        r.status = Status::NumericalFailure, r.message = "non-finite objective";
    if (work_limit) {  // deterministic report: hash covers status, iterations, objective and x bits only
        char hx[32];
        std::snprintf(hx, sizeof hx, "%016llx", (unsigned long long)result_hash(status_name(r.status), r));
        g_json_extra = std::string(", \"work_limit\": ") + std::to_string(work_limit) + ", \"wall_limit_active\": " + (time_limit_given ? "true" : "false") + ", \"work_used\": " +
                       std::to_string(std::min(work_budget().used, work_limit)) + ", \"result_hash\": \"" + hx + "\"";
        std::printf("work_used %ld of %ld result_hash %s\n", std::min(work_budget().used, work_limit), work_limit, hx);
    }
    if (r.status == Status::Optimal && do_presolve) {
        bool ok = true;
        presolve_report(r.x, r.objective, ok);
        if (!ok) r.status = Status::NumericalFailure, r.message = "presolve audit or log write failed";
    }
    double w = wall();
    if (do_presolve) presolve_sanitize(r, orig_md);
    if (json) write_json(json, status_name(r.status), &r, w, r.message);
    if (sol && r.status == Status::Optimal) write_sol(sol, do_presolve ? orig_md : md, r.x);
    std::printf("status %s objective %.12g iterations %ld wall %.3fs %s\n", status_name(r.status), r.objective,
                r.iterations, w, r.message.c_str());
    return exit_code(status_name(r.status));
}
