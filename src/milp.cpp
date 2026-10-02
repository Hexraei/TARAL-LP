// LP-based branch and bound for mixed-integer models. Best-bound node selection with a depth-first
// dive after every branching (the child on the rounding side of the branching variable), most-
// fractional branching, a rounding heuristic at every node, and warm-started node LPs.
// A node is pruned only on a verified LP status: Optimal (by bound) or Infeasible. Any other LP
// outcome keeps the node as unresolved; its inherited bound stays in the reported best bound.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <queue>

#include "taral.hpp"

namespace {
constexpr double kIntTol = 1e-6;   // |x - round(x)| for an integral value
constexpr double kFeasTol = 1e-6;  // row and bound violation, relative to 1 + |bound|
static_assert(kIntTol == kFeasTol, "violation() compares integrality on the feasibility scale");
constexpr double kGapTol = 1e-6;   // optimal when objective - bound <= kGapTol * max(1, |objective|)

double tol_at(double z) { return kGapTol * std::max(1.0, std::abs(z)); }

// Largest violation of rows, column bounds and integrality (all scaled as documented above).
double violation(const Model& md, const std::vector<double>& x) {
    int n = int(md.cols.size()), m = int(md.row_lo.size());
    std::vector<double> act(m, 0.0);
    double worst = 0;
    auto check = [&](double v, double lo, double up) {
        if (std::isfinite(lo)) worst = std::max(worst, (lo - v) / (1 + std::abs(lo)));
        if (std::isfinite(up)) worst = std::max(worst, (v - up) / (1 + std::abs(up)));
    };
    for (int j = 0; j < n; ++j) {
        if (!std::isfinite(x[j])) return kInf;
        for (const Entry& e : md.cols[j]) act[e.index] += e.value * x[j];
        check(x[j], md.col_lo[j], md.col_up[j]);
        // kIntTol == kFeasTol, so the integrality residual is compared on the same scale
        if (md.is_int[j]) worst = std::max(worst, std::abs(x[j] - std::round(x[j])));
    }
    for (int i = 0; i < m; ++i) check(act[i], md.row_lo[i], md.row_up[i]);
    return worst;
}

double objective_of(const Model& md, const std::vector<double>& x) {
    double z = md.obj_const;
    for (size_t j = 0; j < x.size(); ++j) z += md.cost[j] * x[j];
    return z;
}

struct Change {
    int j;
    double lo, up;
};

struct Node {
    double bound;  // valid lower bound on every solution in this node (min form)
    std::vector<Change> changes;  // bound tightenings relative to the root, in branching order
    std::shared_ptr<const std::vector<char>> basis;  // parent's optimal basis (warm start)
    int dir = 0;      // last branching: -1 down, +1 up (variable = changes.back().j), 0 root
    double dist = 0;  // distance the branched variable was moved (fraction)
};

struct Worse {  // priority_queue keeps the smallest bound on top; deeper first on ties
    bool operator()(const Node& a, const Node& b) const {
        return a.bound != b.bound ? a.bound > b.bound : a.changes.size() < b.changes.size();
    }
};

enum class End { Complete, TimeLimit, NodeLimit, RootUnbounded, RootFailed };

// Search on a minimisation model. Integer column bounds are already rounded in root_lo/root_up.
class Search {
public:
    Search(const Model& md, std::vector<double> root_lo, std::vector<double> root_up, double time_limit, long node_limit)
        : md_(md), rlo_(std::move(root_lo)), rup_(std::move(root_up)), limit_(time_limit), node_limit_(node_limit) {}

    End run();
    bool has_inc = false;
    std::vector<double> inc;
    double z = kInf;      // incumbent objective
    double bound = -kInf;  // proven global lower bound after run()
    long nodes = 0, iters = 0, unresolved = 0;
    std::string note;

private:
    const Model& md_;
    std::vector<double> rlo_, rup_;
    double limit_;
    long node_limit_;
    std::chrono::steady_clock::time_point t0_ = std::chrono::steady_clock::now();
    double elapsed() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count(); }

    // Offer an LP point: integer columns rounded first, then as is. Returns true if accepted.
    bool offer(const std::vector<double>& x) {
        std::vector<double> r = x;
        for (size_t j = 0; j < r.size(); ++j)
            if (md_.is_int[j]) r[j] = std::round(r[j]);
        for (const std::vector<double>* c : {static_cast<const std::vector<double>*>(&r), &x}) {
            if (violation(md_, *c) > kFeasTol) continue;
            double v = objective_of(md_, *c);
            if (!has_inc || v < z) has_inc = true, z = v, inc = *c;
            return true;
        }
        return false;
    }
};

End Search::run() {
    std::priority_queue<Node, std::vector<Node>, Worse> open;
    std::vector<double> held;  // bounds of nodes that could not be resolved
    double pruned = kInf;      // smallest bound among nodes pruned by bound
    int n = int(md_.cols.size());
    Node cur{-kInf, {}, nullptr, 0, 0};
    // Pseudocosts: objective gain per unit of fractionality, per variable and direction.
    std::vector<double> pc_sum[2] = {std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
    std::vector<int> pc_cnt[2] = {std::vector<int>(n, 0), std::vector<int>(n, 0)};
    double all_sum[2] = {0, 0};
    long all_cnt[2] = {0, 0};
    auto pseudo = [&](int d, int j) {
        if (pc_cnt[d][j]) return pc_sum[d][j] / pc_cnt[d][j];
        return all_cnt[d] ? all_sum[d] / double(all_cnt[d]) : 1.0;
    };
    bool have_cur = true;
    End end = End::Complete;
    std::vector<double> lo, up;
    for (;;) {
        if (!have_cur) {
            if (open.empty()) break;
            cur = open.top();
            open.pop();
        }
        have_cur = false;
        if (has_inc && cur.bound >= z - tol_at(z)) {
            pruned = std::min(pruned, cur.bound);
            continue;
        }
        if (elapsed() > limit_ || nodes >= node_limit_) {
            end = elapsed() > limit_ ? End::TimeLimit : End::NodeLimit;
            open.push(std::move(cur));
            break;
        }
        ++nodes;
        lo = rlo_, up = rup_;
        for (const Change& c : cur.changes) lo[c.j] = c.lo, up[c.j] = c.up;
        // A child keeps its parent's optimal basis, which stays dual feasible under the tightened
        // bound: re-solve it with the dual simplex. The root and cold retries use the primal.
        Result r = cur.basis ? solve_lp_dual(md_, lo, up, cur.basis.get(), limit_ - elapsed())
                             : solve_lp(md_, lo, up, nullptr, limit_ - elapsed());
        iters += r.iterations;
        if (cur.basis && r.status != Status::Optimal && r.status != Status::Infeasible && r.status != Status::TimeLimit) {
            r = solve_lp(md_, lo, up, nullptr, limit_ - elapsed());  // cold retry
            iters += r.iterations;
        }
        if (r.status == Status::TimeLimit) {
            end = End::TimeLimit;
            open.push(std::move(cur));
            break;
        }
        if (r.status == Status::Infeasible) continue;
        if (r.status != Status::Optimal) {
            if (cur.changes.empty()) {  // root
                note = r.message;
                return r.status == Status::Unbounded ? End::RootUnbounded : End::RootFailed;
            }
            // Unbounded below a bounded root is a numerical artefact; never prune on it.
            ++unresolved;
            held.push_back(cur.bound);
            note = std::string("node LP ") + status_name(r.status) + (r.message.empty() ? "" : ": " + r.message);
            continue;
        }
        double b = std::max(cur.bound, r.objective);
        if (cur.dir) {  // learn from the parent -> child objective change
            int d = cur.dir > 0, j = cur.changes.back().j;
            double g = std::max(0.0, r.objective - cur.bound) / cur.dist;
            pc_sum[d][j] += g, ++pc_cnt[d][j], all_sum[d] += g, ++all_cnt[d];
        }
        offer(r.x);
        if (has_inc && b >= z - tol_at(z)) {
            pruned = std::min(pruned, b);
            continue;
        }
        // Pseudocost product score; uninitialised variables use the average. Before any
        // pseudocost exists every score is f(1-f), i.e. most-fractional branching.
        int bj = -1;
        double best = -1;
        for (int j = 0; j < n; ++j) {
            if (!md_.is_int[j]) continue;
            double f = r.x[j] - std::floor(r.x[j]);
            if (std::min(f, 1 - f) <= kIntTol) continue;
            double sc = std::max(pseudo(0, j) * f, 1e-6) * std::max(pseudo(1, j) * (1 - f), 1e-6);
            if (sc > best) bj = j, best = sc;
        }
        if (bj < 0) {  // integral LP point that failed the feasibility check
            ++unresolved;
            held.push_back(b);
            note = "integral node point rejected by the feasibility check";
            continue;
        }
        auto basis = std::make_shared<const std::vector<char>>(std::move(r.basis));
        double v = r.x[bj], fl = std::floor(v);
        Node down{b, cur.changes, basis, -1, v - fl}, upn{b, std::move(cur.changes), basis, 1, fl + 1 - v};
        down.changes.push_back({bj, lo[bj], fl});
        upn.changes.push_back({bj, fl + 1, up[bj]});
        bool dive_up = v - fl >= 0.5;
        open.push(dive_up ? std::move(down) : std::move(upn));
        cur = dive_up ? std::move(upn) : std::move(down);
        have_cur = true;
    }
    bound = std::min(pruned, has_inc ? z : kInf);
    if (have_cur) bound = std::min(bound, cur.bound);
    if (!open.empty()) bound = std::min(bound, open.top().bound);
    for (double h : held)
        if (!(has_inc && h >= z - tol_at(z))) bound = std::min(bound, h);
        else --unresolved;  // pruned by the final incumbent after all
    return end;
}

}  // namespace

MilpResult solve_milp(const Model& model, double time_limit_s, long node_limit) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    MilpResult res;
    double sign = model.maximize ? -1 : 1;
    Model md = model;  // minimisation form
    md.maximize = false;
    for (double& c : md.cost) c *= sign;
    md.obj_const *= sign;
    md.is_int.resize(md.cols.size(), 0);
    std::vector<double> lo = md.col_lo, up = md.col_up;
    for (size_t j = 0; j < lo.size(); ++j) {
        if (md.is_int[j]) lo[j] = std::ceil(lo[j] - kIntTol), up[j] = std::floor(up[j] + kIntTol);
        if (lo[j] > up[j]) {
            res.status = "infeasible";
            res.best_bound = sign * kInf;
            res.message = "empty integer bound range on column " + md.col_names[j];
            return res;
        }
    }

    Search s(md, lo, up, time_limit_s, node_limit);
    End end = s.run();
    res.nodes = s.nodes, res.lp_iterations = s.iters, res.unresolved_nodes = s.unresolved;
    res.message = s.note;
    if (end == End::RootUnbounded) {
        // For rational data an unbounded relaxation plus one integer-feasible point means the MILP
        // is unbounded (the integer hull has the relaxation's recession cone, Meyer 1974). Look for
        // such a point with a zero objective.
        Model feas = md;
        std::fill(feas.cost.begin(), feas.cost.end(), 0.0);
        Search f(feas, lo, up, time_limit_s - elapsed(), node_limit - s.nodes);
        End fe = f.run();
        res.nodes += f.nodes, res.lp_iterations += f.iters;
        if (f.has_inc && violation(model, f.inc) <= kFeasTol) {
            res.status = "unbounded";
            res.has_solution = true;
            res.x = f.inc;
            res.objective = objective_of(model, f.inc);
            res.best_bound = -sign * kInf;
            res.message = "LP relaxation unbounded and an integer-feasible point exists";
        } else if (fe == End::Complete && f.unresolved == 0) {
            res.status = "infeasible";
            res.best_bound = sign * kInf;
            res.message = "LP relaxation unbounded but no integer-feasible point";
        } else {
            res.status = "unbounded_relaxation";
            res.best_bound = -sign * kInf;
            res.message = "LP relaxation unbounded; integer feasibility undetermined";
        }
        return res;
    }
    if (end == End::RootFailed) {
        res.status = "numerical_failure";
        res.best_bound = -sign * kInf;
        res.message = "root LP failed: " + s.note;
        return res;
    }
    if (s.has_inc) {
        double v = violation(model, s.inc);
        if (v > kFeasTol) {  // cannot happen if offer() is consistent; never report such a point
            res.status = "numerical_failure";
            res.best_bound = sign * s.bound;
            res.message = "incumbent fails the original-model check by " + std::to_string(v);
            return res;
        }
        res.has_solution = true;
        res.x = s.inc;
        res.objective = objective_of(model, s.inc);
    }
    res.best_bound = sign * s.bound;
    if (s.has_inc) res.gap = std::abs(s.z - s.bound) / std::max(1.0, std::abs(s.z));
    bool proven = end == End::Complete && s.unresolved == 0;
    if (end == End::TimeLimit) res.status = "time_limit";
    else if (end == End::NodeLimit) res.status = "node_limit";
    else if (!proven) res.status = "numerical_failure";
    else res.status = s.has_inc ? "optimal" : "infeasible";
    return res;
}
