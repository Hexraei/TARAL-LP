// LP-based branch and bound for mixed-integer models. Best-bound node selection with a depth-first
// dive after every branching (the child on the rounding side of the branching variable), most-
// fractional branching, a rounding heuristic at every node, and warm-started node LPs.
// Every node first runs activity-based bound propagation on its integer columns; the tightenings
// are kept in the node's own change list, so siblings are unaffected.
// A node is pruned on a verified LP status (Optimal by bound, or Infeasible) or when propagation
// proves an INTEGER column's domain empty: ceil(lb - eps) > floor(ub + eps), eps >= 1e-6, with
// bounds that carry a safety margin on the row activities. Such a node can have a feasible LP
// relaxation (an integer bound pair like [0.3, 0.7]) and still hold no integer point, so the
// check that guards the prune is a plain branch and bound over the node (--audit-prop,
// Search::audit_pruned). MilpOptions::no_prop_prune hands those nodes to the LP instead. Any other
// LP outcome keeps the node as unresolved; its inherited bound stays in the reported best bound.
// With an incumbent, reduced-cost fixing tightens the bounds of nonbasic integer columns whose
// reduced cost proves that moving them further cannot beat the incumbent. Like propagation, these
// tightenings belong to the node's subtree only; a box they empty holds no improving point.
// The call refactorizes the basis, so it backs off (skips 1, 3, 7, then 15 later calls) while calls fix
// nothing; a new incumbent resets it. Skipping a fixing is always valid.
// A node whose LP fails (neither Optimal nor Infeasible, after the cold retry) has its widest integer
// column bisected instead of being left unresolved, at most 8 times in a row along one path; the halves
// cover the box and each gets its own LP, so nothing is pruned on the failed solve.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <queue>

#include "taral.hpp"

namespace {
constexpr double kIntTol = 1e-6;   // |x - round(x)| for an integral value
constexpr double kFeasTol = 1e-6;  // row and bound violation, relative to 1 + |bound|
static_assert(kIntTol == kFeasTol, "violation() compares integrality on the feasibility scale");
constexpr double kGapTol = 1e-6;   // optimal when objective - bound <= kGapTol * max(1, |objective|)
constexpr int kRootPasses = 25;    // propagation passes (a pass visits every pending row once)
constexpr int kNodePasses = 10;
constexpr double kPropEps = 1e-6;  // integer bounds round outward by this much (plus a scale term)
constexpr double kPropMax = 1e9;   // a derived bound beyond this is numerical noise: ignored
constexpr double kPropPivot = 1e-9;  // coefficients below this are not used to derive bounds
constexpr double kRcMin = 1e-7;      // never fix on a reduced cost with a smaller magnitude
constexpr double kRcDualTol = 1e-7;  // dual infeasibility (scaled by 1 + max|cost|) that voids a basis for fixing
constexpr int kMaxSplitDepth = 8;    // consecutive bisections of nodes whose LP failed, along one path
constexpr int kRcMaxBackoff = 4;     // after a call that fixes nothing, skip up to 2^4 - 1 later calls

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
    std::vector<Change> changes;  // bound tightenings relative to the root: branchings and propagation
    std::shared_ptr<const std::vector<char>> basis;  // parent's optimal basis (warm start)
    int dir = 0;      // last branching: -1 down, +1 up, 0 root
    double dist = 0;  // distance the branched variable was moved (fraction)
    int bj = -1;      // the branched variable
    int depth = 0;    // number of branchings from the root
    std::vector<int> seeds;  // columns whose bounds changed since the parent was propagated
    int splits = 0;  // bisections in a row on this path (the node's LP failed, so its integer box was halved)
};

struct Worse {  // priority_queue keeps the smallest bound on top; deeper first on ties
    bool operator()(const Node& a, const Node& b) const {
        return a.bound != b.bound ? a.bound > b.bound : a.depth < b.depth;
    }
};

enum class End { Complete, TimeLimit, NodeLimit, RootUnbounded, RootFailed };

// Search on a minimisation model. Integer column bounds are already rounded in root_lo/root_up.
class Search {
public:
    Search(const Model& md, std::vector<double> root_lo, std::vector<double> root_up, double time_limit, long node_limit)
        : md_(md), rlo_(std::move(root_lo)), rup_(std::move(root_up)), limit_(time_limit), node_limit_(node_limit) {
        int m = int(md.row_lo.size());
        rows_.assign(m, {});
        for (size_t j = 0; j < md.cols.size(); ++j)
            for (const Entry& e : md.cols[j]) rows_[e.index].push_back({int(j), e.value});
        queued_.assign(m, 0);
        mark_.assign(md.cols.size(), 0);
    }

    End run();
    bool has_inc = false;
    std::vector<double> inc;
    double z = kInf;      // incumbent objective
    double bound = -kInf;  // proven global lower bound after run()
    long nodes = 0, iters = 0, unresolved = 0;
    long prop_tightened = 0;   // integer bounds tightened by propagation
    long prop_crossed = 0;     // nodes whose propagated integer domain came out empty
    long prop_crossed_lp_infeasible = 0;  // ... of which the LP then confirmed infeasible
    long prop_pruned = 0;      // ... pruned without the LP (all of them unless --no-prop-prune)
    long rc_fixed = 0, rc_skipped = 0;  // reduced-cost tightenings; nodes whose basis was not trusted
    bool use_prop = true;      // false: plain branch and bound (the audit's oracle)
    bool use_rc = true;
    bool prop_prune = true;    // prune an emptied domain directly; false: the LP decides (--no-prop-prune)
    bool audit_prop = false;   // re-check every node whose propagated domain came out empty (debug/test)
    // Audit of those nodes: the LP on the box the node had before propagating, and for the ones the
    // LP finds feasible a plain branch and bound over that box.
    long audit_nodes = 0, audit_lp_infeasible = 0, audit_lp_feasible = 0, audit_lp_other = 0;
    long audit_int_empty = 0, audit_int_feasible = 0, audit_int_undecided = 0, audit_cutoff = 0;
    long audit_rc_checked = 0, audit_rc_bad = 0;  // reduced-cost fixings re-checked / found to cut off a better point
    std::string note;

private:
    const Model& md_;
    std::vector<double> rlo_, rup_;
    double limit_;
    long node_limit_;
    std::chrono::steady_clock::time_point t0_ = std::chrono::steady_clock::now();
    std::vector<std::vector<Entry>> rows_;  // row-wise A: (column, value)
    std::vector<char> queued_, mark_;       // propagation scratch: row pending, column changed
    std::vector<int> work_, next_, changed_;
    int rc_idle_ = 0;       // consecutive reduced-cost calls that fixed nothing (backoff exponent)
    long rc_skip_left_ = 0; // calls still to skip; both reset when the incumbent improves
    double elapsed() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count(); }

    bool propagate(std::vector<double>& lo, std::vector<double>& up, const std::vector<int>* seeds, int max_passes);
    void fix_by_reduced_cost(const Result& r, std::vector<double>& lo, std::vector<double>& up, std::vector<int>& fixed);
    void audit_fixing(std::vector<double> lo, std::vector<double> up, int j, double cut_lo, double cut_up);
    void audit_pruned(const std::vector<double>& lo, const std::vector<double>& up);

    // Offer an LP point: integer columns rounded first, then as is. Returns true if accepted.
    bool offer(const std::vector<double>& x) {
        std::vector<double> r = x;
        for (size_t j = 0; j < r.size(); ++j)
            if (md_.is_int[j]) r[j] = std::round(r[j]);
        for (const std::vector<double>* c : {static_cast<const std::vector<double>*>(&r), &x}) {
            if (violation(md_, *c) > kFeasTol) continue;
            double v = objective_of(md_, *c);
            if (!has_inc || v < z) has_inc = true, z = v, inc = *c, rc_idle_ = 0, rc_skip_left_ = 0;
            return true;
        }
        return false;
    }
};

// Activity-based bound propagation on the integer columns. From each row's activity range (smallest
// and largest value of a'x over the current bounds) every integer column gets the bound the row's
// sides force on it, rounded outward by a small tolerance so no feasible point is ever cut off.
// Only integer bounds change; continuous columns enter the activities but are never tightened.
// seeds == nullptr starts from every row (root), otherwise from the rows of those columns. The columns whose
// bounds changed are left in changed_. Returns false only when an integer column's domain is empty,
// ceil(lb - eps) > floor(ub + eps) with eps >= 1e-6 (the eps also grows with the row's magnitude, so a
// borderline case is never pruned). Crossed bounds of a continuous column cannot occur here: they are
// never tightened, and a row with no integer column is left to the LP.
bool Search::propagate(std::vector<double>& lo, std::vector<double>& up, const std::vector<int>* seeds, int max_passes) {
    changed_.clear();
    work_.clear(), next_.clear();
    auto push = [&](std::vector<int>& to, int i) {
        if (!queued_[i]) queued_[i] = 1, to.push_back(i);
    };
    auto drop = [&] {
        for (int i : work_) queued_[i] = 0;
        for (int i : next_) queued_[i] = 0;
        for (int j : changed_) mark_[j] = 0;
    };
    if (!seeds) {
        for (int i = 0; i < int(rows_.size()); ++i) push(work_, i);
    } else {
        for (int j : *seeds)
            for (const Entry& e : md_.cols[j]) push(work_, e.index);
    }
    for (int pass = 0; pass < max_passes && !work_.empty(); ++pass) {
        for (int i : work_) {
            queued_[i] = 0;
            const double rl = md_.row_lo[i], ru = md_.row_up[i];
            const bool hl = std::isfinite(rl), hu = std::isfinite(ru);
            if (!hl && !hu) continue;
            double minf = 0, maxf = 0;  // finite parts of the smallest and largest activity
            int mini = 0, maxi = 0;     // terms that make them infinite
            for (const Entry& e : rows_[i]) {
                double a = e.value, l = lo[e.index], u = up[e.index];
                double lowv = a > 0 ? l : u, highv = a > 0 ? u : l;
                if (std::isfinite(lowv)) minf += a * lowv; else ++mini;
                if (std::isfinite(highv)) maxf += a * highv; else ++maxi;
            }
            const bool can_up = hu && mini <= 1, can_lo = hl && maxi <= 1;
            if (!can_up && !can_lo) continue;
            for (const Entry& e : rows_[i]) {
                int j = e.index;
                double a = e.value;
                if (!md_.is_int[j] || std::abs(a) < kPropPivot) continue;
                double l = lo[j], u = up[j], nl = l, nu = u;
                if (can_up) {  // a x_j <= ru - (smallest activity of the other terms)
                    double c = a > 0 ? a * l : a * u;
                    bool inf = !std::isfinite(c);
                    if (mini - inf == 0) {
                        double v = (ru - (minf - (inf ? 0.0 : c))) / a;
                        double eps = kPropEps + 1e-9 * (1 + std::abs(ru) + std::abs(minf)) / std::abs(a);
                        if (a > 0) nu = std::min(nu, std::floor(v + eps));
                        else nl = std::max(nl, std::ceil(v - eps));
                    }
                }
                if (can_lo) {  // a x_j >= rl - (largest activity of the other terms)
                    double c = a > 0 ? a * u : a * l;
                    bool inf = !std::isfinite(c);
                    if (maxi - inf == 0) {
                        double v = (rl - (maxf - (inf ? 0.0 : c))) / a;
                        double eps = kPropEps + 1e-9 * (1 + std::abs(rl) + std::abs(maxf)) / std::abs(a);
                        if (a > 0) nl = std::max(nl, std::ceil(v - eps));
                        else nu = std::min(nu, std::floor(v + eps));
                    }
                }
                if (std::abs(nl) > kPropMax) nl = l;
                if (std::abs(nu) > kPropMax) nu = u;
                if (nl > nu) {
                    drop();
                    return false;
                }
                if (nl <= l && nu >= u) continue;
                lo[j] = nl, up[j] = nu;
                if (!mark_[j]) mark_[j] = 1, changed_.push_back(j);
                for (const Entry& r : md_.cols[j])
                    if (r.index != i) push(next_, r.index);
            }
        }
        work_.swap(next_);
        next_.clear();
    }
    drop();
    return true;
}

// Audit of one node whose propagated integer domain came out empty, given the box it had before
// propagating. (1) the LP relaxation of the box: Infeasible agrees with the prune; Optimal only
// means the emptiness came from integer rounding, which an LP cannot see, so it is not a
// disagreement. (2) For every pruned node, LP-infeasible ones included, a plain branch and bound
// (no propagation) over the same box with a zero objective: an integer-feasible point there would
// make the prune wrong (audit_int_feasible,
// must stay 0) unless an incumbent exists and the box was cut by reduced-cost fixing, in which case
// the prune only has to leave no point better than the incumbent: a second plain search with the
// true objective and the incumbent as cutoff must come back empty (audit_cutoff counts those).
void Search::audit_pruned(const std::vector<double>& lo, const std::vector<double>& up) {
    ++audit_nodes;
    Result r = solve_lp(md_, lo, up, nullptr, 60);
    const bool lp_infeasible = r.status == Status::Infeasible;
    if (lp_infeasible) {
        ++audit_lp_infeasible;
    } else if (r.status != Status::Optimal) {
        ++audit_lp_other;
        return;
    } else {
        ++audit_lp_feasible;
    }
    // Every pruned node gets the nested search, including those the LP calls infeasible (there it
    // must find nothing at its root); only the LP-feasible ones are counted as int_empty.
    Model feas = md_;
    std::fill(feas.cost.begin(), feas.cost.end(), 0.0);
    Search sub(feas, lo, up, 60, 2'000'000);
    sub.use_prop = false, sub.use_rc = false;
    End e = sub.run();
    if (lp_infeasible && sub.has_inc) {
        ++audit_int_feasible;
        std::fprintf(stderr, "AUDIT FAILURE: a node pruned by propagation (LP infeasible) holds an integer-feasible point\n");
        return;
    }
    if (lp_infeasible) return;
    if (!sub.has_inc && e == End::Complete && sub.unresolved == 0) {
        ++audit_int_empty;
        return;
    }
    if (!sub.has_inc || !has_inc) {
        if (sub.has_inc) {
            ++audit_int_feasible;
            std::fprintf(stderr, "AUDIT FAILURE: a node pruned by propagation holds an integer-feasible point\n");
        } else {
            ++audit_int_undecided;
        }
        return;
    }
    Search imp(md_, lo, up, 60, 2'000'000);  // integer points exist: none may beat the incumbent
    imp.use_prop = false, imp.use_rc = false;
    imp.has_inc = true, imp.z = z;
    End ie = imp.run();
    if (!imp.inc.empty()) {
        ++audit_int_feasible;
        std::fprintf(stderr, "AUDIT FAILURE: a node pruned by propagation holds a point better than the incumbent\n");
    } else if (ie == End::Complete && imp.unresolved == 0) {
        ++audit_cutoff;
    } else {
        ++audit_int_undecided;
    }
}

// Audit of one reduced-cost fixing: the part of the node's box it removes (column j restricted to
// [cut_lo, cut_up]) is searched by a plain branch and bound with the true objective and the incumbent
// as cutoff. A point better than the incumbent there means the fixing was wrong.
void Search::audit_fixing(std::vector<double> lo, std::vector<double> up, int j, double cut_lo, double cut_up) {
    ++audit_rc_checked;
    lo[j] = cut_lo, up[j] = cut_up;
    Search imp(md_, lo, up, 60, 2'000'000);
    imp.use_prop = false, imp.use_rc = false;
    imp.has_inc = true, imp.z = z;
    imp.run();
    if (!imp.inc.empty()) {
        ++audit_rc_bad;
        std::fprintf(stderr, "AUDIT FAILURE: a reduced-cost fixing cuts off a point better than the incumbent\n");
    }
}

// Reduced-cost fixing. At the optimal basis of a node LP with value zlp, every feasible point x of the
// node has c'x = zlp + sum over nonbasic j of d_j (x_j - xhat_j), each term >= 0 when the basis is dual
// feasible. A solution better than the incumbent z therefore moves a nonbasic integer column j at
// its lower (upper) bound by at most floor((z - zlp) / |d_j|); the bound on the far side is tightened
// to that. The duals come from a fresh factorization of the basis (y = B^-T c_B, d = c - A'y; the
// logical of row i is the column -e_i). The slack in the gap and the 1e-7 floor on |d_j| keep
// numerical noise from cutting off a point that could still win; if any nonbasic reduced cost has the
// wrong sign by more than a tolerance the basis is not trusted and nothing is fixed. The columns
// tightened are appended to `fixed`.
void Search::fix_by_reduced_cost(const Result& r, std::vector<double>& lo, std::vector<double>& up,
                                 std::vector<int>& fixed) {
    const int n = int(md_.cols.size()), m = int(md_.row_lo.size());
    const std::vector<char>& bs = r.basis;
    if (int(bs.size()) != n + m) return;
    std::vector<int> head;
    for (int j = 0; j < n + m; ++j)
        if (bs[j] == 0) head.push_back(j);
    if (int(head.size()) != m) return;
    std::vector<std::vector<Entry>> bcols(m);
    std::vector<double> cb(m), y;
    for (int p = 0; p < m; ++p) {
        int j = head[p];
        if (j < n) {
            bcols[p] = md_.cols[j], cb[p] = md_.cost[j];
        } else {
            bcols[p].push_back({j - n, -1.0}), cb[p] = 0;
        }
    }
    SparseLU lu;
    if (!lu.factor(m, bcols)) {
        ++rc_skipped;
        return;
    }
    lu.btran(cb, y);
    double cmax = 0;
    for (double c : md_.cost) cmax = std::max(cmax, std::abs(c));
    const double dtol = kRcDualTol * (1 + cmax);
    std::vector<double> d(n, 0.0);
    for (int j = 0; j < n + m; ++j) {
        if (bs[j] == 0) continue;
        double dj;
        if (j < n) {
            dj = md_.cost[j];
            for (const Entry& e : md_.cols[j]) dj -= e.value * y[e.index];
            d[j] = dj;
        } else {
            dj = y[j - n];
        }
        double l = j < n ? lo[j] : md_.row_lo[j - n], u = j < n ? up[j] : md_.row_up[j - n];
        if (l == u) continue;  // pinned: the sign of its reduced cost does not matter
        bool bad = bs[j] == 1 ? dj < -dtol : bs[j] == 2 ? dj > dtol : std::abs(dj) > dtol;
        if (bad) {
            ++rc_skipped;
            return;
        }
    }
    const double gap = z - r.objective + tol_at(z);
    for (int j = 0; j < n; ++j) {
        if (!md_.is_int[j] || bs[j] == 0 || bs[j] == 3 || lo[j] == up[j] || std::abs(d[j]) < kRcMin) continue;
        double steps = std::floor(gap / std::abs(d[j]) + kPropEps);
        if (steps > kPropMax) continue;
        if (bs[j] == 1 && d[j] > 0 && lo[j] + steps < up[j]) {
            if (audit_prop) audit_fixing(lo, up, j, lo[j] + steps + 1, up[j]);
            up[j] = lo[j] + steps;
            fixed.push_back(j);
            ++rc_fixed;
        } else if (bs[j] == 2 && d[j] < 0 && up[j] - steps > lo[j]) {
            if (audit_prop) audit_fixing(lo, up, j, lo[j], up[j] - steps - 1);
            lo[j] = up[j] - steps;
            fixed.push_back(j);
            ++rc_fixed;
        }
    }
}

End Search::run() {
    std::priority_queue<Node, std::vector<Node>, Worse> open;
    std::vector<double> held;  // bounds of nodes that could not be resolved
    double pruned = kInf;      // smallest bound among nodes pruned by bound
    int n = int(md_.cols.size());
    Node cur{-kInf, {}, nullptr, 0, 0, -1, 0, {}};
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
        const bool is_root = nodes == 1;
        lo = rlo_, up = rup_;
        for (const Change& c : cur.changes) lo[c.j] = c.lo, up[c.j] = c.up;
        // Propagate this node's domain; the root starts from every row, a child from the rows of
        // its branching variable. The tightenings join the node's own change list (children copy
        // it, siblings do not see them).
        bool crossed = false;
        if (use_prop) {
            if (propagate(lo, up, is_root ? nullptr : &cur.seeds, is_root ? kRootPasses : kNodePasses)) {
                for (int j : changed_) cur.changes.push_back({j, lo[j], up[j]});
                prop_tightened += long(changed_.size());
            } else {
                // Empty integer domain. Restore the box the node had; the audit looks at it, then
                // the node is pruned, or (no_prop_prune) the LP gets to say so.
                ++prop_crossed;
                crossed = true;
                lo = rlo_, up = rup_;
                for (const Change& c : cur.changes) lo[c.j] = c.lo, up[c.j] = c.up;
                if (audit_prop) audit_pruned(lo, up);
                if (prop_prune) {
                    ++prop_pruned;
                    continue;
                }
            }
        }
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
        if (r.status == Status::Infeasible) {
            prop_crossed_lp_infeasible += crossed;
            continue;
        }
        if (r.status != Status::Optimal) {
            if (is_root) {
                note = r.message;
                return r.status == Status::Unbounded ? End::RootUnbounded : End::RootFailed;
            }
            // Unbounded below a bounded root is a numerical artefact; never prune on it. The node's
            // box is split instead: the two halves cover it exactly, each half gets an LP of its own
            // (nothing is pruned without a certified LP status), and the halves inherit this bound.
            int sj = -1;
            double range = 0;
            for (int j = 0; j < n && cur.splits < kMaxSplitDepth; ++j)
                if (md_.is_int[j] && std::isfinite(lo[j]) && std::isfinite(up[j]) && up[j] - lo[j] > range)
                    sj = j, range = up[j] - lo[j];
            if (sj >= 0) {
                double mid = std::floor((lo[sj] + up[sj]) / 2);
                Node a{cur.bound, cur.changes, cur.basis, 0, 0, -1, cur.depth + 1, {sj}, cur.splits + 1}, b = a;
                a.changes.push_back({sj, lo[sj], mid});
                b.changes.push_back({sj, mid + 1, up[sj]});
                open.push(std::move(a));
                open.push(std::move(b));
                continue;
            }
            ++unresolved;
            held.push_back(cur.bound);
            note = std::string("node LP ") + status_name(r.status) + (r.message.empty() ? "" : ": " + r.message);
            continue;
        }
        double b = std::max(cur.bound, r.objective);
        if (cur.dir) {  // learn from the parent -> child objective change
            int d = cur.dir > 0, j = cur.bj;
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
        std::vector<int> seeds{bj};  // the children propagate from the branching and from any fixing
        if (has_inc && use_rc) {
            // Each call refactorizes the basis. When calls keep fixing nothing, back off exponentially;
            // skipping a fixing is always valid, it only forgoes a tightening.
            if (rc_skip_left_ > 0) {
                --rc_skip_left_;
            } else {
                size_t at = seeds.size();
                fix_by_reduced_cost(r, lo, up, seeds);
                for (size_t t = at; t < seeds.size(); ++t) cur.changes.push_back({seeds[t], lo[seeds[t]], up[seeds[t]]});
                if (seeds.size() > at) rc_idle_ = 0;
                else rc_idle_ = std::min(rc_idle_ + 1, kRcMaxBackoff);
                rc_skip_left_ = (1L << rc_idle_) - 1;
            }
        }
        auto basis = std::make_shared<const std::vector<char>>(std::move(r.basis));
        double v = r.x[bj], fl = std::floor(v);
        Node down{b, cur.changes, basis, -1, v - fl, bj, cur.depth + 1, seeds},
            upn{b, std::move(cur.changes), basis, 1, fl + 1 - v, bj, cur.depth + 1, seeds};
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

MilpResult solve_milp(const Model& model, double time_limit_s, long node_limit, MilpOptions opt) {
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
    s.audit_prop = opt.audit_prop, s.prop_prune = !opt.no_prop_prune;
    End end = s.run();
    res.nodes = s.nodes, res.lp_iterations = s.iters, res.unresolved_nodes = s.unresolved;
    res.prop_tightened = s.prop_tightened, res.prop_crossed = s.prop_crossed, res.prop_pruned = s.prop_pruned;
    res.prop_crossed_lp_infeasible = s.prop_crossed_lp_infeasible;
    res.rc_fixed = s.rc_fixed, res.rc_skipped = s.rc_skipped;
    res.audit = {s.audit_nodes, s.audit_lp_infeasible, s.audit_lp_feasible, s.audit_lp_other,
                 s.audit_int_empty, s.audit_int_feasible, s.audit_int_undecided, s.audit_cutoff,
                 s.audit_rc_checked, s.audit_rc_bad};
    res.message = s.note;
    if (end == End::RootUnbounded) {
        // For rational data an unbounded relaxation plus one integer-feasible point means the MILP
        // is unbounded (the integer hull has the relaxation's recession cone, Meyer 1974). Look for
        // such a point with a zero objective.
        Model feas = md;
        std::fill(feas.cost.begin(), feas.cost.end(), 0.0);
        Search f(feas, lo, up, time_limit_s - elapsed(), node_limit - s.nodes);
        f.prop_prune = s.prop_prune, f.use_rc = false;  // the feasibility search has no objective to cut on
        End fe = f.run();
        res.nodes += f.nodes, res.lp_iterations += f.iters;
        res.prop_tightened += f.prop_tightened, res.prop_crossed += f.prop_crossed, res.prop_pruned += f.prop_pruned;
        res.prop_crossed_lp_infeasible += f.prop_crossed_lp_infeasible;
        res.rc_fixed += f.rc_fixed, res.rc_skipped += f.rc_skipped;
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
