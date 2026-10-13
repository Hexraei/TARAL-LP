// Root-node c-MIR / flow-cover style cuts, v2. Single mechanism, root only, opt-in (TARAL_MIR=1).
// Cuts are derived ONLY from original model rows (separate_mir takes the original row count; cut rows are never sources),
// with no tableau or basis access. Every cut is   sum val[k] * x[idx[k]] <= rhs   and is valid for the integer hull.
//
// Arithmetic: all derived quantities use outward-rounded interval arithmetic (src/mir_cuts.cpp). Supported
// assumptions, checked at run time (fp_assumptions_ok): IEEE-754 binary64, round-to-nearest-even, no extended
// precision. Any non-finite or oversized value rejects the cut (fail closed).
#pragma once
#include "taral.hpp"
#include <chrono>

struct MirCut {
    std::vector<int> idx;
    std::vector<double> val;
    double rhs = 0;
    double efficacy = 0;
    int source_row = -1;
    int side = 0; // 0: row <= row_up, 1: -row <= -row_lo
};

struct MirStats {
    long rows_tried = 0, rows_skipped = 0, cuts_found = 0, cuts_rejected_numerics = 0, cuts_rejected_domain = 0,
         vub_used = 0, timeouts = 0, agg_tried = 0, agg_cuts = 0;
    bool fp_ok = true;
};

bool fp_assumptions_ok();

// Minimisation model, integer-rounded root column bounds lo/up, LP point x. Only rows [0, nsrc_rows) are sources.
// Returns cuts violated by x with efficacy >= min_eff, best first, at most max_cuts, no two nearly parallel.
// Stops (returning what it has) at `deadline`; stats->timeouts counts that.
std::vector<MirCut> separate_mir(const Model& md, const std::vector<double>& lo, const std::vector<double>& up,
                                 const std::vector<double>& x, int nsrc_rows, int max_cuts, double min_eff,
                                 std::chrono::steady_clock::time_point deadline, MirStats* stats, int max_agg = 0);

enum class CutCheck { Holds, Violated, Inconclusive };
// Exact (dyadic) evaluation of sum val*x against rhs + tol*(1+|rhs|). Inconclusive when the exact evaluation fails.
CutCheck check_cut_exact(const MirCut& c, const std::vector<double>& x, double tol);
