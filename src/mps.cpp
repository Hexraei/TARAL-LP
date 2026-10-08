// MPS reader (fixed and free format). Tries whitespace tokens first and falls back to
// fixed columns when a line does not parse, mirroring the independent checker's rules.
// Supports OBJSENSE, RANGES, the usual bound types, integer MARKER blocks and quadratic
// objective sections (QUADOBJ, QMATRIX, QSECTION on the objective row). Only the first
// RHS / RANGES / BOUNDS set is used; later named sets are ignored, as the format specifies.
// Malformed input is an error, not a guess: a missing ENDATA, a repeated ROWS or COLUMNS section, a duplicate
// row name, a column listed in two separate blocks and non-decimal numbers (hex, inf, nan) are all rejected.
// Names must be declared before they are used: a BOUNDS, RHS or RANGES entry for a column or row that COLUMNS
// or ROWS never declared is rejected (the format defines them over declared names; HiGHS silently ignores
// the entry, which would hide a typo in the model).
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "taral.hpp"

namespace {

constexpr double kHugeBound = 1e20;

// Every numeric field (coefficient, RHS, range, bound, quadratic term) must be a plain decimal,
// [+-] digits [. digits] [(e|E|d|D) [+-] digits], and finite: "nan", "inf" and literals that overflow (1e400)
// are parse errors, not values that turn into a NaN objective later. The grammar is spelled out because
// strtod alone would also take hex floats ("0x10" reads as 16), which the format does not define.
double number(const std::string& s) {
    size_t i = 0, n = s.size();
    auto digits = [&] {
        size_t a = i;
        while (i < n && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
        return i - a;
    };
    if (i < n && (s[i] == '+' || s[i] == '-')) ++i;
    size_t mantissa = digits();
    if (i < n && s[i] == '.') ++i, mantissa += digits();
    bool ok = mantissa > 0;
    if (ok && i < n && (s[i] == 'e' || s[i] == 'E' || s[i] == 'd' || s[i] == 'D')) {
        ++i;
        if (i < n && (s[i] == '+' || s[i] == '-')) ++i;
        ok = digits() > 0;
    }
    if (!ok || i != n) throw ParseError("bad number '" + s + "'");
    double v;
    if (s.find_first_of("Dd") == std::string::npos) {
        v = std::strtod(s.c_str(), nullptr);
    } else {
        std::string t = s;
        for (char& ch : t)
            if (ch == 'D' || ch == 'd') ch = 'E';
        v = std::strtod(t.c_str(), nullptr);
    }
    if (!std::isfinite(v)) throw ParseError("number '" + s + "' is not finite");
    return v;
}

std::string field(const std::string& line, size_t a, size_t b) {
    if (line.size() <= a) return "";
    std::string s = line.substr(a, b - a);
    size_t i = s.find_first_not_of(" \t"), j = s.find_last_not_of(" \t\r");
    return i == std::string::npos ? "" : s.substr(i, j - i + 1);
}

std::vector<std::string> tokens(const std::string& line, const std::string& section, bool fixed) {
    std::vector<std::string> t;
    if (!fixed) {
        // Same split as `istream >> string` in the C locale: runs of space, \t, \n, \v, \f, \r separate tokens.
        const char* p = line.data();
        const char* e = p + line.size();
        auto ws = [](char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'; };
        while (p < e) {
            while (p < e && ws(*p)) ++p;
            const char* s = p;
            while (p < e && !ws(*p)) ++p;
            if (p > s) t.emplace_back(s, p - s);
        }
        return t;
    }
    auto push = [&](size_t a, size_t b) { t.push_back(field(line, a, b)); };
    if (section == "ROWS") {
        push(1, 3), push(4, 12);
    } else if (section == "BOUNDS") {
        push(1, 3), push(4, 12), push(14, 22);
        if (!field(line, 24, 36).empty()) push(24, 36);
    } else { // COLUMNS, RHS, RANGES
        push(4, 12), push(14, 22), push(24, 36);
        if (!field(line, 39, 47).empty()) push(39, 47), push(49, 61);
    }
    return t;
}

Model parse(const std::string& path, bool fixed, std::chrono::steady_clock::time_point deadline) {
    std::ifstream in(path);
    if (!in) throw ParseError("cannot open " + path);
    Model md;
    std::string section, obj_row, line;
    std::unordered_map<std::string, int> row_id, col_id;
    std::unordered_set<std::string> free_rows; // extra N rows are dropped
    std::vector<char> row_type;
    std::vector<double> rhs;
    std::vector<double> range;
    std::vector<char> has_range;
    bool in_int = false;
    bool saw_endata = false;
    int cur_col = -1;               // column whose entries are being read
    std::vector<Entry> col_entries; // entries of the column being read, in file order
    // Turns the pending entries of column cur_col into its sorted row list: duplicate rows are summed in file
    // order (as the former per-column std::map did) and zero sums are dropped.
    auto flush_column = [&] {
        if (cur_col < 0 || col_entries.empty()) return;
        std::stable_sort(col_entries.begin(), col_entries.end(), [](const Entry& x, const Entry& y) { return x.index < y.index; });
        std::vector<Entry>& out = md.cols[cur_col];
        for (const Entry& e : col_entries) {
            if (!out.empty() && out.back().index == e.index)
                out.back().value += e.value;
            else
                out.push_back(e);
        }
        out.erase(std::remove_if(out.begin(), out.end(), [](const Entry& e) { return e.value == 0; }), out.end());
        col_entries.clear();
    };
    std::unordered_set<std::string> seen_sections; // ROWS and COLUMNS may each appear once
    std::string rhs_set, range_set, bound_set;     // first named set of each kind
    std::map<std::pair<int, int>, double> q;       // (row >= col) -> Q value
    bool q_full = false;                           // QMATRIX/QSECTION list both triangles

    auto first_set = [](std::string& keep, const std::string& name) { // true if name is the set in use
        if (keep.empty()) keep = name;
        return keep == name;
    };
    auto col_of = [&](const std::string& name) {
        auto it = col_id.find(name);
        if (it == col_id.end()) throw ParseError("unknown column '" + name + "'");
        return it->second;
    };

    auto row_of = [&](const std::string& name) -> int {
        auto it = row_id.find(name);
        if (it != row_id.end()) return it->second;
        if (name == obj_row) return -1;
        if (free_rows.count(name)) return -2;
        throw ParseError("unknown row '" + name + "'");
    };

    unsigned long lines_read = 0;
    auto check_deadline = [&] {
        if (std::chrono::steady_clock::now() > deadline) throw ParseTimeLimit("parse exceeded the time limit");
    };
    while (std::getline(in, line)) {
        if ((++lines_read & 4095) == 0) check_deadline();
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.find_first_not_of(" \t") == std::string::npos || line[0] == '*') continue;
        if (line[0] != ' ' && line[0] != '\t') {
            std::istringstream hs(line);
            hs >> section;
            std::string arg;
            hs >> arg;
            if (section == "NAME") {
                md.name = arg;
            } else if (section == "ENDATA") {
                saw_endata = true;
                break;
            } else if (section == "OBJSENSE") {
                if (!arg.empty()) {
                    if (arg == "MAX" || arg == "MAXIMIZE")
                        md.maximize = true;
                    else if (arg != "MIN" && arg != "MINIMIZE")
                        throw ParseError("bad OBJSENSE " + arg);
                }
            } else if (section == "QUADOBJ" || section == "QMATRIX") {
                q_full = section == "QMATRIX";
                section = "QUAD";
            } else if (section == "QSECTION" || section == "QCMATRIX") {
                if (section == "QCMATRIX" || arg != obj_row)
                    throw ParseError("quadratic constraints are not supported (" + section + " " + arg + ")");
                q_full = true;
                section = "QUAD";
            } else if (section != "ROWS" && section != "COLUMNS" && section != "RHS" && section != "RANGES" && section != "BOUNDS") {
                throw ParseError("unsupported section " + section);
            }
            if ((section == "ROWS" || section == "COLUMNS") && !seen_sections.insert(section).second)
                throw ParseError("section " + section + " appears twice");
            continue;
        }
        std::vector<std::string> t = tokens(line, section, fixed);
        if (section == "OBJSENSE") {
            std::string v = t.empty() ? "" : t[0];
            if (v == "MAX" || v == "MAXIMIZE")
                md.maximize = true;
            else if (v != "MIN" && v != "MINIMIZE")
                throw ParseError("bad OBJSENSE " + v);
            continue;
        }
        if (section == "ROWS") {
            // Exactly two fields; more means a name with spaces, which only the fixed-column reading handles.
            if (t.size() != 2 || t[1].empty()) throw ParseError("bad ROWS line");
            char type = t[0].empty() ? '?' : static_cast<char>(std::toupper(t[0][0]));
            if (row_id.count(t[1]) || t[1] == obj_row || free_rows.count(t[1])) throw ParseError("duplicate row name '" + t[1] + "'");
            if (type == 'N') {
                if (obj_row.empty())
                    obj_row = t[1];
                else
                    free_rows.insert(t[1]);
            } else if (type == 'E' || type == 'L' || type == 'G') {
                row_id[t[1]] = static_cast<int>(md.row_names.size());
                md.row_names.push_back(t[1]);
                row_type.push_back(type);
            } else {
                throw ParseError("bad row type " + t[0]);
            }
        } else if (section == "COLUMNS") {
            if (t.size() >= 2 && t[1] == "'MARKER'") {
                bool start = false, end = false;
                for (const std::string& w : t) start |= w == "'INTORG'", end |= w == "'INTEND'";
                if (start == end) throw ParseError("bad MARKER line");
                in_int = start;
                continue;
            }
            if (t.size() != 3 && t.size() != 5) throw ParseError("bad COLUMNS line");
            int j;
            if (cur_col >= 0 && t[0] == md.col_names[cur_col]) { // entries of one column are consecutive: skip the lookup
                j = cur_col;
            } else {
                auto it = col_id.find(t[0]);
                if (it == col_id.end()) {
                    flush_column();
                    j = static_cast<int>(md.col_names.size());
                    col_id[t[0]] = j;
                    md.col_names.push_back(t[0]);
                    md.cost.push_back(0);
                    md.is_int.push_back(in_int);
                    md.cols.emplace_back();
                } else {
                    j = it->second;
                    if (j != cur_col) throw ParseError("column '" + t[0] + "' is listed in two separate blocks");
                }
            }
            cur_col = j;
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                int r = row_of(t[k]);
                double v = number(t[k + 1]);
                if (r == -1)
                    md.cost[j] += v;
                else if (r >= 0)
                    col_entries.push_back({r, v});
            }
        } else if (section == "RHS" || section == "RANGES") {
            size_t start = t.size() % 2 == 1 ? 1 : 0; // optional set name
            if (fixed) start = 1;
            if (start == 1 && !t[0].empty() && !first_set(section == "RHS" ? rhs_set : range_set, t[0])) continue;
            if (rhs.empty()) {
                rhs.assign(md.row_names.size(), 0);
                range.assign(md.row_names.size(), 0);
                has_range.assign(md.row_names.size(), 0);
            }
            for (size_t k = start; k + 1 < t.size(); k += 2) {
                int r = row_of(t[k]);
                double v = number(t[k + 1]);
                if (section == "RHS") {
                    if (r == -1)
                        md.obj_const = -v;
                    else if (r >= 0)
                        rhs[r] = v;
                } else if (r >= 0) {
                    range[r] = v;
                    has_range[r] = 1;
                }
            }
        } else if (section == "BOUNDS") {
            if (md.col_lo.empty()) {
                md.col_lo.assign(md.col_names.size(), 0);
                md.col_up.assign(md.col_names.size(), kInf);
            }
            std::string type = t.empty() ? "" : t[0];
            bool no_value = type == "FR" || type == "MI" || type == "PL" || type == "BV";
            std::string set, col, val;
            if (fixed) {
                set = t.size() > 1 ? t[1] : "";
                col = t.size() > 2 ? t[2] : "";
                val = t.size() > 3 ? t[3] : "";
            } else if (t.size() == 4) {
                set = t[1], col = t[2], val = t[3];
            } else if (t.size() == 3) {
                if (no_value)
                    set = t[1], col = t[2];
                else
                    col = t[1], val = t[2];
            } else if (t.size() == 2) {
                col = t[1];
            } else {
                throw ParseError("bad BOUNDS line");
            }
            if (!set.empty() && !first_set(bound_set, set)) continue;
            auto it = col_id.find(col);
            if (it == col_id.end()) throw ParseError("unknown column '" + col + "' in BOUNDS");
            int j = it->second;
            double v = no_value || val.empty() ? 0 : number(val);
            double &lo = md.col_lo[j], &up = md.col_up[j];
            if (type == "BV" || type == "LI" || type == "UI") md.is_int[j] = 1;
            if (type == "LO" || type == "LI")
                lo = v;
            else if (type == "UP" || type == "UI") {
                up = v;
                if (v < 0 && lo == 0) lo = -kInf;
            } else if (type == "FX")
                lo = up = v;
            else if (type == "FR")
                lo = -kInf, up = kInf;
            else if (type == "MI")
                lo = -kInf;
            else if (type == "PL")
                up = kInf;
            else if (type == "BV")
                lo = 0, up = 1;
            else
                throw ParseError("unsupported bound type " + type);
        } else if (section == "QUAD") {
            if (t.size() != 3) throw ParseError("bad quadratic objective line");
            int a = col_of(t[0]), b = col_of(t[1]);
            double v = number(t[2]);
            if (q_full && a < b) continue; // full listing: keep one triangle
            q[{std::max(a, b), std::min(a, b)}] += v;
        }
    }
    flush_column();
    check_deadline();
    if (!saw_endata) throw ParseError("missing ENDATA");
    for (const auto& [rc, v] : q)
        if (v != 0) md.qobj.push_back({rc.first, rc.second, v});
    if (obj_row.empty()) throw ParseError("no objective row");

    size_t m = md.row_names.size(), n = md.col_names.size();
    if (rhs.empty()) {
        rhs.assign(m, 0);
        range.assign(m, 0);
        has_range.assign(m, 0);
    }
    if (md.col_lo.empty()) {
        md.col_lo.assign(n, 0);
        md.col_up.assign(n, kInf);
    }
    md.row_lo.resize(m);
    md.row_up.resize(m);
    for (size_t i = 0; i < m; ++i) {
        double b = rhs[i], R = range[i];
        double &lo = md.row_lo[i], &up = md.row_up[i];
        if (row_type[i] == 'E') {
            lo = up = b;
            if (has_range[i]) (R >= 0 ? up : lo) = b + R;
        } else if (row_type[i] == 'L') {
            up = b;
            lo = has_range[i] ? b - std::abs(R) : -kInf;
        } else {
            lo = b;
            up = has_range[i] ? b + std::abs(R) : kInf;
        }
    }
    // |bound| >= 1e20 means infinity (the 1e20 / 1e30 sentinels written by CPLEX, Gurobi, HiGHS and others).
    // Only the direction a sentinel can stand for is mapped: a lower bound <= -1e20 and an upper bound >= 1e20.
    // The solvers branch on isinf(), so a sentinel left finite would be pivoted to as a real bound.
    for (std::vector<double>* lo : {&md.col_lo, &md.row_lo})
        for (double& v : *lo)
            if (v <= -kHugeBound) {
                if (std::isfinite(v)) ++md.sentinel_bounds;
                v = -kInf;
            }
    for (std::vector<double>* up : {&md.col_up, &md.row_up})
        for (double& v : *up)
            if (v >= kHugeBound) {
                if (std::isfinite(v)) ++md.sentinel_bounds;
                v = kInf;
            }
    return md;
}

} // namespace

bool Model::has_integers() const {
    for (char c : is_int)
        if (c) return true;
    return false;
}

Model read_mps(const std::string& path, std::chrono::steady_clock::time_point deadline) {
    try {
        return parse(path, false, deadline);
    } catch (const ParseError& free_err) {
        try {
            return parse(path, true, deadline);
        } catch (const ParseError& fixed_err) { // report both readings, not just the fallback's
            if (std::string(free_err.what()) == fixed_err.what()) throw;
            throw ParseError(std::string(free_err.what()) + " (free format); " + fixed_err.what() + " (fixed columns)");
        }
    }
}
