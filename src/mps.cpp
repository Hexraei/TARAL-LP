// MPS reader (fixed and free format). Tries whitespace tokens first and falls back to
// fixed columns when a line does not parse, mirroring the independent checker's rules.
// Supports OBJSENSE, RANGES, the usual bound types, integer MARKER blocks and quadratic
// objective sections (QUADOBJ, QMATRIX, QSECTION on the objective row). Only the first
// RHS / RANGES / BOUNDS set is used; later named sets are ignored, as the format specifies.
// A BOUNDS column that never appears in COLUMNS is read as an empty (zero-column) variable.
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

double number(const std::string& s) {
    std::string t = s;
    for (char& ch : t)
        if (ch == 'D' || ch == 'd') ch = 'E';
    char* end = nullptr;
    double v = std::strtod(t.c_str(), &end);
    if (t.empty() || *end) throw ParseError("bad number '" + s + "'");
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
        std::istringstream in(line);
        for (std::string w; in >> w;) t.push_back(w);
        return t;
    }
    auto push = [&](size_t a, size_t b) { t.push_back(field(line, a, b)); };
    if (section == "ROWS") {
        push(1, 3), push(4, 12);
    } else if (section == "BOUNDS") {
        push(1, 3), push(4, 12), push(14, 22);
        if (!field(line, 24, 36).empty()) push(24, 36);
    } else {  // COLUMNS, RHS, RANGES
        push(4, 12), push(14, 22), push(24, 36);
        if (!field(line, 39, 47).empty()) push(39, 47), push(49, 61);
    }
    return t;
}

Model parse(const std::string& path, bool fixed, bool allow_bound_only_cols) {
    std::ifstream in(path);
    if (!in) throw ParseError("cannot open " + path);
    Model md;
    std::string section, obj_row, line;
    std::unordered_map<std::string, int> row_id, col_id;
    std::unordered_set<std::string> free_rows;  // extra N rows are dropped
    std::vector<char> row_type;
    std::vector<double> rhs;
    std::vector<double> range;
    std::vector<char> has_range;
    std::vector<std::map<int, double>> colmap;
    bool in_int = false;
    std::string rhs_set, range_set, bound_set;  // first named set of each kind
    std::map<std::pair<int, int>, double> q;      // (row >= col) -> Q value
    bool q_full = false;                          // QMATRIX/QSECTION list both triangles

    auto first_set = [](std::string& keep, const std::string& name) {  // true if name is the set in use
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

    while (std::getline(in, line)) {
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
                break;
            } else if (section == "OBJSENSE") {
                if (!arg.empty()) {
                    if (arg == "MAX" || arg == "MAXIMIZE") md.maximize = true;
                    else if (arg != "MIN" && arg != "MINIMIZE") throw ParseError("bad OBJSENSE " + arg);
                }
            } else if (section == "QUADOBJ" || section == "QMATRIX") {
                q_full = section == "QMATRIX";
                section = "QUAD";
            } else if (section == "QSECTION" || section == "QCMATRIX") {
                if (section == "QCMATRIX" || arg != obj_row)
                    throw ParseError("quadratic constraints are not supported (" + section + " " + arg + ")");
                q_full = true;
                section = "QUAD";
            } else if (section != "ROWS" && section != "COLUMNS" && section != "RHS" &&
                       section != "RANGES" && section != "BOUNDS") {
                throw ParseError("unsupported section " + section);
            }
            continue;
        }
        std::vector<std::string> t = tokens(line, section, fixed);
        if (section == "OBJSENSE") {
            std::string v = t.empty() ? "" : t[0];
            if (v == "MAX" || v == "MAXIMIZE") md.maximize = true;
            else if (v != "MIN" && v != "MINIMIZE") throw ParseError("bad OBJSENSE " + v);
            continue;
        }
        if (section == "ROWS") {
            // Exactly two fields; more means a name with spaces, which only the fixed-column reading handles.
            if (t.size() != 2 || t[1].empty()) throw ParseError("bad ROWS line");
            char type = t[0].empty() ? '?' : static_cast<char>(std::toupper(t[0][0]));
            if (type == 'N') {
                if (obj_row.empty()) obj_row = t[1];
                else free_rows.insert(t[1]);
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
            auto it = col_id.find(t[0]);
            int j;
            if (it == col_id.end()) {
                j = static_cast<int>(md.col_names.size());
                col_id[t[0]] = j;
                md.col_names.push_back(t[0]);
                md.cost.push_back(0);
                md.is_int.push_back(in_int);
                colmap.emplace_back();
            } else {
                j = it->second;
            }
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                int r = row_of(t[k]);
                double v = number(t[k + 1]);
                if (r == -1) md.cost[j] += v;
                else if (r >= 0) colmap[j][r] += v;
            }
        } else if (section == "RHS" || section == "RANGES") {
            size_t start = t.size() % 2 == 1 ? 1 : 0;  // optional set name
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
                    if (r == -1) md.obj_const = -v;
                    else if (r >= 0) rhs[r] = v;
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
                if (no_value) set = t[1], col = t[2];
                else col = t[1], val = t[2];
            } else if (t.size() == 2) {
                col = t[1];
            } else {
                throw ParseError("bad BOUNDS line");
            }
            if (!set.empty() && !first_set(bound_set, set)) continue;
            auto it = col_id.find(col);
            if (it == col_id.end()) {
                // Zero-column convention: a BOUNDS column absent from COLUMNS is an empty column
                // (zero cost, no coefficients), appended after the COLUMNS columns in order of first
                // appearance. Only accepted once strict readings have failed (see read_mps), so a
                // misparsed line in the wrong format still falls through to the other reading.
                if (!allow_bound_only_cols || col.empty()) throw ParseError("unknown column '" + col + "' in BOUNDS");
                int nj = static_cast<int>(md.col_names.size());
                col_id[col] = nj;
                md.col_names.push_back(col);
                md.cost.push_back(0);
                md.is_int.push_back(0);
                colmap.emplace_back();
                md.col_lo.push_back(0);
                md.col_up.push_back(kInf);
                it = col_id.find(col);
            }
            int j = it->second;
            double v = no_value || val.empty() ? 0 : number(val);
            double &lo = md.col_lo[j], &up = md.col_up[j];
            if (type == "BV" || type == "LI" || type == "UI") md.is_int[j] = 1;
            if (type == "LO" || type == "LI") lo = v;
            else if (type == "UP" || type == "UI") {
                up = v;
                if (v < 0 && lo == 0) lo = -kInf;
            } else if (type == "FX") lo = up = v;
            else if (type == "FR") lo = -kInf, up = kInf;
            else if (type == "MI") lo = -kInf;
            else if (type == "PL") up = kInf;
            else if (type == "BV") lo = 0, up = 1;
            else throw ParseError("unsupported bound type " + type);
        } else if (section == "QUAD") {
            if (t.size() != 3) throw ParseError("bad quadratic objective line");
            int a = col_of(t[0]), b = col_of(t[1]);
            double v = number(t[2]);
            if (q_full && a < b) continue;  // full listing: keep one triangle
            q[{std::max(a, b), std::min(a, b)}] += v;
        }
    }
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
    md.cols.resize(n);
    for (size_t j = 0; j < n; ++j)
        for (auto [r, v] : colmap[j])
            if (v != 0) md.cols[j].push_back({r, v});
    return md;
}

}  // namespace

bool Model::has_integers() const {
    for (char c : is_int)
        if (c) return true;
    return false;
}

Model read_mps(const std::string& path) {
    // Strict readings first (free, then fixed); only if both fail, retry allowing BOUNDS columns
    // that never appear in COLUMNS.
    try {
        return parse(path, false, false);
    } catch (const ParseError&) {
    }
    try {
        return parse(path, true, false);
    } catch (const ParseError&) {
    }
    try {
        return parse(path, false, true);
    } catch (const ParseError&) {
        return parse(path, true, true);
    }
}
