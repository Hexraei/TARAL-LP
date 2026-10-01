// CLI per the engine contract:  taral MODEL.mps --time-limit S --sol OUT.sol --json OUT.json
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>

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
}  // namespace

int main(int argc, char** argv) {
    const char *model = nullptr, *sol = nullptr, *json = nullptr;
    double limit = 60;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--time-limit") && i + 1 < argc) limit = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--sol") && i + 1 < argc) sol = argv[++i];
        else if (!std::strcmp(argv[i], "--json") && i + 1 < argc) json = argv[++i];
        else model = argv[i];
    }
    if (!model) {
        std::fprintf(stderr, "usage: taral MODEL.mps [--time-limit S] [--sol OUT.sol] [--json OUT.json]\n");
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
        return 0;
    }
    Result r = solve(md, limit - wall());
    double w = wall();
    if (json) write_json(json, status_name(r.status), &r, w, r.message);
    if (sol && r.status == Status::Optimal) {
        if (std::FILE* f = std::fopen(sol, "w")) {
            for (size_t j = 0; j < md.col_names.size(); ++j) std::fprintf(f, "%s %.17g\n", md.col_names[j].c_str(), r.x[j]);
            std::fclose(f);
        }
    }
    std::printf("status %s objective %.12g iterations %ld wall %.3fs %s\n", status_name(r.status), r.objective,
                r.iterations, w, r.message.c_str());
    return 0;
}
