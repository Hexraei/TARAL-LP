#pragma once
#include <string>
#include <vector>
#include "taral.hpp"
namespace miqp_proto {
struct MiqpResult {
    std::string status = "numerical_failure";
    double objective = 0, best_bound = 0;
    long nodes = 0, unresolved = 0;
    bool has_solution = false;
    std::vector<double> x;
    std::string message;
};
MiqpResult solve_miqp(const Model& md, double time_limit, long node_limit);
}
