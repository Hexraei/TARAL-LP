// Minimal C++17 LP interface example. Link solver sources except src/main.cpp.
#include "taral.hpp"
#include <cmath>
#include <iostream>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        Model model = read_mps(argv[1]);
        if (model.has_integers() || !model.qobj.empty()) return 2; // this example is LP-only
        Result result = solve(model, 30.0);
        std::cout << status_name(result.status) << ' ' << result.objective << '\n';
        return result.status == Status::Optimal ? 0 : 1;
    } catch (const ParseError& e) {
        std::cerr << e.what() << '\n';
        return 3;
    }
}
