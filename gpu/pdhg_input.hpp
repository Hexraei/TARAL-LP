// Host-only input gate for the LP-only PDHG prototype (both CPU and GPU backends).
#pragma once
#include <cstdio>

#include "../src/taral.hpp"

inline const char* pdhg_input_error(const Model& model) {
    const bool quadratic = !model.qobj.empty();
    const bool integer = model.has_integers();
    if (quadratic && integer)
        return "PDHG supports continuous LPs only; quadratic objectives and integer variables are unsupported";
    if (quadratic)
        return "PDHG supports continuous LPs only; quadratic objectives are unsupported";
    if (integer)
        return "PDHG supports continuous LPs only; integer variables are unsupported (no LP relaxation is solved)";
    return nullptr;
}

// Call before scaling, device initialization, or writing solution files.
inline bool pdhg_reject_unsupported_input(const Model& model) {
    const char* error = pdhg_input_error(model);
    if (!error) return false;
    std::printf("status unsupported_model\n");
    std::fprintf(stderr, "%s\n", error);
    return true;
}
