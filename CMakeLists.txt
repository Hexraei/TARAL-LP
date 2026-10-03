cmake_minimum_required(VERSION 3.16)
project(TARAL_LP LANGUAGES CXX)

add_executable(taral
    src/crossover.cpp
    src/dual.cpp
    src/explain.cpp
    src/ipm.cpp
    src/kkt_gate.cpp
    src/lu.cpp
    src/main.cpp
    src/milp.cpp
    src/mps.cpp
    src/nonoptimal_certificate.cpp
    src/presolve.cpp
    src/simplex.cpp
)
target_compile_features(taral PRIVATE cxx_std_17)
set_target_properties(taral PROPERTIES CXX_EXTENSIONS OFF)
