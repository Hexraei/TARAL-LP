cmake_minimum_required(VERSION 3.16)
project(TARAL_LP LANGUAGES CXX)

add_executable(taral
    src/main.cpp src/mps.cpp src/simplex.cpp src/dual.cpp
    src/lu.cpp src/milp.cpp src/ipm.cpp)
target_compile_features(taral PRIVATE cxx_std_17)
set_target_properties(taral PROPERTIES CXX_EXTENSIONS OFF)
