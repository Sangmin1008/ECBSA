cmake_minimum_required(VERSION 3.27)
project(ECBSA C)

set(CMAKE_C_STANDARD 11)
include_directories(C:/msys64/gmp-6.3.0/include)
link_directories(C:/msys64/gmp-6.3.0/lib)

find_library(GMP_LIBRARY NAMES gmp gmpxx PATHS C:/msys64/gmp-6.3.0/lib)

if(NOT GMP_LIBRARY)
    message(FATAL_ERROR "GMP library not found!")
endif()


add_executable(ECBSA main.c
        elliptic_curve.c
        elliptic_curve.h
        ecc.c
        ecc.h
        ecbsa.h
        ecbsa.c
)
target_link_libraries(ECBSA gmp)