cmake_minimum_required(VERSION 3.27)
project(ECBSA C)

set(CMAKE_C_STANDARD 11)

add_executable(ECBSA main.c
        elliptic_curve.c
        elliptic_curve.h
        ecc.c
        ecc.h
)
