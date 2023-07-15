cmake_minimum_required(VERSION 3.17)
project(ProvaRtreeV2)

set(CMAKE_CXX_STANDARD 14)

add_executable(prova main.cpp)

find_package(NLopt CONFIG REQUIRED)
target_link_libraries(prova PRIVATE NLopt::nlopt)