cmake_minimum_required(VERSION 3.17)
project(ProvaRtreeV2)

set(CMAKE_CXX_STANDARD 14)

add_executable(main main.cpp)

find_package(NLopt CONFIG REQUIRED)
target_link_libraries(main PRIVATE NLopt::nlopt)