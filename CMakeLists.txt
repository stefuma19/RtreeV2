cmake_minimum_required(VERSION 3.17)
project(ProvaRtreeV2)

set(CMAKE_CXX_STANDARD 14)

add_executable(prova main.cpp)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(" ")
    message("CMAKE IN DEBUG MODE")
    message(" ")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message(" ")
    message("CMAKE IN RELEASE MODE")
    message(" ")
endif ()

find_package(NLopt CONFIG REQUIRED)
target_link_libraries(prova PRIVATE NLopt::nlopt)