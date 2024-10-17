cmake_minimum_required(VERSION 3.13.4)

project(test_project)

find_package(deal.II 9.6.0 REQUIRED
HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)
deal_ii_initialize_cached_variables()

add_executable(test test.cpp)
deal_ii_setup_target(test)
