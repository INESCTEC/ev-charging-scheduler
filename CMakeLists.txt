# CMake minimum version required
cmake_minimum_required(VERSION 3.10)

# Project name
project(MyProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(HIGHS_DIR highs/build/lib/cmake/highs)

find_package(HIGHS REQUIRED)



# Add the executable
add_executable(charging_optimization charging_optimization.cpp)


# Link the HiGHS library
target_link_libraries(charging_optimization highs::highs)