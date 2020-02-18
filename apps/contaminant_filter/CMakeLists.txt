cmake_minimum_required(VERSION 3.0)

project(ContaminantFilter LANGUAGES CXX)

add_executable(ContaminantFilter src/main.cpp)

set_property(TARGET ContaminantFilter PROPERTY CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(OpenMP COMPONENTS CXX REQUIRED)

target_include_directories(ContaminantFilter PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(ContaminantFilter PUBLIC "../../src") 
target_link_libraries(ContaminantFilter ${Boost_LIBRARIES} OpenMP::OpenMP_CXX
  fastxio)