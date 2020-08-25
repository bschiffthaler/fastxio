cmake_minimum_required(VERSION 3.0)

project(TASites LANGUAGES CXX)

add_executable(TASites src/main.cpp)

set_property(TARGET TASites PROPERTY CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(OpenMP COMPONENTS CXX REQUIRED)

target_include_directories(TASites PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(TASites PUBLIC "../../src") 
target_link_libraries(TASites ${Boost_LIBRARIES} OpenMP::OpenMP_CXX
  fastxio)