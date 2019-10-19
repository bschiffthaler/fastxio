cmake_minimum_required(VERSION 3.0)

project(ExtractRegRegions LANGUAGES CXX)

add_executable(ExtractRegRegions src/main.cpp)

set_property(TARGET ExtractRegRegions PROPERTY CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)

target_include_directories(ExtractRegRegions PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(ExtractRegRegions PUBLIC "../../src") 
target_link_libraries(ExtractRegRegions ${Boost_LIBRARIES} fastxio)