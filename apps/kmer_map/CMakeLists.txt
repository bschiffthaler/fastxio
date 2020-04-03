cmake_minimum_required(VERSION 3.0)

project(KmerMap LANGUAGES CXX)

add_executable(KmerMap src/main.cpp)

set_property(TARGET KmerMap PROPERTY CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)

target_include_directories(KmerMap PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(KmerMap PUBLIC "../../src") 
target_link_libraries(KmerMap ${Boost_LIBRARIES} fastxio)