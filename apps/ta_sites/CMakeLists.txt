cmake_minimum_required(VERSION 3.0)

project(TASites LANGUAGES CXX)

if (FASTX_STATIC_APPS)
  set(Boost_USE_STATIC_LIBS ON)
  set(FXLIB fastxioS)
  set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-static ${CMAKE_EXE_LINKER_FLAGS}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
else()
  set(FXLIB fastxio)
endif()

add_executable(TASites src/main.cpp)

set_property(TARGET TASites PROPERTY CXX_STANDARD 11)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(OpenMP COMPONENTS CXX REQUIRED)

target_include_directories(TASites PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(TASites PUBLIC "../../src") 
target_link_libraries(TASites -static ${Boost_LIBRARIES} OpenMP::OpenMP_CXX
  ${FXLIB})