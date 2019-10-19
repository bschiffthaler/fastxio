cmake_minimum_required(VERSION 2.8)
project (fastx_io)

if(NO_CHECKING)
   set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DNO_ERROR_CHECKING=1" )
endif()

if (SANITIZE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
endif()

if(DEBUG_SYM)
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB SOURCES src/*.cpp)
file(GLOB GZSTREAM gzstream/*.cpp)
file(GLOB BZ2STREAM bz2stream/*.cpp)

file(GLOB HEADERS src/*.h)
file(GLOB EXTHEADERS generic_matrix/*.h)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/bz2stream
  ${CMAKE_CURRENT_SOURCE_DIR}/gzstream
  ${CMAKE_CURRENT_SOURCE_DIR}/generic_matrix
  ${CMAKE_CURRENT_SOURCE_DIR}/external/smhasher/src
  ${CMAKE_CURRENT_SOURCE_DIR}/external/libbs/src
  )

add_subdirectory(external/libbs)

add_library(fastxio SHARED ${SOURCES} ${GZSTREAM} ${BZ2STREAM}
            ${CMAKE_CURRENT_SOURCE_DIR}/external/smhasher/src/MurmurHash3.cpp)
add_library(fastxioS STATIC ${SOURCES} ${GZSTREAM} ${BZ2STREAM}
            ${CMAKE_CURRENT_SOURCE_DIR}/external/smhasher/src/MurmurHash3.cpp)

set_target_properties(fastxioS PROPERTIES OUTPUT_NAME fastxio)

set_property(TARGET fastxio PROPERTY CXX_STANDARD 11)
set_property(TARGET fastxioS PROPERTY CXX_STANDARD 11)

target_link_libraries(fastxio z bz2 bs)
target_link_libraries(fastxioS z bz2 bs)

###
## Apps
###

add_subdirectory(apps/extr_prom)

if (PERFTOOLS)
  target_link_libraries(fastxio profiler)
endif()

install(TARGETS fastxio fastxioS
        DESTINATION lib)

install(FILES ${HEADERS} ${EXTHEADERS}
	DESTINATION include/fastxio)
