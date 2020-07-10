cmake_minimum_required (VERSION 3.1)

project(sdf LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# Add march=native?
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  message("Compiler option march=native enabled.")
endif()

project(sdf)

find_package(Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Catch2 QUIET)

set(TheFiles
  src/xtcio.h
  src/xtcio.cpp
  src/trajectory.h
  src/trajectory.cpp
  src/pbc.h
  src/pbc.cpp
  src/utils.h
  src/utils.cpp
  src/sdf.h
  src/sdf.cpp
  src/config.h
  src/config.cpp
  )

add_executable(sdf src/main-sdf.cpp ${TheFiles})
target_link_libraries(sdf Eigen3::Eigen ${CMAKE_THREAD_LIBS_INIT})

set(TestFiles
  test/test-xtcio.cpp
  test/test-pbc.cpp
  test/test-sdf.cpp
  test/testutils.h
  )

IF(Catch2_FOUND)
  add_executable(mdtest test/main.cpp ${TestFiles} ${TheFiles})
  target_link_libraries(mdtest Eigen3::Eigen Catch2::Catch2 ${CMAKE_THREAD_LIBS_INIT})
  target_include_directories(mdtest PRIVATE src)
  include(CTest)
  include(Catch)
  catch_discover_tests(mdtest)
ELSE()
  message(WARNING "Catch2 not found.")
ENDIF()
