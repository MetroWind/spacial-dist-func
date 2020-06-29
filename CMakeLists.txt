cmake_minimum_required (VERSION 3.1)

project(sdf LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(sdf)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Catch2)

set(TheFiles
  src/xtcio.h
  src/xtcio.cpp
  src/trojectory.h
  src/trojectory.cpp
  src/pbc.h
  src/pbc.cpp
  )

add_executable(not_test src/test.cpp ${TheFiles})
target_link_libraries(not_test Eigen3::Eigen)

set(TestFiles
  test/test-xtcio.cpp
  test/test-pbc.cpp
  )

IF(Catch2_FOUND)
  add_executable(mdtest test/main.cpp ${TestFiles} ${TheFiles})
  target_link_libraries(mdtest Eigen3::Eigen Catch2::Catch2)
  target_include_directories(mdtest PRIVATE src)
  include(CTest)
  include(Catch)
  catch_discover_tests(mdtest)
ENDIF()
# add_test(NAME mdtest COMMAND mdtest)
