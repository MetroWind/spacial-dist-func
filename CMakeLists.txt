cmake_minimum_required (VERSION 3.1)

project(sdf LANGUAGES CXX VERSION 1.1.0)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic-errors -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add march=native?
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  message("Compiler option march=native enabled.")
endif()

# Are we running in stupid Ubuntu?
if(STUPID_UBUNTU)
  add_compile_definitions(STUPID_UBUNTU)
  message(WARNING "Ubuntu is stupid. You should consider a distribution with rolling release, like Gentoo and Arch.")
endif()

find_package(Threads REQUIRED)
find_package(Eigen3 3.3 QUIET NO_MODULE)
find_package(Catch2 QUIET)
if(NOT STUPID_UBUNTU)
  find_package(pugixml REQUIRED)
endif()

IF(NOT EIGEN3_INCLUDE_DIR)
  set(LEGACY_EIGEN ON)
  message(WARNING "Eigen >= 3.3 not found. I'll maybe still compile, but you should update.")
  find_package(Eigen3 REQUIRED NO_MODULE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIRS}")
ELSE()
  set(LEGACY_EIGEN OFF)
ENDIF()

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
target_link_libraries(sdf PUBLIC Threads::Threads)
IF(LEGACY_EIGEN)
  target_link_libraries(sdf INTERFACE Eigen3::Eigen)
ELSE()
  target_link_libraries(sdf PUBLIC Eigen3::Eigen)
ENDIF()
if(NOT STUPID_UBUNTU)
  target_link_libraries(sdf PUBLIC pugixml)
endif()

set(TestFiles
  test/test-xtcio.cpp
  test/test-pbc.cpp
  test/test-sdf.cpp
  test/testutils.h
  )

if(NOT STUPID_UBUNTU)
  IF(Catch2_FOUND)
    add_executable(mdtest test/main.cpp ${TestFiles} ${TheFiles})
    target_link_libraries(mdtest PUBLIC Threads::Threads)
    IF(LEGACY_EIGEN)
      target_link_libraries(mdtest INTERFACE Eigen3::Eigen)
    ELSE()
      target_link_libraries(mdtest PUBLIC Eigen3::Eigen)
    ENDIF()
    target_link_libraries(mdtest PUBLIC Catch2::Catch2)
    target_link_libraries(mdtest PUBLIC pugixml)

    target_include_directories(mdtest PRIVATE src)
    include(CTest)
    include(Catch)
    catch_discover_tests(mdtest)
  ELSE()
    message("Catch2 not found.")
  ENDIF()
endif()
