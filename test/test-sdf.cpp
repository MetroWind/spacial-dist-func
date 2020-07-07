#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "testutils.h"
#include "sdf.h"

using v3 = Eigen::Vector3f;

TEST_CASE("Nearest filter")
{
    libmd::Trojectory t;
    t.open("../test/test.xtc", "../test/test.gro");
    t.nextFrame();
    t.close();

    auto Filtered = findNearest("BCDEF", 0.15, t);
    REQUIRE(Filtered.meta().AtomCount == 3);
    CHECK(Filtered.vec("C65").isApprox(v3(4.146, 2.661, 4.503)));
    CHECK(Filtered.vec("C66").isApprox(v3(4.037, 2.473, 4.607)));
    CHECK(Filtered.vec("BCDEF").isApprox(v3(4.145, 2.535, 4.553)));
}
