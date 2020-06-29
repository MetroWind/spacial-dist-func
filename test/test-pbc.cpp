#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "pbc.h"

using v3 = Eigen::Vector3f;

TEST_CASE("PBC")
{
    libmd::RectPbc3d Pbc(4, 8, 10);
    {
        v3 a(0, 0, 9);
        v3 b(0, 0, 11);
        CHECK(Pbc.dist(a, b) == 2.0f);
    }
    {
        v3 a(0, 0, 0);
        v3 b(0, 0, 11);
        CHECK(Pbc.dist(a, b) == 1.0f);
    }
    {
        v3 a(1, 0, 0);
        v3 b(2, 0, 11);
        CHECK(Pbc.dist(a, b) == Approx(std::sqrt(2.0f)));
    }
    {
        v3 a(0, 0, 0);
        v3 b(0, 0, 0);
        CHECK(Pbc.dist(a, b) == 0.0f);
    }
    {
        v3 a(0, 2, 0);
        v3 b(0, 6, 0);
        CHECK(Pbc.dist(a, b) == 4.0f);
    }

}
