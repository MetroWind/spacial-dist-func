#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "testutils.h"
#include "pbc.h"

using v3 = Eigen::Vector3f;

TEST_CASE("PBC dist")
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

TEST_CASE("PBC wrap")
{
    libmd::RectPbc3d Pbc(4, 8, 10);
    {
        const v3 a(0, 0, 9);
        v3 b(0, 0, 11);
        v3 c = b;
        Pbc.wrapVec(a, b);
        CHECK(b == c);
    }
    {
        v3 a(0, 0, 9);
        const v3 b(0, 0, 11);
        v3 c = a;
        Pbc.wrapVec(b, a);
        CHECK(a == c);
    }
    {
        const v3 a(0, 0, 1);
        v3 b(0, 0, 11);
        Pbc.wrapVec(a, b);
        CHECK(a == b);
    }
    {
        v3 a(0, 0, 1);
        const v3 b(0, 0, 11);
        Pbc.wrapVec(b, a);
        CHECK(a == b);
    }

    // Self-wrap
    {
        v3 a = randVec({-100, 100}, {-100, 100}, {-100, 100});
        v3 b = a;
        Pbc.wrapVec(a, a);
        CHECK(a == b);
    }
}

TEST_CASE("PBC wrap forward")
{
    libmd::RectPbc3d Pbc(4, 8, 10);
    const v3 a(0, 0, -4);
    v3 b(0, 0, -12);
    Pbc.wrapVec(a, b);
    CHECK(b == v3(0, 0, -2));
}

TEST_CASE("PBC wrap random vec")
{
    libmd::RectPbc3d Pbc(4, 8, 10);
    for(int i = 0; i < 100; i++)
    {
        v3 a = randVec({-100, 100}, {-100, 100}, {-100, 100});
        v3 b = randVec({-100, 100}, {-100, 100}, {-100, 100});
        float Dist = Pbc.dist(a, b);
        // std::cout << "\nBefore wrap: \n" << a << "\n\n" << b << std::endl;
        Pbc.wrapVec(a, b);
        // std::cout << "\nAfter wrap: \n" << a << "\n\n" << b << std::endl;
        CHECK((a-b).norm() == Approx(Dist).margin(0.0001));
        CHECK(Dist <= Pbc.DiagLength * 0.5);
    }
}
