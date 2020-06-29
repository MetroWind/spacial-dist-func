#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "xtcio.h"
#include "trojectory.h"

TEST_CASE("XTC reading")
{
    libmd::XtcFile f;
    f.open("../test/test.xtc");

    auto MetaHead = f.readFrameMeta();
    REQUIRE(MetaHead.AtomCount == 10);
    REQUIRE(MetaHead.Step == 1000000);
    REQUIRE(MetaHead.Time == Approx(1000.0));
    std::vector<float> data(MetaHead.AtomCount * 3, 0.0f);

    f.readFrame(data.data());
    CHECK(data[0] == Approx(4.249f));
    CHECK(data[1] == Approx(2.67f));
    CHECK(data[2] == Approx(4.389f));

    auto Meta = f.readFrameMeta();
    REQUIRE(Meta.AtomCount == 10);
    REQUIRE(Meta.Step == 1000020);
    REQUIRE(Meta.Time == Approx(1000.02));
    f.readFrame(data.data());
    CHECK(data[0] == Approx(4.26f));
    CHECK(data[1] == Approx(2.669f));
    CHECK(data[2] == Approx(4.396f));

    f.readFrame(data.data());
    CHECK(data[0] == Approx(4.269f));
    CHECK(data[1] == Approx(2.675f));
    CHECK(data[2] == Approx(4.405f));

    REQUIRE(f.eof());
    f.close();
}

TEST_CASE("Trojectory")
{
    libmd::Trojectory t;
    t.open("../test/test.xtc", "../test/test.gro");

    REQUIRE(t.nextFrame());
    CHECK(t.vec(0).isApprox(Eigen::Vector3f(4.249, 2.67, 4.389)));
    CHECK(t.vec("BCDEF").isApprox(Eigen::Vector3f(4.145, 2.535, 4.553)));
    CHECK(t.vec("H11").isApprox(Eigen::Vector3f(4.191, 2.707, 4.304)));

    REQUIRE(t.nextFrame());
    CHECK(t.vec("H11").isApprox(Eigen::Vector3f(4.209, 2.698, 4.304)));

    REQUIRE(t.nextFrame());
    CHECK(t.vec(9).isApprox(Eigen::Vector3f(4.007, 2.533, 4.688)));

    REQUIRE_FALSE(t.nextFrame());
    t.close();
}
