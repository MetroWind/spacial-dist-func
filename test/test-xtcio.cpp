// Copyright 2020 MetroWind <chris.corsair@gmail.com>
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see
// <https://www.gnu.org/licenses/>.

#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "utils.h"
#include "xtcio.h"
#include "trajectory.h"

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

TEST_CASE("Trajectory")
{
    libmd::Trajectory t;
    t.open("../test/test.xtc", "../test/test.gro");

    REQUIRE(t.nextFrame());

    CHECK(t.vec(0).isApprox(Eigen::Vector3f(4.249, 2.67, 4.389)));
    CHECK(t.vec("18+BCDEF").isApprox(Eigen::Vector3f(4.145, 2.535, 4.553)));
    CHECK(t.vec("17+H11").isApprox(Eigen::Vector3f(4.191, 2.707, 4.304)));

    REQUIRE(t.nextFrame());
    CHECK(t.vec("17+H11").isApprox(Eigen::Vector3f(4.209, 2.698, 4.304)));

    REQUIRE(t.nextFrame());
    CHECK(t.vec(9).isApprox(Eigen::Vector3f(4.007, 2.533, 4.688)));

    REQUIRE_FALSE(t.nextFrame());
    t.close();
}

TEST_CASE("Trajectory filter")
{
    libmd::Trajectory t;
    t.open("../test/test.xtc", "../test/test.gro");
    t.nextFrame();
    t.close();

    auto Snap = libmd::filterFrame(
        t,
        [](const libmd::AtomIdentifier& name, size_t _1,
           const typename libmd::V3Map& _2)
        {
            UNUSED(_1); UNUSED(_2);
            return name.toStr() == "18+BCDEF";
        });
    CHECK(Snap.meta().AtomCount == 1);
    CHECK(Snap.data().size() == 3);
    CHECK(Snap.vec("18+BCDEF").isApprox(Eigen::Vector3f(4.145, 2.535, 4.553)));
}
