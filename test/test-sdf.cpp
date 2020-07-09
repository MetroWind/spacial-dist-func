#include <sstream>

#include <catch2/catch.hpp>
#include <Eigen/Dense>

#include "testutils.h"
#include "config.h"
#include "sdf.h"

using v3 = Eigen::Vector3f;

// TEST_CASE("Nearest filter")
// {
//     libmd::Trojectory t;
//     t.open("../test/test.xtc", "../test/test.gro");
//     t.nextFrame();
//     t.close();

//     auto Filtered = findNearest("BCDEF", 0.15, t);
//     REQUIRE(Filtered.meta().AtomCount == 3);
//     CHECK(Filtered.vec("C65").isApprox(v3(4.146, 2.661, 4.503)));
//     CHECK(Filtered.vec("C66").isApprox(v3(4.037, 2.473, 4.607)));
//     CHECK(Filtered.vec("BCDEF").isApprox(v3(4.145, 2.535, 4.553)));
// }

TEST_CASE("Rotation")
{
    for(int i = 0; i < 10; i++)
    {
        v3 v1 = randVec({-10, 10}, {-10, 10}, {-10, 10});
        v3 v2 = randVec({-10, 10}, {-10, 10}, {-10, 10});
        auto Rot = libmd::rotateToAlignX(v1, v2);
        v3 v1Rot = Rot * v1;
        v3 v2Rot = Rot * v2;

        CHECK(v1Rot.isApprox(v3(v1.norm(), 0.0, 0.0)));
        CHECK(v2Rot.norm() == Approx(v2.norm()).margin(0.0001));
        CHECK(v2Rot[2] == Approx(0.0f).margin(0.0001));
        CHECK(v2Rot.dot(v1Rot) == Approx(v1.dot(v2)).margin(0.0001));
        CHECK(v2Rot[1] >= 0.0);
    }
}

TEST_CASE("Config reading 1 run")
{
    std::stringstream ss;
    ss << "test.xtc\n"
       << "test.gro\n"
       << "+++\n"
       << "18+BCDEF\n"
       << "17+O2\n"
       << "17+C65\n"
       << "100\n"
       << "101\n"
       << "17+H10\n"
       << "17+H11\n";

    auto Config = sdf::RuntimeConfig::read(ss);
    CHECK(Config.XtcFile == "test.xtc");
    CHECK(Config.GroFile == "test.gro");
    CHECK(Config.Params.size() == 1);
    CHECK(Config.Params[0].Anchor.toStr() == "18+BCDEF");
    CHECK(Config.Params[0].AtomX.toStr() == "17+O2");
    CHECK(Config.Params[0].AtomXY.toStr() == "17+C65");
    CHECK(Config.Params[0].Distance == 100.0);
    CHECK(Config.Params[0].SliceThickness == 101.0);
    CHECK(Config.Params[0].OtherAtoms.size() == 2);
    CHECK(Config.Params[0].OtherAtoms.find("17+H10") !=
          std::end(Config.Params[0].OtherAtoms));
    CHECK(Config.Params[0].OtherAtoms.find("17+H11") !=
          std::end(Config.Params[0].OtherAtoms));
}

TEST_CASE("Config reading 2 runs")
{
    std::stringstream ss;
    ss << "test.xtc\n"
       << "test.gro\n"
       << "+++\n"
       << "18+BCDEF\n"
       << "17+O2\n"
       << "17+C65\n"
       << "100\n"
       << "101\n"
       << "17+H10\n"
       << "17+H11\n"
       << "+++\n"
       << "17+O2\n"
       << "18+BCDEF\n"
       << "17+C65\n"
       << "99\n"
       << "98\n"
       << "17+H12\n"
       << "17+H13";

    auto Config = sdf::RuntimeConfig::read(ss);
    CHECK(Config.XtcFile == "test.xtc");
    CHECK(Config.GroFile == "test.gro");
    CHECK(Config.Params.size() == 2);
    CHECK(Config.Params[0].Anchor.toStr() == "18+BCDEF");
    CHECK(Config.Params[0].AtomX.toStr() == "17+O2");
    CHECK(Config.Params[0].AtomXY.toStr() == "17+C65");
    CHECK(Config.Params[0].Distance == 100.0);
    CHECK(Config.Params[0].SliceThickness == 101.0);
    CHECK(Config.Params[0].OtherAtoms.size() == 2);
    CHECK(Config.Params[0].OtherAtoms.find("17+H10") !=
          std::end(Config.Params[0].OtherAtoms));
    CHECK(Config.Params[0].OtherAtoms.find("17+H11") !=
          std::end(Config.Params[0].OtherAtoms));

    CHECK(Config.Params[1].Anchor.toStr() == "17+O2");
    CHECK(Config.Params[1].AtomX.toStr() == "18+BCDEF");
    CHECK(Config.Params[1].AtomXY.toStr() == "17+C65");
    CHECK(Config.Params[1].Distance == 99.0);
    CHECK(Config.Params[1].SliceThickness == 98.0);
    CHECK(Config.Params[1].OtherAtoms.size() == 2);
    CHECK(Config.Params[1].OtherAtoms.find("17+H12") !=
          std::end(Config.Params[1].OtherAtoms));
    CHECK(Config.Params[1].OtherAtoms.find("17+H13") !=
          std::end(Config.Params[1].OtherAtoms));
}
