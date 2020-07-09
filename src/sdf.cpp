#include "sdf.h"
#include "utils.h"

using namespace libmd;

namespace libmd
{
    Eigen::Matrix3f rotateToAlignX(const VecRefType& vec,
                                   const VecRefType& to_xy)
    {
        Eigen::Matrix3f Rot1;
        {
            Eigen::Vector3f Axis;
            Axis[0] = 0.0f;
            Axis[1] = vec[2];
            Axis[2] = -vec[1];
            Axis /= Axis.norm();

            const float Cos = vec[0] / vec.norm();
            const float Sin = std::sqrt(1.0f - Cos*Cos);

            Rot1(0, 0) = Cos;
            Rot1(0, 1) = -Axis[2] * Sin;
            Rot1(0, 2) = Axis[1] * Sin;
            Rot1(1, 0) = Axis[2] * Sin;
            Rot1(1, 1) = Cos + Axis[1] * Axis[1] * (1.0 - Cos);
            Rot1(1, 2) = Axis[1] * Axis[2] * (1.0 - Cos);
            Rot1(2, 0) = -Axis[1] * Sin;
            Rot1(2, 1) = Axis[2] * Axis[1] * (1.0 - Cos);
            Rot1(2, 2) = Cos + Axis[2] * Axis[2] * (1 - Cos);
        }

        const Eigen::Vector3f Rotated = Rot1 * to_xy;
        Eigen::Matrix3f Rot2 = Eigen::Matrix3f::Zero();
        {
            Eigen::Vector3f VecProj = Rotated;
            float Cos = Rotated[1] / std::sqrt(Rotated[1] * Rotated[1] +
                                                     Rotated[2] * Rotated[2]);
            float Sin = std::sqrt(1.0f - Cos*Cos);
            if(Rotated[2] > 0)
            {
                Sin = -Sin;
            }

            Rot2(0, 0) = 1.0;
            Rot2(1, 1) = Cos;
            Rot2(1, 2) = -Sin;
            Rot2(2, 1) = Sin;
            Rot2(2, 2) = Cos;
        }

        return Rot2 * Rot1;
    }
} // namespace libmd

namespace sdf
{
    using namespace libmd;
    TrojectorySnapshot prepareFrame(
        const Parameters& params, const Trojectory& frame)
    {
        const auto& AtomX = frame.vec(params.AtomX);
        const auto& BoxDim = frame.meta().BoxDim;
        const RectPbc3d Pbc(BoxDim[0][0], BoxDim[1][1], BoxDim[2][2]);

        // Filter by distance
        auto Frame = filterFrame(
            frame,
            [&](const libmd::AtomIdentifier& name, size_t _, const auto& vec)
            {
                UNUSED(_);
                if(name == params.Anchor || name == params.AtomX ||
                   name == params.AtomXY)
                {
                    return true;
                }

                // if(params.OtherAtoms.find(name) == std::end(params.OtherAtoms))
                // {
                //     return false;
                // }

                return Pbc.dist(AtomX, vec) < params.Distance;
            });

        wrapFrame(params.AtomX, Frame);
        Eigen::Vector3f ShiftBy = -(Frame.vec(params.Anchor));
        shiftFrame(ShiftBy, Frame);
        auto Rot = rotateToAlignX(Frame.vec(params.AtomX),
                                  Frame.vec(params.AtomXY));
        rotateFrame(Rot, Frame);
        // std::cout << Frame.debugString() << std::endl;

        // Ditch the anchors, x atoms, and xy atoms, and take a slice
        // at the XY plane.
        float HalfThickness = params.SliceThickness * 0.5;
        return filterFrame(
            Frame,
            [&](const libmd::AtomIdentifier& id, auto _1, const auto& pos)
            {
                UNUSED(_1);
                return !(id.Name == params.Anchor.Name ||
                         id.Name == params.AtomX.Name ||
                         id.Name == params.AtomXY.Name) &&
                    std::fabs(pos[2]) <= HalfThickness;
            });
    }

    void run(const RuntimeConfig& config)
    {
        libmd::Trojectory t;
        t.open(config.XtcFile, config.GroFile);
        t.nextFrame();

        std::cout << prepareFrame(config.Params[0], t).debugString() << std::endl;
    }

} // namespace libmd
