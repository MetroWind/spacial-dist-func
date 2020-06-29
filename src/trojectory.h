#ifndef SDF_TROJECTORY_H
#define SDF_TROJECTORY_H

#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "xtcio.h"

namespace libmd
{
    class Trojectory
    {
    public:
        void open(const std::string& xtc_path, const std::string& gro_path);
        // Return false if EOF is reached.
        bool nextFrame();

        Eigen::Map<Eigen::Vector3f>& vec(const std::string& atom_name)
        {
            return Vecs[AtomNamesReverse[atom_name]];
        }

        const Eigen::Map<Eigen::Vector3f>&
        vec(const std::string& atom_name) const
        {
            return Vecs[AtomNamesReverse.at(atom_name)];
        }

        Eigen::Map<Eigen::Vector3f>& vec(size_t i)
        {
            return Vecs[i];
        }

        const Eigen::Map<Eigen::Vector3f>& vec(size_t i) const
        {
            return Vecs[i];
        }

        const XtcFile::FrameMeta& meta() const { return Meta; }

        void close();

    private:
        std::vector<std::string> AtomNames;
        std::unordered_map<std::string, size_t> AtomNamesReverse;
        XtcFile f;
        XtcFile::FrameMeta Meta;

        // Just a somewhat unrelated FYI: the memory layout of an
        // array of N Eigen Vector3f’s is not nessesarily the same as
        // an array of 3N floats, because of alignment.
        std::vector<float> Data;
        std::vector<Eigen::Map<Eigen::Vector3f>> Vecs;
    };

} // namespace libmd

#endif
