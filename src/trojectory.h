#ifndef SDF_TROJECTORY_H
#define SDF_TROJECTORY_H

#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include <Eigen/Dense>

#include "xtcio.h"

namespace libmd
{
    // A snapshot of a frame of a trojectory. This is only
    // constructable by using Trojectory::filter(), or copying.
    //
    // Has the same interface with Trojectory, minus the file-related
    // part.
    class TrojectorySnapshot
    {
    public:
        TrojectorySnapshot() = delete;
        TrojectorySnapshot(const TrojectorySnapshot& from) = default;
        TrojectorySnapshot& operator=(const TrojectorySnapshot& from)
        = default;

        V3Map& vec(const std::string& atom_name)
        {
            return Vecs[AtomNamesReverse[atom_name]];
        }

        const V3Map& vec(const std::string& atom_name) const
        {
            return Vecs[AtomNamesReverse.at(atom_name)];
        }

        V3Map& vec(size_t i)
        {
            return Vecs[i];
        }

        const V3Map& vec(size_t i) const
        {
            return Vecs[i];
        }

        const XtcFile::FrameMeta& meta() const { return Meta; }
        bool hasAtom(const std::string& name)
        {
            return AtomNamesReverse.find(name) != std::end(AtomNamesReverse);
        }

        size_t size() const { return Vecs.size(); }

        const std::vector<float>& data() const
        {
            return Data;
        }

        std::string debugString() const;

    private:
        TrojectorySnapshot(size_t size)
                : AtomNames(size), Data(size*3) {}

        std::vector<std::string> AtomNames;
        std::unordered_map<std::string, size_t> AtomNamesReverse;
        XtcFile::FrameMeta Meta;
        std::vector<float> Data;
        std::vector<V3Map> Vecs;

        friend class Trojectory;
    };

    // This wraps a XtcFile class and provides high level access
    // throught vec(). If I’m using C++20 I would make a concept out
    // of this and the snapshot type.
    //
    // This class does not care about boxes and boundary conditions.
    class Trojectory
    {
    public:
        void open(const std::string& xtc_path, const std::string& gro_path);
        // Return false if EOF is reached.
        bool nextFrame();

        V3Map& vec(const std::string& atom_name)
        {
            return Vecs[AtomNamesReverse[atom_name]];
        }

        const V3Map& vec(const std::string& atom_name) const
        {
            return Vecs[AtomNamesReverse.at(atom_name)];
        }

        V3Map& vec(size_t i)
        {
            return Vecs[i];
        }

        const V3Map& vec(size_t i) const
        {
            return Vecs[i];
        }

        const XtcFile::FrameMeta& meta() const { return Meta; }

        const std::unordered_map<std::string, size_t>& atoms() const
        {
            return AtomNamesReverse;
        }

        bool hasAtom(const std::string& name)
        {
            return AtomNamesReverse.find(name) != std::end(AtomNamesReverse);
        }

        size_t index(const std::string& name)
        {
            auto Found = AtomNamesReverse.find(name);
            if(Found == std::end(AtomNamesReverse))
            {
                throw std::out_of_range(std::string("unknown atom: ") + name);
            }
            else
            {
                return Found->second;
            }
        }

        size_t size() const { return Vecs.size(); }

        const std::vector<float>& data() const
        {
            return Data;
        }

        void close();

        // Filter the atoms by “func”. The argument “func” is a
        // callable that takes 3 arguments: the name of an atom, the
        // index of the atom, and the vector of the atom, and returns
        // a bool. For each atom, “func” is called once. The atom
        // survives if it returns true, and is filtered out otherwise.
        // As an example, the signature of “func” may be
        //
        //   bool func(const std::string&, size_t, const V3Map&)
        //
        // This does not change meta(), but it does change size().
        template <class FilterFunc>
        TrojectorySnapshot filter(FilterFunc func) const;

        std::string debugString() const;

    private:
        std::vector<std::string> AtomNames;
        std::unordered_map<std::string, size_t> AtomNamesReverse;
        XtcFile f;
        XtcFile::FrameMeta Meta;

        // Just a somewhat unrelated FYI: the memory layout of an
        // array of N Eigen Vector3f’s is not nessesarily the same as
        // an array of 3N floats, because of alignment.
        std::vector<float> Data;
        std::vector<V3Map> Vecs;
    };

    template <class FilterFunc>
    TrojectorySnapshot Trojectory :: filter(FilterFunc func) const
    {
        std::vector<size_t> Passed;
        for(size_t i = 0; i < Vecs.size(); i++)
        {
            if(func(AtomNames[i], i, Vecs[i]))
            {
                Passed.push_back(i);
            }
        }

        TrojectorySnapshot Snap(Passed.size());
        Snap.Meta = Meta;
        Snap.Meta.AtomCount = Passed.size();

        for(size_t i = 0; i < Passed.size(); i++)
        {
            Snap.AtomNames[i] = AtomNames[Passed[i]];
        }

        for(size_t i = 0; i < Passed.size(); i++)
        {
            const size_t NewIdx = Passed[i];
            Snap.Data[i*3] = Data[NewIdx*3];
            Snap.Data[i*3+1] = Data[NewIdx*3+1];
            Snap.Data[i*3+2] = Data[NewIdx*3+2];
        }

        // Rebuild AtomNamesReverse and Vecs.
        Snap.AtomNamesReverse.clear();
        for(size_t i = 0; i < Snap.AtomNames.size(); i++)
        {
            Snap.AtomNamesReverse[Snap.AtomNames[i]] = i;
        }

        Snap.Vecs.clear();
        for(size_t i = 0; i < Snap.AtomNames.size(); i++)
        {
            Snap.Vecs.emplace_back(&(Snap.Data[i*3]));
        }
        return Snap;
    }

} // namespace libmd

#endif
