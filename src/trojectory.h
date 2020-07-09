#ifndef SDF_TROJECTORY_H
#define SDF_TROJECTORY_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include <Eigen/Dense>

#include "xtcio.h"
#include "utils.h"

namespace libmd
{
    struct AtomIdentifier
    {
        AtomIdentifier() = default;
        AtomIdentifier(int s1, const std::string& s2)
                : Res(s1), Name(s2) {}
        AtomIdentifier(const std::string& s);
        AtomIdentifier(const char s[]);

        AtomIdentifier(const AtomIdentifier&) = default;
        AtomIdentifier& operator=(const AtomIdentifier&) = default;

        bool operator==(const AtomIdentifier& rhs) const
        {
            return Res == rhs.Res && Name == rhs.Name;
        }

        std::string toStr() const
        {
            return std::to_string(Res) + "+" + Name;
        }

        int Res;
        std::string Name;
    };

    inline std::ostream& operator<< (std::ostream& stream, const AtomIdentifier& id)
    {
        stream << id.toStr();
        return stream;
    }
}

// Hash for AtomIdentifier
namespace std
{
    template<> struct hash<libmd::AtomIdentifier>
    {
        typedef libmd::AtomIdentifier argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& s) const noexcept
        {
            result_type const h1 ( std::hash<int>{}(s.Res) );
            result_type const h2 ( std::hash<std::string>{}(s.Name) );
            return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };
}

namespace libmd
{
    // A snapshot of a frame of a trojectory. This is only
    // constructable by using filterFrame(), or copying.
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

        V3Map& vec(const AtomIdentifier& atom_name)
        {
            // It is dangerous to use AtomNamesReverse[atom_name]
            // here. Because if atom_name does not exist, it would
            // create a new element in the map and initilize it with
            // 0!!! And this would return the 0th vector!!!
            return Vecs[AtomNamesReverse.at(atom_name)];
        }

        const V3Map& vec(const AtomIdentifier& atom_name) const
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
        bool hasAtom(const AtomIdentifier& name)
        {
            return AtomNamesReverse.find(name) != std::end(AtomNamesReverse);
        }

        size_t size() const { return Vecs.size(); }

        const std::vector<float>& data() const
        {
            return Data;
        }

        const std::vector<V3Map>& vecs() const
        {
            return Vecs;
        }

        std::string debugString() const;

    private:
        TrojectorySnapshot(size_t size)
                : AtomNames(size), Data(size*3) {}

        std::vector<AtomIdentifier> AtomNames;
        std::unordered_map<AtomIdentifier, size_t> AtomNamesReverse;
        XtcFile::FrameMeta Meta;
        std::vector<float> Data;
        std::vector<V3Map> Vecs;

        template <class FrameType, class FilterFunc> friend
        TrojectorySnapshot filterFrame(const FrameType& frame, FilterFunc func);
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

        V3Map& vec(const AtomIdentifier& atom_name)
        {
            // It is dangerous to use AtomNamesReverse[atom_name]
            // here. Because if atom_name does not exist, it would
            // create a new element in the map and initilize it with
            // 0!!! And this would return the 0th vector!!!
            return Vecs[AtomNamesReverse.at(atom_name)];
        }

        const V3Map& vec(const AtomIdentifier& atom_name) const
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

        const std::unordered_map<AtomIdentifier, size_t>& atoms() const
        {
            return AtomNamesReverse;
        }

        bool hasAtom(const AtomIdentifier& name)
        {
            return AtomNamesReverse.find(name) != std::end(AtomNamesReverse);
        }

        size_t index(const AtomIdentifier& name)
        {
            return AtomNamesReverse.at(name);
        }

        size_t size() const { return Vecs.size(); }

        const std::vector<float>& data() const
        {
            return Data;
        }

        const std::vector<V3Map>& vecs() const
        {
            return Vecs;
        }

        void close();
        std::string debugString() const;

    private:
        std::vector<AtomIdentifier> AtomNames;
        std::unordered_map<AtomIdentifier, size_t> AtomNamesReverse;
        XtcFile f;
        XtcFile::FrameMeta Meta;

        // Just a somewhat unrelated FYI: the memory layout of an
        // array of N Eigen Vector3f’s is not nessesarily the same as
        // an array of 3N floats, because of alignment.
        std::vector<float> Data;
        std::vector<V3Map> Vecs;

        template <class FrameType, class FilterFunc> friend
        TrojectorySnapshot filterFrame(const FrameType& frame, FilterFunc func);
    };

    // Filter the atoms by “func”. The argument “func” is a callable
    // that takes 3 arguments: the name of an atom, the index of the
    // atom, and the vector of the atom, and returns a bool. For each
    // atom, “func” is called once. The atom survives if it returns
    // true, and is filtered out otherwise. As an example, the
    // signature of “func” may be
    //
    //   bool func(const AtomIdentifier&, size_t, const V3Map&)
    //
    // This does not change meta(), but it does change size().
    template <class FrameType, class FilterFunc>
    TrojectorySnapshot filterFrame(const FrameType& frame, FilterFunc func)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");

        std::vector<size_t> Passed;
        for(size_t i = 0; i < frame.Vecs.size(); i++)
        {
            if(func(frame.AtomNames[i], i, frame.Vecs[i]))
            {
                Passed.push_back(i);
            }
        }

        TrojectorySnapshot Snap(Passed.size());
        Snap.Meta = frame.Meta;
        Snap.Meta.AtomCount = Passed.size();

        for(size_t i = 0; i < Passed.size(); i++)
        {
            Snap.AtomNames[i] = frame.AtomNames[Passed[i]];
        }

        for(size_t i = 0; i < Passed.size(); i++)
        {
            const size_t NewIdx = Passed[i];
            Snap.Data[i*3] = frame.Data[NewIdx*3];
            Snap.Data[i*3+1] = frame.Data[NewIdx*3+1];
            Snap.Data[i*3+2] = frame.Data[NewIdx*3+2];
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
