// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.1.0 (2016/01/25)

#ifndef ORIENTED_BOX_HPP
#define ORIENTED_BOX_HPP

// A box has center C, axis directions U[i], and extents e[i].  The set
// {U[0],...,U[N-1]} is orthonormal, which means the vectors are
// unit-length and mutually perpendicular.  The extents are nonnegative;
// zero is allowed, meaning the box is degenerate in the corresponding
// direction.  A point X is represented in box coordinates by
// X = C + y[0]*U[0] + y[1]*U[1].  This point is inside or on the
// box whenever |y[i]| <= e[i] for all i.

template <int N>
class OrientedBox {
public:
    // Construction and destruction.  The default constructor sets the center
    // to (0,...,0), axis d to Vector<N,Real>::Unit(d), and extent d to +1.
    OrientedBox();
    OrientedBox(Vector const &inCenter, Mat const &inAxis, Vector const &inExtent);

    // Public member access.  It is required that extent[i] >= 0.
    Vector center;
    Mat axis;
    Vector extent;

public:
    // Comparisons to support sorted containers.
    bool operator==(OrientedBox const& box) const;
    bool operator!=(OrientedBox const& box) const;
    bool operator< (OrientedBox const& box) const;
    bool operator<=(OrientedBox const& box) const;
    bool operator> (OrientedBox const& box) const;
    bool operator>=(OrientedBox const& box) const;
};

// Template aliases for convenience.
using OrientedBox2 = OrientedBox<2>;

using OrientedBox3 = OrientedBox<3>;


template <int N>
OrientedBox<N>::OrientedBox() {
    center = Vector::Zero(N);
    for (int i = 0; i < N; ++i) {
        //axis[i].MakeUnit(i);
        axis.row(i) = Vector::Zero(N);
        axis(i,i) = 1;
        extent[i] = 1.0;
    }
}

template <int N>
OrientedBox<N>::OrientedBox(Vector const &inCenter, Mat const &inAxis, Vector const &inExtent)
    : center(inCenter), axis(inAxis), extent(inExtent) {
}

template <int N>
bool OrientedBox<N>::operator==(OrientedBox const &box) const {
    return center == box.center && axis == box.axis && extent == box.extent;
}

template <int N>
bool OrientedBox<N>::operator!=(OrientedBox const &box) const {
    return !operator==(box);
}

template <int N>
bool OrientedBox<N>::operator<(OrientedBox const& box) const {
    if (center < box.center) {
        return true;
    }

    if (center > box.center) {
        return false;
    }

    if (axis < box.axis) {
        return true;
    }

    if (axis > box.axis) {
        return false;
    }

    return extent < box.extent;
}

template <int N>
bool OrientedBox<N>::operator<=(OrientedBox const& box) const {
    return operator<(box) || operator==(box);
}

template <int N>
bool OrientedBox<N>::operator>(OrientedBox const& box) const {
    return !operator<=(box);
}

template <int N>
bool OrientedBox<N>::operator>=(OrientedBox const& box) const {
    return !operator<(box);
}

#endif // ORIENTED_BOX_HPP