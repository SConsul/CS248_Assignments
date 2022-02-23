
#include "../geometry/spline.h"
#include "debug.h"

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

    // TODO (Animation): Task 1a
    // Given time in [0,1] compute the cubic spline coefficients and use them to compute
    // the interpolated value at time 'time' based on the positions & tangents

    // Note that Spline is parameterized on type T, which allows us to create splines over
    // any type that supports the * and + operators.
    float t = time, t2 = time*time, t3 = time*time*time;
    float h_00 = 2*t3-3*t2+1;
    float h_10 = t3-2*t2+t;
    float h_01 = -2*t3+3*t2;
    float h_11 = t3-t2;
    return h_00*position0 + h_10*tangent0 + h_01*position1 + h_11*tangent1;
}

template<typename T> T Spline<T>::at(float time) const {

    // TODO (Animation): Task 1b

    // Given a time, find the nearest positions & tangent values
    // defined by the control point map.

    // Transform them for use with cubic_unit_spline

    // Be wary of edge cases! What if time is before the first knot,
    // before the second knot, etc...

    return cubic_unit_spline(0.0f, T(), T(), T(), T());
}
