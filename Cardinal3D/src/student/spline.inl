
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

    if(!this->any()){  // No control points
        return T(); 
    }
    else if(this->control_points.size() == 1){  // Only one knot. Just return its value
        return this->control_points.begin()->second;
    }
    else if(this->has(time)){  // Current time is a knot time
        return this->control_points.at(time);
    }
    else if (time < this->control_points.begin()->first) {  // Current time is less than smallest knot time
        return this->control_points.begin()->second;
    }
    else if (time > this->control_points.rbegin()->first) { // Current time is greater than largest knot time
        return this->control_points.rbegin()->second;
    }
    else{
        auto lowestGreaterThanTime = this->control_points.upper_bound(time);
        auto greatestLowerThanTime = std::prev(lowestGreaterThanTime);
        
        float t[4];
        t[1] = greatestLowerThanTime->first; 
        t[2] = lowestGreaterThanTime->first;
        T k[4];
        k[1] = greatestLowerThanTime->second;
        k[2] = lowestGreaterThanTime->second;

        if(greatestLowerThanTime == this->control_points.begin()){ //  No t0 exists. Mirroring
            t[0] = t[1] - (t[2] - t[1]);
            k[0] = k[1] - (k[2] - k[1]); 
        }
        else{
            t[0] = std::prev(greatestLowerThanTime)->first;
            k[0] = std::prev(greatestLowerThanTime)->second;
        }

        if(lowestGreaterThanTime == std::prev(this->control_points.end())){ //  No t3 exists. Mirroring
            t[3] = t[2] + (t[2] - t[1]);
            k[3] = k[2] + (k[2] - k[1]);
        }  
        else{
            t[3] = std::next(lowestGreaterThanTime)->first;
            k[3] = std::next(lowestGreaterThanTime)->second;
        }

        // Transform them for use with cubic_unit_spline
        float offset = t[0], normFactor = t[3]-t[0]; // Convert times to lie between 0 and 1
        for (size_t i = 0; i < sizeof(t) / sizeof(t[0]); i++) {
            t[i] = (t[i] - offset) / normFactor;
        }
        float timeNorm = (time - offset) / normFactor;

        T tangent1 = (k[2] - k[0]) / (t[2] - t[0]);
        T tangent2 = (k[3] - k[1]) / (t[3] - t[1]);

        return cubic_unit_spline(timeNorm, k[1], k[2], tangent1, tangent2);
    }

    // Be wary of edge cases! What if time is before the first knot,
    // before the second knot, etc...
}
