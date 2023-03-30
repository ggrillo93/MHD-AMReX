#include "limiters.H"
#include <algorithm>

double minbee(double r) {
    if (r <= 0) {
        return 0;
    }
    else if (r > 1) {
        return std::min(1., 2./(1+r));
    }
    else {
        return r;
    }
}

double forcezero(double r) {
    return 0.0;
}

double superbee(double r) {
    if (r <= 0) {
        return 0;
    }
    else if (r > 0 && r <= 0.5) {
        return 2*r;
    }
    else if (r > 1) {
        return std::min(std::min(r, 2./(1.+r)), 2.);
    }
    else {
        return 1;
    }
}

double vanleer(double r) {
    if (r <= 0) {
        return 0;
    }
    else {
        return std::min(2.*r/(1. + r), 2./(1. + r));
    }
}

double vanalbada(double r) {
    if (r <= 0) {
        return 0;
    }
    else {
        return std::min(r*(1 + r)/(1 + r*r), 2./(1 + r));
    }
}