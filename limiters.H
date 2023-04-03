#include <functional>

extern std::function<double(double)> limFunc; // slope limiter for MUSCL-Hancock

double minbee(double);
double forcezero(double);
double superbee(double);
double vanleer(double);
double vanalbada(double);