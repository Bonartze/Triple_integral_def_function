#include "integral.h"

double f(double x, double y, double z) {
    return x * x * x * x + y * y * y * y + z * z * z * z;
}

double integral_3_f(int n1_begin, int n1_end, int n2_begin, int n2_end, int n3_begin, int n3_end) {
    double I = 0.0;
    double dx = 0.001, dy = 0.004, dz = 0.001; //integration steps
    for (int k = n1_begin; k < n1_end; k++)
        for (int j = n2_begin; j < n2_end; j++)
            for (int i = n3_begin; i < n3_end; i++)
                I += f(i * dx, j * dy, k * dz) * dx * dy * dz;
    return I;
}
