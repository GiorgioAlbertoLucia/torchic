#pragma once

#include <cmath>

/**
 *  Bethe-Bloch formula for the energy loss of a particle in a medium.
 *  Reference: https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/Detectors/TPC/include/DataFormatsTPC/BetheBlochAleph.h
 *
 *  Args:
 *    bg (double): beta*gamma -> pass it as pTPC / mass
 *    kp1 (double): constant
 *    kp2 (double): constant
 *    kp3 (double): constant
 *    kp4 (double): constant
 *    kp5 (double): constant
*/
double BetheBlochAleph(double bg, double kp1, double kp2, double kp3, double kp4, double kp5)
{
    double beta = bg / sqrt(1. + bg * bg);
    double aa = pow(beta, kp4);
    double bb = pow((1 / bg), kp5);
    bb = log(kp3 + bb);
    return (kp2 - aa - bb) * kp1 / aa;

}

double BetheBloch(double * xs, double * params)
{
    double bg = xs[0];
    double kp1 = params[0];
    double kp2 = params[1];
    double kp3 = params[2];
    double kp4 = params[3];
    double kp5 = params[4];

    return BetheBlochAleph(bg, kp1, kp2, kp3, kp4, kp5);
}