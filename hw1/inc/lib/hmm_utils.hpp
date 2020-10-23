#ifndef __INCLUDE_SOLVER__
#define __INCLUDE_SOLVER__
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <assert.h>
#include "hmm.hpp"

namespace HMM {
    using namespace std;
    double hmmForward( HMM *hmm, int *observe, int T, double *alpha );
    double hmmBackward( HMM *hmm, int *observe, int T, double *beta );
    double hmmViterbi( HMM *hmm, int *observe, int T, int *result );
    void baumWelchSolver ( HMM *hmm, int *observe, int T, double *accumInitial, double *accumTranstition, double *accumObservation );
}

#endif
