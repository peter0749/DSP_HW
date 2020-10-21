#include "hmm.h"
#include "lib/hmm_utils.h"

#define EPSILON 1e-40
#define INDEX2(A,N,I,J) A[(I)*(N)+(J)]
#define INDEX3(A,N,M,I,J,K) A[(K)+(M)*((J)+(N)*(I))]

double hmmForward( HMM *hmm, int *observe, int T, double *alpha )
{
    int s_n = hmm->state_num; //  number of states
    char return_alpha = ( alpha == NULL ) ? 0 : 1;
    if ( !return_alpha )
        alpha = malloc( sizeof(double) * s_n * T ); 
    //  assert( alpha != NULL );
    //  Initialize alpha
    for ( int i = 0; i < s_n; ++i ){
        INDEX2(alpha, s_n, 0, i) = hmm->initial[i] * hmm->observation[i][observe[0]];
    }
    double prob = 0.0;
    for ( int t = 1; t < T; ++t ){
        for ( int i = 0; i < s_n; ++i ){
            prob = 0.0;
            for ( int j = 0; j < s_n; ++j ){
                prob += INDEX2(alpha, s_n, t-1, j) * hmm->transition[j][i];
            }
            INDEX2(alpha, s_n, t, i) = prob * hmm->observation[i][observe[t]];
        }
    }
    prob = 0.0;
    for ( int i = 0; i < s_n; ++i )
        prob += INDEX2(alpha, s_n, T-1, i);
    if ( !return_alpha )
        free(alpha);
    return prob;
}


double hmmBackward( HMM *hmm, int *observe, int T, double *beta )
{
    int s_n = hmm->state_num; //  number of states
    char return_beta = ( beta == NULL ) ? 0 : 1;
    if ( !return_beta )
        beta = malloc( sizeof(double) * s_n * T ); 
    //  assert( beta != NULL );
    for ( int i = 0; i < s_n; ++i )
        INDEX2(beta, s_n, T-1, i) = 1.0;
    double prob = 0.0;
    for ( int t = T-2; t >= 0; --t ){
        for ( int i = 0; i < s_n; ++i ){
            prob = 0.0;
            for ( int j = 0; j < s_n; ++j ){
                prob += hmm->transition[i][j] * hmm->observation[j][observe[t+1]] * INDEX2(beta, s_n, t+1, j);
            }
            INDEX2(beta, s_n, t, i) = prob;
        }
    }
    if ( !return_beta )
        free(beta);
    return 0.0;
}


double hmmViterbi( HMM *hmm, int *observe, int T, int *result )
{
    int s_n = hmm->state_num;
    double *delta = NULL;
    int *backtrack = NULL;
    delta = malloc( sizeof(double) * s_n * 2 ); 
    backtrack = malloc( sizeof(int) * s_n * T ); 
    //  assert( result != NULL );
    //  assert( delta != NULL );
    //  assert( backtrack != NULL );
    for ( int i = 0; i < s_n; ++i )
        delta[i] = hmm->initial[i] * hmm->observation[i][observe[0]];
    double prob = -2e9;
    for ( int t = 1; t < T; ++t ){
        for ( int i = 0; i < s_n; ++i ){
            prob = -2e9;
            for ( int j = 0; j < s_n; ++j ){
                double prob_new = INDEX2(delta, s_n, (t-1)&1, j) * hmm->transition[j][i];
                if (prob_new > prob) {
                    prob = prob_new;
                    INDEX2(backtrack, s_n, t, i) = j;
                }
            }
            INDEX2(delta, s_n, t&1, i) = prob * hmm->observation[i][observe[t]];
        }
    }
    prob = -2e9;
    for ( int i = 0; i < s_n; ++i )
        if ( INDEX2(delta, s_n, (T-1)&1, i) > prob ){
            prob = INDEX2(delta, s_n, (T-1)&1, i);
            if (result != NULL)
                result[T-1] = i;
        }
    if (result != NULL) {
        for ( int t = T-1; t > 0; --t )
            result[t-1] = INDEX2(backtrack, s_n, t, result[t]);
    }
    free(delta);
    free(backtrack);
    return prob;
}


void baumWelchSolver ( HMM *hmm, int *observe, int T, double *accumInitial, double *accumTranstition, double *accumObservation ){
    int s_n = hmm->state_num; //  number of states
    double *alpha = NULL;
    double *beta = NULL;
    double *gamma = NULL;
    double *ps = NULL; // temporary memory
    double *K = NULL;
    alpha = malloc( sizeof(double) * s_n * T );   
    beta  = malloc( sizeof(double) * s_n * T );   
    gamma  = malloc( sizeof(double) * s_n * T );   
    K = malloc( sizeof(double) * s_n * s_n * T );   
    ps = malloc(sizeof(double)*s_n);
    //  assert (alpha != NULL);
    //  assert (beta != NULL);
    //  assert (gamma != NULL);
    //  assert (K != NULL);
    //  assert(ps != NULL);

    hmmForward(hmm, observe, T, alpha);
    hmmBackward(hmm, observe, T, beta);

    double prob = 0.0;
    for ( int t = 0; t < T; ++t ){
        prob = 0.0;
        for ( int i = 0; i < s_n; ++i )
            prob += INDEX2(alpha, s_n, t, i) * INDEX2(beta, s_n, t, i);
        if ( prob < EPSILON ) prob = EPSILON;
        for ( int i = 0; i < s_n; ++i ){
            INDEX2(gamma, s_n, t, i) = INDEX2(alpha, s_n, t, i) * INDEX2(beta, s_n, t, i) / prob;
        }
    }

    for ( int t = 1; t < T; ++t ){
        double prob = 0.0;
        for ( int i = 0; i < s_n; ++i )
            for (int j = 0; j < s_n; ++j)
                prob += INDEX2(alpha, s_n, t-1, i) * \
                        INDEX2(beta, s_n, t, j) * \
                        hmm->transition[i][j] * \
                        hmm->observation[j][observe[t]];
        if ( prob < EPSILON ) prob = EPSILON;
        for (int i = 0; i < s_n; ++i)
            for ( int j = 0; j < s_n; ++j )
                INDEX3(K, s_n, s_n, t-1, i, j) = INDEX2(alpha, s_n, t-1, i) * \
                                                 INDEX2(beta, s_n, t, j) * \
                                                 hmm->transition[i][j] * \
                                                 hmm->observation[j][observe[t]] / prob;
    }

    // Update parameters
    for ( int i = 0; i < s_n; ++i ){
        // Transition
        double p = 0.0, q = 0.0;
        for ( int t = 0; t < T-1; ++t ){
            p += INDEX2(gamma, s_n, t, i);
        }
        if (p < EPSILON) p = EPSILON;
        for ( int j = 0; j < s_n; ++j ){
            q = 0.0;
            for ( int t = 0; t < T-1; ++t )
                q += INDEX3(K, s_n, s_n, t, i, j);
            INDEX2(accumTranstition, s_n, i, j) += q / p;
        }

        // Observation
        q = 0.0;
        memset(ps, 0x00, sizeof(double)*s_n);
        for ( int t = 0; t < T; ++t ){
            double w = INDEX2(gamma, s_n, t, i);
            ps[observe[t]] += w;
            q += w;
        }
        if (q < EPSILON) q = EPSILON;
        for ( int j = 0; j < s_n; ++j ) {
            INDEX2(accumObservation, s_n, i, j) += ps[j] / q;
        }

        // Initial
        accumInitial[i] += INDEX2(gamma, s_n, 0, i);
    }
    
    free(alpha);
    free(beta);
    free(gamma);
    free(K);
    free(ps);
}
