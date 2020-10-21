#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "hmm.hpp"
#include "lib/hmm_utils.hpp"

int main(const int argc, const char **argv)
{
    using namespace std;
	assert(argc>4);
    HMM::HMM hmm;
	int n_iter = atoi(argv[1]);
	const char *data_path = argv[3];
	const char *save_path = argv[4];
    char buffer[MAX_LINE];
    vector<vector<int> > sequence;
    HMM::loadHMM( &hmm, argv[2] );
    double accumInit[MAX_STATE] = {0};
    double accumTrans[MAX_STATE*MAX_STATE] = {0};
    double accumObserv[MAX_STATE*MAX_STATE] = {0};
    int n_state = hmm.state_num;
    // int n_observe = hmm->observ_num;
    FILE *fpSequence = HMM::open_or_die(data_path, "r");
    while ( fscanf(fpSequence, "%s", buffer) > 0 ){
        int T = strlen(buffer);
        vector<int> new_seq(T, 0);
        for (int i=0; i<T; ++i) {
            new_seq[i] = buffer[i]-'A';
        }
        sequence.push_back(new_seq);
    }
    fclose(fpSequence);
    fpSequence = NULL;
    for (int i=0; i<n_iter; ++i) {
        random_shuffle(sequence.begin(), sequence.end());
        memset(accumInit, 0x00, n_state*sizeof(double));
        memset(accumTrans, 0x00, n_state*n_state*sizeof(double));
        memset(accumObserv, 0x00, n_state*n_state*sizeof(double));
        for (int j=0; j<sequence.size(); ++j){
            HMM::baumWelchSolver(&hmm, sequence[j].data(), sequence[j].size(), accumInit, accumTrans, accumObserv);
        }
        for (int j=0; j<n_state; ++j) {
            hmm.initial[j] = accumInit[j] / (double)sequence.size();
            for (int k = 0; k<n_state; ++k) {
                hmm.transition[j][k] = accumTrans[j*n_state+k] / (double)sequence.size();
                hmm.observation[j][k] = accumObserv[j*n_state+k] / (double)sequence.size();
            }
        }
        HMM::dumpHMM(stdout, &hmm);
        fprintf(stdout, "Iteration: %d / %d\n", i+1, n_iter);
    }
    FILE *fpOutputHMM = HMM::open_or_die(save_path, "w");
    HMM::dumpHMM(fpOutputHMM, &hmm);
    fclose(fpOutputHMM);
    fpOutputHMM = NULL;
	return 0;
}
