#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include "hmm.hpp"
#include "lib/hmm_utils.hpp"

int main(const int argc, const char **argv)
{
    using namespace std;
	assert(argc>3);
	const char *model_list_path = argv[1];
	const char *seq_path = argv[2];
	const char *output_path = argv[3];
    char buffer[MAX_LINE];
    vector<vector<int> > sequence;
    HMM::HMM *hmm = (HMM::HMM*)malloc(sizeof(HMM::HMM)*5);
    int hmm_n = HMM::load_models( model_list_path, hmm, 5 );
    FILE *fpSequence = HMM::open_or_die(seq_path, "r");
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
    FILE *fpOutput = HMM::open_or_die(output_path, "w");
    for (int i=0; i<sequence.size(); ++i) {
        int best_m = 0;
        double best = -2e9;
        // vector<int> best_result;
        for (int m=0; m<hmm_n; ++m) {
            // vector<int> result(sequence[i].size(), 0);
            double l = HMM::hmmViterbi(&hmm[m], sequence[i].data(), sequence[i].size(), NULL);
            if (l>best) {
                best = l;
                best_m = m;
                // best_result = result;
            }
        }
        /*
        for (auto v : best_result) printf("%c", (char)(v+'A'));
        printf("\n");
        for (auto v : sequence[i]) printf("%c", (char)(v+'A'));
        printf("\n");
        */
        fprintf(fpOutput, "%s %e\n", hmm[best_m].model_name, best);
    }
    fclose(fpOutput);
    fpOutput = NULL;
    free(hmm);
	return 0;
}
