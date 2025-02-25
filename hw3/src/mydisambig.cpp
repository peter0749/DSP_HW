#include <omp.h>
#include <iostream>
#include <deque>
#include <vector>
#include <list>
#include <string>
#include <cstdio>
#include <cfloat>
#include <utility>
#include <unordered_map>
#include <map>
#include <queue>
#include <algorithm>
// #include <assert.h>
#include "Ngram.h"

using namespace std;

// Hash function for unordered_map< pair<T1, T2> >
typedef struct hash_pair_struct {
    template <class T1, class T2>
    size_t operator() (const pair<T1, T2> &p) const {
        size_t lhs = hash<T1>{}(p.first);
        size_t rhs = hash<T2>{}(p.second);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        return lhs;
    }
    template <class T1, class T2, class T3>
    size_t operator() (const pair<pair<T1, T2>, T3> &p) const {
        size_t lhs = hash<T1>{}(p.first.first);
        size_t rhs = hash<T2>{}(p.first.second);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        rhs = hash<T3>{}(p.second);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        return lhs;
    }
    template <class T1, class T2, class T3>
    size_t operator() (const pair<T1, pair<T2, T3> > &p) const {
        size_t lhs = hash<T2>{}(p.second.first);
        size_t rhs = hash<T3>{}(p.second.second);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        rhs = lhs;
        lhs = hash<T1>{}(p.first);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        return lhs;
    }
} hash_pair;

list<VocabIndex> viterbi_bigram(
        const vector<string> &seq,
        const unordered_map<string, vector<string> > &mapping,
        Vocab &voc,
        Ngram &lm) {
    unordered_map< pair<int, VocabIndex>, pair<VocabIndex, double>, hash_pair > backtrack;
    VocabIndex best_qi_at_T = 0;
    double best_at_T = -DBL_MAX;
    int T = seq.size();
    // Unigram
    vector<string> s_0 = (mapping.count(seq[0])==0) ? vector<string>(1, seq[0]) : mapping.at(seq[0]);
    for (auto s_qi : s_0) {
        VocabIndex qi = voc.getIndex(s_qi.c_str());
        if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
        VocabIndex c0[] = {Vocab_None};
        VocabIndex c1[] = {voc.getIndex(Vocab_SentStart), Vocab_None};
        double logP = max(lm.wordProb(qi, c0), lm.wordProb(qi, c1) + lm.wordProb(voc.getIndex(Vocab_SentStart), c0));
        pair<int, VocabIndex> key(0, qi);
        backtrack[key] = pair<VocabIndex, double>(0, logP); // [0, qi] -> logP
        if (T == 1) {
            if (logP > best_at_T) {
                best_at_T = logP;
                best_qi_at_T = qi;
            }
        }
    }
    // Bigram
    deque<vector<string> > s_n;
    s_n.push_back({}); // s_-1
    s_n.push_back(s_0); // s_0
    for (int t=1; t<T; ++t) {
        s_n.pop_front(); // discard s_-1
        s_n.push_back((mapping.count(seq[t]) == 0) ? vector<string>(1,seq[t]) : mapping.at(seq[t])); // s_1
        for (auto s_qi : s_n[1]) {
            VocabIndex qi = voc.getIndex(s_qi.c_str());
            if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
            VocabIndex best_qj = 0;
            double best_p = -DBL_MAX;
            for (auto s_qj : s_n[0]) {
                VocabIndex qj = voc.getIndex(s_qj.c_str());
                if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
                pair<int, VocabIndex> key(t-1, qj);
                double log_delta_tm1_qj = backtrack[key].second;
                VocabIndex c0[] = {qj, Vocab_None};
                double logP_qi_qj = lm.wordProb(qi, c0);
                double score = logP_qi_qj + log_delta_tm1_qj;
                if (score > best_p) {
                    best_qj = qj;
                    best_p = score;
                }
            }
            // delta_t_qi == best_p
            backtrack[pair<int, VocabIndex>(t, qi)] = pair<VocabIndex, double>(best_qj, best_p);
            if (t == T-1) {
                if (best_p > best_at_T) {
                    best_qi_at_T = qi;
                    best_at_T = best_p;
                }
            }
        }
    }
    s_n.clear();
    list<VocabIndex> best_path;
    best_path.push_front(best_qi_at_T); // actually in indiex T-1
    for (int t=T-1; t>=1; --t) {
        VocabIndex best_qi_prev = backtrack[pair<int, VocabIndex>(t, best_path.front())].first;
        best_path.push_front(best_qi_prev);
    }
    return best_path;
}

list<VocabIndex> viterbi_trigram(
        const vector<string> &seq,
        const unordered_map<string, vector<string> > &mapping,
        Vocab &voc,
        Ngram &lm) {
    const int n_beams = 200;
    // backtracking: (t, i, j) -> (best_k, value)
    unordered_map< pair<int, pair<VocabIndex, VocabIndex> >, pair<VocabIndex, double>, hash_pair > backtrack;
    // map< pair<int, pair<VocabIndex, VocabIndex> >, pair<VocabIndex, double> > backtrack;
    pair<VocabIndex, VocabIndex> best_i_j_at_T(0,0);
    double best_at_T = -DBL_MAX;
    int T = seq.size();
    if (T<3) {
        // fprintf(stderr, "warning: String length < 3. Automatically switch to bigram mode.\n");
        return viterbi_bigram(seq, mapping, voc, lm);
    }
    //Bigram
    vector<string> s_1 = (mapping.count(seq[1]) == 0) ? vector<string>(1,seq[1]) : mapping.at(seq[1]);
    vector<string> s_0 = (mapping.count(seq[0]) == 0) ? vector<string>(1,seq[0]) : mapping.at(seq[0]);
    for (auto s_qi : s_1) {
        VocabIndex qi = voc.getIndex(s_qi.c_str());
        if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
        for (auto s_qj : s_0) {
            VocabIndex qj = voc.getIndex(s_qj.c_str());
            if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
            // delta_2(q_i,q_j) = max_{qk} { P(q_i|q_j) delta_1(q_j, q_k) }
            //                  = max_{qk} { P(q_i|q_j) P(q_j) }
            //                  = P(q_i|q_j) P(q_j)
            VocabIndex c0[] = {qj, Vocab_None};
            VocabIndex c1[] = {Vocab_None};
            VocabIndex c2[] = {voc.getIndex(Vocab_SentStart), Vocab_None};
            double log_delta = lm.wordProb(qi, c0) + max(lm.wordProb(qj, c1), lm.wordProb(qj, c2) + lm.wordProb(voc.getIndex(Vocab_SentStart), c1));
            pair<int, pair<VocabIndex, VocabIndex> > key = {1, {qi, qj}};
            backtrack[key] = pair<VocabIndex, double>(0, log_delta); // -1: This is the first element.
        }
    }
    // Trigram
    deque<vector<string> > s_n;
    s_n.push_back({}); // s_-1
    s_n.push_back(s_0); // s_0
    s_n.push_back(s_1); // s_1
    for (int t=2; t<T; ++t) {
        s_n.pop_front(); // discard s_-1
        s_n.push_back((mapping.count(seq[t]) == 0) ? vector<string>(1,seq[t]) : mapping.at(seq[t])); // s_2
        // cerr << s_n[0].size() << " " << s_n[1].size() << " " << s_n[2].size() << endl;
        for (auto s_qi : s_n[2]) {
            VocabIndex qi = voc.getIndex(s_qi.c_str());
            if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
            for (auto s_qj : s_n[1]) {
                VocabIndex qj = voc.getIndex(s_qj.c_str());
                if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
                VocabIndex best_qk = 0;
                double best_p = -DBL_MAX;
                for (auto s_qk : s_n[0]) {
                    VocabIndex qk = voc.getIndex(s_qk.c_str());
                    if (qk == Vocab_None) qk = voc.getIndex(Vocab_Unknown);
                    pair<int, pair<VocabIndex, VocabIndex> > key = {t-1, {qj, qk}};
                    double log_delta_tm1_qj_qk = backtrack[key].second;
                    VocabIndex c0[] = {qj, qk, Vocab_None};
                    double logP_qi_qj_qk = lm.wordProb(qi, c0);
                    double score = logP_qi_qj_qk + log_delta_tm1_qj_qk;
                    if (score > best_p) {
                        best_qk = qk;
                        best_p = score;
                    }
                }
                pair<int, pair<VocabIndex, VocabIndex> > key = {t, {qi, qj}};
                backtrack[key] = pair<VocabIndex, double>(best_qk, best_p);
                if (t == T-1) {
                    if (best_p > best_at_T) {
                        best_i_j_at_T = pair<VocabIndex, VocabIndex>(qi, qj);
                        best_at_T = best_p;
                    }
                }
            }
        }
        // Beam Pruning
        priority_queue<pair<double,VocabIndex>, vector<pair<double,VocabIndex> > > beam; // find pruned s_n[2]
        for (auto s_qi : s_n[2]) {
            VocabIndex qi = voc.getIndex(s_qi.c_str());
            if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
            VocabIndex best_j = 0;
            double best_at_i = -DBL_MAX;
            for (auto s_qj : s_n[1]) {
                VocabIndex qj = voc.getIndex(s_qj.c_str());
                if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
                pair<int, pair<VocabIndex, VocabIndex> > key = {t, {qi, qj}};
                double prob = backtrack[key].second;
                if (prob > best_at_i) {
                    best_j = qj;
                    best_at_i = prob;
                }
            }
            if (beam.size() < n_beams) beam.push({-best_at_i, qi}); // initialize top-k
            else if (best_at_i > -beam.top().first) { // update top-kth
                beam.pop();
                beam.push({-best_at_i, qi});
            }
        }
        s_n[2].clear(); // prune s_n[2]
        while (!beam.empty()) {
            VocabIndex qi = beam.top().second;
            beam.pop();
            string s_qi(voc.getWord(qi));
            s_n[2].push_back(s_qi);
        }
    }
    s_n.clear();
    list<VocabIndex> best_path;
    best_path.push_front(best_i_j_at_T.first); // T-1
    best_path.push_front(best_i_j_at_T.second); // T-2
    VocabIndex curr_t2 = best_i_j_at_T.first;
    VocabIndex curr_t1 = best_i_j_at_T.second;
    for (int t=T-1; t>=2; --t) {
        pair<int, pair<VocabIndex, VocabIndex> > key = {t, {curr_t2, curr_t1}};
        VocabIndex curr_t0 = backtrack[key].first;
        // Update path
        best_path.push_front(curr_t0);
        // Move to next
        curr_t2 = curr_t1;
        curr_t1 = curr_t0;
    }
    
    return best_path;
}

int main(int argc, char *argv[]) {
    const char *input_file = argv[1];
    const char *mapping_file = argv[2];
    const char *lm_file = argv[3];
    const char *output_file = argv[4];
    char *buffer = NULL;
    const char *tok = " \n\t", *ptr = NULL;
    int ngram = 0;
    Vocab voc;
    {
        File lmFile(lm_file, "r");
        char *temp = NULL;
        while(temp = lmFile.getline()) {
            char *pos = strstr(temp, "-grams");
            if (pos != NULL && pos - temp > 1) {
                char numbers[16];
                size_t len = pos-temp-1;
                strncpy(numbers, temp+1, len);
                numbers[len] = '\0';
                ngram = max(ngram, atoi(numbers));
            }
        }
        lmFile.close();
    }
    // The viterbi function pointer
    list<VocabIndex> (*viterbi)(
        const vector<string>&,
        const unordered_map<string, vector<string> > &,
        Vocab &,
        Ngram &);
    switch (ngram) {
        case 2:
            viterbi = &viterbi_bigram;
            break;
        case 3:
            viterbi = &viterbi_trigram;
            break;
        default:
            fprintf(stderr, "Error: Only supports bigram and trigram!\n");
            exit(2);
    }
    Ngram lm(voc, ngram);
    {
        File lmFile(lm_file, "r");
        lm.read(lmFile);
        lmFile.close();
    }
    unordered_map<string, vector<string> > mapping;
    {
        File mapping_fp(mapping_file, "r");
        while((buffer = mapping_fp.getline()) != NULL) {
            ptr = strtok(buffer, tok);
            // assert(ptr != NULL);
            string key(ptr);
            mapping[key] = vector<string>();
            ptr = strtok(NULL, tok);
            while(ptr != NULL) {
                mapping[key].push_back(string(ptr));
                ptr = strtok(NULL, tok);
            }
        }
        mapping_fp.close();
    }
    VocabIndex unk_ind = voc.getIndex(Vocab_Unknown);
    File input_text(input_file, "r");
    File output_text(output_file, "w");
    size_t num_threads = omp_get_max_threads();
    size_t pool_size = 4;
    vector<vector<string> > input_strings;
    vector<list<VocabIndex> > output_sequence(num_threads * pool_size, list<VocabIndex>() );
    while((buffer = input_text.getline()) != NULL) {
        ptr = strtok(buffer, tok);
        vector<string> s;
        while (ptr != NULL) {
            s.push_back(string(ptr));
            ptr = strtok(NULL, tok);
        }
        s.push_back(string(Vocab_SentEnd));
        input_strings.push_back(s);
        if (input_strings.size() >= num_threads * pool_size) {
            #pragma omp parallel for schedule(dynamic, 1) shared(input_strings, output_sequence, mapping, voc, lm)
            for (int i=0; i<input_strings.size(); ++i) {
                output_sequence[i] = viterbi(input_strings[i], mapping, voc, lm);
            }
            for (int i=0; i<input_strings.size(); ++i) {
                int t = 0;
                fprintf(output_text, "%s", Vocab_SentStart);
                for (auto j : output_sequence[i]) {
                    fprintf(output_text, " %s", (j == unk_ind) ? input_strings[i][t].c_str() : voc.getWord(j));
                    ++t;
                }
                fprintf(output_text, "\n");
                fflush(output_text);
            }
            input_strings.clear();
        }
        
    }
    if (input_strings.size() > 0) {
        #pragma omp parallel for schedule(dynamic, 1) shared(input_strings, output_sequence, mapping, voc, lm)
        for (int i=0; i<input_strings.size(); ++i) {
            output_sequence[i] = viterbi(input_strings[i], mapping, voc, lm);
        }
        for (int i=0; i<input_strings.size(); ++i) {
            int t = 0;
            fprintf(output_text, "%s", Vocab_SentStart);
            for (auto j : output_sequence[i]) {
                fprintf(output_text, " %s", (j == unk_ind) ? input_strings[i][t].c_str() : voc.getWord(j));
                ++t;
            }
            fprintf(output_text, "\n");
            fflush(output_text);
        }
        input_strings.clear();
    }
    output_sequence.clear();
    output_text.close();
    input_text.close();
    return 0;
}
