#include <omp.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <cstdio>
#include <cfloat>
#include <utility>
#include <unordered_map>
#include <map>
#include <assert.h>
#include "Ngram.h"

using namespace std;

// Hash function for unordered_map< pair<T1, T2> >
typedef struct hash_pair_struct {
    template <class T1, class T2>
    size_t operator() (const pair<T1, T2> &p) const {
        return hash<T1>{}(p.first) ^ hash<T2>{}(p.second);
    }
    template <class T1, class T2, class T3>
    size_t operator() (const pair<pair<T1, T2>, T3> &p) const {
        return hash<T1>{}(p.first.first) ^ hash<T2>{}(p.first.second) ^ hash<T3>{}(p.second);
    }
    template <class T1, class T2, class T3>
    size_t operator() (const pair<T1, pair<T2, T3> > &p) const {
        return hash<T1>{}(p.first) ^ hash<T2>{}(p.second.first) ^ hash<T3>{}(p.second.second);
    }
} hash_pair;

list<VocabIndex> viterbi_bigram(
        const vector<string> &seq,
        const unordered_map<string, vector<string> > &mapping,
        Vocab &voc,
        Ngram &lm) {
    unordered_map< pair<unsigned int, VocabIndex>, pair<VocabIndex, double>, hash_pair > backtrack;
    VocabIndex best_qi_at_T = 0;
    double best_at_T = -DBL_MAX;
    unsigned int T = seq.size();
    // Unigram
    for (string s_qi : mapping.at(seq[0])) {
        VocabIndex qi = voc.getIndex(s_qi.c_str());
        if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
        VocabIndex context[] = {Vocab_None};
        double logP = lm.wordProb(qi, context);
        pair<unsigned int, VocabIndex> key(0, qi);
        backtrack[key] = pair<VocabIndex, double>(-1, logP); // [0, qi] -> logP
        if (T == 1) {
            if (logP > best_at_T) {
                best_at_T = logP;
                best_qi_at_T = qi;
            }
        }
    }
    // Bigram
    for (unsigned int t=1; t<T; ++t) {
        for (string s_qi : mapping.at(seq[t])) {
            VocabIndex qi = voc.getIndex(s_qi.c_str());
            if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
            VocabIndex best_qj = 0;
            double best_p = -DBL_MAX;
            for (string s_qj : mapping.at(seq[t-1])) {
                VocabIndex qj = voc.getIndex(s_qj.c_str());
                if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
                pair<unsigned int, VocabIndex> key(t-1, qj);
                double log_delta_tm1_qj = backtrack[key].second;
                VocabIndex context[] = {qj, Vocab_None};
                double logP_qi_qj = lm.wordProb(qi, context);
                double score = logP_qi_qj + log_delta_tm1_qj;
                if (score > best_p) {
                    best_qj = qj;
                    best_p = score;
                }
            }
            // delta_t_qi == best_p
            backtrack[pair<unsigned int, VocabIndex>(t, qi)] = pair<VocabIndex, double>(best_qj, best_p);
            if (t == T-1) {
                if (best_p > best_at_T) {
                    best_qi_at_T = qi;
                    best_at_T = best_p;
                }
            }
        }
    }
    list<VocabIndex> best_path;
    best_path.push_front(best_qi_at_T); // actually in indiex T-1
    for (int t=T-1; t>=1; --t) {
        VocabIndex best_qi_prev = backtrack[pair<unsigned int, VocabIndex>(t, best_path.front())].first;
        best_path.push_front(best_qi_prev);
    }
    return best_path;
}

list<VocabIndex> viterbi_trigram(
        const vector<string> &seq,
        const unordered_map<string, vector<string> > &mapping,
        Vocab &voc,
        Ngram &lm) {
    // backtracking: (t, i, j) -> (best_k, value)
    unordered_map< pair<unsigned int, pair<VocabIndex, VocabIndex> >, pair<VocabIndex, double>, hash_pair > backtrack;
    // map< pair<unsigned int, pair<VocabIndex, VocabIndex> >, pair<VocabIndex, double> > backtrack;
    pair<VocabIndex, VocabIndex> best_i_j_at_T(0,0);
    double best_at_T = -DBL_MAX;
    unsigned int T = seq.size();
    if (T<3) {
        fprintf(stderr, "warning: String length < 3. Automatically switch to bigram mode.\n");
        return viterbi_bigram(seq, mapping, voc, lm);
    }
    //Bigram
    for (string s_qi : mapping.at(seq[1])) {
        VocabIndex qi = voc.getIndex(s_qi.c_str());
        if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
        for (string s_qj : mapping.at(seq[0])) {
            VocabIndex qj = voc.getIndex(s_qj.c_str());
            if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
            // delta_2(q_i,q_j) = max_{qk} { P(q_i|q_j) delta_1(q_j, q_k) }
            //                  = max_{qk} { P(q_i|q_j) P(q_j) }
            //                  = P(q_i|q_j) P(q_j)
            VocabIndex c0[] = {Vocab_None};
            VocabIndex c1[] = {qj, Vocab_None};
            double log_delta = lm.wordProb(qi, c1) + lm.wordProb(qj, c0);
            pair<unsigned int, pair<VocabIndex, VocabIndex> > key = {1, {qi, qj}};
            backtrack[key] = pair<VocabIndex, double>(-1, log_delta); // -1: This is the first element.
        }
    }
    // Trigram
    for (unsigned int t=2; t<T; ++t) {
        for (string s_qi : mapping.at(seq[t])) {
            VocabIndex qi = voc.getIndex(s_qi.c_str());
            if (qi == Vocab_None) qi = voc.getIndex(Vocab_Unknown);
            for (string s_qj : mapping.at(seq[t-1])) {
                VocabIndex qj = voc.getIndex(s_qj.c_str());
                if (qj == Vocab_None) qj = voc.getIndex(Vocab_Unknown);
                VocabIndex best_qk = 0;
                double best_p = -DBL_MAX;
                for (string s_qk : mapping.at(seq[t-2])) {
                    VocabIndex qk = voc.getIndex(s_qk.c_str());
                    if (qk == Vocab_None) qk = voc.getIndex(Vocab_Unknown);
                    pair<unsigned int, pair<VocabIndex, VocabIndex> > key = {t-1, {qj, qk}};
                    double log_delta_tm1_qj_qk = backtrack[key].second;
                    VocabIndex context[] = {qj, qk, Vocab_None};
                    double logP_qi_qj_qk = lm.wordProb(qi, context);
                    double score = logP_qi_qj_qk + log_delta_tm1_qj_qk;
                    if (score > best_p) {
                        best_qk = qk;
                        best_p = score;
                    }
                }
                pair<unsigned int, pair<VocabIndex, VocabIndex> > key = {t, {qi, qj}};
                backtrack[key] = pair<VocabIndex, double>(best_qk, best_p);
                if (t == T-1) {
                    if (best_p > best_at_T) {
                        best_i_j_at_T = pair<VocabIndex, VocabIndex>(qi, qj);
                        best_at_T = best_p;
                    }
                }
            }
        }
    }
    list<VocabIndex> best_path;
    best_path.push_front(best_i_j_at_T.first); // T-1
    best_path.push_front(best_i_j_at_T.second); // T-2
    VocabIndex curr_t2 = best_i_j_at_T.first;
    VocabIndex curr_t1 = best_i_j_at_T.second;
    for (int t=T-1; t>=2; --t) {
        pair<unsigned int, pair<VocabIndex, VocabIndex> > key = {t, {curr_t2, curr_t1}};
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
    cout << ngram << endl;
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
            assert(ptr != NULL);
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
    File input_text(input_file, "r");
    File output_text(output_file, "w");
    size_t num_threads = omp_get_max_threads();
    size_t pool_size = 2;
    vector<vector<string> > input_strings;
    vector<list<VocabIndex> > output_sequence(num_threads * pool_size, list<VocabIndex>() );
    while((buffer = input_text.getline()) != NULL) {
        ptr = strtok(buffer, tok);
        vector<string> s;
        while (ptr != NULL) {
            s.push_back(string(ptr));
            ptr = strtok(NULL, tok);
        }
        input_strings.push_back(s);
        if (input_strings.size() >= num_threads * pool_size) {
            #pragma omp parallel for schedule(dynamic) shared(input_strings, output_sequence, mapping, voc, lm)
            for (int i=0; i<input_strings.size(); ++i) {
                output_sequence[i] = viterbi(input_strings[i], mapping, voc, lm);
            }
            for (int i=0; i<input_strings.size(); ++i) {
                fprintf(output_text, "<s>");
                for (auto j : output_sequence[i]) {
                    fprintf(output_text, " %s", voc.getWord(j));
                }
                fprintf(output_text, " </s>\n");
                fflush(output_text);
            }
            input_strings.clear();
        }
        
    }
    if (input_strings.size() > 0) {
        #pragma omp parallel for schedule(dynamic) shared(input_strings, output_sequence, mapping, voc, lm)
        for (int i=0; i<input_strings.size(); ++i) {
            output_sequence[i] = viterbi(input_strings[i], mapping, voc, lm);
        }
        for (int i=0; i<input_strings.size(); ++i) {
            fprintf(output_text, "<s>");
            for (auto j : output_sequence[i]) {
                fprintf(output_text, " %s", voc.getWord(j));
            }
            fprintf(output_text, " </s>\n");
            fflush(output_text);
        }
        input_strings.clear();
    }
    output_sequence.clear();
    output_text.close();
    input_text.close();
    return 0;
}
