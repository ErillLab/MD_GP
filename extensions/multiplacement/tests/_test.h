#include "_organism.h"
/*
int min(int n, int k);
unsigned long long bin(unsigned long long n, unsigned long long k);
float norm_cdf(float x, float mu, float sigma);
float norm_pf(float x, float mu, float sigma);
float get_numerator(int dna_length, int distance, float mu, float sigma);
float get_denominator(int d, int N, int L);
float get_score(float *arr, int dna_length, int effective_length, int num_rec, int gap_size, int max_length, int curr_conn, bool is_precomputed);
int get_forward_offset(int index, int cols[], int num_rec);
int get_reverse_offset(int index, int cols[], int num_rec);
int max_index(float *arr, int size);
*/

bool test_min();
bool test_bin();
bool test_norm_cdf();
bool test_norm_pf();
bool test_get_numerator();
bool test_get_denominator();
bool test_get_score();
bool test_get_forward_offset();
bool test_get_reverse_offset(); 

float random_scores[100] = {0.00};
float random_mus[4] = {0.00, 0.00, 0.00, 0.00};
float random_sigmas[4] = {0.00, 0.00, 0.00, 0.00};
int random_sizes[5] = {4, 3, 2, 4, 5};

