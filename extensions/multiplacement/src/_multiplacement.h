#ifndef _MULTIPLACEMENT_H
#define _MULTIPLACEMENT_H
/*
#define PY_SSIZE_T_CLEAN
#include <python3.11/Python.h>
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "_aux.h"

void traceback(int num_rec, int len_seq, float* con_matrices, float* rec_score_matrix, int num_alignments, float* rec_alignments, int* con_alignments, float* rec_scores, float* con_scores, int* con_lengths, int max_seq_len, int effective_len, bool is_precomputed);

void fill_traceback_matrix(float *score_matrix, int num_alignments, float *gapMatrix, int *cols, int num_rec, int len_seq, float *con_scores, float *rec_scores, int *con_lengths, int max_length, bool is_precomputed);

void fill_matrix(const char seq[], int len_seq, float pssm[], int cols[], int num_rec, float score_matrix[], int num_alignments);

#endif
