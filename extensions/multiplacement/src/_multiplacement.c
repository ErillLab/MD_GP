#include "_multiplacement.h"

void traceback(int num_rec, int len_seq, float* con_matrices, float* rec_score_matrix, int num_alignments, float* rec_alignments, int* con_alignments, float* rec_scores, float* con_scores, int* con_lengths, int max_seq_len, int effective_len, bool is_precomputed){
  // finding gap lengths to ref orma alignments
  // starting from the index of the greatest score in alignments
  // trace back of gap alignments is conducted by subtracting the value
  // at the max score index from the index and going up a row
  // this is repeated until we know all the gap lengths in the alignment
  int index = max_index(rec_alignments, num_alignments);
  con_lengths[num_rec - 1] = con_alignments[(num_rec - 2) * num_alignments + index];

  // the cumulative best score is written into the last index of the scores
  // array
  rec_scores[num_rec] = rec_alignments[index];

  //the start and end of this for loop are weird because the first and last
  //indices are filled manually, so these bounds only encompass the middle
  //scores
  for (int i = num_rec - 3; i >= 0; i--) {
    con_lengths[i + 1] = con_alignments[i * num_alignments + index - con_lengths[i + 2]];
    index -= con_lengths[i + 2];
  }
  con_lengths[0] = index - con_lengths[1];

  int gapOffset = 0;
  for (int i = 0; i < num_rec - 1; i++) {
    //con_scores[i] = get_score(con_matrices[i * max_seq_len + con_lengths[i + 1]], effective_len, num_rec, con_lengths[i + 1]);
    con_scores[i] = get_score(con_matrices, len_seq, effective_len, num_rec, con_lengths[i + 1], max_seq_len, i, is_precomputed);
  }

  // scores for each PSSM are filled by iterating over the score matrix
  // and using the appropriate gap lengths as a cumulative offset
  for (int i = 0; i < num_rec; i++) {
    gapOffset += con_lengths[i];
    rec_scores[i] = rec_score_matrix[i * num_alignments + gapOffset];
  }

}

void fill_traceback_matrix(float *score_matrix, int num_alignments, float *gapMatrix, int *cols, int num_rec, int len_seq, float *con_scores, float *rec_scores, int *con_lengths, int max_length, bool is_precomputed) {
  int gap_length = 0;
  int effective_length = len_seq;
  int sum_of_lengths = 0;

  for (int i = 0; i < num_rec; i++){
    sum_of_lengths += cols[i];
  }

  effective_length -= sum_of_lengths;
  //  number of total alignments by number of pssms
  //  first index in each column holds current max score for that index
  //  all other indices hold gap lengths that got that alignment
  float *alignments = calloc(num_alignments, sizeof(*score_matrix));
  int *gap_alignments = calloc(num_alignments * (num_rec - 1), sizeof(*con_lengths));
  float *temp_max_scores = calloc(num_alignments, sizeof(*score_matrix));
  int *temp_gap_lengths = calloc(num_alignments, sizeof(*con_lengths));
  float temp_gap_score = 0.0;

  // start with first row as our current max
  for (int i = 0; i < num_alignments; i++) {
    alignments[i] = score_matrix[i];
    temp_max_scores[i] = score_matrix[i];
  }

  // for each connector (number of pssms - 1) populate alignments with
  // current maximum score for that index

  // overview: the algorithm iterates over the next row in the scoresMatrix
  // (first row is our starting cumulative score), for each index, we compare
  // the scores obtained by summing the PSSM score at that index with each
  // previous cumulative alignment score for k <= j, when k == j, the gap length
  // is 0
  for (int i = 1; i < num_rec; i++) {

    for (int j = 0; j < num_alignments; j++) {

      for (int k = 0; k <= j; k++) {

        // every column before or equal to the current column is a valid
        // alignment. Scores for each alignment
        // gap_length = difference between column normalized j and k
        gap_length = j - k;

        //compute on fly
        //temp_gap_score = get_score(gapMatrix[(i - 1) * max_length + gap_length], effective_length, num_rec, gap_length);
        temp_gap_score = get_score(gapMatrix, len_seq, effective_length, num_rec, gap_length, max_length, i - 1, is_precomputed);
        if (k == 0) {
          temp_max_scores[j] = alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j];
          temp_gap_lengths[j] = gap_length;
        }else{
          if (temp_max_scores[j] < alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j]) {
              temp_max_scores[j] = alignments[k] + temp_gap_score + score_matrix[i * num_alignments + j];
              temp_gap_lengths[j] = gap_length;
          }
        }
      }
    }

    // we must reset our temp arrays so that they will be overwritten
    // when doing the comparison of scores
    for (int l = 0; l < num_alignments; l++) {
      alignments[l] = temp_max_scores[l];
      gap_alignments[(i - 1) * num_alignments + l] = temp_gap_lengths[l];
      temp_max_scores[l] = -INFINITY;
      temp_gap_lengths[l] = 0;
    }
  }
  free(temp_max_scores);
  free(temp_gap_lengths);

  traceback(num_rec, len_seq, gapMatrix, score_matrix, num_alignments, alignments, gap_alignments, rec_scores, con_scores, con_lengths, max_length, effective_length, is_precomputed);
  free(alignments);
  free(gap_alignments);
  
}

void fill_matrix(const char seq[], int len_seq, float pssm[], int cols[], int num_rec, float score_matrix[], int num_alignments) {
  // length of the seq by number of pssms

  // printf("last in fill_matrix\n");
  float score = 0;
  int forward_offset = 0;
  int reverse_offset = 0;

  // pre computes alignments of each pssm at each possible position
  // i = current recognizer
  // j = starting position on seq for computing score
  // k = current column in recognizer

  for (int i = 0; i < num_rec; i++) {
    forward_offset = get_forward_offset(i, cols, num_rec);
    reverse_offset = get_reverse_offset(i, cols, num_rec);

    for (int j = forward_offset; j < len_seq - reverse_offset; j++) {
      score = 0;
      for (int k = 0; k < cols[i]; k++) {
        switch (seq[j + k]) {
        case 'A':
        case 'a':
          score += pssm[(forward_offset + k) * 4 + 0];
          break;
        case 'G':
        case 'g':
          score += pssm[(forward_offset + k) * 4 + 1];
          break;
        case 'C':
        case 'c':
          score += pssm[(forward_offset + k) * 4 + 2];
          break;
        case 'T':
        case 't':
          score += pssm[(forward_offset + k) * 4 + 3];
          break;
        }
      }
      score_matrix[(i * num_alignments) + j - forward_offset] = score;
    }
  }
}
