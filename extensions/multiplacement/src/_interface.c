#include "_interface.h"

static int matrix_converter(PyObject *object, void *address) {
  const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
  char datatype;
  Py_buffer *view = address;

  if (object == NULL)
    goto exit;
  if (PyObject_GetBuffer(object, view, flags) == -1) {
    PyErr_SetString(PyExc_RuntimeError,
                    "position-weight matrix is not an array");
    return 0;
  }
  datatype = view->format[0];
  switch (datatype) {
  case '@':
  case '=':
  case '<':
  case '>':
  case '!':
    datatype = view->format[1];
    break;
  default:
    break;
  }
  return Py_CLEANUP_SUPPORTED;

exit:
  PyBuffer_Release(view);
  return 0;
}

static char calculate__doc__[] =
    "    calculate(sequence, recognizers, recognizer_lengths, connectors, "
    "recognizer_scores, connector_scores, connector_legnths)\n"
    "\n"
    "This function computes optimal placement for a \n"
    "transcription factor composed of PSSM recogniers\n"
    "and variable length connectors.\n";

static PyObject *py_calculate(PyObject *self, PyObject *args,
                              PyObject *keywords) {
  const char *seq;
  static char *kwlist[] = {
      "sequence",   "rec_matrices", "rec_lengths", "con_matrices",
      "rec_scores", "con_scores",   "con_lengths", "max_length", NULL};
  Py_ssize_t len_seq;
  Py_ssize_t num_rec;
  int max_length;
  PyObject *result = Py_None;
  Py_buffer rec_matrices;
  Py_buffer con_matrices;
  Py_buffer rec_lengths;
  Py_buffer rec_scores;
  Py_buffer con_scores;
  Py_buffer con_lengths;

  rec_matrices.obj = NULL;
  con_matrices.obj = NULL;
  rec_lengths.obj = NULL;
  rec_scores.obj = NULL;
  con_scores.obj = NULL;
  con_lengths.obj = NULL;
  if (!PyArg_ParseTupleAndKeywords(
          args, keywords, "y#O&O&O&O&O&O&i", kwlist, &seq, &len_seq,
          matrix_converter, &rec_matrices, matrix_converter, &rec_lengths,
          matrix_converter, &con_matrices, matrix_converter, &rec_scores,
          matrix_converter, &con_scores, matrix_converter, &con_lengths, &max_length))
    return NULL;
  // sequence:     DNA sequence used for placement
  // rec_matrices: one dimensional flattened representation of the
  //               scoring matrices for all recognizers
  // rec_lengths:  length of each recognizer (the number of columns)
  // con_matrices: one dimensional flattened representation of pre-computed
  //               scores for each connector for each gap length
  // rec_scores:   buffer used to store our calculated scores for recognizers
  // con_scores:   buffer used to score our calculated scores for connectors
  // con_lengths:  buffer used to store the length of each connector for our
  // placment
  num_rec = rec_lengths.shape[0];
  float *rec_matrices_ptr = rec_matrices.buf;
  float *con_matrices_ptr = con_matrices.buf;
  float *con_scores_ptr = con_scores.buf;
  float *rec_scores_ptr = rec_scores.buf;
  int *rec_lengths_ptr = rec_lengths.buf;
  int *con_lengths_ptr = con_lengths.buf;
  bool is_precomputed = true;
  // getting the number of alignments is needed for calculating the size
  // of the array we will use to store the scores for each recognizer alignment
  int forward_offset = get_forward_offset(0, rec_lengths_ptr, num_rec);
  int reverse_offset = get_reverse_offset(0, rec_lengths_ptr, num_rec);
  int num_alignments = len_seq - forward_offset - reverse_offset;
  float *score_matrix =
      PyMem_Calloc(num_alignments * num_rec, sizeof(*rec_matrices_ptr));
  fill_matrix(seq, len_seq, rec_matrices_ptr, rec_lengths_ptr, num_rec,
              score_matrix, num_alignments);

  if (con_matrices.shape[0] == (num_rec - 1) * 2)
    is_precomputed = false;
  // traceback function breaks when the number of recognizers is less than
  // two since it opperates on the assumption of having at least one connector
  if (num_rec == 1) {
    con_lengths_ptr[0] =
        max_index(score_matrix, len_seq - forward_offset - reverse_offset);
    rec_scores_ptr[0] = score_matrix[con_lengths_ptr[0]];
    con_scores_ptr[0] = 0.00;
  } else {
    fill_traceback_matrix(score_matrix, num_alignments, con_matrices_ptr, rec_lengths_ptr,
              num_rec, len_seq, con_scores_ptr, rec_scores_ptr,
              con_lengths_ptr, max_length, is_precomputed);
  }

  PyMem_Free(score_matrix);
  Py_INCREF(Py_None);
  result = Py_None;
  matrix_converter(NULL, &rec_matrices);
  matrix_converter(NULL, &rec_lengths);
  matrix_converter(NULL, &con_matrices);
  matrix_converter(NULL, &rec_scores);
  matrix_converter(NULL, &con_scores);
  matrix_converter(NULL, &con_lengths);
  return result;
}

static struct PyMethodDef methods[] = {
    {
        "calculate",
        (PyCFunction)py_calculate,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR(calculate__doc__),
    },
    {NULL, NULL, 0, NULL} // sentinel
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_multiplacement",
    PyDoc_STR("Fast calculations involving multiple connected PSSMs"),
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL};

PyObject *PyInit__multiplacement(void) { return PyModule_Create(&moduledef); }