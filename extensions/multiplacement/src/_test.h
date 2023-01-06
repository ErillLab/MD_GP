#ifndef _TEST_H
#define _TEST_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "_aux.h"

struct _rec
{
  float* scores;
  int size;
};

struct _con
{
  float mu;
  float sigma;
};

struct _org
{
  struct _rec* recognizers;
  struct _con* connectors;
  int num_rec;
  int num_con;
};

struct _seq
{
  char* sequence;
  int length;
};

#endif
