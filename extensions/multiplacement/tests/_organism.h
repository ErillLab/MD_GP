#ifndef _TEST_H
#define _TEST_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jansson.h>
#include "../src/_multiplacement.h"

typedef struct _rec _rec;
struct _rec
{
  float* scores;
  int size;
};

typedef struct _con _con;
struct _con
{
  float mu;
  float sigma;
};

typedef struct _org _org;
struct _org
{
  _rec* recognizers;
  _con* connectors;
  int num_rec;
  int num_con;
};

typedef struct _seq _seq;
struct _seq
{
  char* sequence;
  int length;
};

_org import_org(char* filename);

void set_con(_con *connector, float mu, float sigma);
void print_con(_con *connector);

void set_rec(_rec *recongnizer, int size, float* scores);
void set_rec_size(_rec *recognizer, int size);
void set_rec_scores(_rec *recognizer, float* scores);
void print_rec(_rec *recongnizer);

void set_org(_org *organism, int num_rec, int* sizes, float* scores, float* mus, float* sigmas);
void print_org(_org *organism);

void print_seq(_seq sequence);
#endif
