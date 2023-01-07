#include "_test.h"

bool test_min(){
  int a = 5;
  int b = 6;
  int c = 3;

  if (min(a,b) != a)
    return false;
  if (min(b,a) != a)
    return false;
  if (min(c,c) != c)
    return false;
  return true;
}

bool test_bin(){
  unsigned long long a = 200;
  unsigned long long b = 4;

  if ((int)bin(a,b) != 64684950)
    return false;
  if ((int)bin(a,a) != 1)
    return false;
  if ((int)bin(b,a) != 0)
    return false;
  return true;
}

bool test_norm_cdf(){
  float a     = 0.5;
  float b     = 2.5;
  float mu    = 4;
  float sigma = 2;
  float test1 = norm_cdf(a, mu, sigma);
  float test2 = norm_cdf(b, mu, sigma);

  if (test1 <= 0.038 || test1 >= 0.042)
    return false;
  if (test2 <= 0.224 || test2 >= 0.228)
    return false;
  return true;

}

bool test_norm_pf(){
  float a     = 0.5;
  float b     = 2.5;
  float mu    = 4;
  float sigma = 2;
  float test1 = norm_pf(a, mu, sigma);
  float test2 = norm_pf(b, mu, sigma);

  if (test1 <= 0.042 || test1 >= 0.045)
    return false;
  if (test2 <= 0.148 || test2 >= 0.152)
    return false;
  return true;
}

bool test_get_numerator(){

}

bool test_get_denominator(){

}

bool test_get_score(){

}

bool test_get_forward_offset(){

}

bool test_get_reverse_offset(){

}


int main (int argc, char *argv[])
{
  printf("running tests...\n============================\n");
  if (test_bin() == true)
    printf("combinations:         PASSED\n");
  else
    printf("combinations:         FAILED\n");

  if (test_min() == true)
    printf("minimum:              PASSED\n");
  else
    printf("minimum:              FAILED\n");

  if (test_norm_cdf() == true)
    printf("norm cdf:             PASSED\n");
  else
    printf("norm cdf:             FAILED\n");

  if (test_norm_pf() == true)
    printf("norm pf:              PASSED\n");
  else
    printf("norm pf:              FAILED\n");

  return 0;
}

