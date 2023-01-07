#include "_organism.h"

void set_rec(_rec *recognizer, int size, float* scores)
{
  set_rec_size(recognizer, size);
  recognizer->scores = calloc(recognizer->size * 4, sizeof(float));
  set_rec_scores(recognizer, scores);
}

void set_rec_size(_rec *recognizer, int size)
{
  recognizer->size = size;
}

void set_rec_scores(_rec *recognizer, float* scores)
{
  for (int i = 0; i < recognizer->size; i++){
    for (int j = 0; j < 4; j++){
      *(recognizer->scores + (j + i * 4)) = *(scores + (j + i * 4));
    }
  }
}

void print_rec(_rec *recognizer)
{
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < recognizer->size; j++){
      printf("%1.2f ", *(recognizer->scores + (j * 4 + i)));
    }
    printf("\n");
  }

}

void set_con(_con *connector, float mu, float sigma)
{
  connector->mu = mu;
  connector->sigma = sigma;
}

void print_con(_con *connector)
{
  printf("mu = %2.2f\n", connector->mu);
  printf("sigma = %2.2f\n", connector->sigma);
}

void set_org(_org *organism, int num_rec, int* sizes, float* scores, float* mus, float* sigmas)
{
  organism->num_rec = num_rec;
  organism->num_con = num_rec - 1;
  organism->recognizers = calloc(num_rec, sizeof(_rec));
  organism->connectors = calloc((num_rec - 1), sizeof(_con));

  for (int i = 0; i < organism->num_rec; i++){
    set_rec(organism->recognizers + i, *(sizes + i), scores + (i * sizes[i] * 4));
  }

  for (int i = 0; i < organism->num_con; i++){
    set_con(organism->connectors, *(mus + i), *(sigmas + i));
  }

}

void print_org(_org *organism)
{
  for (int i = 0; i < organism->num_rec; i++){
    printf("\nRecognizer |%i|\n", i);
    print_rec(organism->recognizers + i);
    if (i < organism->num_rec - 1){
      printf("\nConnector |%i|\n", i);
      print_con(organism->connectors + i);
    }
  }
}
