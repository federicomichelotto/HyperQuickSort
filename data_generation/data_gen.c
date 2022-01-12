#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <string.h>

int compare(const void *a, const void *b)
{
  return (*(u_int64_t *)a - *(u_int64_t *)b);
}

void random_int_vector(u_int64_t *v, u_int64_t n)
{
  for (u_int64_t i = 0; i < n; ++i)
  {
    v[i] = (u_int64_t)rand();
  }
}

// sorted array with 0.5*n*(1-perc) random swaps
void partially_sorted_int_vector(u_int64_t *v, u_int64_t n, float perc)
{
  random_int_vector(v, n);
  qsort(v, n, sizeof(u_int64_t *), compare);
  int count = (n - n*perc)*0.5;
  for (u_int64_t i = 0; i < count; ++i)
  {
    u_int64_t j = (u_int64_t)rand() % n;
    u_int64_t k = (u_int64_t)rand() % n;
    u_int64_t temp = v[j];
    v[j] = v[k];
    v[k] = temp;
  }
}

int main(int argc, char **argv)
{

  if (argc == 1)
  {
    printf("Input error: \n");
    printf("first input: [-r | -ps perc] (r: random; ps: generate partially sorted instances with size(1-perc) random swaps) \n");
    printf("second input: size \n");
    printf("third input: number of instances to generate\n");
    return -1;
  }

  // set the seed for the random number generator
  struct timeval time;
  gettimeofday(&time, NULL);
  srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

  char *option = malloc(100 * sizeof(char));
  strcpy(option, argv[1]);

  if (strncmp(option, "-r", 2) == 0) //random
  {
    if (argc != 4)
    {
      printf("Input error: \nfirst input: [-r | -ps perc] (r: random; ps: generate partially sorted instances with size(1-perc) random swaps \nsecond input: size \nthird input: number of instances to generate\n");
      return -1;
    }

    char *size_string = malloc(100 * sizeof(char));
    strcpy(size_string, argv[2]);

    char *n_string = malloc(100 * sizeof(char));
    strcpy(n_string, argv[3]);

    u_int64_t size = strtoll(size_string, NULL, 10);
    u_int64_t n = strtoll(n_string, NULL, 10);

    char *folder_name = malloc(100 * sizeof(char));
    sprintf(folder_name, "random_%ld", size);
    mkdir(folder_name, 0777);
    for (u_int64_t i = 0; i < n; i++)
    {
      char *filename = malloc(100 * sizeof(char));
      sprintf(filename, "%s/random_%ld", folder_name, i);
      u_int64_t *v = malloc(size * sizeof(u_int64_t));

      random_int_vector(v, size);

      FILE *fp = fopen(filename, "w+");
      fprintf(fp, "size:%ld\n", size);
      for (u_int64_t j = 0; j < size; ++j)
      {
        fprintf(fp, "%ld\n", v[j]);
      }
      fclose(fp);
      free(v);
    }
    free(size_string);
    free(n_string);
  }
  else if (strncmp(option, "-ps", 3) == 0) //partially_sorted
  {
    if (argc != 5)
    {
      printf("Input error: \nfirst input: [-r | -ps perc] (r: random; ps: generate partially sorted instances with size(1-perc) random swaps \nsecond input: size \nthird input: number of instances to generate\n");
      return -1;
    }
    char *perc_string = malloc(100 * sizeof(char));
    strcpy(perc_string, argv[2]);

    char *size_string = malloc(100 * sizeof(char));
    strcpy(size_string, argv[3]);

    char *n_string = malloc(100 * sizeof(char));
    strcpy(n_string, argv[4]);

    float perc = strtof(perc_string, NULL);
    u_int64_t size = strtoll(size_string, NULL, 10);
    u_int64_t n = strtoll(n_string, NULL, 10);

    if (perc < 0 || perc > 1)
    {
      printf("Error: perc value must be between 0 and 1.\n\n");
      printf("Input requirements: \nfirst input: [-r | -ps perc] (r: random; ps: generate partially sorted instances with size(1-perc) random swaps \nsecond input: size \nthird input: number of instances to generate\n");
      return -1;
    }

    //0.5*n(1-perc) swaps applied to a sorted random sequence
    char *folder_name = malloc(100 * sizeof(char));
    sprintf(folder_name, "partial_%s_%ld", perc_string + 2, size);
    mkdir(folder_name, 0777);

    for (u_int64_t i = 0; i < n; i++)
    {
      char *filename = malloc(100 * sizeof(char));
      sprintf(filename, "%s/partially_sorted_%ld", folder_name, i);

      u_int64_t *v = malloc(size * sizeof(u_int64_t));
      partially_sorted_int_vector(v, size, perc);

      FILE *fp = fopen(filename, "w+");
      fprintf(fp, "size:%ld\n", size);
      for (u_int64_t j = 0; j < size; ++j)
      {
        fprintf(fp, "%ld\n", v[j]);
      }
      fclose(fp);
      free(v);
    }
    free(perc_string);
    free(size_string);
    free(n_string);
  }
  else
  {
    printf("Input error: \nfirst input: [-r | -ps perc] (r: random; ps: generate partially sorted instances with size(1-perc) random swaps \nsecond input: size \nthird input: number of instances to generate\n");
    return -1;
  }

  free(option);
  return 0;
}
