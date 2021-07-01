#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

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
void sorted_int_vector(u_int64_t *v, u_int64_t n)
{
  for (u_int64_t i = 0; i < n; ++i)
  {
    v[i] = i;
  }
}

// sorted array with n(1-perc) random swaps
void partially_sorted_int_vector(u_int64_t *v, u_int64_t n, float perc)
{
  random_int_vector(v, n);
  qsort(v, n, sizeof(u_int64_t *), compare);
  u_int64_t count = (u_int64_t)n * (1 - perc);
  for (int i = 0; i < count; ++i)
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
  u_int64_t n = 1 << 22;
  // random
  mkdir("random", 0777);
  for (int i = 0; i < 20; i++)
  {
    char *filename = malloc(100 * sizeof(char));
    sprintf(filename, "random/random_%d", i);
    u_int64_t *v = malloc((n + 1) * sizeof(u_int64_t));
    srand(time(NULL));
    random_int_vector(v, n);

    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "size:%ld\n", n);
    for (u_int64_t j = 1; j < n + 1; ++j)
    {
      fprintf(fp, "%ld\n", v[j]);
    }
    fclose(fp);
    free(v);
  }

  // partially_sorted
  // 0.05*n swaps applied to a sorted random sequence
  mkdir("partial", 0777);
  for (int i = 0; i < 20; i++)
  {
    char *filename = malloc(100 * sizeof(char));
    sprintf(filename, "partial/partially_sorted_95_%d", i);

    u_int64_t *v = malloc((n + 1) * sizeof(u_int64_t));
    srand(time(NULL));
    partially_sorted_int_vector(v, n, 0.95);

    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "size:%ld\n", n);
    for (u_int64_t j = 1; j < n + 1; ++j)
    {
      fprintf(fp, "%ld\n", v[j]);
    }
    fclose(fp);
    free(v);
  }

  // sorted
  char *filename = malloc(100 * sizeof(char));
  sprintf(filename, "sorted/sorted");
  mkdir("sorted", 0777);
  u_int64_t *v = malloc((n + 1) * sizeof(u_int64_t));
  sorted_int_vector(v, n);
  FILE *fp = fopen(filename, "w+");
  fprintf(fp, "size:%ld\n", n);
  for (u_int64_t j = 1; j < n + 1; ++j)
  {
    fprintf(fp, "%ld\n", v[j]);
  }
  fclose(fp);
  free(v);

  return 0;
}
