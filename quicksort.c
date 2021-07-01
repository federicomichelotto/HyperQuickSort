#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

void random_int_vector(u_int64_t *v, u_int64_t n)
{
  for (u_int64_t i = 0; i < n; ++i)
  {
    v[i] = (u_int64_t)rand();
  }
}

u_int64_t partition(u_int64_t *seq, u_int64_t low, u_int64_t high, u_int64_t pivot)
{
  // create the two partitions
  u_int64_t i = low;
  u_int64_t j = high;
  while (1)
  {
    while (seq[i] < pivot)
      i++;

    while (seq[j] > pivot)
      j--;

    if (i >= j)
      break;
    if (seq[i] == seq[j]) // seq[i] == pivot && seq[j] == pivot
    {
      i++;
      j--;
      continue;
    }
    //printf("swap\n");
    // swap
    // after the swap i and j can be increased and decreased j by one
    u_int64_t temp = seq[j];
    seq[j--] = seq[i];
    seq[i++] = temp;
  }
  // DEBUG
  // printf("sep = %ld\n", j);
  if (j == high)
    return high - 1;
  return j;
}

void seq_quicksort(u_int64_t *seq, u_int64_t low, u_int64_t high)
{
  // base step
  if (high == low)
    return;
  if (high == low + 1)
  {
    if (seq[low] > seq[high])
    {
      u_int64_t temp = seq[high];
      seq[high] = seq[low];
      seq[low] = temp;
    }
    return;
  }
  u_int64_t pivot;
  u_int64_t a = seq[low];
  u_int64_t b = seq[(u_int64_t)((high + low) / 2)];
  u_int64_t c = seq[high];
  if ((a > b) ^ (a > c))
    pivot = a;
  else if ((b < a) ^ (b < c))
    pivot = b;
  else
    pivot = c;
  u_int64_t p = partition(seq, low, high, pivot);
  seq_quicksort(seq, low, p);
  seq_quicksort(seq, p + 1, high);
}

int main(int argc, char **argv)
{

  double t0, t1;
  t0 = MPI_Wtime();

  /** P MUST BE A POWER OF TWO**/
  MPI_Init(&argc, &argv);
  int rank, P;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &P); // number of processors

  if (argc != 2)
  {
    printf("Error: input path required.");
    return -1;
  }

  u_int64_t n;
  u_int64_t *v = NULL;
  // only the roots fills the initial array v

  char *filename = malloc(100 * sizeof(char));
  strcpy(filename, argv[1]);

  if (rank == 0) // parse the input file
  {

    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
      printf("Error: Could not open the file\n");
    char line[30];
    fgets(line, sizeof(line), fp);
    line[strcspn(line, " \n\r\t")] = '\0'; // removing trailing white spaces
    if (strncmp(line, "size:", 5) == 0)
      n = strtoll(&line[5], NULL, 10);
    else
    {
      printf("Error: size not found in input file.\n");
      return -1;
    }
    if (!((n != 0) && ((n & (n - 1)) == 0))) // check if n is a power of two
    { 
      printf("Error: the size of the input list must be a power of 2.\n");
      return -1;
    }
    v = malloc(n * sizeof(u_int64_t));
    u_int64_t i = 0;
    while (fgets(line, sizeof(line), fp) != NULL)
    {
      line[strcspn(line, " \n\r\t")] = '\0'; // removing trailing white spaces
      if (strncmp(line, "\0", 1) == 0)
        break;
      v[i++] = strtoll(line, NULL, 10);
    }
    if (i!=n){
      printf("Error: mismatch between the size declared and the numbers of values detected.\n");
      return -1;
    }
    fclose(fp);
  }
  
  // broadcast the size of the list
  MPI_Bcast(&n, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  // each process prepare the buffer
  // we assume for convenvenience that P divides n
  u_int64_t list_size = n / P;
  u_int64_t *part_v = malloc(list_size * sizeof(u_int64_t));

  u_int64_t *list_medians;

  int iter = 1;
  u_int64_t dim = P / iter; // dimension of each subgroup of processors

  // scatter portions of the array
  MPI_Scatter(v, list_size, MPI_UINT64_T, part_v, list_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  seq_quicksort(part_v, 0, list_size - 1);

  // if (list_size < 2)
  //   printf("*** (%ld)[%ld] list_size < 2 ***\n", iter, rank);
  // printf("(%ld)[%ld] median: %ld\n", rank, iter, part_v[list_size / 2]);

  for (iter = 1; iter <= log2(P); iter++)
  {
    dim = P / iter;

    // Split the communicator based on the group size
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, (u_int64_t)rank / dim, rank, &row_comm);
    int col_rank, row_size;
    MPI_Comm_rank(row_comm, &col_rank);
    MPI_Comm_size(row_comm, &row_size);

    if (col_rank == 0)
    {
      // printf("### dim = %ld, rank = %ld , col_rank = %ld\n", dim, rank, col_rank);
      list_medians = malloc(dim * sizeof(u_int64_t));
    }

    // printf("(%d)[%d] part_v: ", rank, iter);
    // for (int i = 0; i < list_size; ++i)
    // {
    //   printf("%ld, ", part_v[i]);
    // }
    // printf("\n");

    // printf("(%ld)[%ld] send median %ld to %ld \n", rank, iter, part_v[list_size / 2], rank);

    MPI_Gather(&part_v[list_size / 2], 1, MPI_UINT64_T, list_medians, 1, MPI_UINT64_T, 0, row_comm);

    if (col_rank == 0)
    {
      // printf("### rank: %ld\n", rank);
      seq_quicksort(list_medians, 0, dim - 1);
      // printf("(%ld)[%ld] list_medians(ordered): %ld", rank, iter, list_medians[0]);

      // for (u_int64_t i = 1; i < dim; i++)
      // {
      //   printf(",%ld", list_medians[i]);
      // }
      // printf("; median of medians: %ld\n", list_medians[dim / 2]);
    }

    u_int64_t pivot = -1;
    // broadcast the median of medians to all processors that belong to the same group
    if (col_rank == 0) // root processors send the pivot to others processors
    {
      pivot = list_medians[dim / 2];
      for (u_int64_t i = 1; i < dim; i++)
      {
        MPI_Send(&pivot, 1, MPI_UINT64_T, i, 0, row_comm);
      }
    }
    else
    { // non root processors receive the pivot
      MPI_Status status_rcv;
      MPI_Recv(&pivot, 1, MPI_UINT64_T, 0, 0, row_comm, &status_rcv);
    }
    // printf("(%ld)[%ld] pivot = %ld\n", rank, iter, pivot);

    // count how many elements there are <= pivot
    u_int64_t low_size = 0;
    // printf("(%ld)[%ld] low list: ", rank, iter);
    while (low_size < list_size && part_v[low_size] <= pivot)
    {
      // printf("%ld,", part_v[low_size]);
      low_size++;
    }
    // printf("\n");
    u_int64_t high_size = list_size - low_size;
    // printf("(%ld)[%ld] low_size = %ld, hi_size = %ld\n", rank, iter, low_size, high_size);

    u_int64_t new_size, count_rcv;
    u_int64_t *rcv_buffer;
    if (col_rank < dim / 2)
    {
      MPI_Request request_send;
      MPI_Status status_rcv;
      // send the size of the high list to P[rank+dim/2]
      MPI_Send(&high_size, 1, MPI_UINT64_T, col_rank + dim / 2, 0, row_comm);
      // send the high list to P[rank+P/2]
      MPI_Isend(&part_v[low_size], high_size, MPI_UINT64_T, col_rank + dim / 2, 1, row_comm, &request_send);

      // receive the size of the low list of P[rank+dim/2]
      MPI_Recv(&count_rcv, 1, MPI_UINT64_T, col_rank + dim / 2, 0, row_comm, &status_rcv);
      // printf("(%ld)[%ld] count_rcv = %ld\n", rank, iter, count_rcv);
      new_size = low_size + count_rcv; // P of the new list

      rcv_buffer = malloc(count_rcv * sizeof(u_int64_t));
      // receive the low list from P[rank+dim/2]
      MPI_Recv(rcv_buffer, count_rcv, MPI_UINT64_T, col_rank + dim / 2, 1, row_comm, &status_rcv);

      // printf("(%ld)[%ld] rcv_buffer: ", rank, iter);
      // for (u_int64_t i = 0; i < count_rcv; ++i)
      //   printf("%ld,", rcv_buffer[i]);
    }
    else
    {
      MPI_Request request_send;
      MPI_Status status_rcv;
      // send the size of the low list to P[rank-dim/2]
      MPI_Send(&low_size, 1, MPI_UINT64_T, col_rank - dim / 2, 0, row_comm);

      //if (index > 0)

      // send the low list to P[rank-dim/2]
      MPI_Isend(part_v, low_size, MPI_UINT64_T, col_rank - dim / 2, 1, row_comm, &request_send);
      // receive the size of the high list of P[rank-dim/2]
      MPI_Recv(&count_rcv, 1, MPI_UINT64_T, col_rank - dim / 2, 0, row_comm, &status_rcv);
      // printf("(%ld)[%ld] count_rcv = %ld\n", rank, iter, count_rcv);

      new_size = high_size + count_rcv; // P of the new list
      rcv_buffer = malloc(count_rcv * sizeof(u_int64_t));

      // receive the high list from P[rank-P/2]
      MPI_Recv(rcv_buffer, count_rcv, MPI_UINT64_T, col_rank - P / 2, 1, row_comm, &status_rcv);

      // printf("(%ld)[%ld] rcv_buffer:", rank, iter);
      // for (u_int64_t i = 0; i < count_rcv; ++i)
      //   printf("%ld,", rcv_buffer[i]);
      // printf("\n");
    }

    // merge the local list with the received one

    u_int64_t *new_list = malloc(new_size * sizeof(u_int64_t));
    u_int64_t start = 0; // index where to start in the local list
    u_int64_t end = low_size;
    if (col_rank >= dim / 2)
    {
      start = low_size;
      end = high_size;
    }
    u_int64_t i = 0, j = 0; // indices of the two lists
    for (u_int64_t index = 0; index < new_size; ++index)
    {
      if (i == end) // reached the end of the local list
      {
        new_list[index] = rcv_buffer[j++];
        continue;
      }
      if (j == count_rcv) // reached the end of the received list
      {
        new_list[index] = part_v[start + i++];
        continue;
      }
      if (part_v[start + i] <= rcv_buffer[j])
        new_list[index] = part_v[start + i++];
      else
        new_list[index] = rcv_buffer[j++];
    }
    // printf("(%d)[%d] new_size = %ld\n", rank, iter, new_size);

    // printf("(%d)[%d] new_list: ", rank, iter);
    // for (u_int64_t i = 0; i < new_size; ++i)
    //   printf("%ld,", new_list[i]);
    // printf("\n");

    free(part_v);
    free(rcv_buffer);
    if (col_rank == 0)
      free(list_medians);
    part_v = new_list;
    list_size = new_size;
    printf("(%d)[%d] list size = %ld\n", rank, iter, list_size);
  }

  // save result to output file [filename]_sorted
  if (rank == 0)
  {
    strncat(filename, "_sorted", 8);
    FILE *fp = fopen(filename, "w+");
    for (u_int64_t i = 0; i < list_size; ++i)
    {
      fprintf(fp, "%ld\n", part_v[i]);
      //printf("%ld, ", part_v[i]);
    }
    // retrieve the sorted elements from processor 1 to P-1
    for (u_int64_t i = 1; i < P; i++)
    {
      u_int64_t rcv_size = 0;
      MPI_Status status_rcv;
      MPI_Recv(&rcv_size, 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &status_rcv);
      u_int64_t *rcv_buff = malloc(rcv_size * sizeof(u_int64_t));
      MPI_Recv(rcv_buff, rcv_size, MPI_UINT64_T, i, 1, MPI_COMM_WORLD, &status_rcv);
      for (u_int64_t j = 0; j < rcv_size; ++j)
      {
        fprintf(fp, "%ld\n", rcv_buff[j]);
        //printf("%ld, ", rcv_buff[j]);
      }
      free(rcv_buff);
    }
    fclose(fp);
  }
  else
  { // send list to the root processor
    MPI_Send(&list_size, 1, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
    MPI_Send(part_v, list_size, MPI_UINT64_T, 0, 1, MPI_COMM_WORLD);
  }

  free(filename);
  if (rank == 0)
    free(v);
  else
    free(part_v);
  MPI_Finalize();

  t1 = MPI_Wtime();
  printf("Time elapsed: %f\n", t1 - t0);
  return 0;
}
