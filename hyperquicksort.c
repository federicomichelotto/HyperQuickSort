#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

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
		// swap
		// after the swap i can be increased by one and j can be decreased by one
		u_int64_t temp = seq[j];
		seq[j--] = seq[i];
		seq[i++] = temp;
	}
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
	// PIVOT selected as the median of a,b and c
	// ^ operator: Bitwise exclusive OR operator
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
	/** P MUST BE A POWER OF TWO**/
	MPI_Init(&argc, &argv);
	int rank, P;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P); // number of processors

	double t0, t1;

	if (rank == 0)
	{
		if (argc != 3 && argc != 4)
		{
			printf("Error: it is required the input path, and the csv filename where to store the time elapsed (and the optional flag -save to store print to file the ordered list).");
			return -1;
		}
		if (!((P != 0) && ((P & (P - 1)) == 0))) // check if P is a power of two
		{
			printf("Error: the number of processors must be a power of 2.\n");
			return -1;
		}
	}

	u_int64_t n, N; // n: length list, N = n + (P - n mod P)
	u_int64_t *v = NULL;
	// only the roots fills the initial array v

	char *filename;
	char *csv_filename;
	int save_flag = 0;

	if (rank == 0) // parse the input file
	{
		filename = malloc(100 * sizeof(char));
		csv_filename = malloc(100 * sizeof(char));
		strcpy(filename, argv[1]);
		strcpy(csv_filename, argv[2]);
		if (argc == 4)
		{
			if (strncmp(argv[3], "-save", 5) == 0)
				save_flag = 1;
		}

		FILE *fp = fopen(filename, "r");
		if (fp == NULL)
			printf("Error: Could not open the file\n");
		char line[30];
		if (fgets(line, sizeof(line), fp) == NULL)
		{
			printf("Error during fgets\n");
		}
		if (strncmp(line, "size:", 5) == 0)
		{
			n = strtoll(&line[5], NULL, 10);
		}
		else
		{
			printf("Error: the first line must be \"size:{n}\" where n is the size of the list. \nThen, each line must contain one non-negative number.\n");
			return -1;
		}

		N = n;
		if (n % P)
			N += (P - n % P); // 0 padding the list in order to have a list size multiple of P

		v = calloc(N, sizeof(u_int64_t));
		for (u_int64_t i = 0; i < n; i++)
		{
			char line[30];
			if (fgets(line, sizeof(line), fp) == NULL)
			{
				printf("Error during fgets\n");
			}
			v[i] = strtoll(line, NULL, 10);
		}
		fclose(fp);
		t0 = MPI_Wtime(); // start counting time (ignore parsing times)
	}

	// broadcast the size of the list
	MPI_Bcast(&N, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

	// each process prepare the buffer
	u_int64_t list_size = N / P;
	u_int64_t *part_v = malloc(list_size * sizeof(u_int64_t));

	// scatter portions of the array
	MPI_Scatter(v, list_size, MPI_UINT64_T, part_v, list_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);

	seq_quicksort(part_v, 0, list_size - 1);

	for (int iter = 0; iter < log2(P); iter++)
	{
		u_int64_t *list_medians;

		u_int64_t dim = P >> iter; // dimension of each subgroup of processors

		// Split the communicator based on the group size
		MPI_Comm row_comm;
		MPI_Comm_split(MPI_COMM_WORLD, rank / dim, rank, &row_comm);
		int col_rank, row_size;
		MPI_Comm_rank(row_comm, &col_rank);
		MPI_Comm_size(row_comm, &row_size);

		if (col_rank == 0)
		{
			list_medians = malloc(dim * sizeof(u_int64_t));
		}

		// all processors send to their master node its median value
		u_int64_t median = 0;
		if (list_size)
			median = part_v[list_size / 2];
		MPI_Gather(&median, 1, MPI_UINT64_T, list_medians, 1, MPI_UINT64_T, 0, row_comm);

		// each master node sort its list of median values
		if (col_rank == 0)
		{
			seq_quicksort(list_medians, 0, dim - 1);
		}

		// broadcast the pivot (median of medians) to all processors that belong to the same group
		u_int64_t pivot;
		if (col_rank == 0) // root processors send the pivot to others processors
		{
			// ignore medians with value 0
			u_int64_t k;
			for (k = 0; k < dim; k++)
			{
				if (list_medians[k])
					break;
			}
			pivot = list_medians[k + (dim - k) / 2];
		}
		MPI_Bcast(&pivot, 1, MPI_UINT64_T, 0, row_comm);

		// each processor counts how many elements there are <= pivot (part_v is sorted => binary search)
		u_int64_t low_size = 0, high_size = 0;

		if (list_size)
		{
			u_int64_t first = 0, last = list_size - 1;
			u_int64_t idx = last / 2;
			while (first < last)
			{
				if (pivot < part_v[idx])
				{
					if (idx == 0)
						break;
					last = idx - 1;
				}
				else
					first = idx + 1;
				idx = first + (last - first) / 2;
			}
			// Three cases:
			// 1) last < first : part_v[0:first] <= pivot -> low_size = first
			// 2) first = last: must be checked if part_v[first] <= pivot
			// 3) first = 0, and last = 1 with part_v[0] > pivot -> low_size = 0 = first
			low_size = first;
			if (first == last && part_v[first] <= pivot) // case 2)
				low_size++;
		}
		high_size = list_size - low_size;

		// exchange low lists and high lists
		u_int64_t new_size, count_rcv;
		u_int64_t *rcv_buffer;
		if (col_rank < dim / 2)
		{
			MPI_Request request_send;
			MPI_Status status_rcv;
			// send the size of the high list to P[rank+dim/2]
			MPI_Send(&high_size, 1, MPI_UINT64_T, col_rank + dim / 2, 0, row_comm);
			// send the high list to P[rank+P/2]
			if (high_size)
				MPI_Isend(&part_v[low_size], high_size, MPI_UINT64_T, col_rank + dim / 2, 1, row_comm, &request_send);

			// receive the size of the low list of P[rank+dim/2]
			MPI_Recv(&count_rcv, 1, MPI_UINT64_T, col_rank + dim / 2, 0, row_comm, &status_rcv);
			new_size = low_size + count_rcv; // P of the new list

			rcv_buffer = malloc(count_rcv * sizeof(u_int64_t));
			// receive the low list from P[rank+dim/2]
			if (count_rcv)
				MPI_Recv(rcv_buffer, count_rcv, MPI_UINT64_T, col_rank + dim / 2, 1, row_comm, &status_rcv);
		}
		else
		{
			MPI_Request request_send;
			MPI_Status status_rcv;
			// send the size of the low list to P[rank-dim/2]
			MPI_Send(&low_size, 1, MPI_UINT64_T, col_rank - dim / 2, 0, row_comm);

			// send the low list to P[rank-dim/2]
			if (low_size)
				MPI_Isend(part_v, low_size, MPI_UINT64_T, col_rank - dim / 2, 1, row_comm, &request_send);
			// receive the size of the high list of P[rank-dim/2]
			MPI_Recv(&count_rcv, 1, MPI_UINT64_T, col_rank - dim / 2, 0, row_comm, &status_rcv);

			new_size = high_size + count_rcv; // P of the new list
			rcv_buffer = malloc(count_rcv * sizeof(u_int64_t));

			// receive the high list from P[rank-P/2]
			if (count_rcv)
				MPI_Recv(rcv_buffer, count_rcv, MPI_UINT64_T, col_rank - dim / 2, 1, row_comm, &status_rcv);
		}

		// merge the local list with the received one
		u_int64_t *new_list = malloc(new_size * sizeof(u_int64_t)); // if (list_size == 0 && count_rcv == 0) -> new_list = NULL (ok)

		// merge the two lists
		u_int64_t i = 0, j = 0;			 // indices of the two lists
		u_int64_t list_size2 = low_size; // size of the (low/high) list in current processor
		if (col_rank >= dim / 2)
		{
			i = low_size;
			list_size2 = list_size;
		}

		for (u_int64_t index = 0; index < new_size; ++index)
		{
			if (i == list_size2) // copy all the remaining elements of rcv_buffer
			{
				new_list[index] = rcv_buffer[j++];
				continue;
			}
			if (j == count_rcv) // copy all the remaining elements of part_v
			{
				new_list[index] = part_v[i++];
				continue;
			}
			if (part_v[i] <= rcv_buffer[j])
				new_list[index] = part_v[i++];
			else
				new_list[index] = rcv_buffer[j++];
		}

		// free memory
		free(part_v);
		free(rcv_buffer);
		if (col_rank == 0)
			free(list_medians);

		// free MPI communicator object
		MPI_Comm_free(&row_comm);

		// let part_v point the new_list of size new_size
		part_v = new_list;
		list_size = new_size;
	}

	// save ordered sublist asynchronously
	u_int64_t **ordered_sublists = malloc(P * sizeof(u_int64_t *));
	u_int64_t *size_sublists = malloc(P * sizeof(u_int64_t));

	MPI_Request requests[P - 1];

	if (rank == 0)
	{
		// ignore padding (n % P zeros)
		ordered_sublists[0] = &part_v[P - n % P];
		size_sublists[0] = list_size - (P - n % P);

		//retrieve the sorted elements from processor 1 to P-1 asynchronously
		for (u_int64_t i = 1; i < P; i++)
		{
			MPI_Status status_rcv;
			MPI_Request request;
			MPI_Recv(&size_sublists[i], 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &status_rcv);
			ordered_sublists[i] = malloc(size_sublists[i] * sizeof(u_int64_t));
			MPI_Irecv(ordered_sublists[i], size_sublists[i], MPI_UINT64_T, i, 1, MPI_COMM_WORLD, &requests[i - 1]);
		}
		MPI_Waitall(P - 1, requests, MPI_STATUSES_IGNORE);
	}
	else
	{
		// send list to the root processor
		MPI_Send(&list_size, 1, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
		MPI_Send(part_v, list_size, MPI_UINT64_T, 0, 1, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		t1 = MPI_Wtime(); // stop counting time (ignore write to file times)

		char *filename_sorted = malloc(100 * sizeof(char));
		strcpy(filename_sorted, filename);
		strncat(filename_sorted, "_sorted", 8);

		if (save_flag) // save to file ordered list
		{
			u_int64_t prev = 0;
			FILE *fp = fopen(filename_sorted, "w+");
			for (u_int64_t i = 0; i < P; ++i)
			{
				for (u_int64_t j = 0; j < size_sublists[i]; ++j)
				{
					fprintf(fp, "%ld\n", ordered_sublists[i][j]);
				}
			}
			fclose(fp);
		}

		// write result to file
		FILE *res_csv = fopen(csv_filename, "a");
		fprintf(res_csv, "%s, %f\n", filename, t1 - t0);
		fclose(res_csv);

		free(filename_sorted);
	}

	free(part_v);
	if (rank == 0)
	{
		free(v);
		for (u_int64_t i = 1; i < P; ++i) // ordered_sublists[0] = part_v (already freed)
		{
			free(ordered_sublists[i]);
		}
		free(ordered_sublists);
		free(size_sublists);
		free(filename);
		free(csv_filename);
		printf("Input list sorted!\n");
		printf("Time required to sort the list (excluding I/O times): %f seconds\n", t1 - t0);
	}

	MPI_Finalize();

	return 0;
}
