#include "mpi.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_INT 300024
#define ARRAY_SIZE 30000

void print_int_array(int *array, int size) {
  printf("[");
  for (int i = 0; i < size - 1; i++) {
    printf("%d, ", array[i]);
  }
  printf("%d]\n", array[size - 1]);
}

void simple_sort(int *array, int size) {
  int sorted = 0;
  for (int i = 0; i < size - 1; i++) {
    sorted = 0;
    for (int j = 0; j < size - i - 1; j++) {
      if (array[j] > array[j + 1]) {
        int tmp = array[j];
        array[j] = array[j + 1];
        array[j + 1] = tmp;
        sorted = 1;
      }
    }
    if (sorted == 0)
      return;
  }
}
/*
int *DiagonalIntersection(int *merge_matrix, int diag_number, int size,
                          int rank) {
  int *coord = malloc(sizeof(int) * 2);
  if (diag_number + 1 == ARRAY_SIZE * 2) {
    coord[0] = ARRAY_SIZE - 1;
    coord[1] = ARRAY_SIZE - 1;
    return coord;
  }
  int chunk_size = ARRAY_SIZE / size;
  int i = 0;
  int a = diag_number;
  int b = 0;
  if (diag_number >= ARRAY_SIZE) {
    a = ARRAY_SIZE - 1;
    b = diag_number - (ARRAY_SIZE - 1);
  }

  while (merge_matrix[(a - i) * ARRAY_SIZE + b + i] == 1) {
    if (((a - (i + 1)) < 0) || ((b + i + 1) >= ARRAY_SIZE) ||
        (merge_matrix[(a - (i + 1)) * ARRAY_SIZE + b + i + 1] == 0)) {
      break;
    } else {
      i++;
    }
  }
  coord[0] = a - i;
  coord[1] = b + i;
  return coord;
}*/

int main(int argc, char *argv[]) {

  printf("main... \n");

  int rank, size;
  int merge_matrix[ARRAY_SIZE * ARRAY_SIZE];
  int array1 = malloc(sizeof(int) * ARRAY_SIZE);
  int array2 = malloc(sizeof(int) * ARRAY_SIZE);
  int result = malloc(sizeof(int) * ARRAY_SIZE * 2);
  //
  // MPI_INIT
  //
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm custom_comm;
  MPI_Group group_world;
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);
  MPI_Comm_create(MPI_COMM_WORLD, group_world, &custom_comm);

  int chunk_size = ARRAY_SIZE / size;
  int local_array[chunk_size];
  //
  // SEQUENTIAL INITIALIZATION...
  //
  if (rank == 0) {
    // printf("filling arrays ... \n");

    for (int i = 0; i < ARRAY_SIZE; i++) {
      array1[i] = rand() % MAX_INT;
      array2[i] = rand() % MAX_INT;
    }

    // print_int_array(array1, ARRAY_SIZE);
    // print_int_array(array2, ARRAY_SIZE);

    simple_sort(array1, ARRAY_SIZE);
    simple_sort(array2, ARRAY_SIZE);
    /**DEBUG PRINTING
    printf("arrays filled \nsorting...\n");
    print_int_array(array1, ARRAY_SIZE);
    print_int_array(array2, ARRAY_SIZE);
    */
  }

  //
  // PARALLEL WORK
  //
  MPI_Barrier(custom_comm);
  double t1, t2;
  t1 = MPI_Wtime();

  // rank 0 broadcast array2 and scatter array1 to custom_comm
  MPI_Bcast(array2, ARRAY_SIZE, MPI_INT, 0, custom_comm);
  MPI_Scatter(array1, chunk_size, MPI_INT, local_array, chunk_size, MPI_INT, 0,
              custom_comm);

  // filling merge_matrix
  for (int i = 0; i < chunk_size; i++) {
    for (int j = 0; j < ARRAY_SIZE; j++) {
      int offset = (rank * chunk_size + i) * ARRAY_SIZE + j;
      if (local_array[i] >= array2[j]) {
        merge_matrix[offset] = 1;
      } else {
        merge_matrix[offset] = 0;
      }
    }
  }

  if (rank == 0) {
    t2 = MPI_Wtime();
    printf("before Gather in %f\n", t2 - t1);
  }
  MPI_Gather(&merge_matrix[rank * ARRAY_SIZE * chunk_size],
             ARRAY_SIZE * chunk_size, MPI_INT, merge_matrix,
             ARRAY_SIZE * chunk_size, MPI_INT, 0, custom_comm);
  if (rank == 0) {
    t2 = MPI_Wtime();
    printf("Gather in %f\n", t2 - t1);
  }
  /** DEBUG PRINTING
  if (rank == 0) {
    printf("\0nPRINTING MERGE_MATRIX \n");
    for (int i = 0; i < ARRAY_SIZE; i++) {
      for (int j = 0; j < ARRAY_SIZE; j++) {
        int offset = i * ARRAY_SIZE + j;
        printf("%d ", merge_matrix[offset]);
      }
      printf("\n");
    }
  }*/

  int diag_number = ((rank + 1) * (ARRAY_SIZE + ARRAY_SIZE) / size) - 1;
  int length = (ARRAY_SIZE + ARRAY_SIZE) / size;

  // Compute Merge Path
  int path[ARRAY_SIZE * 2 - 1];
  int localpath[ARRAY_SIZE / 2];
  if (rank == 0) {
    int index = 0;
    int nb1 = 0;
    int nb2 = 0;
    if (merge_matrix[0] == 1) {
      path[index] = array2[0];
      nb2++;
    } else {
      path[index] = array1[0];
      nb1++;
    }
    index++;
    int i = 0;
    int j = 0;
    while (index < ARRAY_SIZE * 2) {
      if (j == ARRAY_SIZE - 1) {
        j = 0;
        i++;
      }
      if (merge_matrix[i * ARRAY_SIZE + j + 1] == 0) {
        path[index] = array1[nb1];
        index++;
        nb1++;
        i++;
      } else {
        path[index] = array2[nb2];
        index++;
        nb2++;
        j++;
      }
    }
    // print_int_array(path, ARRAY_SIZE * 2);
  }

  if (rank == 0) {
    t2 = MPI_Wtime();
    printf("Done in %f\n", t2 - t1);
  }
  MPI_Finalize();
  return 0;
}
