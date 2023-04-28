#include <stdio.h>
#include "integral.h"
#include <mpi/mpi.h>
#include <pthread.h>

#define N 12
double I_inter = 0.0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
typedef struct ags {
    int n;
    int n_t;
    int k;
    int p;
    int l_p;
    int f_p;
} arguments;

void first_sigma_sharing(int, int, int, int, int);

void *second_sigma_sharing(void *);

int main(int argc, char *argv[]) {
    int rank, size;
    int n1 = 300, n2 = 200, n3 = 400;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    first_sigma_sharing(rank, size, n1, n2, n3);
    MPI_Finalize();
    return 0;
}

void first_sigma_sharing(int rank, int size, int n1, int n2, int n3) {
    if (rank != 0) {
        int k = rank - 1;
        int first_process = k * n2 / size;
        int last_process = (n2 * (k + 1)) / size + 1;
        pthread_t threads[N];
        arguments args[N];
        for (int i = 0; i < N; i++) {
            args[i].n = N;
            args[i].k = i;
            args[i].l_p = last_process;
            args[i].f_p = first_process;
            args[i].p = n2;
            args[i].n_t = n3;
        }
        for (int i = 0; i < N; i++) {
            pthread_create(&threads[i], 0, second_sigma_sharing, &args[i]);
        }
        double inter_result = 0.0;
        double *temp;
        for (int i = 0; i < N; i++) {
            pthread_join(threads[i], (void **) &temp);
            inter_result += *temp;
        }
        MPI_Reduce(&inter_result, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        double result;
        if (size == 1) {
            result = integral_3_f(0, n1, 0, n2, 0, n3);
        } else
            MPI_Reduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        fprintf(stderr, "%.3f\n", result);
    }
}

void *second_sigma_sharing(void *args) {
    arguments *a = args;
    pthread_mutex_lock(&mutex);
    int first = (a->p * a->k) / a->n;
    int last = (a->p * (a->k + 1)) / a->n;
    I_inter += integral_3_f(a->f_p, a->l_p, first, last, 0, a->n_t);
    pthread_mutex_unlock(&mutex);
    return &I_inter;
}