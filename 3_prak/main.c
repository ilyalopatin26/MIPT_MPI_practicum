#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <time.h>
#include <semaphore.h>

#define DEBUG 0

sem_t mutex;

typedef struct Thread_Information
{
    pthread_t thread_id;
    double sum;
    double a, b, c, d;
    int num;
    int N;
    double *S;
} ThreadInfo, *pThreadInfo;

int InRegion(double xi, double eta)
{
    return (xi <= M_PI) && (xi >= 0) && (eta <= sin(xi)) && (eta >= 0);
}

void *thread_function(void *arg)
{
    pThreadInfo thread = (pThreadInfo)arg;
    unsigned int xi_k = thread->num, eta_k = thread->num;
    int N = thread->N;
    double a = thread->a;
    double b = thread->b;
    double c = thread->c;
    double d = thread->d;
    int i = 0;
    double xi, eta;
    double sum = 0;
    for (i = 0; i < N; i++)
    {
        xi = (double)rand_r(&eta_k) / RAND_MAX * (b - a) + a;  // float in range a to b
        eta = (double)rand_r(&eta_k) / RAND_MAX * (d - c) + c; // float in range -1 to 1
        if (InRegion(xi, eta) == 1)
        {
            sum += xi * eta;
        }
    }
    sem_wait(&mutex);
    *(thread->S) += sum;
    sem_post(&mutex);
    thread->sum = sum;
    return (void *)thread;
}

int main(int argc, char *argv[])
{
    time_t t;
    double a = 0;
    double b = M_PI;
    double c = 0;
    double d = 1;
    int N = 1e9;
    if (argc != 2)
        return 1;
    int p = atoi(argv[1]);
    int N_per = N / p;
    sem_init(&mutex, 0, 1);
    pThreadInfo *Thread_Info = (pThreadInfo *)malloc(sizeof(pThreadInfo) * p);
    double S;
    int i;
    for (i = 0; i < p; ++i)
    {
        pThreadInfo thread = (pThreadInfo)malloc(sizeof(ThreadInfo));
        Thread_Info[i] = thread;
    }

    struct timespec start, end; // требует использования ключа –lrt  !
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start);

    for (i = 0; i < p; ++i)
    {
        pThreadInfo thread = Thread_Info[i];
        thread->N = N_per;
        thread->a = a;
        thread->b = b;
        thread->c = c;
        thread->d = d;
        thread->num = i;
        thread->S = &S; //= S + i;
        int e = pthread_create(&(Thread_Info[i]->thread_id), NULL, thread_function, Thread_Info[i]);
        if (e == -1)
        {
            printf("Can't create thread :(\n");
            return 0;
        }
    }

    // wait all threads
    for (i = 0; i < p; ++i)
    {
        void *ret;
        pthread_join(Thread_Info[i]->thread_id, &ret);
#if DEBUG
        printf(" %d returned \n", ((pThreadInfo)ret)->num);
#endif
    }
    clock_gettime(CLOCK_REALTIME, &end);
    double parallel_time = end.tv_sec - start.tv_sec;
    parallel_time += (end.tv_nsec - start.tv_nsec) / (double)1e9;

    sem_destroy(&mutex);
    for (i = 0; i < p; ++i)
    {
        free(Thread_Info[i]);
    }
#if DEBUG
    printf("Answer is %f\n", S * (b - a) * (d - c) / N);
    printf("Good answer is %f\n", M_PI * M_PI / 8);
#else
    // by ne process
    int Nout = 0;
    double xi = 0.5, eta = 0.5;
    double sum = 0;

    clock_gettime(CLOCK_REALTIME, &start);
    for (i = 0; i < N; i++)
    {
        xi = (double)rand() / RAND_MAX * (b - a) + a;  // float in range a to b
        eta = (double)rand() / RAND_MAX * (d - c) + c; // float in range -1 to 1
        if (InRegion(xi, eta) == 1)         
        {
            sum += xi * eta;
            Nout += 1;
        }
    }
    clock_gettime(CLOCK_REALTIME, &end);
    double conseq_time = end.tv_sec - start.tv_sec;
    conseq_time += (end.tv_nsec - start.tv_nsec) / (double)1e9;

    printf("%d %f\n", p, conseq_time / parallel_time);
#endif
    return 0;
}