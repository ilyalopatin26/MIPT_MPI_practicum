#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <unistd.h>

long double a = 0.0 ,
            b = 1.0;

long int N = 1e6;

double target_function ( double x ) 
{
    return 4 / ( 1 + x * x );
}

long double consistently () 
{
    
    long double summ = 0.0,
                delta = ( b- a ) / N;
    for ( long long i = 0; i < N; ++i ) 
    {
        long double idx = a + i * delta;
        long double ds = delta * ( target_function(idx + delta) + target_function( idx ) ) / 2;
        summ += ds; 
    }
    return summ;
}


int main (int argc, char* argv[])
{
    MPI_Status status;                                            


    // переменные интегрирования
    long double delta = (b-a) / N ,  loc_summ  =  0, global_summ = 0;
    long double* interval =  new long double[2];
    long double* work_interval =  new long double[2];
    
    int myrank, size, error;
    long double begin, end ; // время

    error = MPI_Init (&argc, &argv);
    if ( error != 0 ) {
        printf("error MPI_Init\n");
        exit(0);
    }

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &size );
    MPI_Barrier( MPI_COMM_WORLD );

    if ( myrank == 0 ) 
    {
        begin = MPI_Wtime();
        int k = N / size;
        int tail = N % size;
        long double step = k * delta; 

        for ( int i = 1; i  < size; i++ ) 
        {
            if ( tail == 0 ) 
            {
                interval[0] = a + step * i ;
                interval[1] = interval[0] + step ;
            }
            else 
            {
                if ( i < tail ) 
                {
                    interval[0] =  a + i * ( step + delta ) ; 
                    interval[1] = interval[0] + ( step + delta ) ;
                }
                else 
                {
                    interval[0] = a + tail * ( step + delta) + ( i -tail ) * step;
                    interval[1] = interval[0] + ( step ) ;
                }
            }
            MPI_Send ( &interval[0], 2, MPI_LONG_DOUBLE, i, 0, MPI_COMM_WORLD );
        }
        work_interval[0] = a ;
        work_interval[1] = work_interval[0 ] + step;
        if ( tail != 0 )
            work_interval[1] += delta;   
    }
    
    if ( myrank != 0 ) 
        MPI_Recv( &work_interval[0] , 2 , MPI_LONG_DOUBLE , 0 , 0 , MPI_COMM_WORLD, &status );
    
    // здесь все процессы знают свою обласль интегрирования


    for ( long double idx = work_interval[0]; idx < work_interval[1]; idx += delta ) 
        loc_summ += delta * ( target_function(idx + delta) + target_function( idx ) ) / 2;
    
    if ( myrank != 0)
        MPI_Send ( &loc_summ, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD );


    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0 ) {
        end = MPI_Wtime(); 
        global_summ = loc_summ;
        printf("0: %Lf;\n", loc_summ );
        for ( int n_pros = 1; n_pros < size; ++n_pros ) 
        {
            MPI_Recv( &loc_summ , 1 , MPI_LONG_DOUBLE , n_pros , 0 , MPI_COMM_WORLD, &status );
                printf("%d: %Lf;\n", n_pros ,  loc_summ );
                global_summ += loc_summ;
        }
        printf("result: %Lf\n", global_summ);
        long double cons = consistently();
        printf("consistently: %Lf\n", cons);
        printf("time: %Lf\n",  end-begin);      
    } 

    error=MPI_Finalize();
    
    return 0;
}