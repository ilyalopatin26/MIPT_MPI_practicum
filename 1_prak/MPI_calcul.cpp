#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

#include "parser.h"

double target_function ( double x ) 
{
    return 4 / ( 1 + x * x );
}

long double a = 0.0 , b = 1.0;

int main (int argc, char* argv[])
{
    MPI_Status status;                     

    bool flag_res = true, flag_segment = false;  // -значения по умолчанию
    long long N = 1e3;                           //
    
    std::vector<std::string> Arg(argc-1);
    for (int i = 1; i < argc; ++i )
        Arg[i-1] = argv[i];

    bool er =  parser( Arg, flag_res, flag_segment, N); 


    long double delta = (b-a) / N , begin = 0.0, end = 0.0, loc_summ  =  0, global_summ = 0;
    int myrank, size, error, i;

    error = MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &size );
    MPI_Barrier( MPI_COMM_WORLD );

    if ( myrank == 0 )
        begin = MPI_Wtime();

    
    for ( i = myrank; i < N; i += size)
    {
        long double idx = delta * i + a;
        long double s = delta * ( target_function (idx) + target_function (idx + delta)   )/2;
        loc_summ += s;
    }
    
    
     MPI_Reduce( &loc_summ, &global_summ, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0)
        end = MPI_Wtime();
    
    if ( flag_segment  && ( myrank != 0 ) )
        MPI_Send( &loc_summ , 1 , MPI_LONG_DOUBLE , 0  , 0 , MPI_COMM_WORLD );

    
    if (myrank == 0)
    {
        
        if ( flag_res ) {
            printf("np=%d; Time=%Lf; PI=%Lf\n", size, end-begin, global_summ); 
            if ( flag_segment ) 
            {
                printf("0: %Lf ", loc_summ );
                for ( int n_pros = 1; n_pros < size; ++n_pros ) 
                {
                    MPI_Recv( &loc_summ , 1 , MPI_LONG_DOUBLE , n_pros , 0 , MPI_COMM_WORLD, &status );
                    printf("%d: %Lf ", n_pros ,  loc_summ );
                }
                printf(" \n ");
            }
        }
        else
            printf("%Lf\n",  end-begin);            
    }
    
    error=MPI_Finalize();
    
    return 0;
}