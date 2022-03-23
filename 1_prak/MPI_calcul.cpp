#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <unistd.h>

#include "parser.h"

long double a = 0.0 , b = 1.0;

double target_function ( double x ) 
{
    return 4 / ( 1 + x * x );
}

long double consistently (const long long &N) 
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
    // -значения по умолчанию
    bool flag_res = true, 
        flag_segment = false, 
        flag_consistently =false ;  
    long long N = 1e3 ;                          
    
    //обработка флагов
    std::vector<std::string> Arg(argc-1);
    for (int i = 1; i < argc; ++i )
        Arg[i-1] = argv[i];
    bool er =  parser( Arg, flag_res, flag_segment, flag_consistently ,N); 
    if ( !er ) {
        printf("error parser\n");
        exit(0);
    }
    // переменные интегрирования
    long double delta = (b-a) / N , begin = 0.0, end = 0.0, loc_summ  =  0, global_summ = 0;
    int myrank, size, error;

    error = MPI_Init (&argc, &argv);
    if ( error != 0 ) {
        printf("error MPI_Init\n");
        exit(0);
    }

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &size );
    MPI_Barrier( MPI_COMM_WORLD );

    if ( flag_segment && ( size != 1 ) ) // задержка для корректной печати
        sleep(0.1);

    if ( myrank == 0 )
        begin = MPI_Wtime();

    for ( int i = myrank; i < N; i += size )
    {
        long double idx = delta * i + a ;
        long double s = delta * (  target_function ( idx ) + target_function ( idx + delta ) ) / 2 ;
        loc_summ += s ;
    }
    
    if ( flag_segment  && ( myrank != 0 ) )
        MPI_Send( &loc_summ , 1 , MPI_LONG_DOUBLE , 0  , 0 , MPI_COMM_WORLD );


    MPI_Reduce( &loc_summ, &global_summ, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0 )
        end = MPI_Wtime();
    


    if ( myrank == 0 )
    {        
        if ( flag_res ) {
            if ( flag_consistently ) {
                long double res = consistently( N );
                printf("consistently: %Lf\n", res);
            }
            printf("np=%d; Time=%Lf; Result=%Lf\n", size, end-begin, global_summ); 
            if ( flag_segment ) 
            {
                printf("partial sum:\n");
                printf("0: %Lf;\n", loc_summ );
                for ( int n_pros = 1; n_pros < size; ++n_pros ) 
                {
                    MPI_Recv( &loc_summ , 1 , MPI_LONG_DOUBLE , n_pros , 0 , MPI_COMM_WORLD, &status );
                    printf("%d: %Lf;\n", n_pros ,  loc_summ );
                }
                printf("\n");
            }
        }
        else
            printf("%Lf\n",  end-begin);            
    }
    
    error=MPI_Finalize();
    
    return 0;
}