#include <stdio.h>
#include<iostream>
#include <mpi.h>
#include <vector>
#include<string>


double target_function ( double x ) 
{
    return 4 / ( 1 + x * x );
}

void parser (std::vector<std::string> &A,  long long & N ) {
    std::string number;
    for ( auto it: A ) 
    {
        if ( it[0] == '-' && it[1] == 'N' ) 
        {
            for (int i = 2; i < it.size(); ++i )
                number += it[i];  
            N = atoi(number.c_str());
            return;
        }
    }
}


int main ( int argc, char* argv[] ) 
{
    std::ios_base::sync_with_stdio (false);
    std::cin.tie( nullptr );
    std::cout.tie( nullptr );

    int rank, size;
    double start_time, end_time;

    size_t i;

    
    long long N = 1000 ; // значение по умолчанию

    
    
    std::vector<std::string> A(argc-1);
    for (int i = 1; i < argc; ++i ){
        A[i-1] = argv[i];
    }
    

    parser(A, N);
    long double delta =  1.0 / N ;

    long double loc_summ = 0.0;
    long double global_summ = 0.0;

    int error = MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    
    if ( rank == 0 )
        start_time = MPI_Wtime();

    for ( i = rank; i < N ; i += size ) 
    {
        long double idx = delta * i;
        long double s = delta * ( target_function (idx) + target_function (idx + delta)   )/2;
        loc_summ += s;
    }
    // сумматор
    MPI_Reduce(&loc_summ, &global_summ, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    


    if ( rank==0 ) // вывод
    {
        end_time = MPI_Wtime();   
        std::cout.setf(std::ios::fixed);  // вывод в фиксированном формате 
        std::cout.precision(9); 

        //std::cout<<global_summ<<'\n';
        std::cout<< end_time- start_time<<'\n' ;
    }
    
    error=MPI_Finalize();
    
    return 0;

} 