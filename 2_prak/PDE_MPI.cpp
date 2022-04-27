#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <cmath>

#include "parser.h"

//default
long double h = 0.02 ,
            k = 1.0 , 
            l = 1.0 ,
            u0 = 1.0 ,
            dt = 0.0002 ,
            t_target = 0.1 ,
            pi = 3.14159265358979323846 ; 

#include "exact_sol.h"

long double exch_mess ( 
    const int &k, const long double* left, const long double* right, const int &ph ) 
{
    MPI_Status status;
    long double* get_value = new long double[1];
    
    if ( k % 2 == 0 )
        if ( ph == 1 )
            MPI_Send ( &right[0], 1, MPI_LONG_DOUBLE , k+1 , 0, MPI_COMM_WORLD );
        else
            MPI_Recv( &get_value[0] , 1 , MPI_LONG_DOUBLE , k-1 , 0 , MPI_COMM_WORLD, &status );
    else
        if ( ph == 1 )
            MPI_Recv( &get_value[0] , 1 , MPI_LONG_DOUBLE , k-1 , 0 , MPI_COMM_WORLD, &status );
        else
            MPI_Send ( &right[0], 1, MPI_LONG_DOUBLE , k+1 , 0, MPI_COMM_WORLD );

    if ( k % 2 == 0 )
        if ( ph == 1 )
            MPI_Recv( &get_value[0] , 1 , MPI_LONG_DOUBLE , k+1 , 0 , MPI_COMM_WORLD, &status );    
        else
            MPI_Send ( &left[0], 1, MPI_LONG_DOUBLE , k-1 , 0, MPI_COMM_WORLD );
    else
        if ( ph == 1 )
            MPI_Send ( &left[0], 1, MPI_LONG_DOUBLE , k-1 , 0, MPI_COMM_WORLD );
        else
            MPI_Recv( &get_value[0] , 1 , MPI_LONG_DOUBLE , k+1 , 0 , MPI_COMM_WORLD, &status );

    long double ans = get_value[0];
    delete [] get_value;
    return ans;
}

int main ( int argc, char* argv[] )
{
    MPI_Status status;   
    int myrank, size, error;
    long double beginT, endT;

    bool flag_time = false;
    
    long int Nnode = 51 ;

    //обработка флагов
    std::vector<std::string> Arg(argc-1);
    for (int i = 1; i < argc; ++i )
        Arg[i-1] = argv[i];
    bool er =  parser( Arg, flag_time ); 
    if ( !er ) {
        std::cerr<<"error parser" << '\n';
        MPI_Finalize();
        return 0 ;
    }

    std::ifstream IN ("./IN.txt");
    std::string skip_l ,
                skip_k ,
                skip_u ,
                skip_T ,
                skip_N ,
                skip_dt ;
    
    if (IN.is_open())
    {
        IN >> l >> skip_l >> k >> skip_k >> u0 >> skip_u ;
        IN >> t_target >> skip_T >> Nnode >> skip_N >> dt >> skip_dt ;
    }
    else 
    {
        std::cerr<<"Can not open IN.txt"<<'\n';
        MPI_Finalize();
        return 0 ;
    }
    IN.close(); 

    long double h = l / ( Nnode - 1 );

    long double kr = (dt * k) / ( h * h );
    std::vector<long double> ans;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &size );
    MPI_Barrier( MPI_COMM_WORLD );
    if ( myrank == 0 )
        beginT = MPI_Wtime();

    if (  dt > ( h * h ) / ( 2 * k ) ) 
    {
        dt =  ( h * h ) / ( 2 * k ) ;
        kr = (dt * k) / ( h * h );
        if ( myrank == 0 && !flag_time )
            printf("the Courant condition is not satisfied \n");
    }


    long double t_curr = 0;

    long int my_quant = Nnode / size;
    long int rem = Nnode % size;
    if (  (rem > 0) && ( myrank < rem ) )
        ++my_quant; 

    long double* U = new long double [my_quant ];
    for ( size_t i = 0; i < my_quant; ++i ) 
        U[i] = u0;
    long double* U_n = new long double [my_quant];

    if ( myrank == 0 )
        U[0] = 0;
    if ( myrank == size - 1 )
        U[my_quant - 1] = 0;    

    while ( t_curr < t_target  ) 
    {
        long double l, r;
        
        if ( myrank % 2 != 0 )
            l = exch_mess( myrank, &U[0], &U[my_quant-1] , 1 );
        if ( ( myrank % 2 == 0 ) && ( myrank != size - 1 )  )
            r = exch_mess( myrank, &U[0], &U[my_quant-1] , 1 );
        
        MPI_Barrier( MPI_COMM_WORLD );
    
        if (  ( myrank % 2 == 0 )  && ( myrank != 0 )  )
            l = exch_mess( myrank, &U[0], &U[my_quant-1] , 2);
        if ( (myrank % 2 != 0) && (  myrank + 1 < size )  )
            r = exch_mess( myrank, &U[0], &U[my_quant-1] , 2 );    

        for ( int i = 1 ; i < my_quant - 1 ; i++ )
            U_n[i] = kr * ( U[i+1] + U[i-1] - 2 * U[i] ) ;

        if ( myrank != 0 )
            U_n[0] = kr * ( U[1] + l - 2 * U[0] ) ;
        else
            U_n[0] = 0 ;   

        if ( myrank != size - 1 )    
            U_n[ my_quant-1 ] = kr * ( U[my_quant-2] + r - 2 * U[ my_quant-1 ]  ) ;
        else
            U_n[ my_quant-1 ] = 0;
        
        for ( int i = 0; i < my_quant; ++i)
            U[i] += U_n[i];

        t_curr += dt;
        MPI_Barrier( MPI_COMM_WORLD );
    } 

    if ( myrank > 0  )
        MPI_Send( &U[0] , my_quant  , MPI_LONG_DOUBLE , 0 , 0 , MPI_COMM_WORLD );

    if ( myrank == 0) {
        for ( int i = 0; i < my_quant ; ++i )
            ans.push_back( U[i] );
        for ( int n = 1; n < size ; ++n ){
            int quant_k = Nnode / size ;
            if (  (rem > 0) && ( n < rem ) )
                ++quant_k;
            long double* mes = new long double[quant_k];     
            MPI_Recv( &mes[0] , quant_k , MPI_LONG_DOUBLE , n , 0 , MPI_COMM_WORLD, &status );
            for ( int j = 0; j < quant_k ; ++ j)
                ans.push_back (mes[j]);   
            delete [] mes; 
        }
        endT = MPI_Wtime(); 
       
       if ( !flag_time ) {
            for ( int i = 0; i <= 50 ; i += 5 ) {
                long double temp = get_exact_sol ( i * h, t_target, 1e-7 );
                printf("%Lf, %Lf \n", ans[i], temp );
            }
            printf("time: %Lf \n", (endT - beginT) ); 
        }
        else 
            printf("%Lf \n", (endT - beginT) );
    }    
    
    delete [] U;
    delete [] U_n; 

    MPI_Finalize();

    return 0;
}