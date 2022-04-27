unsigned int get_limit_sum ( const long double &eps ) 
{
    long double temp = pow(u0, 1.0 / 6 ) / ( pow( k * t_target, 1.0 / 2 ) * pow (eps, 1.0 / 6 ) );
    temp *=  (l / ( pi * pow( 3*pi , 1.0 / 6 ) ) );
    temp = ( 1 + temp ) / 2;
    long double m = ( l * sqrt( 5 / (k* t_target ) ) / pi - 1 ) / 2.0 ;
    if ( temp >= m )
        return ceil( temp );
    else
        return ceil ( m );    
}

long double get_exact_sol ( 
    const long double &x, const long double &t , const long double &eps ) 
{
    if (   x <= 1e-8 ||  l - x <= 1e-8   )
        return 0;
    long double sum = 0; 
    unsigned int lim_sum =  get_limit_sum ( eps );
    for ( unsigned int m = 0; m <= lim_sum ; ++m ) {
        long double a_m = exp( - k * pow(pi, 2) * pow(2*m + 1 , 2) * t / ( l * l )  ) / ( 2 * m +1 );
        sum +=  a_m  * sin ( pi * x * ( 2 * m +1 ) / l );
    }  
    return sum * 4 * u0 / ( pi );
}