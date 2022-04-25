bool parser (std::vector<std::string> &Arg, bool &flag_time ) {
    for ( auto it : Arg) {
        if ( it[0] == '-' ) {
            if ( it.size() == 1 ) {
                std::cerr<<"no flag after -"<<'\n';
                return false;
            }
            if ( it[1] =='t' && it.size() == 2  )
                flag_time = true;
            else {
                std::cerr<<"unknow flag: "<< it <<'\n';
                return false;
            }
        }
        else {
            std::cerr<<"not - before flag" <<'\n'; 
            return false;
        }
    }
    return true;
}