
void print_unknow_flag ( const std::string &flag ) {
    std::cerr<<"unknow flag: "<<flag<<'\n';
    return;
}

bool parser (std::vector<std::string> &Arg, bool &flag_res, bool &flag_segment, long long &N ) {
    for ( auto it : Arg) {
        if ( it[0] == '-' ) {
            if ( it.size() == 1 ) {
                std::cerr<<"нет флага после -"<<'\n';
                return false;
            }
            switch ( it[1] )
            {
            case 't':
                if ( it.size() == 2) {
                    flag_res = false;
                    break;
                }
                else 
                {
                    print_unknow_flag ( it );
                    return false;
                    break; 
                }
            case 's':
                if ( it.size() == 2) {
                    flag_segment = true;
                    break;
                }
                else 
                {
                    print_unknow_flag ( it );
                    return false;
                    break; 
                }
            case 'N': {
                unsigned int idx = 2;
                if ( it[2] == 'e' )
                    ++idx;
                std::string meaning;
                for ( size_t j = idx; j < it.size(); ++j )
                    meaning += it[j];
                long long meaning_ll = atoi( meaning.c_str() );
                N = 1;
                if ( it[2] == 'e' )
                    for ( size_t j = 1; j <= meaning_ll ; ++j )
                        N *= 10;
                else
                    N = meaning_ll;
                break;     
                }
            case 'n':
                if ( it[2] == 'p' && it.size() == 3 )
                    break;
                else 
                {
                    print_unknow_flag ( it );
                    return false;
                    break;
                }                        
            default:
                print_unknow_flag ( it );
                return false;
                break;
            }
        }
    }
    return true;
}