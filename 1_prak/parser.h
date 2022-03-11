
bool parser (std::vector<std::string> &Arg, bool &flag_res, bool &flag_segment, long long &N ) {
    for ( auto it : Arg) {
        if ( it[0] == '-' ) {
            if ( it.size() == 1 ) {
                std::cout<<"нет флага после -"<<'\n';
                return false;
            }
            switch ( it[1] )
            {
            case 't':
                flag_res = false;
                break;
            case 's':
                flag_segment = true;
                break;
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
                break;                       
            default:
                std::cout<< "unknow flag \n";
                return false;
                break;
            }
        }
    }
    return true;
}