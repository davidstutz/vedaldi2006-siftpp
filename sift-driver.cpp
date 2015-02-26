// file:        sift-driver.cpp
// author:      Andrea Vedaldi
// description: SIFT command line utility implementation

// AUTORIGTHS

#include<sift.hpp>

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

extern "C" {
#include<getopt.h>
}

using namespace std ;

int
main(int argc, char** argv)
{
  int first    = -1 ;
  int octaves  = -1 ;
  int levels   = 3 ;
  int nodescr  = 0 ;
  int noorient = 0 ;
  int savegss  = 0 ;
  int verbose  = 0 ;
  std::string prefix("") ;

  static struct option longopts[] = {
    { "verbose",         no_argument,            NULL,              'v' },
    { "help",            no_argument,            NULL,              'h' },
    { "output",          required_argument,      NULL,              'o' },
    { "first-octave",    required_argument,      NULL,              'f' },
    { "octaves",         required_argument,      NULL,              'O' },
    { "levels",          required_argument,      NULL,              'S' },
    { "no-descriptors",  no_argument,            &nodescr,          1   },
    { "no-orientations", no_argument,            &noorient,         1   },
    { "save-gss",        no_argument,            &savegss,          1   },
    { NULL,              0,                      NULL,              0   }
  };

  int ch ;
  try {
    while ( (ch = getopt_long(argc, argv, "hvo:O:S:", longopts, NULL)) != -1) {
      switch (ch) {

      case '?' :
        VL_THROW("Invalid option '" << argv[optind-1] << "'.") ;
        break;
        
      case ':' :
        VL_THROW("Missing option argument for '" << argv[optind-1] << "'.") ;
        break;
        
      case 'h' :
        std::cout
          << argv[0] << " [--verbose] [--help] [--output OUTPREFIX]" << endl
          << "   [--save-gss] [--no-descriptors] [--no-orientations] " << endl
          << "   [--levels NLEVELS] [--octaves NOCTAVES] [--first-octave FIRSTOCTAVE] " << endl
          << "   IMAGE [IMAGE2 ...]" << endl ;
        return 0 ;

      case 'v' :
        verbose = 1 ;
        break ;
        
      case 'f': // first octave
        {
          std::istringstream iss(optarg) ;
          iss >> first ;
          if( iss.fail() )
            VL_THROW("Invalid argument '" << optarg << "'.") ;
        }
        break ;
        
      case 'O' : // octaves
        {
          std::istringstream iss(optarg) ;
          iss >> octaves ;
          if( iss.fail() )
            VL_THROW("Invalid argument '" << optarg << "'.") ;
          if( octaves < 1 ) {
            VL_THROW("Number of octaves cannot be smaller than one."); 
          }
        }
        break ;
        
      case 'S' : // levels
        {
          std::istringstream iss(optarg) ;
          iss >> levels ;
          if( iss.fail() )
            VL_THROW("Invalid argument '" << optarg << "'.") ;
          if( levels < 1 ) {
            VL_THROW("Number of levels cannot be smaller than one.") ;
          }
        }      
        break ;

      case 'o' : // output
        {
          prefix = std::string(optarg) ;
          break ;
        }
        
      case 0 : // all other options
        break ;
        
      default:
        assert(false) ;
      }
    }
    
    argc -= optind;
    argv += optind;

    if( argc == 0 ) VL_THROW("No input image specfied.") ;
    if( prefix.size() != 0 && argc > 1 ) {
      VL_THROW("--output cannot be used with multiple images.") ;
    }
  }
  catch( VL::Exception const & e ) {
    std::cerr<<"error: "<<e.msg<<std::endl ;
    exit(1) ;
  }  

  while( argc > 0 ) {
    VL::PgmBuffer buffer ;

    std::string name(argv[0]) ;
    
    // get prefix
    if( prefix.size() == 0 ) {
      int i ;
      std::size_t const not_found = std::numeric_limits<std::size_t>::max() - 1 ;
      if(  ( (i = name.rfind(".pgm",not_found)) != not_found ) ||
           ( (i = name.rfind(".PGM",not_found)) != not_found )  ) {
        prefix = std::string(name, 0, i) ;
      } else {
        prefix = name ;
      }
    }

    // ---------------------------------------------------------------
    //                                                  Load PGM image
    // ---------------------------------------------------------------    
    try {          
      ifstream in( name.c_str() ) ; 
      if(! in.good()) VL_THROW("Could not open '"<< name <<"'.") ;      
      extractPgm( in, buffer ) ;
    }    
    catch( VL::Exception const& e ) {
      cerr<<e.msg<<endl ;
      return 1 ;
    }
    
    if( verbose )
      cout << "Loaded PGM image ("
           << buffer.width<<" x "<<buffer.height<<")."
           << endl ;

    // ---------------------------------------------------------------
    //                                            Gaussian scale space
    // ---------------------------------------------------------------    
    if( verbose ) 
      cout << "Computing Gaussian scale space ... " << std::flush ;

    int         O      = octaves ;    
    int const   S      = levels ;
    int const   omin   = first ;
    float const sigman = .5 ;
    float const sigma0 = 1.6 * powf(2.0f, 1.0f / S) ;

    // autoselect of number of octaves?
    if( O < 1 ) {
      O = 
        std::max
        ( int(floor(log2(std::min(buffer.width,buffer.height))) - omin -4),
          1 ) ;
    }
    
    VL::Sift sift( buffer.data, buffer.width, buffer.height, 
                   sigman, sigma0,
                   O, S,
                   omin, -1, S+1 ) ;

    if(verbose) 
      cout << "(" << O << " octaves starting from " << omin <<", "
           << S << " levels) done." << endl ;

    // ---------------------------------------------------------------
    //                                       Save Gaussian scale space
    // ---------------------------------------------------------------    

    if( savegss ) {
      if(verbose) 
        cout << "Saving Gaussian scale space." << endl ;

      for(int o = omin ; o < omin + O ; ++o) {
        for(int s = 0 ; s < S ; ++s ) {

          std::ostringstream suffix ;
          suffix<<'.'<<o<<'.'<<s<<".pgm" ;
          std::string outname = prefix + suffix.str() ;

          if(verbose) 
            cout << "Saving octave " <<o
                 << " level "<<s
                 << " to '"<< outname << "' ... " << std::flush ;
          
          ofstream out( outname.c_str() ) ;          
          VL::insertPgm( out,
                         sift.getLevel(o,s),
                         sift.getOctaveWidth(o),
                         sift.getOctaveHeight(o) ) ;
          
          if(verbose) cout << "done." <<endl ;
        }
      }
    }
    
    // ---------------------------------------------------------------
    //                                                   SIFT detector
    // ---------------------------------------------------------------    
    
    if(verbose) 
      cout << "Running SIFT detector ... " << std::flush ;

    sift.detectKeypoints() ;

    if(verbose) 
      cout << "done (" 
           << sift.keypointsEnd() - sift.keypointsBegin() 
           << " keypoints)." << endl ;

    // ---------------------------------------------------------------
    //                               SIFT orientations and descriptors
    // ---------------------------------------------------------------    

    if( verbose ) {
      if( ! noorient &&   nodescr) cout << "Computing keypoint orientations" ;
      if(   noorient && ! nodescr) cout << "Computing keypoint descriptors" ;
      if( ! noorient && ! nodescr) cout << "Computing orientations and descriptors" ;
      if(   noorient &&   nodescr) cout << "Finalizing" ; 
    }

    {
      std::string outname = prefix + ".key" ;
      ofstream out( outname.c_str() ) ;
      
      if(verbose)
        std::cout << " (saving keypoints to '" << outname << "') ... " 
                  << std::flush ;
      
      for( VL::Sift::KeypointsConstIter iter = sift.keypointsBegin() ;
           iter != sift.keypointsEnd() ;
           ++iter ) {
        
        // detect orientations
        float angles [4] ;
        int nangles ;
        if( ! noorient ) {
          nangles = sift.computeKeypointOrientations(angles, *iter) ;
        } else {
          nangles = 1;
          angles[0] = 0.0f ;
        }
        
        // compute descriptors
        for(int a = 0 ; a < nangles ; ++a) {
          
          out << iter->x << " "
              << iter->y << " "
              << iter->sigma ;
          
          if( ! noorient ) {
            out << " " << angles[a] ;
          }
          
          if( ! nodescr ) {
            float descr [128] ;
            sift.computeKeypointDescriptor(descr, *iter, angles[a]) ;
            for(int i = 0 ; i < 128 ; ++i) {
              out << " " << (int)( 512.0f * descr[i] ) ;
            }
          }

          // next keypoint
          out << endl ;
        } 
      } 
      
        
      if(verbose)
        std::cout << "done." << std::endl ;
    }
    
    argc-- ;
    argv++ ;
  } // next image

  return 0 ;
}
