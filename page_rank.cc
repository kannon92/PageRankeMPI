#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <ctime>
#include <tuple>
#include "SparsePower.h"
#include "mpi.h"
#include <cstring>
using namespace std;

std::vector<std::pair< int, int> >parse_file(std::string file, int& number_of_nodes);

int main(int argc, char* argv[])
{
    MPI::Init(argc, argv);
    std::string file_name;
    bool debug = false;
    bool stupid_pagerank = false;
    int max_iter = 100;

    std::vector<std::pair<int, int> > parsed_file;
    int number_of_nodes;
    for(int i = 1; i < argc; i++)
    {
        if ( !( std::strcmp(argv[i] , "-file_name") ) )
        {
            file_name = argv[i+1];
            parsed_file = parse_file(file_name, number_of_nodes);

        }
        else if( !(std::strcmp(argv[i], "-debug") ) )
        {
            debug = true;
        }
        else if( !(std::strcmp(argv[i], "-maxiter")))
        {
            max_iter = atoi(argv[i+1]);
        }
        else if( !(std::strcmp(argv[i], "-stupid_pagerank")))
        {
            stupid_pagerank = true;
        }
    }
    SparsePower page_rank(parsed_file, number_of_nodes);
    page_rank.debug_mode(debug);
    page_rank.set_max_iterations(max_iter);
    double start, end;
    if(MPI::COMM_WORLD.Get_rank() == 0 && debug)
    {
        page_rank.display_csr_information();
	start = MPI::Wtime();
        page_rank.solve_page_rank_serial();
	end = MPI::Wtime();
	    cout << "Serial Sparse Page Rank takes " << end - start << " s " << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPI::COMM_WORLD.Get_rank() == 0){start = MPI::Wtime();}
    page_rank.solve_page_rank_parallel();
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPI::COMM_WORLD.Get_rank() == 0){end   = MPI::Wtime();}
    
    if(MPI::COMM_WORLD.Get_rank() == 0)
    {
	    cout << "Parallel Sparse Page Rank takes " << end - start << " s." << endl;
    }

    MPI::Finalize();
}
std::vector<std::pair<int, int> > parse_file(std::string file, int& number_of_nodes)
{
    std::vector<std::pair<int, int> > return_pair;
    std::ifstream file_stream;
    file_stream.open(file);
    std::string file_line;
    std::pair<int, int> ij;
    file_stream >> number_of_nodes;
    while(std::getline(file_stream, file_line))
    {
        int nodeone;
        int nodetwo;
        file_stream >> nodeone;
        file_stream >> nodetwo;
        ij = std::make_pair(nodeone, nodetwo);
        if(file_stream.eof() ) break;
        return_pair.push_back(ij);
    }
    file_stream.close();
    return return_pair;
    

}
