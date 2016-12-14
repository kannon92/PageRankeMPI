#include <vector>
#include <cstdlib>
class SparsePower
{
private: 
    // Sparse representation for Adjaceny matrix
    std::vector<double> non_zeroH_;
    std::vector<size_t> col_indH_;
    std::vector<size_t> row_ptrH_;

    // Sparse representation for matrix with no connections
    std::vector<double> non_zeroA_;
    std::vector<size_t> col_indA_;
    std::vector<size_t> row_ptrA_;
    int matrix_size_;
    int max_iter_ = 100;
    double norm_ = 0.0;
    /// MPI Information
    int nproc_;
    bool debug_ = false;
    bool print_csr_ = false;
    std::vector<double> SerialSolution_;
    std::vector<double> ParallelSolution_;
    /// Take the parsed_data and create data_value, row_ptr, and col_index
    void convert_to_csr(const std::vector<std::pair<int, int> >);
    /// Print the final solutions
    void print_solution(const std::vector<double>&);
    
public:
    SparsePower(const std::vector<std::pair<int, int> >&, int& size_of_matrix);
    std::vector<double> non_zeroH(){return non_zeroH_;}
    std::vector<size_t> row_ptrH(){return row_ptrH_;} 
    std::vector<size_t> col_indH(){return col_indH_;} 
    std::vector<double> non_zeroA(){return non_zeroA_;}
    std::vector<size_t> row_ptrA(){return row_ptrA_;} 
    std::vector<size_t> col_indA(){return col_indA_;} 
    void display_csr_information();
    void solve_page_rank_serial();
    void solve_page_rank_parallel();
    void solve_page_rank_stupid();
    void debug_mode(bool debug){debug_ = debug;}
    void print_csr(){print_csr_ = true;}
    std::vector<double> SerialSolution(){return SerialSolution_;}
    std::vector<double> ParallelSolution(){return ParallelSolution_;}
    void set_max_iterations(int max_iter){max_iter_ = max_iter;}
};
