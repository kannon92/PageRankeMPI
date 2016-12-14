#include "SparsePower.h"
#include <iostream>
#include "mpi.h"
#include <cmath>
#include <algorithm>
using namespace std;

SparsePower::SparsePower(const std::vector<std::pair<int, int> >& parsed_data, int& size_of_matrix)
{
    matrix_size_ = size_of_matrix;
    MPI_Bcast(&matrix_size_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(MPI::COMM_WORLD.Get_rank() == 0)
    {
        convert_to_csr(parsed_data);
    }
}
void SparsePower::convert_to_csr(const std::vector<std::pair<int, int> > parsed_data)
{
    std::vector<std::pair<int, int> > parsed_data_modified(parsed_data);
    std::sort(parsed_data_modified.begin(), parsed_data_modified.end(), [](const std::pair<int, int>& left, const std::pair<int, int>& right) {
        return left.second < right.second;
    });

    std::vector<size_t> col_index_H;
    std::vector<size_t> row_index_H;
    std::vector<size_t> row_ptr_H(matrix_size_ + 1);
    ///The first row of parsed data is the column_index in csr
    std::vector<double> non_zero_H(parsed_data_modified.size());
    for(size_t i = 0; i < parsed_data.size(); i++)
    {
        col_index_H.push_back(parsed_data_modified[i].first);
        row_index_H.push_back(parsed_data_modified[i].second);
    }
    for(size_t i = 0; i < parsed_data.size(); i++)
    {
        int value_column = std::count(col_index_H.begin(), col_index_H.end(), col_index_H[i]);
        non_zero_H[i] = 1.0 /(double) value_column;
    }
    for(int i = 0; i < matrix_size_; i++)
    {
        int value_row    = std::count(row_index_H.begin(), row_index_H.end(), i);
        if(value_row > 0)
        {
            if(i == 0)
            {
                row_ptr_H[0] = 0;
                row_ptr_H[1] = value_row;
            }
            else
            {
                row_ptr_H[i+1] = row_ptr_H[i] + value_row;
            }

        }
        else{
            row_ptr_H[i + 1] = row_ptr_H[i];
        }
            
    }

    std::vector<size_t> non_zero_A_where;
    for(int i = 0; i < matrix_size_; i++){
        if (std::find(col_index_H.begin(), col_index_H.end(), i) == col_index_H.end())
        {
            non_zero_A_where.push_back(i);
        }
    }
    std::vector<double> non_zero_A ( non_zero_A_where.size() * matrix_size_ );
    std::vector<size_t> col_index_A( non_zero_A_where.size() * matrix_size_ );
    std::vector<size_t> row_ptr_A(matrix_size_ + 1);
    for(size_t i = 0; i < non_zero_A_where.size() * matrix_size_; i++)
    {
        non_zero_A[i] = 1.0  / matrix_size_;
        col_index_A[i] = non_zero_A_where[i % non_zero_A_where.size()];

    }
    row_ptr_A[0] = 0;
    for(int i = 1; i < matrix_size_ + 1; i++){
        row_ptr_A[i] = row_ptr_A[i-1] + non_zero_A_where.size();
    }

    non_zeroA_ = non_zero_A;
    col_indA_  = col_index_A;
    row_ptrA_  = row_ptr_A;

    col_indH_   = col_index_H;
    non_zeroH_  = non_zero_H;
    row_ptrH_   = row_ptr_H;
        
    display_csr_information();
}
void SparsePower::solve_page_rank_serial()
{
    std::vector<double> x_0(matrix_size_);
    std::vector<double> x_i(matrix_size_);
    
    for(int i = 0; i < matrix_size_; i++)
    {
        x_0[i] = 1.0 / matrix_size_;
        x_i[i] = 0.15 * 1.0 / matrix_size_;
    }

    ///Form H * x
    for(int iter = 0; iter < max_iter_; iter++)
    {
        for (int i = 0; i < matrix_size_; i = i + 1)
        {  
            for (size_t k = row_ptrH_[i]; k < row_ptrH_[i+1]; k = k + 1)
            {  
                x_i[i] = x_i[i] + 0.85 * non_zeroH_[k] * x_0[col_indH_[k]];
             
            }  
            for (size_t k = row_ptrA_[i]; k < row_ptrA_[i+1]; k = k + 1)
            {  
                x_i[i] = x_i[i] + 0.85 * non_zeroA_[k] * x_0[col_indA_[k]];

            }  
        }
        double sum = 0.0;
        for(int i = 0; i < matrix_size_; i++)
        {
            sum += (x_i[i] - x_0[i]) * (x_i[i] - x_0[i]);
        }
        double rms = std::sqrt(sum);
        for(int i = 0; i < matrix_size_; i++)
        {
            x_0[i] = x_i[i];
            x_i[i] = 0.15 * 1.0 / matrix_size_;
        }
        printf("\n iter: %d  rms: %8.8f", iter, rms);
        if(rms < 1e-4)
        {
            //cout << "PageRank converged" << endl;
            break;
        }
        
    }
    //cout <<"Serial ";
    print_solution(x_0);
    //SerialSolution_ = x_0;
    //delete x_0;
    //delete x_i;

}
void SparsePower::display_csr_information()
{
    if(MPI::COMM_WORLD.Get_rank() == 0 && debug_ == true)
    {
    //    cout << "Column IndexH:" << endl;
    //    for(auto& col_ind : col_indH_)
    //    {
    //       cout << col_ind << endl;
    //    }
    //    cout << "Row PtrH:" << endl;
    //    for(auto& row_ptr : row_ptrH_)
    //    {
    //       cout << row_ptr << endl;
    //    }
    //    cout << "NonZeroH:" << endl;
    //    for(auto& non_zero : non_zeroH_)
    //    {
    //       cout << non_zero << endl;
    //    }
    //    cout << "Column IndexA:" << endl;
    //    for(auto& col_ind : col_indA_)
    //    {
    //       cout << col_ind << endl;
    //    }
    //    cout << "Row PtrA:" << endl;
    //    for(auto& row_ptr : row_ptrA_)
    //    {
    //       cout << row_ptr << endl;
    //    }
    //    cout << "NonZeroA:" << endl;
    //    for(auto& non_zero : non_zeroA_)
    //    {
    //       cout << non_zero << endl;
    //    }
    }
}
void SparsePower::solve_page_rank_parallel()
{
    std::vector<double> x_0(matrix_size_);
    int iproc = MPI::COMM_WORLD.Get_rank();
    int nproc = MPI::COMM_WORLD.Get_size();
    int non_zeroH_size, row_ptrH_size, col_indH_size, non_zeroA_size, col_indA_size, row_ptrA_size;
    if(iproc == 0)
    {
        
        non_zeroH_size = non_zeroH_.size();
        row_ptrH_size = row_ptrH_.size() ;
        col_indH_size  = col_indH_.size() ;
        non_zeroA_size = non_zeroA_.size();
        col_indA_size  = col_indA_.size() ;
        row_ptrA_size  = row_ptrA_.size() ;
    }
    MPI_Bcast(&non_zeroA_size ,  1,  MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col_indA_size  ,  1,   MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&row_ptrA_size  ,  1,   MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col_indH_size  ,  1,   MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&non_zeroH_size ,  1,  MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&row_ptrH_size  ,  1,   MPI_INT, 0, MPI_COMM_WORLD);
    if(iproc != 0 )
    {
        row_ptrH_.resize(row_ptrH_size);
        row_ptrA_.resize(row_ptrA_size);
        col_indH_.resize(col_indH_size);
        col_indA_.resize(col_indA_size);
        non_zeroH_.resize(non_zeroH_size);
        non_zeroA_.resize(non_zeroA_size);
        
    }
    MPI_Bcast(&non_zeroA_[0], non_zeroA_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col_indA_[0], col_indA_size,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&row_ptrA_[0], row_ptrA_size,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&col_indH_[0], col_indH_size,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&non_zeroH_[0], non_zeroH_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&row_ptrH_[0], row_ptrH_size,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(int i = 0; i < matrix_size_; i++)
    {
        x_0[i] = 1.0 / matrix_size_;
    }
    MPI_Bcast(&x_0[0], matrix_size_,  MPI_DOUBLE, 0, MPI_COMM_WORLD);

    ///Form H * x
    int mystart,myend;
    mystart = (matrix_size_ / nproc) * iproc;
    std::vector<int> offset (nproc);
    std::vector<int> offset_length (nproc);
    if (matrix_size_ % nproc > iproc){
        mystart += iproc;
        myend = mystart + (matrix_size_ / nproc) + 1;
    }else{
        mystart += matrix_size_ % nproc;
        myend = mystart + (matrix_size_ / nproc);
    }
    int offsetval = myend - mystart;
    MPI_Allgather(&offsetval, 1, MPI_INT, &offset_length[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&(mystart), 1, MPI_INT, &offset[0], 1,MPI_INT, MPI_COMM_WORLD);


    std::vector<double> local_x_i(matrix_size_, 0.0);
    std::vector<double> x_i(matrix_size_, 0.0);
    MPI_Barrier(MPI_COMM_WORLD);
    for(int iter = 0; iter < max_iter_; iter++)
    {
        for(int i = mystart; i < myend; i++)
            local_x_i[i - mystart] = 0.15 / matrix_size_;

        for (int i = mystart; i < myend; i++)
        {
            for (size_t k = row_ptrH_[i]; k < row_ptrH_[i+1]; k = k + 1)
            {
               local_x_i[i - mystart] +=  0.85 * non_zeroH_[k] *  x_0[col_indH_[k]];
            }
        
            for (size_t k = row_ptrA_[i]; k < row_ptrA_[i+1]; k = k + 1)
            {
                
               local_x_i[i - mystart]   +=  0.85 * non_zeroA_[k] * x_0[col_indA_[k]];

            }
        }
        if(iproc == 0)
        {
            for(int i = 0; i < matrix_size_; i++)
               x_i[i] = x_0[i];
        }

        MPI_Allgatherv(&local_x_i[0], offset_length[iproc], MPI_DOUBLE, &x_0[0],&offset_length[0], &offset[0], MPI_DOUBLE, MPI_COMM_WORLD);
        double sum = 0.0;
        double rms = 0.0;
        if(iproc == 0)
        {
        for(int i = 0; i < matrix_size_; i++)
            sum += (x_i[i] - x_0[i]) * (x_i[i] - x_0[i]);
        }
        MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        rms = std::sqrt(sum);
        MPI_Bcast(&rms, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(rms < 1e-4)
        {
            
           if(iproc == 0)
           {
           }
           break; 
        }
    }

    if(MPI::COMM_WORLD.Get_rank() == 0)
    {
                cout << "Page Rank converged" << endl;
                cout << "Parallel " << endl;
                print_solution(x_0);
    }
    //ParallelSolution_ = x_0;
}
void SparsePower::print_solution(const std::vector<double>& final_solution)
{
    //cout << "Solution \n\n" << endl;
    std::vector<std::pair<double, int> > x_0_i;
    for(int i = 0; i < final_solution.size(); i++)
    {
        x_0_i.push_back(std::make_pair(final_solution[i], i));

    }
    std::sort(x_0_i.rbegin(), x_0_i.rend());
    int max_size = (x_0_i.size() < 100 ? x_0_i.size() : 100);
    for(int i = 0; i < max_size; i++)
    {
        printf("%4d\t%8.6f\n", x_0_i[i].second, x_0_i[i].first);
    }

}
