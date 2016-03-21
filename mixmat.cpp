#include "mixmat.hpp"

#include <fstream>
#include <string>

int main(void)
{
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> powSpecType;
    // read the file
    std::vector<double> W;//(4097,double(1));

    //std::ifstream infile("../unitPower.dat");
    //std::ifstream infile("/home/sbalan/Downloads/cl_mask_1024.dat");
    std::ifstream infile("/home/sbalan/Downloads/cl_mask_8192.dat");

    if(!infile.is_open())
    {
        std::cout<<"Input file cannot be opened"<<std::endl;
        return 0;
    }

    std::string content;

    while(infile >> content)
    {
        W.push_back( std::stod(content) );
    }
    infile.close();


    powSpecType WEll(W.size());
    for(size_t i=0;i<(size_t)WEll.rows();++i)
    {
        WEll(i) = W[i];
        //std::cout<<i<<"\t"<<WEll(i)<<std::endl;
    }

    /*
    MixingMatrix<double> mm(WEll,TT);
    mm.computeMixingMatrixWithRecursion();
    mm.save(std::string("mixmat_test.dat"));
    */

    MixingMatrix<double> mm(WEll,EE);
    mm.computeMixingMatrixWithRecursion();
    //mm.save(std::string("mixmat_test_EE.dat"));
    mm.save(std::string("mixmat_test_EE_8192.dat"));

    return 0;
}
