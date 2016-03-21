#ifndef MIXMAT_HPP
#define MIXMAT_HPP

#include <vector>
#include <cstddef>
#include <iostream>
#include <wignerSymbols.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <exception>

#include <Eigen/Core>
#include <Eigen/Dense>

enum mixMatType
{
    TT,
    EE,
};

template<typename _realScalarType>
class MixingMatrix
{
public:
    /**
     * \typedef _realScalarType realScalarType
     * \brief real floating point type
     */
    typedef _realScalarType realScalarType;

    /**
     * \typedef typename Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> mixingMatrixType
     * \brief mixing matrix type
     */
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> mixingMatrixType;


    /**
     * \typedef typename mixingMatrixType::Index indexType
     * \brief integral type
     */
    typedef typename mixingMatrixType::Index indexType;


    /**
     * \typedef Eigen::LLT<mixingMatrixType> LLTType;
     * \brief Cholesky decompostion type
     */
    typedef Eigen::LLT<mixingMatrixType> LLTType;

    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> powSpecType;


    MixingMatrix(powSpecType const & W,mixMatType const mt)
    :mW(W)
    {
        const indexType lMax = W.rows()-1;
        const indexType matDim(lMax+1);
        mM = mixingMatrixType::Zero(matDim,matDim);
        mLMax = lMax;
        mMT = mt;
        std::cout<<"lMax = "<<lMax<<std::endl;
    }

    void computeMixingMatrixWithRecursion()
    {
        if(mMT == TT)
        {
            std::cout<<"Computing TT mixing matrix"<<std::endl;

            const indexType matDim(mLMax+1);
            mM = mixingMatrixType::Zero(matDim,matDim);

            for(indexType l1=0;l1<=mLMax;++l1)
            {
                for(indexType l2=0;l2<=mLMax;++l2)
                {
                    std::vector<double> wig3j =
                        WignerSymbols::wigner3j((double)l1,(double)l2,
                            (double)0,(double)0,(double)0);

                    indexType l3Min = std::max(std::abs(l1-l2),indexType(0));
                    indexType l3max = std::min( l3Min + indexType(wig3j.size()-1),mLMax ) ;

                    for(indexType l3=l3Min;l3<=l3max;++l3)
                    {
                        if( selectWig3j( l3,l1,l2,indexType(0),indexType(0),indexType(0) )  )
                        {
                            mM(l1,l2) += realScalarType(2*l3+1)*mW(l3)*wig3j[l3 - l3Min]*wig3j[l3 - l3Min];
                        }
                    }

                    mM(l1,l2) *= realScalarType(2*l2+1)/(4*M_PI);
                }
            }
        }
        else if(mMT == EE)
        {
            std::cout<<"Computing EE mixing matrix"<<std::endl;

            const indexType matDim(mLMax+1);
            mM = mixingMatrixType::Zero(matDim,matDim);

            for(indexType l1=0;l1<=mLMax;++l1)
            {
                for(indexType l2=0;l2<=mLMax;++l2)
                {
                    std::vector<double> wig3j =
                        WignerSymbols::wigner3j((double)l1,(double)l2,
                            (double)0,(double)2,(double)-2);

                    indexType l3Min = std::max(std::abs(l1-l2),indexType(0));
                    indexType l3max = std::min( l3Min + indexType(wig3j.size()-1),mLMax ) ;

                    for(indexType l3=l3Min;l3<=l3max;++l3)
                    {
                        if( selectWig3j( l3,l1,l2,indexType(0),indexType(2),indexType(-2) ) && (l1+l2+l3)%2 == 0 )
                        {
                            mM(l1,l2) += realScalarType(2*l3+1)*realScalarType(2)*mW(l3)*wig3j[l3 - l3Min]*wig3j[l3 - l3Min];
                        }
                    }

                    mM(l1,l2) *= realScalarType(2*l2+1)/(8*M_PI);
                }
            }

        }

    }

    bool selectWig3j(indexType const l1,indexType const l2,indexType const l3,
        indexType const m1,indexType const m2,indexType const m3)
    {
        bool select(true);

        select = (
            std::abs(m1+m2+m3) == indexType(0)
            && l3 >= std::abs(l1-l2)
            && l3 <= l1+l2
            && std::abs(m1) <= l1
            && std::abs(m2) <= l2
            && std::abs(m3) <= l3
        );

        return ( select );
    }

    void save(std::string const fileName)
    {
        writeMatrixToFile(fileName,mM,10,",");
    }


    void applyMixingMatrix(powSpecType & C)
    {
        C = mM*C;
    }

    void applyMixingMatrixInv(powSpecType & C)
    {
        LLTType lltOfM;
        lltOfM.compute(mM);
        assert(lltOfM.info() == Eigen::Success);
        C = lltOfM.solve(C);
    }

private:

    void writeMatrixToFile(std::string const & fileName,
        mixingMatrixType const& mat,unsigned int precision,
        std::string const& separation)
    {
        std::ofstream outFile;
        outFile.open(fileName.c_str(),std::ios::trunc);

        if(outFile.is_open())
        {
            outFile<<std::scientific;
            for(indexType i=0;i<mat.rows();++i)
            {
                for(indexType j=0;j<mat.cols()-1;++j)
                {
                    outFile<<std::setprecision(precision)<<mat(i,j)<<separation;
                }
                outFile<<mat(i,(mat.cols()-1) )<<std::endl;
            }
            outFile.close();
        }
        else
        {
            std::string msg = std::string("Error in opening the file ") +
                fileName;
            throw std::runtime_error(msg);
        }
    }

    powSpecType mW;
    mixingMatrixType mM;
    indexType mLMax;
    mixMatType mMT;

};

#endif //MIXMAT_HPP
