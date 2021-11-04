//
// Created by laura on 9/26/18.
//

#ifndef AFQMCLAB_REALMATERIALSOC_H
#define AFQMCLAB_REALMATERIALSOC_H


#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/common/common.h"
#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/afqmc/blocks/oneBodyOperator/hop/include/hop.h"
#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/afqmc/blocks/oneBodyOperator/logHop/include/logHop.h"
#include "realMaterialSOCInteract.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/common/common.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/afqmc/blocks/oneBodyOperator/hop/include/hop.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/afqmc/blocks/oneBodyOperator/logHop/include/logHop.h"

//Details about t, K, Kp are in my note.
//t is the non-interacting matrix.
//K is t - choleskyVecs*choleskyVecs.
//Kp is K + choleskyBg*choleskyVecs.

#ifdef MPI_HAO
class RealMaterialSOC;
void MPIBcast(RealMaterialSOC &buffer, int root=0,  const MPI_Comm& comm=MPI_COMM_WORLD);
#endif

class RealMaterialSOC
{
private:
    size_t L, N, Nup, Ndn; //L=NMO
    size_t eigenNum;
    std::complex<double> SumEngShf;
    tensor_hao::TensorHao<std::complex<double>,2> t, K;
    tensor_hao::TensorHao<double,1> eigenValue;
    tensor_hao::TensorHao<double,3> phoSum, phoMinus;
    tensor_hao::TensorHao<double,4> P_upup, P_dndn, P_updn;
    tensor_hao::TensorHao<double,1> phoSumBg, phoMinusBg;

    size_t KpEigenStatus; //0: void, 1: calculated Kp, 2: calculated Kp Eigens.
    tensor_hao::TensorHao<std::complex<double>,2> Kp;
    tensor_hao::TensorHao<double,1> KpEigenValue;
    tensor_hao::TensorHao<std::complex<double>,2> KpEigenVector;

public:
    RealMaterialSOC();
    RealMaterialSOC(const std::string &filename);
    ~RealMaterialSOC();

    size_t getL() const;
    size_t getN() const;
    size_t getNup() const;
    size_t getNdn() const;
    size_t getEigenNumber() const;
    const tensor_hao::TensorHao<std::complex<double>,2> &getT() const;
    const tensor_hao::TensorHao<std::complex<double>,2> &getK() const;
    const tensor_hao::TensorHao<double,4> &getPupup() const;
    const tensor_hao::TensorHao<double,4> &getPdndn() const;
    const tensor_hao::TensorHao<double,4> &getPupdn() const;
    //const tensor_hao::TensorHao<double,3> &getphoUpSum() const;
    //const tensor_hao::TensorHao<double,3> &getphoUpMinus() const;
    //const tensor_hao::TensorHao<double,3> &getphoDnSum() const;
    //const tensor_hao::TensorHao<double,3> &getphoDnMinus() const;
    const tensor_hao::TensorHao<double,3> &getphoSum() const;
    const tensor_hao::TensorHao<double,3> &getphoMinus() const;
    const tensor_hao::TensorHao<double,1> &getphoSumBg() const;
    const tensor_hao::TensorHao<double,1> &getphoMinusBg() const;
    size_t getKpEigenStatus() const;
    const tensor_hao::TensorHao<std::complex<double>,2> &getKp() const;
    const tensor_hao::TensorHao<double,1> &getKpEigenValue() const;
    const tensor_hao::TensorHao<std::complex<double>,2> &getKpEigenVector() const;

    void read(const std::string &filename);
    void checkpotential(const std::string &filename) const;
    void checkkinetic(const std::string &filename) const;
    void write(const std::string &filename) const;
    void readFCIDUMP();
    void readFCIDUMPSOC();
#ifdef MPI_HAO
    friend void MPIBcast(RealMaterialSOC &buffer, int root,  const MPI_Comm& comm);
#endif

    void writeBackGround(const std::string &filename) const;
    void updateBackGround(const tensor_hao::TensorHao<double,1> &backgroundOne, const tensor_hao::TensorHao<double,1> &backgroundTwo);
    void updateBackGround(tensor_hao::TensorHao<double,1> &&backgroundOne, tensor_hao::TensorHao<double,1> &&backgroundTwo);

    Hop returnExpMinusAlphaK(double alpha);
    LogHop returnLogExpMinusAlphaK(double alpha);
    RealMaterialSOCInteract returnExpMinusAlphaV(double alpha);

    //void setpho();
    void setKp();
    void setKpEigenValueAndVector();

    double getMemory() const;

private:
    RealMaterialSOC(const RealMaterialSOC& x);
    RealMaterialSOC & operator  = (const RealMaterialSOC& x);
};

#endif //AFQMCLAB_REALMATERIALSOC_H