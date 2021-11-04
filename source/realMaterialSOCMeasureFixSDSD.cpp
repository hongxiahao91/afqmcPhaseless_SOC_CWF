//

#include "../include/realMaterialSOCMeasureFixSDSD.h"
#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/afqmc/utilities/manipulateMCData/include/writeThreadSum.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/afqmc/utilities/manipulateMCData/include/writeThreadSum.h"
using namespace std;
using namespace tensor_hao;

RealMaterialSOCMeasureFixSDSD::RealMaterialSOCMeasureFixSDSD()
{
    initModelWalkerNullptr();
    reSet();
}

RealMaterialSOCMeasureFixSDSD::RealMaterialSOCMeasureFixSDSD(const RealMaterialSOC &realMaterialSOC_,
                                                                                     const SD &walkerLeft_)
{
    setModelWalker(realMaterialSOC_, walkerLeft_);
    reSet();
}

RealMaterialSOCMeasureFixSDSD::~RealMaterialSOCMeasureFixSDSD()
{

}

void RealMaterialSOCMeasureFixSDSD::initModelWalkerNullptr()
{
    realMaterialSOC = nullptr;
    walkerLeft = nullptr;
}

void RealMaterialSOCMeasureFixSDSD::setModelWalker(const RealMaterialSOC &realMaterialSOC_, const SD &walkerLeft_)
{
    if( 2*realMaterialSOC_.getL() != walkerLeft_.getL() ) {cout<<"Model L does not consistent with walker L!"<<endl; exit(1);}
    if( realMaterialSOC_.getN() != walkerLeft_.getN() ) {cout<<"Model N does not consistent with walker N!"<<endl; exit(1);}

    realMaterialSOC = &realMaterialSOC_;
    walkerLeft = &walkerLeft_;
    initWfDaggerT();
}

void RealMaterialSOCMeasureFixSDSD::reSet()
{
    complex<double> zero(0,0);
    den = zero;
    TNum = zero;
    HNum = zero;
    phoSumBgNum = zero; phoMinusBgNum = zero;
    NupNum = zero; NdnNum = zero; SplusNum = zero; SminusNum = zero;
    NupTotNum = zero; NdnTotNum = zero; SplusTotNum = zero; SminusTotNum = zero;
}

complex<double> RealMaterialSOCMeasureFixSDSD::returnEnergy()
{
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;
    if( MPIRank() == 0 ) energy = Htot/denTot;
    MPIBcast(energy);
    return energy;
}

TensorHao<double,1> RealMaterialSOCMeasureFixSDSD::returnPhoSumBg()
{
    size_t eigenNum = realMaterialSOC->getEigenNumber();

    TensorHao<complex<double>, 1> phoSumBgTot( eigenNum );
    MPISum( eigenNum, phoSumBgNum.data(), phoSumBgTot.data() );

    complex<double> denTot = MPISum(den);

    TensorHao<double,1> phoSumBg( eigenNum );
    for(size_t i = 0; i < eigenNum; ++i)
    {
        phoSumBg(i) = ( phoSumBgTot(i)/denTot ).real();
    }
    MPIBcast(phoSumBg);

    return phoSumBg;
}

TensorHao<double,1> RealMaterialSOCMeasureFixSDSD::returnPhoMinusBg()
{
    size_t eigenNum = realMaterialSOC->getEigenNumber();

    TensorHao<complex<double>, 1> phoMinusBgTot(eigenNum);

    MPISum( eigenNum, phoMinusBgNum.data(), phoMinusBgTot.data() );

    complex<double> denTot = MPISum(den);

    TensorHao<double,1> phoMinusBg(eigenNum);
    for(size_t i = 0; i < eigenNum; ++i)
    {
        phoMinusBg(i) = ( phoMinusBgTot(i)/denTot ).real();
    }

    MPIBcast(phoMinusBg);

    return phoMinusBg;
}

void RealMaterialSOCMeasureFixSDSD::addMeasurement(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    checkWalkerLeft(sdsdOperation);

    den += denIncrement;

    addEnergy(sdsdOperation, denIncrement);
    addPhoSumBg(sdsdOperation, denIncrement);
    addPhoMinusBg(sdsdOperation, denIncrement);
    addNupNdn(sdsdOperation, denIncrement);
    addSplusSminus(sdsdOperation, denIncrement);
}

RealMaterialSOCInteractForce RealMaterialSOCMeasureFixSDSD::getForce(const RealMaterialSOCInteract &realMaterialSOCInteract,
                                                                     SDSDOperation &sdsdOperation,
                                                                      double cap)
{
    checkWalkerLeft(sdsdOperation);

    size_t eigenNum = realMaterialSOC->getEigenNumber();
    const TensorHao<double,1> & currentPhoSumBg = realMaterialSOC->getphoSumBg();
    const TensorHao<double,1> & currentPhoMinusBg = realMaterialSOC->getphoMinusBg();
    const TensorHao<complex<double>,1> minusGamma = realMaterialSOCInteract.getMinusGamma();
    const TensorHao<complex<double>,1> plusGamma = realMaterialSOCInteract.getPlusGamma();
    TensorHao<complex<double>, 1> phoSumBg = calculatePhoSumBg(sdsdOperation);
    TensorHao<complex<double>, 1> phoMinusBg = calculatePhoMinusBg(sdsdOperation);

    RealMaterialSOCInteractForce force(eigenNum);
    complex<double> plusForce, minusForce;
    for(size_t i = 0; i < eigenNum; ++i)
    {
        plusForce = (phoSumBg(i)-currentPhoSumBg(i)) * minusGamma(i);
        if( abs(plusForce) > cap ) force.PlusForce(i) = plusForce*cap/abs(plusForce);
        else force.PlusForce(i) = plusForce;

        minusForce = (phoMinusBg(i)-currentPhoMinusBg(i)) * plusGamma(i);
        if( abs(minusForce) > cap ) force.MinusForce(i) = minusForce*cap/abs(minusForce);
        else force.MinusForce(i) = minusForce;
    }

    return force;
}

void RealMaterialSOCMeasureFixSDSD::write() const
{
    writeThreadSum(den, "den.dat", ios::app);
    writeThreadSum(TNum, "TNum.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);
    writeThreadSum(phoSumBgNum.size(), phoSumBgNum.data(), "phoSumBgNum.dat", ios::app);
    writeThreadSum(phoMinusBgNum.size(), phoMinusBgNum.data(), "phoMinusBgNum.dat", ios::app);

    writeThreadSum(NupNum.size(),    NupNum.data(),    "NupNum.dat",    ios::app);
    writeThreadSum(NdnNum.size(),    NdnNum.data(),    "NdnNum.dat",    ios::app);
    writeThreadSum(SplusNum.size(),  SplusNum.data(),  "SplusNum.dat",  ios::app);
    writeThreadSum(SminusNum.size(), SminusNum.data(), "SminusNum.dat", ios::app);

    writeThreadSum(NupTotNum,   "NupTotNum.dat",   ios::app);
    writeThreadSum(NdnTotNum,   "NdnTotNum.dat",   ios::app);
    writeThreadSum(SplusTotNum, "SplusTotNum.dat", ios::app);
    writeThreadSum(SminusTotNum, "SminusTotNum.dat", ios::app);
}

double RealMaterialSOCMeasureFixSDSD::getMemory() const
{
    double mem(0.0);
    mem += 8.0*2;
    mem += 16.0*3; //den, TNum, HNum
    mem += phoSumBgNum.getMemory()+phoMinusBgNum.getMemory();
    mem += NupNum.getMemory()+NdnNum.getMemory()+SplusNum.getMemory()+SminusNum.getMemory();
    mem += 16.0*4;
    mem += wfDaggerT.getMemory();

    return mem;
}

RealMaterialSOCMeasureFixSDSD::RealMaterialSOCMeasureFixSDSD(const RealMaterialSOCMeasureFixSDSD &x)
{

}

RealMaterialSOCMeasureFixSDSD & RealMaterialSOCMeasureFixSDSD::operator=(const RealMaterialSOCMeasureFixSDSD &x)
{
    return *this;
}

void RealMaterialSOCMeasureFixSDSD::initWfDaggerT()
{
    size_t L = walkerLeft->getL(); size_t N = walkerLeft->getN();

    wfDaggerT.resize(N, L);

    BL_NAME(gmm)(walkerLeft->getWf(), realMaterialSOC->getT(), wfDaggerT, 'C');
}

//void RealMaterialSOCMeasureFixSDSD::initWfDaggerCholeskyVecs()
//{
//    size_t L = walkerLeft->getL(); size_t N = walkerLeft->getN(); size_t Ndn = walkerLeft->getNdn();
//    size_t choleskyNumber = realMaterialSOC->getCholeskyNumber();
//
//    wfUpDaggerCholeskyVecs.resize(Nup, L, choleskyNumber); wfDnDaggerCholeskyVecs.resize(Ndn, L, choleskyNumber);
//
//    const TensorHao<double,3> & choleskyVecs = realMaterialSOC->getCholeskyVecs();
//    TensorHao<complex<double>, 2> choleskyVecComplex(L,L);
//    TensorHaoRef<complex<double>, 2> wfDaggerCholeskyVec;
//    for(size_t k = 0; k < choleskyNumber ; ++k)
//    {
//        for(size_t i = 0; i < L; ++i)
//        {
//            for(size_t j = 0; j < L; ++j) choleskyVecComplex(j,i) = choleskyVecs(j,i,k);
//        }
//        wfDaggerCholeskyVec=wfUpDaggerCholeskyVecs[k];
//        BL_NAME(gmm)(walkerLeft->getWfUp(), choleskyVecComplex, wfDaggerCholeskyVec, 'C');
//        wfDaggerCholeskyVec=wfDnDaggerCholeskyVecs[k];
//        BL_NAME(gmm)(walkerLeft->getWfDn(), choleskyVecComplex, wfDaggerCholeskyVec, 'C');
//    }
//}

void RealMaterialSOCMeasureFixSDSD::checkWalkerLeft(const SDSDOperation &sdSDOperation)
{
    if( walkerLeft != sdSDOperation.getWalkerLeft() )
    {
        cout<<"Error!!! realMaterialSOCMeasureFixSDSD only accept SDSDOperation with fixed SD!"<<endl;
        exit(1);
    }
}

void RealMaterialSOCMeasureFixSDSD::addEnergy(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{

    complex<double> Tenergy=calculateTenergy(sdsdOperation);
    complex<double> Venergy=calculateVenergy(sdsdOperation);

    complex<double> Henergy(0.0);

    Henergy += Venergy;
    cout<<"Venergy: "<< Henergy<<endl;
    Henergy += Tenergy;
    cout<<"Tenergy: "<< Tenergy<<endl;


    TNum += ( Tenergy * denIncrement );
    HNum += ( Henergy * denIncrement );
}

/*
void RealMaterialSOCMeasureFixSDSD::addEnergy(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    size_t choleskyNumber = realMaterialSOC->getCholeskyNumber();

    if( choleskyBgNum.rank(0) != choleskyNumber ) { choleskyBgNum.resize(choleskyNumber); choleskyBgNum = complex<double>(0,0); }
    if( choleskyExNum.rank(0) != choleskyNumber ) { choleskyExNum.resize(choleskyNumber); choleskyExNum = complex<double>(0,0); }

    complex<double> Tenergy=calculateTenergy(sdsdOperation);
    TensorHao<complex<double>, 1> choleskyBg = calculateCholeskyBg(sdsdOperation);
    TensorHao<complex<double>, 1> choleskyEx = calculateCholeskyEx(sdsdOperation);

    complex<double> Henergy(0.0);
    for(size_t i = 0; i < choleskyNumber; ++i)
    {
        Henergy += ( choleskyBg(i)*choleskyBg(i) -choleskyEx(i) );
    }
    Henergy *= 0.5;
    cout<<"Henergy: "<< Henergy<<endl;
    Henergy += Tenergy;
    cout<<"Tenergy: "<< Tenergy<<endl;


    TNum += ( Tenergy * denIncrement );
    choleskyBgNum += ( choleskyBg * denIncrement );
    choleskyExNum += ( choleskyEx * denIncrement );
    HNum += ( Henergy * denIncrement );
}
*/

complex<double> RealMaterialSOCMeasureFixSDSD::calculateTenergy(SDSDOperation &sdSDOperation)
{
    size_t L = walkerLeft->getL(); size_t N = walkerLeft->getN();
    const TensorHao<complex<double>, 2> &theta_T =  sdSDOperation.returnTheta_T();


    complex<double> Tenergy(0,0);
    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < N; ++j) Tenergy += wfDaggerT(j,i) * theta_T(j,i);
    }
    return Tenergy;
}

complex<double> RealMaterialSOCMeasureFixSDSD::calculateVenergy(SDSDOperation &sdsdOperation)
{
    size_t L = realMaterialSOC->getL();
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    const TensorHao< double, 4> &P_upup = realMaterialSOC->getPupup();
    const TensorHao< double, 4> &P_updn = realMaterialSOC->getPupdn();
    const TensorHao< double, 4> &P_dndn = realMaterialSOC->getPdndn();

    complex<double> Venergy(0,0);
    for(size_t i = 0; i < L; ++i){
        for(size_t l = 0; l < L; ++l){
            for(size_t j = 0; j < L; ++j){
                for(size_t k = 0; k < L; ++k){
                    Venergy += 0.5*P_upup(i, l, j, k)*( greenMatrix(i,   l  )*greenMatrix(j,   k  ) - greenMatrix(i,   k  )*greenMatrix(j,   l  ));
                    Venergy += 0.5*P_dndn(i, l, j, k)*( greenMatrix(i+L, l+L)*greenMatrix(j+L, k+L) - greenMatrix(i+L, k+L)*greenMatrix(j+L, l+L));
                    Venergy += 0.5*P_updn(i, l, j, k)*( greenMatrix(i,   l  )*greenMatrix(j+L, k+L) - greenMatrix(i,   k+L)*greenMatrix(j+L, l  ));
                    Venergy += 0.5*P_updn(i, l, j, k)*( greenMatrix(i+L, l+L)*greenMatrix(j,   k  ) - greenMatrix(i+L, k  )*greenMatrix(j,   l+L));
                }
            }
        }
    }
    return Venergy;
}


TensorHao<complex<double>, 1> RealMaterialSOCMeasureFixSDSD::calculatePhoSumBg(SDSDOperation &sdsdOperation)
{
    size_t L = realMaterialSOC->getL();  //NMO
    size_t eigenNum = realMaterialSOC->getEigenNumber();
    TensorHao<complex<double>, 1> phoSumBg(eigenNum);
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();
    const TensorHao<double,3> &phoSum = realMaterialSOC->getphoSum();

    for(size_t k = 0; k < eigenNum; k++)
    {
        phoSumBg(k)=0.0;

        for(size_t j = 0; j < L; j++)
        {
            for(size_t i = 0; i < L; i++)
            {
                phoSumBg(k) += phoSum(i,j,k) * greenMatrix(i, j);
                phoSumBg(k) += phoSum(i,j,k) * greenMatrix(i+L, j+L);
            }
        }
    }

    return phoSumBg;
}

TensorHao<complex<double>, 1> RealMaterialSOCMeasureFixSDSD::calculatePhoMinusBg(SDSDOperation &sdsdOperation)
{
    size_t L = realMaterialSOC->getL();  //NMO
    size_t eigenNum = realMaterialSOC->getEigenNumber();
    TensorHao<complex<double>, 1> phoMinusBg(eigenNum);
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();
    const TensorHao<double,3> &phoMinus = realMaterialSOC->getphoMinus();

    for(size_t k = 0; k < eigenNum; k++)
    {
        phoMinusBg(k)=0.0;

        for(size_t j = 0; j < L; j++)
        {
            for(size_t i = 0; i < L; i++)
            {
                phoMinusBg(k) += phoMinus(i,j,k) * greenMatrix(i, j);
                phoMinusBg(k) += phoMinus(i,j,k) * greenMatrix(i+L, j+L);
            }
        }
    }

    return phoMinusBg;
}

void RealMaterialSOCMeasureFixSDSD::addPhoSumBg(SDSDOperation &sdsdOperation, complex<double> denIncrement) {
    size_t L = realMaterialSOC->getL();
    size_t eigenNum = realMaterialSOC->getEigenNumber();

    if (phoSumBgNum.rank(0) != eigenNum) {
        phoSumBgNum.resize(eigenNum);
        phoSumBgNum = complex<double>(0, 0);
    }

    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    for (size_t k = 0; k < eigenNum; ++k) {
        for (size_t i = 0; i < L; ++i) {
            for (size_t j = 0; j < L; ++j) {
                phoSumBgNum(k) += (greenMatrix(j, i)  + greenMatrix(j + L, i + L) ) * denIncrement;
            }
        }
    }
}

void RealMaterialSOCMeasureFixSDSD::addPhoMinusBg(SDSDOperation &sdsdOperation, complex<double> denIncrement) {
    size_t L = realMaterialSOC->getL();
    size_t eigenNum = realMaterialSOC->getEigenNumber();

    if (phoMinusBgNum.rank(0) != eigenNum) {
        phoMinusBgNum.resize(eigenNum);
        phoMinusBgNum = complex<double>(0, 0);
    }

    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    for (size_t k = 0; k < eigenNum; ++k) {
        for (size_t i = 0; i < L; ++i) {
            for (size_t j = 0; j < L; ++j) {
                phoMinusBgNum(k) += (greenMatrix(j, i)  + greenMatrix(j + L, i + L) ) * denIncrement;
            }
        }
    }
}

void RealMaterialSOCMeasureFixSDSD::addNupNdn(SDSDOperation &sdsdOperation, std::complex<double> denIncrement)
{
    size_t L = realMaterialSOC->getL();
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    if( NupNum.rank(0) != L ) { NupNum.resize(L); NupNum = complex<double>(0,0); }
    if( NdnNum.rank(0) != L ) { NdnNum.resize(L); NdnNum = complex<double>(0,0); }

    for (size_t i = 0; i < L; ++i)
    {
        NupNum(i) += greenMatrix(i, i) * denIncrement;
        NdnNum(i) += greenMatrix(i+L, i+L) * denIncrement;
    }

    complex<double> diagUpAll(0,0), diagDnAll(0,0);
    for(size_t i = 0; i < L ; ++i)
    {
        diagUpAll += greenMatrix(i, i);
        diagDnAll += greenMatrix(i+L, i+L);
    }

    NupTotNum += diagUpAll * denIncrement;
    NdnTotNum += diagDnAll * denIncrement;
}

void RealMaterialSOCMeasureFixSDSD::addSplusSminus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement)
{
    size_t L = realMaterialSOC->getL();
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    if( SplusNum.rank(0) != L )  { SplusNum.resize(L);  SplusNum = complex<double>(0,0); }
    if( SminusNum.rank(0) != L ) { SminusNum.resize(L); SminusNum = complex<double>(0,0); }

    for (size_t i = 0; i < L; ++i)
    {
        SplusNum(i)  += greenMatrix(i, i+L) * denIncrement;
        SminusNum(i) += greenMatrix(i+L, i) * denIncrement;
    }

    complex<double> SplusAll(0,0), SminusAll(0,0);
    for(size_t i = 0; i < L; ++i)
    {
        SplusAll +=  greenMatrix(i, i+L);
        SminusAll += greenMatrix(i+L, i);
    }

    SplusTotNum  += SplusAll * denIncrement;
    SminusTotNum += SminusAll * denIncrement;
}
