//
// Created by Hongxia Hao on 10/3/18.
//

#ifndef AFQMCLAB_REALMATERIALSOCMEASUREFIXSDSD_H
#define AFQMCLAB_REALMATERIALSOCMEASUREFIXSDSD_H

#include "realMaterialSOC.h"
#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/afqmc/blocks/walkerWalkerOperation/SDSDOperation/include/SDSDOperation.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/afqmc/blocks/walkerWalkerOperation/SDSDOperation/include/SDSDOperation.h"

class RealMaterialSOCMeasureFixSDSD
{
 private:
    const RealMaterialSOC *realMaterialSOC; //L= NMO
    const SD  *walkerLeft; //L= SMO

    std::complex<double> den;
    std::complex<double> TNum, HNum;
    tensor_hao::TensorHao<std::complex<double>, 1> phoSumBgNum, phoMinusBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> NupNum, NdnNum, SplusNum, SminusNum;
    std::complex<double> NupTotNum, NdnTotNum, SplusTotNum, SminusTotNum;

    tensor_hao::TensorHao<std::complex<double>,2> wfDaggerT;

 public:
    RealMaterialSOCMeasureFixSDSD();
    RealMaterialSOCMeasureFixSDSD(const RealMaterialSOC& realMaterialSOC_, const SD &walkerLeft_);
    ~RealMaterialSOCMeasureFixSDSD();

    void initModelWalkerNullptr();
    void setModelWalker(const RealMaterialSOC& realMaterialSOC_, const SD &walkerLeft_);
    void reSet();
    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnPhoSumBg();
    tensor_hao::TensorHao<double, 1> returnPhoMinusBg();
    void addMeasurement(SDSDOperation &sdSDOperation, std::complex<double> denIncrement);
    RealMaterialSOCInteractForce getForce(const RealMaterialSOCInteract &realMaterialSOCInteract, SDSDOperation &sdSDOperation, double cap=1.0);

    void write() const;
    double getMemory() const;
 private:
    RealMaterialSOCMeasureFixSDSD(const RealMaterialSOCMeasureFixSDSD& x);
    RealMaterialSOCMeasureFixSDSD & operator  = (const RealMaterialSOCMeasureFixSDSD& x);

    void initWfDaggerT();
    void checkWalkerLeft(const SDSDOperation &sdSDOperation);

    void addEnergy(SDSDOperation &sdSDOperation, std::complex<double> denIncrement);
    std::complex<double> calculateTenergy(SDSDOperation &sdSDOperation);
    std::complex<double> calculateVenergy(SDSDOperation &sdsdOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculatePhoSumBg(SDSDOperation &sdSDOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculatePhoMinusBg(SDSDOperation &sdSDOperation);
    void addPhoSumBg(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addPhoMinusBg(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addNupNdn(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addSplusSminus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
};

#endif //AFQMCLAB_REALMATERIALMOLECULEFIXEDSD2SSD2IS_H