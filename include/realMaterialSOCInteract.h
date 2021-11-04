//
//

#ifndef AFQMCLAB_REALMATERIALSOCINTERACT_H
#define AFQMCLAB_REALMATERIALSOCINTERACT_H

#include "realMaterialSOCInteractAux.h"
#include "realMaterialSOCInteractForce.h"
#include "realMaterialSOCInteractSample.h"

class RealMaterialSOCInteract
{
 private:
    double dt;
    const tensor_hao::TensorHao<double, 1> *eigenValue;
    tensor_hao::TensorHao<double,1> dtEigenValue;
    size_t eigenNum;
    const tensor_hao::TensorHao<double, 3> *phoSum, *phoMinus;
    tensor_hao::TensorHao<std::complex<double>,1> plusGamma, minusGamma;
    const tensor_hao::TensorHao<double, 1> *phoSumBg;
    const tensor_hao::TensorHao<double, 1> *phoMinusBg;
    tensor_hao::TensorHao<std::complex<double>, 3> sqrtMinusGammaPhoSum, sqrtPlusGammaPhoMinus; /* 2L*2L*eigenNum */


 public:
    RealMaterialSOCInteract();
    RealMaterialSOCInteract(double dt,
                 const tensor_hao::TensorHao<double, 1> &eigenValue,
                 size_t eigenNum,
                 const tensor_hao::TensorHao<double, 3> &phoSum, /* 2L*2L*eigenNum */
                 const tensor_hao::TensorHao<double, 3> &phoMinus,
                 const tensor_hao::TensorHao<double, 1> &phoSumBg,
                 const tensor_hao::TensorHao<double, 1> &phoMinusBg);
    RealMaterialSOCInteract(const RealMaterialSOCInteract& x);
    RealMaterialSOCInteract(RealMaterialSOCInteract&& x);
    ~RealMaterialSOCInteract();

    RealMaterialSOCInteract & operator  = (const RealMaterialSOCInteract& x);
    RealMaterialSOCInteract & operator  = (RealMaterialSOCInteract&& x);


    //const tensor_hao::TensorHao<std::complex<double>, 3> &getSqrtMinusGammaPhoSum() const; /* 2L*2L*eigenNum  */
    //const tensor_hao::TensorHao<double, 1> *getPlusBg() const;
    //const tensor_hao::TensorHao<double, 1> *getMinusBg() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getPlusGamma() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getMinusGamma() const;

    std::complex<double> calculateAuxForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force);
    RealMaterialSOCInteractForce readForce(const std::string &filename) const;
    RealMaterialSOCInteractAux sampleAuxFromForce(const RealMaterialSOCInteractForce &force) const;
    std::complex<double> logProbOfAuxFromForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force) const;
    RealMaterialSOCInteractSample getTwoBodySampleFromAux(const RealMaterialSOCInteractAux &aux) const;
    RealMaterialSOCInteractSample getTwoBodySampleFromAuxForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force) const;

    double getMemory() const;
 private:
    void copy_deep(const RealMaterialSOCInteract &x);
    void move_deep(RealMaterialSOCInteract &x);

    void setPlusGamma();
    void setMinusGamma();

    void initialSqrtMinusGammaPhoSum(const tensor_hao::TensorHao<double, 3> &phoSum);
    void initialSqrtPlusGammaPhoMinus(const tensor_hao::TensorHao<double, 3> &phoMinus);
    void setTwoBodySampleMatrix(RealMaterialSOCInteractSample &realMaterialSOCInteractSample, const RealMaterialSOCInteractAux &aux) const;
};

#endif //AFQMCLAB_REALMATERIALSOCINTERACT_H
