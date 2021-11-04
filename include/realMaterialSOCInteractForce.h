//

#ifndef AFQMCLAB_REALMATERIALSOCINTERACTFORCE_H
#define AFQMCLAB_REALMATERIALSOCINTERACTFORCE_H

#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/common/tensorHao/include/tensor_all.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/common/tensorHao/include/tensor_all.h"

class RealMaterialSOCInteractForce
{
public:
    tensor_hao::TensorHao<std::complex<double>, 1> PlusForce;
    tensor_hao::TensorHao<std::complex<double>, 1> MinusForce;

    RealMaterialSOCInteractForce();
    RealMaterialSOCInteractForce(size_t eigenNum);
    RealMaterialSOCInteractForce(const RealMaterialSOCInteractForce& x);
    RealMaterialSOCInteractForce(RealMaterialSOCInteractForce&& x);
    ~RealMaterialSOCInteractForce();

    RealMaterialSOCInteractForce & operator  = (const RealMaterialSOCInteractForce& x);
    RealMaterialSOCInteractForce & operator  = (RealMaterialSOCInteractForce&& x);

    size_t getEigenNumber() const;
    double getMemory() const;

private:
    void copy_deep(const RealMaterialSOCInteractForce &x);
    void move_deep(RealMaterialSOCInteractForce &x);
};

#endif //AFQMCLAB_REALMATERIALSOCINTERACTFORCE_H
