//

#ifndef AFQMCLAB_REALMATERIALSOCINTERACTAUX_H
#define AFQMCLAB_REALMATERIALSOCINTERACTAUX_H

#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/common/tensorHao/include/tensor_all.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/common/tensorHao/include/tensor_all.h"

class RealMaterialSOCInteractAux
{
public:
    tensor_hao::TensorHao<std::complex<double>, 1> PlusAux;
    tensor_hao::TensorHao<std::complex<double>, 1> MinusAux;

    RealMaterialSOCInteractAux();
    RealMaterialSOCInteractAux(size_t eigenNum);
    RealMaterialSOCInteractAux(const RealMaterialSOCInteractAux& x);
    RealMaterialSOCInteractAux(RealMaterialSOCInteractAux&& x);
    ~RealMaterialSOCInteractAux();

    RealMaterialSOCInteractAux & operator  = (const RealMaterialSOCInteractAux& x);
    RealMaterialSOCInteractAux & operator  = (RealMaterialSOCInteractAux&& x);

    size_t getEigenNumber() const;
    double getMemory() const;

private:
    void copy_deep(const RealMaterialSOCInteractAux &x);
    void move_deep(RealMaterialSOCInteractAux &x);
};

#endif //AFQMCLAB_REALMATERIALSOCINTERACTAUX_H




