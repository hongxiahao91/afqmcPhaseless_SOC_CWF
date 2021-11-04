//

#ifndef AFQMCLAB_REALMATERIALSOCSDOPERATION_H
#define AFQMCLAB_REALMATERIALSOCSDOPERATION_H

#include "realMaterialSOC.h"
#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/afqmc/blocks/walker/SD/include/SD.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/afqmc/blocks/walker/SD/include/SD.h"

void fillWalkerRandomly(SD &walker, const RealMaterialSOC &model);
void fillWalkerFromModel(SD &walker, RealMaterialSOC &model);

#endif //AFQMCLAB_REALMATERIALSOCSDSOPERATION_H
