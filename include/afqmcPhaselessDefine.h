//

#ifndef AFQMCLAB_AFQMCPHASELESSDEFINE_H
#define AFQMCLAB_AFQMCPHASELESSDEFINE_H

#include "afqmclab.h"
#include "realMaterialSOCAll.h"

typedef Hop OneBody;

typedef RealMaterialSOCInteract       TwoBody;
typedef RealMaterialSOCInteractAux    TwoBodyAux;
typedef RealMaterialSOCInteractForce  TwoBodyForce;
typedef RealMaterialSOCInteractSample TwoBodySample;

typedef RealMaterialSOC Model;
typedef RealMaterialSOCMeasureFixSDSD ModelMeasureMixed;

//typedef SD2s   WalkerLeft;
//typedef SD2is  WalkerRight;
//typedef Hop2isSD2isOperation OneBodyWalkerRightOperation;
//typedef CholeskyRealSampleSD2isOperation TwoBodySampleWalkerRightOperation;
//typedef SD2sSD2isOperation WalkerWalkerOperation;
//typedef RealMaterialMoleculeMeasureFixedSD2sSD2is ModelMeasureMixed;

//typedef MDCas2s WalkerLeft;
//typedef SD2is  WalkerRight;
//typedef Hop2isSD2isOperation OneBodyWalkerRightOperation;
//typedef CholeskyRealSampleSD2isOperation TwoBodySampleWalkerRightOperation;
//typedef MDCas2sSD2isOperation WalkerWalkerOperation;
//typedef RealMaterialMoleculeMeasureFixedMDCas2sSD2is ModelMeasureMixed;

typedef SD WalkerLeft;
typedef SD  WalkerRight;
typedef HopSDOperation OneBodyWalkerRightOperation;
typedef LogHopSDOperation TwoBodySampleWalkerRightOperation;
typedef SDSDOperation WalkerWalkerOperation;


#endif //AFQMCLAB_AFQMCPHASELESSDEFINE_H
