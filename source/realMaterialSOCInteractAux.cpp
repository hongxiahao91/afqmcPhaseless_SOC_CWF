//
// Created by laura on 2019-09-30.
//

#include "../include/realMaterialSOCInteractAux.h"

using namespace std;

RealMaterialSOCInteractAux::RealMaterialSOCInteractAux() { }

RealMaterialSOCInteractAux::RealMaterialSOCInteractAux(size_t eigenNum)
{
    PlusAux.resize(eigenNum);
    MinusAux.resize(eigenNum);
}

RealMaterialSOCInteractAux::RealMaterialSOCInteractAux(const RealMaterialSOCInteractAux &x) { copy_deep(x); }

RealMaterialSOCInteractAux::RealMaterialSOCInteractAux(RealMaterialSOCInteractAux&&x) { move_deep(x); }

RealMaterialSOCInteractAux::~RealMaterialSOCInteractAux() { }

RealMaterialSOCInteractAux &RealMaterialSOCInteractAux::operator=(const RealMaterialSOCInteractAux &x)  { copy_deep(x); return *this; }

RealMaterialSOCInteractAux &RealMaterialSOCInteractAux::operator=(RealMaterialSOCInteractAux &&x) { move_deep(x); return *this; }

size_t RealMaterialSOCInteractAux::getEigenNumber() const { return PlusAux.size(); }

double RealMaterialSOCInteractAux::getMemory() const
{
    double mem(0.0);
    mem += PlusAux.getMemory();
    mem += MinusAux.getMemory();
    return mem;
}

void RealMaterialSOCInteractAux::copy_deep(const RealMaterialSOCInteractAux &x)
{
    PlusAux = x.PlusAux;
    MinusAux = x.MinusAux;
}

void RealMaterialSOCInteractAux::move_deep(RealMaterialSOCInteractAux &x)
{
    PlusAux = move( x.PlusAux );
    MinusAux = move( x.MinusAux );
}