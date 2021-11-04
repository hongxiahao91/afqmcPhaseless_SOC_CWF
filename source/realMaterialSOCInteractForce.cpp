//
// Created by laura on 2019-09-30.
//

#include "../include/realMaterialSOCInteractForce.h"

using namespace std;

RealMaterialSOCInteractForce::RealMaterialSOCInteractForce()  { }

RealMaterialSOCInteractForce::RealMaterialSOCInteractForce(size_t eigenNum)
{
    PlusForce.resize(eigenNum);
    MinusForce.resize(eigenNum);
}

RealMaterialSOCInteractForce::RealMaterialSOCInteractForce(const RealMaterialSOCInteractForce &x) { copy_deep(x); }

RealMaterialSOCInteractForce::RealMaterialSOCInteractForce(RealMaterialSOCInteractForce &&x) { move_deep(x); }

RealMaterialSOCInteractForce::~RealMaterialSOCInteractForce() { }

RealMaterialSOCInteractForce &RealMaterialSOCInteractForce::operator=(const RealMaterialSOCInteractForce &x)  { copy_deep(x); return *this; }

RealMaterialSOCInteractForce &RealMaterialSOCInteractForce::operator=(RealMaterialSOCInteractForce &&x) { move_deep(x); return *this; }

size_t RealMaterialSOCInteractForce::getEigenNumber() const  { return PlusForce.size(); }

double RealMaterialSOCInteractForce::getMemory() const
{
    double mem(0.0);
    mem += PlusForce.getMemory();
    mem += MinusForce.getMemory();
    return mem;
}

void RealMaterialSOCInteractForce::copy_deep(const RealMaterialSOCInteractForce &x)
{
    PlusForce = x.PlusForce;
    MinusForce = x.MinusForce;
}

void RealMaterialSOCInteractForce::move_deep(RealMaterialSOCInteractForce &x)
{
    PlusForce = move( x.PlusForce );
    MinusForce = move( x.MinusForce );
}
