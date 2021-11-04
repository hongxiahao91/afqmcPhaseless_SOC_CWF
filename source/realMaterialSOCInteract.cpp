//

#include "/Users/laolabaobei/Research/Clion/AFQMCLAB_bitbucket/afqmclab/include/common/common.h"
#include "../include/realMaterialSOCInteract.h"
//#include "/users/hh8/script/AFQMCLAB-bitbucket-ccv3/afqmclab/include/common/common.h"

using namespace std;
using namespace tensor_hao;

#define pi 3.14159265358979324

RealMaterialSOCInteract::RealMaterialSOCInteract(){}

RealMaterialSOCInteract::RealMaterialSOCInteract(double dt,
                           const TensorHao<double, 1> &eigenValue,
                           size_t eigenNum,
                           const TensorHao<double, 3> &phoSum,
                           const TensorHao<double, 3> &phoMinus,
                           const TensorHao<double, 1> &phoSumBg,
                           const TensorHao<double, 1> &phoMinusBg)
{
    RealMaterialSOCInteract::dt = dt;
    RealMaterialSOCInteract::eigenValue = &eigenValue;
    dtEigenValue = dt*eigenValue;
    RealMaterialSOCInteract::eigenNum = eigenNum;
    RealMaterialSOCInteract::phoSum = &phoSum;
    RealMaterialSOCInteract::phoSumBg = &phoSumBg;
    RealMaterialSOCInteract::phoMinusBg = &phoMinusBg;

    setPlusGamma();
    setMinusGamma();

    initialSqrtMinusGammaPhoSum(phoSum);
    initialSqrtPlusGammaPhoMinus(phoMinus);
}

RealMaterialSOCInteract::RealMaterialSOCInteract(const RealMaterialSOCInteract &x) { copy_deep(x); }

RealMaterialSOCInteract::RealMaterialSOCInteract(RealMaterialSOCInteract &&x) { move_deep(x); }

RealMaterialSOCInteract::~RealMaterialSOCInteract() { }

RealMaterialSOCInteract &RealMaterialSOCInteract::operator = (const RealMaterialSOCInteract &x) { copy_deep(x); return *this; }

RealMaterialSOCInteract &RealMaterialSOCInteract::operator = (RealMaterialSOCInteract &&x) { move_deep(x); return *this; }


//const TensorHao<complex<double>, 3> &RealMaterialSOCInteract::getSqrtMinusDtCholeskyVecs() const { return sqrtMinusDtCholeskyVecs; }

//const TensorHao<double, 1> *RealMaterialSOCInteract::getPlusBg() const { return plusBg; }

//const TensorHao<double, 1> *RealMaterialSOCInteract::getMinusBg() const { return minusBg; }

const TensorHao<complex<double>, 1> &RealMaterialSOCInteract::getPlusGamma() const { return plusGamma; }

const TensorHao<complex<double>, 1> &RealMaterialSOCInteract::getMinusGamma() const { return minusGamma; }

std::complex<double> RealMaterialSOCInteract::calculateAuxForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force)
{
    complex<double> auxForce(0,0);

    for(size_t i = 0; i < eigenNum; ++i)
    {
        auxForce += aux.PlusAux(i) * force.PlusForce(i);
        auxForce += aux.MinusAux(i) * force.MinusForce(i);
    }

    return auxForce;
}

RealMaterialSOCInteractForce RealMaterialSOCInteract::readForce(const std::string &filename) const
{
    RealMaterialSOCInteractForce force(eigenNum);

    if( !checkFile(filename) )
    {
        force.PlusForce = 0.0;
        force.MinusForce = 0.0;
    }
    else
    {
        size_t eigenNum_in;

        ifstream file;
        file.open(filename, ios::in);
        if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

        readFile( eigenNum_in, file );
        if( eigenNum_in!=eigenNum )
        {
            cout<<"Error!!! Input force has different numberOfKana "
                <<" "<<eigenNum<<" "<<eigenNum_in<<" !"<<endl;
            exit(1);
        }

        readFile( eigenNum, force.PlusForce.data(), file );
        readFile( eigenNum, force.MinusForce.data(), file );

        file.close();
    }

    return force;
}

RealMaterialSOCInteractAux RealMaterialSOCInteract::sampleAuxFromForce(const RealMaterialSOCInteractForce &force) const
{
    if( eigenNum != force.getEigenNumber() ) { cout<<"Error!!! Force size does not consistent with eigenNumber! "<<force.getEigenNumber()<<endl; exit(1); }

    RealMaterialSOCInteractAux aux( eigenNum );
    for(size_t i = 0; i < eigenNum; ++i)
    {
        aux.PlusAux(i) = gaussianHao() + force.PlusForce(i);
        aux.MinusAux(i) = gaussianHao() + force.MinusForce(i);
    }

    return aux;
}

complex<double> RealMaterialSOCInteract::logProbOfAuxFromForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force) const
{
    if( eigenNum != aux.getEigenNumber() ) { cout<<"Error!!! Aux size does not consistent with eigenNumber! "<<aux.getEigenNumber()<<endl; exit(1); }
    if( eigenNum != force.getEigenNumber() ) { cout<<"Error!!! Force size does not consistent with eigenNumber! "<<force.getEigenNumber()<<endl; exit(1); }

    // Product: 1/sqrt(2.0*Pi) * Exp( -(x-force)^2 / 2.0 )
    complex<double> logProb(0,0), auxMinusForceSquarePlus(0,0), auxMinusForceSquareMinus(0,0), tmpPlus, tmpMinus;
    for(size_t i = 0; i < eigenNum; ++i)
    {
        tmpPlus = aux.PlusAux(i)-force.PlusForce(i);
        auxMinusForceSquarePlus += tmpPlus*tmpPlus;
        tmpMinus = aux.MinusAux(i)-force.MinusForce(i);
        auxMinusForceSquareMinus += tmpMinus*tmpMinus;
    }
    logProb = -0.5*log(2.0*pi)*eigenNum - 0.5*auxMinusForceSquarePlus -0.5*log(2.0*pi)*eigenNum - 0.5*auxMinusForceSquareMinus;

    return logProb;
}

RealMaterialSOCInteractSample RealMaterialSOCInteract::getTwoBodySampleFromAux(const RealMaterialSOCInteractAux &aux) const
{
    if( eigenNum != aux.getEigenNumber() ) { cout<<"Error!!! Aux size does not consistent with eigenNumber! "<<aux.getEigenNumber()<<endl; exit(1); }

    RealMaterialSOCInteractSample twoBodySample( sqrtMinusGammaPhoSum.rank(0) );

    // Product: 1/sqrt(2.0*Pi) * Exp( -x^2 / 2.0 ) * Exp( - minusGamma *x*B ), minusGamma = sqrt(-0.25*dt*Dr)
    // Product: 1/sqrt(2.0*Pi) * Exp( -x^2 / 2.0 ) * Exp( - plusGamma *x*B ), plusGamma = sqrt(0.25*dt*Dr)

    complex<double> aux2SumPlus(0,0), aux2SumMinus(0,0), auxBSumPlus(0,0), auxBSumMinus(0,0);
    for(size_t i = 0; i < eigenNum; ++i)
    {
        aux2SumPlus += aux.PlusAux(i)*aux.PlusAux(i);
        auxBSumPlus += aux.PlusAux(i) * minusGamma(i) * phoSumBg->operator()(i);
        aux2SumMinus += aux.MinusAux(i)*aux.MinusAux(i);
        auxBSumMinus += aux.MinusAux(i) * plusGamma(i) * phoMinusBg->operator()(i);
    }
    twoBodySample.logw = -0.5*log(2.0*pi)*eigenNum - 0.5*aux2SumPlus - auxBSumPlus-0.5*log(2.0*pi)*eigenNum - 0.5*aux2SumMinus - auxBSumMinus;

    setTwoBodySampleMatrix(twoBodySample, aux);

    return twoBodySample;
}

RealMaterialSOCInteractSample RealMaterialSOCInteract::getTwoBodySampleFromAuxForce(const RealMaterialSOCInteractAux &aux, const RealMaterialSOCInteractForce &force) const {
    if (eigenNum != aux.getEigenNumber()) {
        cout << "Error!!! Aux size does not consistent with choleskyNumber! " << aux.getEigenNumber() << endl;
        exit(1);
    }
    if (eigenNum != force.getEigenNumber()) {
        cout << "Error!!! Force size does not consistent with choleskyNumber! " << force.getEigenNumber() << endl;
        exit(1);
    }

    // Product: Exp( force^2/2 - x*force ) Exp( -minusGamma* x*B )
    // Product: Exp( force^2/2 - x*force ) Exp( -plusGamma* x*B )
    RealMaterialSOCInteractSample twoBodySample = getTwoBodySampleFromAux(aux);

    twoBodySample.logw = twoBodySample.logw - logProbOfAuxFromForce(aux, force);

    return twoBodySample;
}


double RealMaterialSOCInteract::getMemory() const
{
    double mem(0.0);
    mem += 8.0;
    mem += 8.0; //pointer
    mem += dtEigenValue.getMemory();
    mem += 8.0;
    mem += 8.0*2;
    mem += plusGamma.getMemory()+minusGamma.getMemory();
    mem += 8.0*2;
    mem += sqrtMinusGammaPhoSum.getMemory()+sqrtPlusGammaPhoMinus.getMemory();
    return 8.0+16.0+8.0+sqrtMinusGammaPhoSum.getMemory()+sqrtPlusGammaPhoMinus.getMemory()+8.0;
}

void RealMaterialSOCInteract::copy_deep(const RealMaterialSOCInteract &x)
{
    dt = x.dt;
    eigenValue = x.eigenValue;
    dtEigenValue = x.dtEigenValue;
    eigenNum = x.eigenNum;
    phoSum = x.phoSum;
    phoMinus = x.phoMinus;
    plusGamma = x.plusGamma;
    minusGamma = x.minusGamma;
    phoSumBg = x.phoSumBg;
    phoMinusBg = x.phoMinusBg;
    sqrtMinusGammaPhoSum = x.sqrtMinusGammaPhoSum;
    sqrtPlusGammaPhoMinus = x.sqrtPlusGammaPhoMinus;
}

void RealMaterialSOCInteract::move_deep(RealMaterialSOCInteract &x)
{
    dt = move( x.dt );
    eigenValue = move( x.eigenValue );
    dtEigenValue = move( x.dtEigenValue );
    eigenNum = move( x.eigenNum );
    phoSum = move( x.phoSum );
    phoMinus = move( x.phoMinus );
    plusGamma = move( x.plusGamma );
    minusGamma = move( x.minusGamma );
    phoSumBg = move( x.phoSumBg );
    phoMinusBg = move( x.phoMinusBg );
    sqrtMinusGammaPhoSum = move( x.sqrtMinusGammaPhoSum );
    sqrtPlusGammaPhoMinus = move( x.sqrtPlusGammaPhoMinus );
}

void RealMaterialSOCInteract::setPlusGamma()
{
    plusGamma.resize(eigenNum);

    //sqrt(0.25*dt*Dr)
    for(size_t i = 0; i < eigenNum; ++i) plusGamma(i) = sqrt( complex<double>(0.25, 0.0) * dtEigenValue(i) );

}

void RealMaterialSOCInteract::setMinusGamma()
{
    minusGamma.resize(eigenNum);

    //sqrt(-0.25*dt*Dr)
    for(size_t i = 0; i < eigenNum; ++i) minusGamma(i) = sqrt( -complex<double>(0.25, 0.0) * dtEigenValue(i) );

}

void RealMaterialSOCInteract::initialSqrtMinusGammaPhoSum(const tensor_hao::TensorHao<double, 3> &phoSum)
{
    //sqrt(-0.25*Dt*Dr)
    sqrtMinusGammaPhoSum.resize( phoSum.getRank() );
    complex<double> *p0 = sqrtMinusGammaPhoSum.data();
    const double *p1 = phoSum.data();

    //There is a cubic scaling, can be faster by OpenMP
    for(size_t i = 0; i < eigenNum; ++i) p0[i] = p1[i] * minusGamma(i);
}

void RealMaterialSOCInteract::initialSqrtPlusGammaPhoMinus(const tensor_hao::TensorHao<double, 3> &phoMinus)
{
    //sqrt(0.25*Dt*Dr)
    sqrtPlusGammaPhoMinus.resize( phoMinus.getRank() );
    complex<double> *p0 = sqrtPlusGammaPhoMinus.data();
    const double *p1 = phoMinus.data();

    //There is a cubic scaling, can be faster by OpenMP
    for(size_t i = 0; i < eigenNum; ++i) p0[i] = p1[i] * plusGamma(i);
}

//TODO: check setTwoBodySampleMatrix function
void RealMaterialSOCInteract::setTwoBodySampleMatrix(RealMaterialSOCInteractSample &twoBodySample, const RealMaterialSOCInteractAux &aux) const
{
    //TensorHao<complex<double>, 2> &matrix = twoBodySample.matrix;
    //cout<<"matrix.size()"<<matrix.size()<<endl;
    //matrix = complex<double>(0.0, 0.0);

    //Calculate aux * sqrtMinusDt * choleskyVecs
    size_t L = sqrtMinusGammaPhoSum.rank(0);  /*  2L  */
    size_t L2 = L * L;
    //TensorHao<complex<double>, 1> vecsAux(L2);
    TensorHaoRef<complex<double>, 1> vecsAux(L2);
    vecsAux.point( twoBodySample.matrix.data() );
    TensorHaoRef<complex<double>, 2> vecsPlus(L2, eigenNum);
    vecsPlus.point( const_cast<complex<double>*> ( sqrtMinusGammaPhoSum.data() ) );
    BL_NAME(gemv)(vecsPlus, aux.PlusAux, vecsAux);  // vecsAux = 1*vecsPlus*PlusAux + 1* vecsAux

    TensorHao<complex<double>, 1> vecsAuxMinus(L2);
    TensorHaoRef<complex<double>, 2> vecsMinus(L2, eigenNum);
    vecsMinus.point( const_cast<complex<double>*> ( sqrtPlusGammaPhoMinus.data() ) );
    BL_NAME(gemv)(vecsMinus, aux.MinusAux, vecsAux, 'N', 1.0, 1.0);  // vecsAux = 1*vecsMinus*MinusAux + 1* vecsAux


    //vecsAux2.point( realMaterialSOCInteractSample.matrix.data() );
    //vecsAux.resize(L,L);
    //cout<<"matrix in twobodysamplematrix: "<<matrix<<endl;

}

