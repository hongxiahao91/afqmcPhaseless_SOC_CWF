//
// Created by Hongxia Hao on 9/26/18.
//

#include "../include/realMaterialSOCSDOperation.h"

using namespace std;
using namespace tensor_hao;

void fillWalkerRandomly(SD &walker, const RealMaterialSOC &model)
{
    size_t L = model.getL();
    size_t N = model.getN();
    walker.resize(2*L, N);

    walker.randomFill();
}

void fillWalkerFromModel(SD &walker, RealMaterialSOC &model)
{
    size_t L = model.getL();
    size_t N = model.getN();
    walker.resize(2*L, N);
    //cout<<"Success at this point!" <<endl;
    model.setKpEigenValueAndVector();
    //cout<<"Success after setKpEigenValueAndVector"<<endl;
    const TensorHao< complex<double>, 2 > &KpEigenVector = model.getKpEigenVector();
    //cout<<"KpEigenVector at fill WalkerFromModel: "<<KpEigenVector<<endl;
    TensorHao<complex<double>,2> &wf = walker.wfRef();

    copy( KpEigenVector.data(), KpEigenVector.data()+2*L*N, wf.data() );
    //cout<<"wf: "<<wf.data()<<endl;
}