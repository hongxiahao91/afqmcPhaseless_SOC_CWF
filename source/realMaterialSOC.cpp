//
// Created by laura on 9/26/18.
//

#include "../include/realMaterialSOC.h"

using namespace std;
using namespace H5;
using namespace tensor_hao;

RealMaterialSOC::RealMaterialSOC():L(0), N(0), Nup(0), Ndn(0), eigenNum(0),KpEigenStatus(0) { }

RealMaterialSOC::RealMaterialSOC(const string &filename) { read(filename); }

RealMaterialSOC::~RealMaterialSOC() { }

size_t RealMaterialSOC::getL() const { return L; }

size_t RealMaterialSOC::getN() const { return N; }

size_t RealMaterialSOC::getNup() const { return Nup; }

size_t RealMaterialSOC::getNdn() const { return Ndn; }

size_t RealMaterialSOC::getEigenNumber() const { return eigenNum; }

const TensorHao<complex<double>,2> &RealMaterialSOC::getT() const { return t; }

const TensorHao<complex<double>,2> &RealMaterialSOC::getK() const { return K; }

const TensorHao<double,4> &RealMaterialSOC::getPupup() const { return P_upup; }

const TensorHao<double,4> &RealMaterialSOC::getPdndn() const { return P_dndn; }

const TensorHao<double,4> &RealMaterialSOC::getPupdn() const { return P_updn; }

//const TensorHao<double,3> &RealMaterialSOC::getphoUpSum() const { return phoUpSum; }

//const TensorHao<double,3> &RealMaterialSOC::getphoUpMinus() const { return phoUpMinus; }

//const TensorHao<double,3> &RealMaterialSOC::getphoDnSum() const { return phoDnSum; }

//const TensorHao<double,3> &RealMaterialSOC::getphoDnMinus() const { return phoDnMinus; }

const TensorHao<double,3> &RealMaterialSOC::getphoSum() const { return phoSum; }

const TensorHao<double,3> &RealMaterialSOC::getphoMinus() const { return phoMinus; }

const TensorHao<double,1> &RealMaterialSOC::getphoSumBg() const { return phoSumBg; }

const TensorHao<double,1> &RealMaterialSOC::getphoMinusBg() const { return phoMinusBg; }

size_t RealMaterialSOC::getKpEigenStatus() const { return KpEigenStatus; }

const TensorHao<complex<double>,2> &RealMaterialSOC::getKp() const { return Kp; }

const TensorHao<double,1> &RealMaterialSOC::getKpEigenValue() const { return KpEigenValue; }

const TensorHao<complex<double>, 2> &RealMaterialSOC::getKpEigenVector() const { return KpEigenVector; }

void RealMaterialSOC::read(const string &filename)
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    readFile(L, file); //NMO
    readFile(N, file);
    readFile(Nup, file);
    readFile(Ndn, file);
    readFile(eigenNum, file);
    readFile(SumEngShf, file);
    t.resize(2*L, 2*L); readFile( t.size(),  t.data(),  file);
    K.resize(2*L, 2*L); readFile( K.size(),  K.data(),  file);
    eigenValue.resize(eigenNum);readFile(eigenValue.size(),eigenValue.data(),file);
    phoSum.resize(2*L, 2*L, eigenNum); readFile(phoSum.size(), phoSum.data(), file);
    phoMinus.resize(2*L, 2*L, eigenNum); readFile(phoMinus.size(), phoMinus.data(), file);
    P_upup.resize(L, L, L, L); readFile(P_upup.size(), P_upup.data(), file);
    P_dndn.resize(L, L, L, L); readFile(P_dndn.size(), P_dndn.data(), file);
    P_updn.resize(L, L, L, L); readFile(P_updn.size(), P_updn.data(), file);
    phoSumBg.resize(eigenNum); readFile(phoSumBg.size(), phoSumBg.data(), file);
    phoMinusBg.resize(eigenNum); readFile(phoMinusBg.size(), phoMinusBg.data(), file);
    file.close();

    cout<<"SumEngShf: "<<SumEngShf<<endl;

    //cout<<K<<endl;
    KpEigenStatus = 0;
    Kp.resize(0,0);
    KpEigenValue.resize( static_cast<size_t>(0) );
    KpEigenVector.resize( 0, 0 );
}

void RealMaterialSOC::checkpotential(const string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    for(size_t i = 0; i < L; ++i){
        for(size_t l = 0; l < L; ++l){
            for(size_t j = 0; j < L; ++j){
                for(size_t k = 0; k < L; ++k){
                    file<<i+1<<"\t"<<l+1<<"\t"<<j+1<<"\t"<<k+1<<"\t"<<P_upup(i,l,j,k)<<"\n";
                }
            }
        }
    }
}

void RealMaterialSOC::checkkinetic(const string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    for(size_t i = 0; i < 2*L; ++i){
        for(size_t l = 0; l < 2*L; ++l){
                    file<<i+1<<"\t"<<l+1<<"\t"<<t(i,l)<<"\n";
        }
    }
}

//void RealMaterialSOC::write(const string &filename) const
//{
//    H5File file(filename, H5F_ACC_TRUNC);
//
//    writeFile( L, file, "L" );
//    writeFile( N, file, "N" );
//    writeFile( Nup, file, "Nup" );
//    writeFile( Ndn, file, "Ndn" );
//    writeFile( choleskyNumber, file, "choleskyNumber" );
//    writeFile( t.size(),  t.data(),  file, "t" );
//    writeFile( K.size(),  K.data(),  file, "K" );
//    writeFile( choleskyVecs.size(), choleskyVecs.data(), file, "choleskyVecs" );
//    writeFile( choleskyBg.size(), choleskyBg.data(), file, "choleskyBg" );
//
//    file.close();
//}

void RealMaterialSOC::write(const string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}


    writeFile( L, file);
    writeFile( N, file);
    writeFile( Nup, file);
    writeFile( Ndn, file);
    writeFile( eigenNum, file);
    writeFile(SumEngShf, file);
    writeFile( t.size(),  t.data(),  file);
    writeFile( K.size(),  K.data(),  file);
    writeFile( eigenValue.size(),  eigenValue.data(),  file);
    writeFile( phoSum.size(), phoSum.data(), file);
    writeFile( phoMinus.size(), phoMinus.data(), file);
    writeFile( P_upup.size(), P_upup.data(), file);
    writeFile( P_dndn.size(), P_dndn.data(), file);
    writeFile( P_updn.size(), P_updn.data(), file);
    writeFile( phoSumBg.size(), phoSumBg.data(), file);
    writeFile( phoMinusBg.size(), phoMinusBg.data(), file);

    file.close();
}

void RealMaterialSOC::readFCIDUMP()  //real orbital, 8 fold symmetry, 5 count
{
    double elementReal;
    int i, j, k, l;
    int count = 0;

    t.resize(2*L,2*L);
    for(size_t i = 0; i < 2*L; ++i){
        for(size_t l = 0; l < 2*L; ++l){
            t(l, i) = 0;
        }
    }
    P_upup.resize(L,L,L,L);
    P_dndn.resize(L,L,L,L);
    P_updn.resize(L,L,L,L);
    for(size_t i = 0; i < L; ++i) {
        for (size_t l = 0; l < L; ++l) {
            for (size_t j = 0; j < L; ++j) {
                for (size_t k = 0; k < L; ++k) {
                    P_upup(i, l, j, k) = 0;
                    P_dndn(i, l, j, k) = 0;
                    P_updn(i, l, j, k) = 0;
                }
            }
        }
    }
    tensor_hao::TensorHao<std::complex<double>,1> EngShf;
    EngShf.resize(6);
    for(size_t i = 0; i < 6; ++i) {EngShf(i)=0;}


    FILE *fcidump = fopen("./FCIDUMP", "r");
    char ignore[1024];
    int norb, nelec;

    /*      1a  Read norb and nelec */
    fscanf(fcidump, "%*[^0-9]%d%*[^0-9]%d%*s\n", &norb, &nelec); // scan first line
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);

    SumEngShf = 0.0;

    while (fscanf(fcidump, "%lf %i %i %i %i\n", &elementReal, &i, &l, &j, &k) == 5){
        if (i == 0 && j == 0 && k == 0 && l == 0){
            EngShf(count) = complex<double>(elementReal, 0);
            SumEngShf += EngShf(count);
            count++;
            continue;
        }

        switch (count){
            case 0:  // up-up 8-fold symmetry
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_upup(i, l, j, k) = elementReal;
                P_upup(k, j, l, i) = elementReal;
                P_upup(j, k, i, l) = elementReal;
                P_upup(l, i, k, j) = elementReal;
                P_upup(i, l, k, j) = elementReal;
                P_upup(j, k, l, i) = elementReal;
                P_upup(k, j, i, l) = elementReal;
                P_upup(l, i, j, k) = elementReal;
                cout << P_upup(i, l, j, k) << endl;
                break;
            case 1: // dn-dn 8-fold symmetry
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_dndn(i, l, j, k) = elementReal;
                P_dndn(k, j, l, i) = elementReal;
                P_dndn(j, k, i, l) = elementReal;
                P_dndn(l, i, k, j) = elementReal;
                P_dndn(i, l, k, j) = elementReal;
                P_dndn(j, k, l, i) = elementReal;
                P_dndn(k, j, i, l) = elementReal;
                P_dndn(l, i, j, k) = elementReal;
                break;
            case 2: //up-dn 4-fold symmetry, we don't have exchange particle symmetry anymore.
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_updn(i, l, j, k) = elementReal;
                P_updn(l, i, k, j) = elementReal; //dumy indice
                P_updn(l, i, j, k) = elementReal;
                P_updn(i, l, k, j) = elementReal;
                break;
            case 3:
                i -= 1; l -= 1;
                t(i, l) = complex<double>(elementReal, 0);
                t(l, i) = complex<double>(elementReal, 0);
                break;
            case 4:
                i -= 1; l -= 1;
                t(i+L, l+L) = complex<double>(elementReal, 0);
                t(l+L, i+L) = complex<double>(elementReal, 0);
                break;
        }
    }

    fclose(fcidump);
    cout <<"SumEnf: "<<SumEngShf<<endl;
    cout <<"P_upup(11,2,11,2): "<<P_upup(11,2,11,2)<<endl;
}


void RealMaterialSOC::readFCIDUMPSOC()  //real orbital, 8 fold symmetry, 5 count
{
    double elementReal, elementImg;
    int i, j, k, l;
    int count = 0;

    t.resize(2*L,2*L);
    for(size_t i = 0; i < 2*L; ++i){
        for(size_t l = 0; l < 2*L; ++l){
            t(l, i) = 0;
        }
    }
    P_upup.resize(L,L,L,L);
    P_dndn.resize(L,L,L,L);
    P_updn.resize(L,L,L,L);
    for(size_t i = 0; i < L; ++i) {
        for (size_t l = 0; l < L; ++l) {
            for (size_t j = 0; j < L; ++j) {
                for (size_t k = 0; k < L; ++k) {
                    P_upup(i, l, j, k) = 0;
                    P_dndn(i, l, j, k) = 0;
                    P_updn(i, l, j, k) = 0;
                }
            }
        }
    }
    tensor_hao::TensorHao<std::complex<double>,1> EngShf;
    EngShf.resize(10);
    for(size_t i = 0; i < 10; ++i) {EngShf(i)=0;}


    FILE *fcidump = fopen("./FCIDUMP", "r");
    char ignore[1024];
    int norb, nelec;

    /*      1a  Read norb and nelec */
    fscanf(fcidump, "%*[^0-9]%d%*[^0-9]%d%*s\n", &norb, &nelec); // scan first line
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);
    fgets(ignore, sizeof(ignore), fcidump);

    SumEngShf = 0.0;

    while (fscanf(fcidump, "%lf %lf %i %i %i %i\n", &elementReal, &elementImg, &i, &l, &j, &k) == 6){
        if (i == 0 && j == 0 && k == 0 && l == 0){
            EngShf(count) = complex<double>(elementReal, elementImg);
            SumEngShf += EngShf(count);
            count++;
            continue;
        }

        switch (count){
            case 0:  // up-up 8-fold symmetry
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_upup(i, l, j, k) = elementReal;
                cout << P_upup(i, l, j, k) << endl;
                break;
            case 1: // dn-dn 8-fold symmetry
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_dndn(i, l, j, k) = elementReal;
                break;
            case 2: //up-dn 4-fold symmetry, we don't have exchange particle symmetry anymore.
                i -= 1; l -= 1; j -= 1; k-= 1;
                P_updn(i, l, j, k) = elementReal;
                break;
            case 3:
                i -= 1; l -= 1;
                t(i, l) = complex<double>(elementReal, elementImg);
                break;
            case 4:
                i -= 1; l -= 1;
                t(i+L, l+L) = complex<double>(elementReal, elementImg);
                break;
            case 5:
                i -= 1; l -= 1;
                t(i,   l+L) = complex<double>(elementReal, elementImg);
                t(i+L, l  ) = complex<double>(elementReal, elementImg);
                break;
        }
    }

    fclose(fcidump);
    cout <<"SumEnf: "<<SumEngShf<<endl;
    cout <<"P_upup(11,2,11,2): "<<P_upup(11,2,11,2)<<endl;
}

#ifdef MPI_HAO
void MPIBcast(RealMaterialSOC &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast(buffer.L, root, comm);
    MPIBcast(buffer.N, root, comm);
    MPIBcast(buffer.Nup, root, comm);
    MPIBcast(buffer.Ndn, root, comm);
    MPIBcast(buffer.eigenNum, root, comm);
    MPIBcast(buffer.SumEngShf, root, comm);
    MPIBcast(buffer.t, root, comm);
    MPIBcast(buffer.K, root, comm);
    MPIBcast(buffer.eigenValue, root, comm);
    MPIBcast(buffer.phoSum, root, comm);
    MPIBcast(buffer.phoMinus, root, comm);
    MPIBcast(buffer.P_upup, root, comm);
    MPIBcast(buffer.P_dndn, root, comm);
    MPIBcast(buffer.P_updn, root, comm);
    MPIBcast(buffer.phoSumBg, root, comm);
    MPIBcast(buffer.phoMinusBg, root, comm);
    MPIBcast(buffer.KpEigenStatus, root, comm);
    MPIBcast(buffer.Kp, root, comm);
    MPIBcast(buffer.KpEigenValue, root, comm);
    MPIBcast(buffer.KpEigenVector, root, comm);
}
#endif

//void RealMaterialSOC::writeBackGround(const string &filename) const
//{
//    H5File file(filename, H5F_ACC_RDWR);
//    writeFile( choleskyBg.size(), choleskyBg.data(), file, "choleskyBg");
//    file.close();
//}

void RealMaterialSOC::writeBackGround(const string &filename) const
{
    write(filename);
}

void RealMaterialSOC::updateBackGround(const TensorHao<double, 1> &backgroundOne, const TensorHao<double, 1> &backgroundTwo)
{
    if( backgroundOne.size() != eigenNum ) {cout<<"Error!!! BackgroundOne size is not eigenNumber!"<<endl; exit(1);}
    if( backgroundTwo.size() != eigenNum ) {cout<<"Error!!! BackgroundTwo size is not eigenNumber!"<<endl; exit(1);}
    KpEigenStatus = 0;
    phoSumBg = backgroundOne;
    phoMinusBg = backgroundTwo;

}

void RealMaterialSOC::updateBackGround(TensorHao<double, 1> &&backgroundOne, TensorHao<double, 1> &&backgroundTwo)
{
    if( backgroundOne.size() != eigenNum ) {cout<<"Error!!! BackgroundOne size is not eigenNumber!"<<endl; exit(1);}
    if( backgroundTwo.size() != eigenNum ) {cout<<"Error!!! BackgroundTwo size is not eigenNumber!"<<endl; exit(1);}
    KpEigenStatus = 0;
    phoSumBg = move(backgroundOne);
    phoMinusBg = move(backgroundTwo);
}

Hop RealMaterialSOC::returnExpMinusAlphaK(double alpha)
{
    setKpEigenValueAndVector();

    Hop hop(2*L);
    double bg2Plus(0.0), bg2Minus(0.0);
    for(size_t i = 0; i < eigenNum; ++i)
    {
        bg2Plus  -= eigenValue(i) * phoSumBg(i) * phoSumBg(i);
        bg2Minus += eigenValue(i) * phoMinusBg(i) * phoMinusBg(i);
    }
    hop.logw = -alpha*0.125*(bg2Plus+bg2Minus);
    cout<<"returnExpMinusAlphaK.logw: "<<hop.logw<<endl;

    //TensorHao<complex<double>,2> matrix(2*L,2*L);
    BL_NAME(gmm)( KpEigenVector, dMultiMatrix( exp(-alpha*KpEigenValue), trans(KpEigenVector) ), hop.matrix );

    return hop;
}

LogHop RealMaterialSOC::returnLogExpMinusAlphaK(double alpha)
{
    setKp();

    LogHop logHop(2*L);
    double bg2Plus(0.0), bg2Minus(0.0);
    for(size_t i = 0; i < eigenNum; ++i)
    {
        bg2Plus  -= eigenValue(i) * phoSumBg(i) * phoSumBg(i);
        bg2Minus += eigenValue(i) * phoMinusBg(i) * phoMinusBg(i);
    }
    logHop.logw = -alpha*0.125*(bg2Plus+bg2Minus);

    logHop.matrix = complex<double>(-alpha, 0) * Kp;
    return logHop;
}

RealMaterialSOCInteract RealMaterialSOC::returnExpMinusAlphaV(double alpha)
{
    return RealMaterialSOCInteract(alpha, eigenValue, eigenNum, phoSum, phoMinus, phoSumBg, phoMinusBg);
}

/*
void RealMaterialSOC::setpho()
{
    phoSum.resize(2*L,2*L,eigenNum);
    phoMinus.resize(2*L,2*L,eigenNum);

    for (size_t k =0; k < eigenNum; k++)
    {
        for (size_t i = 0; i < L; i++)
        {
            for (size_t j = 0; j< L; j++)
            {
                phoSum(i,   j,   k) = phoUpSum(i,j,k);
                phoSum(i+L, j+L, k) = phoDnSum(i,j,k);
                phoMinus(i,   j ,  k) = phoUpMinus(i,j,k);
                phoMinus(i+L, j+L, k) = phoDnMinus(i,j,k);

            }
        }
    }
}
*/

void RealMaterialSOC::setKp()
{
    if( KpEigenStatus >=1 ) return;

    size_t L2 = 2*L;

    Kp = K;
    //cout<<"read Initial K: "<<Kp<<endl;

    TensorHao<double,1> phoSumBgNew(eigenNum);
    TensorHao<double,1> phoMinusBgNew(eigenNum);

    for (size_t i = 0; i < eigenNum; i++)
    {
        phoSumBgNew(i) = 0.25*eigenValue(i)*phoSumBg(i);
        phoMinusBgNew(i) = -0.25*eigenValue(i)*phoMinusBg(i);

    }

    TensorHaoRef<double,2> vecsSum(L2*L2, eigenNum); vecsSum.point( phoSum.data() );
    TensorHaoRef<double,2> vecsMinus(L2*L2, eigenNum); vecsMinus.point( phoMinus.data() );
    //TensorHao<double,2> vecs(L2, choleskyNumber); choleskyVecs.resize(L2, choleskyNumber);
    //TensorHaoRef<double,1> vecsSumBg(L2*L2); vecsSumBg.point( Kp.data() );
    //TensorHaoRef<double,1> vecsMinusBg(L2*L2); vecsMinusBg.point( Kp.data() );

    TensorHao<double,1> vecsBg(L2*L2);

    BL_NAME(gemv)(vecsSum, phoSumBgNew, vecsBg, 'N', 1.0, 1.0);  // vecsBg = 1*vecsSum*phoSumBgNew + 1* vecBg
    BL_NAME(gemv)(vecsMinus, phoMinusBgNew, vecsBg, 'N', 1.0, 1.0); //vecsBg = 1*vecsMinus*phoMinusBgNew + 1*vecBg
    cout<<"vecsSumBg: "<<vecsBg<<endl;

    for(size_t i = 0; i < L2; i++)
    {
        for(size_t j = 0; j < L2; j++)
        {
            Kp(i,  j) +=vecsBg(i+j*L);
        }
    }
    KpEigenStatus=1;
    //checkSymmetry(Kp, 1e-8);
    //checkHermitian(Kp, 1e-8);
    //cout<<"Kp: "<<Kp <<endl;
}

void RealMaterialSOC::setKpEigenValueAndVector()
{
    if( KpEigenStatus >=2 ) return;

    setKp();
    //cout<<"Success after setKp, Kp is: "<<Kp<<endl;
    KpEigenVector = Kp;
    KpEigenValue.resize(2*L);
    BL_NAME(eigen)(KpEigenVector, KpEigenValue);

    cout<<"KpEigenValue: "<<KpEigenValue<<endl;
    KpEigenStatus = 2;
}

double RealMaterialSOC::getMemory() const
{
    double mem(0.0);

    mem += 8.0*5;
    mem += t.getMemory() + K.getMemory();
    mem += P_upup.getMemory() + P_dndn.getMemory() + P_updn.getMemory();
    mem += eigenValue.getMemory();
    mem += 16.0;
    mem += phoSum.getMemory() + phoMinus.getMemory();
    mem += phoSumBg.getMemory() + phoMinusBg.getMemory();
    mem += 8.0;
    mem += Kp.getMemory();
    mem += KpEigenValue.getMemory();
    mem += KpEigenVector.getMemory();

    return mem;
}

RealMaterialSOC::RealMaterialSOC(const RealMaterialSOC &x)  { }

RealMaterialSOC &RealMaterialSOC::operator=(const RealMaterialSOC &x) { return *this; }
