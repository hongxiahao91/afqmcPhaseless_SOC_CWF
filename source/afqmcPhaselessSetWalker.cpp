//

#include "../include/afqmcPhaseless.h"

using namespace std;
using namespace tensor_hao;

void AfqmcPhaseless::initialPhiT()
{
    if(method.initialPhiTFlag == "setFromModel")
    {
        if( MPIRank()==0 ) fillWalkerFromModel(phiT, model);
        if( MPIRank()==0 ) phiT.write("phiT.dat");
        MPIBcast(phiT);
    }
    else if(method.initialPhiTFlag == "setRandomly")
    {
        if( MPIRank()==0 ) fillWalkerRandomly(phiT, model);
        if( MPIRank()==0 ) phiT.write("phiT.dat");
        MPIBcast(phiT);
    }
    else if(method.initialPhiTFlag == "readFromFile")
    {
        if( MPIRank()==0 ) phiT.read("phiT.dat");
        MPIBcast(phiT);
    }
    else
    {
        cout<<"Error!!! Do not recognize initialPhiTFlag!"<<endl;
        exit(1);
    }

    cout<<"Initiate PhiT Success!" <<endl;
}

void AfqmcPhaseless::initialWalker()
{
    walker.resize(method.walkerSizePerThread);
    walkerIsAlive.resize(method.walkerSizePerThread);

    if(method.initialWalkerFlag == "setFromModel")
    {
        if( MPIRank()==0 ) fillWalkerFromModel(walker[0], model);
        if( MPIRank()==0 ) checkOverlap( walker[0] );
        if( MPIRank()==0 ) walker[0].write("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "setRandomly")
    {
        if( MPIRank()==0 ) fillWalkerRandomly(walker[0], model);
        if( MPIRank()==0 ) checkOverlap( walker[0] );
        if( MPIRank()==0 ) walker[0].write("phi.dat");
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "sampleFromPhiT")
    {
        setWalkerFromPhiT(walker, walkerIsAlive, phiT);
    }
    else if(method.initialWalkerFlag == "readFromFile")
    {
        if( MPIRank()==0 ) walker[0].read("phi.dat");
        if( MPIRank()==0 ) checkOverlap( walker[0] );
        MPIBcast(walker[0]);

        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            walker[i] = walker[0];
            walkerIsAlive[i] = true;
        }
    }
    else if(method.initialWalkerFlag == "readAllWalkers")
    {
        string filename;
        int baseNumber = MPIRank() * method.walkerSizePerThread;
        for(int i = 0; i < method.walkerSizePerThread; ++i)
        {
            filename = "./walkers/phi_" + to_string(i+baseNumber) +".dat";
            if( checkFile(filename) ) { walker[i].read(filename); walkerIsAlive[i]=true; }
            else { fillWalkerRandomly(walker[i], model);  walkerIsAlive[i]=false; }
        }
    }
    else
    {
        cout<<"Error!!! Do not recognize initialWalkerFlag!"<<endl;
        exit(1);
    }

    initialMgsAndPopControl();
    cout<<"Initiate Walker Success!" <<endl;
}

void AfqmcPhaseless::writeWalkers()
{
    MPIBarrier();
    if( MPIRank() == 0 )
    {
        int flag;
        flag = system("rm -rf walkers");   if(flag != 0) cout<<"WARNING!!! system command does not exit properly!"<<endl;
        flag = system("mkdir -p walkers"); if(flag != 0) cout<<"WARNING!!! system command does not exit properly!"<<endl;
    }
    MPIBarrier();

    string filename;
    int baseNumber = MPIRank() * method.walkerSizePerThread;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            filename = "./walkers/phi_" + to_string(i+baseNumber) +".dat";
            walker[i].write(filename);
        }
    }
}

void AfqmcPhaseless::checkOverlap(WalkerRight &oneWalker)
{
    walkerWalkerOperation.set(phiT, oneWalker);
    complex<double> overlap = exp( walkerWalkerOperation.returnLogOverlap() );

    cout << "Overlp between trial wave function and walker is " << overlap << endl;

    if( abs(overlap)<1e-12 )
    {
        cout<<"Warning!!!!! walker has tiny overlap with trial wave function!"<<endl;
        cout<<"QMC results might be unreliable!!!!!"<<endl;
        cout<<"\n"<<endl;
    }
}

void AfqmcPhaseless::initialMgsAndPopControl()
{
    modifyGM();
    popControl(1.01);
}
