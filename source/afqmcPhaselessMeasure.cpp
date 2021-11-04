//

#include "../include/afqmcPhaseless.h"

using namespace std;
using namespace tensor_hao;

void AfqmcPhaseless::addPhiTMeasure()
{
    complex<double> logOverlap;

    walkerWalkerOperation.set(phiT, phiT);
    logOverlap = walkerWalkerOperation.returnLogOverlap();

    mixedMeasure.addMeasurement( walkerWalkerOperation, exp(logOverlap) );
}

void AfqmcPhaseless::addPureMeasure()
{
    complex<double> logOverlap;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            walkerWalkerOperation.set(phiT, walker[i]);
            logOverlap = walkerWalkerOperation.returnLogOverlap();
            cout<<"weight: "<< exp(logOverlap) <<endl;
            mixedMeasure.addMeasurement( walkerWalkerOperation, exp(logOverlap) );
        }
    }
}

void AfqmcPhaseless::addMixedMeasurement()
{
    complex<double> logOverlap;
    double logDiff;
    WalkerRight walkerTemp;

    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walker[i], walkerTemp);
            walkerWalkerOperation.set(phiT, walkerTemp);

            logOverlap = walkerWalkerOperation.returnLogOverlap();

            //Cap weight real part
            logDiff = ( logOverlap - logOverlapBackup[i] ).real();
            if( (logDiff-method.logEnergyCap)>1e-12 )
            {
                cout<<"Cap weight in measurement: "
                    <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                    <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                logOverlap+= ( method.logEnergyCap - logDiff );
            }
            cout<<"weight: "<< exp(logOverlap)<<endl;
            mixedMeasure.addMeasurement( walkerWalkerOperation, exp(logOverlap) );
        }
    }
}

void AfqmcPhaseless::writeAndResetMeasurement()
{
    mixedMeasure.write();
    mixedMeasure.reSet();
}

void AfqmcPhaseless::adjustETThenResetMeasurement()
{
    method.ET = ( mixedMeasure.returnEnergy() ).real();

    if( MPIRank()==0 )
    {
        cout<<"\nAdjust trial energy: "<<method.ET<<endl;
    }

    mixedMeasure.reSet();
}

void AfqmcPhaseless::adjustETAndBackGroundThenResetMeasurement()
{
    method.ET = ( mixedMeasure.returnEnergy() ).real();
    model.updateBackGround( mixedMeasure.returnPhoSumBg(), mixedMeasure.returnPhoMinusBg());

    if( MPIRank()==0 )
    {
        cout<<"\nAdjust trial energy: "<<method.ET<<endl;
        cout<<"Adjust phoSum background: "<<model.getphoSumBg()<<"\n"<<endl;
        cout<<"Adjust phoMinus background: "<<model.getphoMinusBg()<<"\n"<<endl;
    }

    mixedMeasure.reSet();
}
