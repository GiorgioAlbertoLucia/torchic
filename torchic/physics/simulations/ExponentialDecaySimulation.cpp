#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <Riostream.h>

TH1F RunExponentialDecaySimulation(double tau, int nevents, int nbins, double totT, int seed)
{
    /*
        Function to generate a time distribution of a decay process

        Args:
            tau (double): The half-life of the decay
            nevents (int): The number of events
            nbins (int): The number of bins
            totT (double): The total time
            seed (int): The seed for the random number generator
        Returns:
            TH1F: The time distribution histogram
    */
    double alfa=1/tau;
    gRandom->SetSeed(seed);
    double delt= totT/nbins;
    TH1F *decayhist = new TH1F("decayhist","Decay",nbins+1,-delt,totT+delt/2);
    double prob=alfa*delt;
    decayhist->Fill(0.,nevents);
    for(double time=delt; time<totT+delt/2; time+=delt)
    {
        int ndec=0;
        for(int i=0;i<nevents;i++)
        {
            if(gRandom->Rndm()<prob)
            {
                ndec++;
            }
        }
        nevents-=ndec;
        decayhist->Fill(time,nevents);
    }
    return *decayhist;
}