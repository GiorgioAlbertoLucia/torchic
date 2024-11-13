#include <Riostream.h>
#include <string>
#include <fstream>

#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <Math/GenVector/Boost.h>
#include <TSystem.h>
#include <TF1.h>
#include <TString.h>

void RunTwoBodyDecaySimulation(const char * MotherFileName,const char * MotherHistName,const char * inputTimeFileName,const char * inputTimeHistName,const char * OutputFileName,float motherMass,float firstDaugtherMass,float secondDaughterMass,float detectorLimit=0.,const int nEvents=1000000,int seed=0)
{
    /*
    *   This function is used to run the simulation of the two body decay of the mother particle
    *   Arguments:
    *   MotherFileName: The name of the file containing the mother particle pt distribution
    *   MotherHistName: The name of the histogram containing the mother particle pt distribution
    *   inputTimeFileName: The name of the file containing the decay time distribution
    *   inputTimeHistName: The name of the histogram containing the decay time distribution
    *   OutputFileName: The name of the file to save the output, it will contain the pt distribution of daughters and mother
    *   motherMass: The mass of the mother particle (GeV/c^2)
    *   firstDaugtherMass: The mass of the first daughter particle (GeV/c^2)
    *   secondDaughterMass: The mass of the second daughter particle (GeV/c^2)
    *   detectorLimit: The limit of the detector  (cm) (optional)
    *   nEvents: The number of events to generate
    *   seed: The seed for the random number generator
    */
    gSystem->AddLinkedLibs("-L$LD_LIBRARY_PATH -lGenVector");
    TFile *MotherFile = new TFile(MotherFileName,"READ");
    TFile *inputTimeFile = new TFile(inputTimeFileName,"READ");
    TFile *OutputFile = new TFile(OutputFileName,"RECREATE");

    TH1F *hTime = (TH1F*)inputTimeFile->Get(inputTimeHistName);
    TH1F *hMotherPt = (TH1F*)MotherFile->Get(MotherHistName);

    float maxMotherPt = hMotherPt->GetMaximum();

    TH1F *hPtFirstDaughter = new TH1F("hPtFirstDaughter","Pt of the first daughter",100,0,MaxMotherPt);
    TH1F *hPtSecondDaughter = new TH1F("hPtSecondDaughter","Pt of the second daughter",100,0,MaxMotherPt);
    TH1F *hPtMotherOutput = new TH1F("hPtMother","Pt of the mother",100,0,MaxMotherPt);

    gRandom->SetSeed(seed);

    for(int ievent=0;ievent<nEvents;i++)
    {
        float decayTime = hTime->GetRandom();
        float motherPt = hMotherPt->GetRandom();
        float motherEta=2*gRandom->Uniform()-1; //Uniform distribution between -1 and 1
        float motherPhi=2*TMath::Pi()*gRandom->Uniform(); //Uniform distribution between 0 and 2pi

        TLorentzVector motherTLV;
        motherTLV.SetPtEtaPhiM(motherPt,motherEta,motherPhi,motherMass);
        float beta = motherTLV.Beta();
        float betax=beta*cos(motherTLV.Phi())*sin(motherTLV.Theta());
        float betay=beta*sin(motherTLV.Phi())*sin(motherTLV.Theta());
        float betaz=beta*cos(motherTLV.Theta());

        ROOT::Math::Boost boost(betax,betay,betaz);

        float energyFirstDaughterCM = (motherMass*motherMass+firstDaugtherMass*firstDaugtherMass-secondDaughterMass*secondDaughterMass)/(2*motherMass); //Energy of the first daughter in the Center of Mass frame
        float energySecondDaughterCM = (motherMass*motherMass+secondDaugtherMass*secondDaughterMass-firstDaughterMass*firstDaughterMass)/(2*motherMass); //Energy of the second daughter in the Center of Mass frame
        float momentumFirstDaughterCM = TMath::Sqrt(energyFirstDaughterCM*energyFirstDaughterCM-firstDaugtherMass*firstDaugtherMass); //Momentum of the first daughter in the Center of Mass frame
        float momentumSecondDaughterCM = TMath::Sqrt(energySecondDaughterCM*energySecondDaughterCM-secondDaugtherMass*secondDaugtherMass); //Momentum of the second daughter in the Center of Mass frame

        float thetaCM = TMath::ACos(1-2*gRandom->Uniform()); //Uniform distribution between 0 and 2pi
        float phiCM = 2*TMath::Pi()*gRandom->Uniform(); //Uniform distribution between 0 and 2pi

        TLorentzVector firstDaughterTLV;
        firstDaughterTLV.SetPxPyPzE(momentumFirstDaughterCM*sin(thetaCM)*cos(phiCM),momentumFirstDaughterCM*sin(thetaCM)*sin(phiCM),momentumFirstDaughterCM*cos(thetaCM),energyFirstDaughterCM);
        TLorentzVector secondDaughterTLV;
        secondDaughterTLV.SetPxPyPzE(-momentumSecondDaughterCM*sin(thetaCM)*cos(phiCM),-momentumSecondDaughterCM*sin(thetaCM)*sin(phiCM),-momentumSecondDaughterCM*cos(thetaCM),energySecondDaughterCM); //The minus sign is because the second daughter is going in the opposite direction for the momentum conservation

        firstDaughterTLV = boost(firstDaughterTLV);
        secondDaughterTLV = boost(secondDaughterTLV);

        float firstDaughterBeta = firstDaughterTLV.Beta();
        float secondDaughterBeta = secondDaughterTLV.Beta();

        if(detectorLimit>0)
        {
            float decayLengthFirstDaughter = firstDaughterBeta*decayTime*TMath::C()*100; //cm
            float decayLengthSecondDaughter = secondDaughterBeta*decayTime*TMath::C()*100; //cm
            if(decayLengthFirstDaughter >= detectorLimit) hPtFirstDaughter->Fill(firstDaughterTLV.Pt());
            if(decayLengthSecondDaughter >= detectorLimit) hPtSecondDaughter->Fill(secondDaughterTLV.Pt());
        }
        else
        {
            hPtFirstDaughter->Fill(firstDaughterTLV.Pt());
            hPtSecondDaughter->Fill(secondDaughterTLV.Pt());
        }
        hPtMotherOutput->Fill(motherPt);
    }
    hPtFirstDaughter->Write();
    hPtSecondDaughter->Write();
    hPtMotherOutput->Write();
    OutputFile->Close();
}
