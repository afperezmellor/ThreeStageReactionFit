//==============================================================================================================================
//===================================          3-Stage multiple fit analysis     ===============================================
// This script do a simultaneous fit analysis of a  general three stage kinetic scheme. It has been developed to extract the
// kinetic parameters reported in the article entitled: Determination of kinetic properties in unimolecular dissociation of
// complex systems from graph-theory based analysis of an ensemble of reactive trajectories. 
// Please for further details please read the  article. (doi: )
//
//  The installation of ROOT can be done easily through the following steps: https://iscinumpy.gitlab.io/post/root-conda/
//
// Author: Ariel Francis Perez Mellor
// Last Update: 2021-05-22
//
// input files:   00prob-SP        =  Probability of the starting point state as a function of time. <time> <Prob> <ePprob> <#events>     
//                01prob-INT       =  Probability of the intermediate state as a function of time.
//                02prob-PF        =  Probability of the primary fragmentation state.
//    
//
//   parameters
//   p[0] = kis
//   p[1] = ksi       
//   p[2] = kpi     
//   p[3] = kip    
//   p[4] = ksp
//   p[5] = kps
//   p[6] = t0   (it is not included in the article, it is a time shift ) 
//
// ouput files: MULTIFIT-RESULTS_FIT    = display all the rate constants values as well as the statistical errors 
//              MULTIFIT-RESULTS_LINES  = display the population of the states obtained from the fit
//              MultipleFit.pdf         = display the figures resulting from the fit    
//
//==============================================================================================================================
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include <TStopwatch.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TVectorT.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TPave.h"

#include "ThreeStageFunctions.h"



using namespace std;

//  here we define the shared parameters
        
    int iparREACT[7] = {                        
        0,        //    kis
        1,        //    ksi
        2,        //    kpi 
        3,        //    kip  
        4,        //    ksp  
        5,        //    kps
        6,        //    t0 
        };    
        
    int iparINT[7] = {                        
        0,        //    kis
        1,        //    ksi
        2,        //    kpi 
        3,        //    kip  
        4,        //    ksp  
        5,        //    kps
        6,        //    t0 
        };    
        
    int iparNO_REACT[7] = {                        
        0,        //    kis
        1,        //    ksi
        2,        //    kpi 
        3,        //    kip  
        4,        //    ksp  
        5,        //    kps
        6,        //    t0 
        };                

//  here we construct the global chi structure                        
    struct GlobalChi2 { GlobalChi2(  ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2, ROOT::Math::IMultiGenFunction & f3) 
                                    : fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}
                
        double operator() (const double *par) const {
        
            double p1[7];
            for (int i = 0; i < 7; ++i) p1[i] = par[ iparREACT[i] ];
            
            double p2[7];
            for (int i = 0; i < 7; ++i) p2[i] = par[ iparINT[i] ];
            
            double p3[7];
            for (int i = 0; i < 7; ++i) p3[i] = par[ iparNO_REACT[i] ];
                    
            return (*fChi2_1)(p1) + (*fChi2_2)(p2)+(*fChi2_3)(p3) ;
        }

        const  ROOT::Math::IMultiGenFunction * fChi2_1;
        const  ROOT::Math::IMultiGenFunction * fChi2_2;
        const  ROOT::Math::IMultiGenFunction * fChi2_3;
        
                        
        };


void ThreeStageReactionFit()
{

    TStopwatch timer;
    timer.Start();    

// loading input data 
    char REACT[]="02prob-PF";
    double cTIME, PROB, ePROB, tmp;
    
    TVectorD PROB_REACT, ePROB_REACT, TIME, eTIME;
    ifstream fin_REACT(REACT);

    int j=0;
    while(fin_REACT >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_REACT.ResizeTo(j+1);
            ePROB_REACT.ResizeTo(j+1);
            TIME.ResizeTo(j+1);
            eTIME.ResizeTo(j+1);
                        
            PROB_REACT[j]=PROB;
            ePROB_REACT[j]=ePROB;
            TIME[j]=cTIME;
            eTIME[j]=0.0;
            
            j++;
    }
    fin_REACT.close();    fin_REACT.clear();    
    

    char INT[]="01prob-INT";
    TVectorD PROB_INT, ePROB_INT;
    ifstream fin_INT(INT);

    j=0;
    while(fin_INT >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_INT.ResizeTo(j+1);
            ePROB_INT.ResizeTo(j+1);
                        
            PROB_INT[j]=PROB;
            ePROB_INT[j]=ePROB;
        
            j++;
    }
    fin_INT.close();    fin_INT.clear();    

    char NO_REACT[]="00prob-SP";

    TVectorD PROB_NO_REACT, ePROB_NO_REACT, TIME_SP, eTIME_SP;
    ifstream fin_NO_REACT(NO_REACT);    
    
    j=0;
    while(fin_NO_REACT >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_NO_REACT.ResizeTo(j+1);
            ePROB_NO_REACT.ResizeTo(j+1);
            TIME_SP.ResizeTo(j+1);
            eTIME_SP.ResizeTo(j+1);
                        
            PROB_NO_REACT[j]=PROB;
            ePROB_NO_REACT[j]=ePROB;
            TIME_SP[j]=cTIME;
            eTIME_SP[j]=0.0;
        
            j++;
    }
    fin_NO_REACT.close();    fin_NO_REACT.clear();    
    

//  starting the fitting
    gROOT->SetStyle("Plain");
        
        
        if(gROOT->FindObject("DATA_PROB_REACT")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_REACT");
        }
        
        TGraphErrors *DATA_PROB_REACT=new TGraphErrors( TIME, PROB_REACT, eTIME, ePROB_REACT);

        DATA_PROB_REACT->SetLineColor(1);
        DATA_PROB_REACT->SetMarkerStyle(20);
        DATA_PROB_REACT->SetMarkerColor(1);
        DATA_PROB_REACT->SetMarkerSize(0.8);
//______________________________________________________________________________________________           
        if(gROOT->FindObject("DATA_PROB_INT")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_INT");
        }
        
        TGraphErrors *DATA_PROB_INT=new TGraphErrors( TIME, PROB_INT, eTIME, ePROB_INT);

        DATA_PROB_INT->SetLineColor(1);
        DATA_PROB_INT->SetMarkerStyle(20);
        DATA_PROB_INT->SetMarkerColor(1);
        DATA_PROB_INT->SetMarkerSize(0.8);
//______________________________________________________________________________________________
        
        if(gROOT->FindObject("DATA_PROB_NO_REACT")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_NO_REACT");
        }        
        
        TGraphErrors *DATA_PROB_NO_REACT=new TGraphErrors( TIME_SP, PROB_NO_REACT, eTIME_SP, ePROB_NO_REACT);

        DATA_PROB_NO_REACT->SetLineColor(1);
        DATA_PROB_NO_REACT->SetMarkerStyle(20);
        DATA_PROB_NO_REACT->SetMarkerColor(1);
        DATA_PROB_NO_REACT->SetMarkerSize(0.8);
//______________________________________________________________________________________________       
        if(gROOT->FindObject("FIT_REACT")!=0)
        {
           gROOT->ProcessLine("delete FIT_REACT");
        }
        TF1 *FIT_REACT=new TF1("FIT_REACT", ADJ_REACT_GEN, TIME.Min(), TIME.Max(), 7);
//______________________________________________________________________________________________             
        if(gROOT->FindObject("FIT_INT")!=0)
        {
           gROOT->ProcessLine("delete FIT_INT");
        }
        
        TF1 *FIT_INT=new TF1("FIT_INT", ADJ_INT_GEN, TIME.Min(), TIME.Max(), 7);
//______________________________________________________________________________________________
        if(gROOT->FindObject("FIT_NO_REACT")!=0)
        {
           gROOT->ProcessLine("delete FIT_NO_REACT");
        }
        
        TF1 *FIT_NO_REACT=new TF1("FIT_NO_REACT", ADJ_NO_REACT_GEN, TIME_SP.Min(), TIME_SP.Max(), 7);    
//______________________________________________________________________________________________        
        ROOT::Math::WrappedMultiTF1 wFIT_REACT(*FIT_REACT, 1);
        ROOT::Math::WrappedMultiTF1 wFIT_INT(*FIT_INT, 1);
        ROOT::Math::WrappedMultiTF1 wFIT_NO_REACT(*FIT_NO_REACT, 1);
    
        ROOT::Fit::DataOptions opt;

        ROOT::Fit::DataRange rangeFIT_REACT;
        rangeFIT_REACT.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_REACT(opt, rangeFIT_REACT);
        ROOT::Fit::FillData(dataFIT_REACT, DATA_PROB_REACT);
//______________________________________________________________________________________________      
        ROOT::Fit::DataRange rangeFIT_INT;
        rangeFIT_INT.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_INT(opt, rangeFIT_INT);
        ROOT::Fit::FillData(dataFIT_INT, DATA_PROB_INT);
//______________________________________________________________________________________________      
        ROOT::Fit::DataRange rangeFIT_NO_REACT;
        rangeFIT_NO_REACT.SetRange(TIME_SP.Min(), TIME_SP.Max());
        ROOT::Fit::BinData dataFIT_NO_REACT(opt, rangeFIT_NO_REACT);
        ROOT::Fit::FillData(dataFIT_NO_REACT, DATA_PROB_NO_REACT);
//______________________________________________________________________________________________        
        ROOT::Fit::Chi2Function chi2_REACT(dataFIT_REACT, wFIT_REACT);
        ROOT::Fit::Chi2Function chi2_INT(dataFIT_INT, wFIT_INT);
        ROOT::Fit::Chi2Function chi2_NO_REACT(dataFIT_NO_REACT, wFIT_NO_REACT);
        
        GlobalChi2 globalChi2(chi2_REACT, chi2_INT, chi2_NO_REACT);

        ROOT::Fit::Fitter fitter;
        
        const int Npar =7;
        double par0[Npar] = {    

        //   initializing the parameters   
        0.000266756,              //    kis       0
        2.1277e-05    ,           //    ksi       1
        0.000189036,              //    kpi       2
        0        ,                //    kip       3
        0,                        //    ksp       4
        0.00,                     //    kps       5
        0,                        //    t0        6
        
        };        

//  parameter settings 
    fitter.Config().SetParamsSettings(7, par0);
    

//    fitter.Config().ParSettings(0).Fix();
//    fitter.Config().ParSettings(1).Fix();
//    fitter.Config().ParSettings(2).Fix();
      fitter.Config().ParSettings(3).Fix();
      fitter.Config().ParSettings(4).Fix();
//    fitter.Config().ParSettings(5).Fix();    
      fitter.Config().ParSettings(6).Fix();
   
    
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit","Migrad");
    
// fiting FCN function directly (specify optionally data size and flag to indicate that is a chi2 fit) 
    fitter.FitFCN(7, globalChi2, 0, dataFIT_REACT.Size()+dataFIT_INT.Size()+dataFIT_NO_REACT.Size(), true);
    
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

// plotting    
    TCanvas *PROBABILITIES = new TCanvas("PROBABILITIES","PROBABILITIES", 0,0,900,900);
        
    PROBABILITIES->Divide(1,3);
    PROBABILITIES->cd(1);
        
        DATA_PROB_REACT->SetTitle("state PF");
        DATA_PROB_REACT->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_REACT->GetXaxis()->CenterTitle();
        DATA_PROB_REACT->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_REACT->GetYaxis()->CenterTitle();
        DATA_PROB_REACT->Draw("apz");    
        
        FIT_REACT->SetFitResult( result, iparREACT);
        FIT_REACT->SetRange(rangeFIT_REACT().first, rangeFIT_REACT().second);

        const int o_REACT=2000;
        double Y_REACT[o_REACT], X_REACT[o_REACT];
        for(int h_REACT=0; h_REACT < o_REACT; h_REACT++)
        {
            X_REACT[h_REACT]=TIME.Min()+h_REACT*(TIME.Max()-TIME.Min())/double(o_REACT-1);
            
            Y_REACT[h_REACT]=FIT_REACT->Eval(X_REACT[h_REACT]);
        
        }    
    
        TGraph *gr_REACT=new TGraph(o_REACT, X_REACT, Y_REACT);
        
        gr_REACT->SetLineColor(4);
        gr_REACT->SetLineWidth(3);
        gr_REACT->Draw("L");    
        
        
//______________________________________________________________________________________________       
        PROBABILITIES->cd(2);
        
        DATA_PROB_INT->SetTitle("state INT");
        DATA_PROB_INT->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_INT->GetXaxis()->CenterTitle();
        DATA_PROB_INT->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_INT->GetYaxis()->CenterTitle();
        DATA_PROB_INT->Draw("apz");    

        FIT_INT->SetFitResult( result, iparINT);
        FIT_INT->SetRange(rangeFIT_INT().first, rangeFIT_INT().second);

        const int o_INT=2000;
        double Y_INT[o_INT], X_INT[o_INT];
        for(int h_INT=0; h_INT<o_INT; h_INT++)
        {
            X_INT[h_INT]=TIME.Min()+h_INT*(TIME.Max()-TIME.Min())/double(o_INT-1);
            
            Y_INT[h_INT]=FIT_INT->Eval(X_INT[h_INT]);
        
        }    
    
        TGraph *gr_INT=new TGraph(o_INT, X_INT, Y_INT);
        
        gr_INT->SetLineColor(4);
        gr_INT->SetLineWidth(3);
        gr_INT->Draw("L");
//______________________________________________________________________________________________
        
        PROBABILITIES->cd(3);
        
        DATA_PROB_NO_REACT->SetTitle("state SP");
        DATA_PROB_NO_REACT->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_NO_REACT->GetXaxis()->CenterTitle();
        DATA_PROB_NO_REACT->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_NO_REACT->GetYaxis()->CenterTitle();
        DATA_PROB_NO_REACT->Draw("apz");    

        FIT_NO_REACT->SetFitResult( result, iparNO_REACT);
        FIT_NO_REACT->SetRange(rangeFIT_NO_REACT().first, rangeFIT_NO_REACT().second);

        const int o_NO_REACT=2000;
        double Y_NO_REACT[o_NO_REACT], X_NO_REACT[o_NO_REACT];
        for(int h_NO_REACT=0; h_NO_REACT<o_NO_REACT; h_NO_REACT++)
        {
            X_NO_REACT[h_NO_REACT]=TIME_SP.Min()+h_NO_REACT*(TIME_SP.Max()-TIME_SP.Min())/double(o_NO_REACT-1);
            
            Y_NO_REACT[h_NO_REACT]=FIT_NO_REACT->Eval(X_NO_REACT[h_NO_REACT]);
        
        }    
    
        TGraph *gr_NO_REACT=new TGraph(o_NO_REACT, X_NO_REACT, Y_NO_REACT);
        
        gr_NO_REACT->SetLineColor(4);
        gr_NO_REACT->SetLineWidth(3);
        gr_NO_REACT->Draw("L");
        
        
        const int o=2000;
        ofstream Results_lines("MULTIFIT-RESULTS_LINES_INT-REACT");
        
        Results_lines << " TIME (fs) " <<  "\t"  << " satate-INT "  <<  "\t"  << " state-PF " <<  endl;
        
        for(int h=0; h < o; h++)
        {
            
            Results_lines << X_REACT[h] << "\t"   << Y_INT[h]  <<  "\t"  << Y_REACT[h] <<  endl;
    
        }
        
        
        ofstream Results_lines_SP("MULTIFIT-RESULTS_LINES_NO-REACT");
        
        Results_lines_SP << " TIME (fs) " <<  "\t"  << " state-SP "  <<  endl;
        
        for(int h=0; h < o; h++)
        {
            
            Results_lines_SP << X_NO_REACT[h] << "\t"   << Y_NO_REACT[h]  <<  endl;
    
        }
        
        PROBABILITIES -> Print("MultipleFit.pdf");  

//  results and printing          
        ofstream Results("MULTIFIT-RESULTS_FIT");
        
        double a01=FIT_REACT->GetParameter(0);
        double ea01=FIT_REACT->GetParError(0);
        
        double b01=FIT_REACT->GetParameter(1);
        double eb01=FIT_REACT->GetParError(1);
        
        double c01=FIT_REACT->GetParameter(2);
        double ec01=FIT_REACT->GetParError(2);
        
        double e01=FIT_REACT->GetParameter(3);
        double ee01=FIT_REACT->GetParError(3);
        
        double f01=FIT_REACT->GetParameter(4);
        double ef01=FIT_REACT->GetParError(4);
        
        double g01=FIT_REACT->GetParameter(5);
        double eg01=FIT_REACT->GetParError(5);
        
        double tao=FIT_REACT->GetParameter(6);
        double etao=FIT_REACT->GetParError(6);
        
        double chi_red=FIT_REACT->GetChisquare()/FIT_REACT->GetNDF();
        
        double eChi_a01=ea01;
        double eChi_b01=eb01;
        double eChi_c01=ec01;
        double eChi_e01=ee01;
        double eChi_f01=ef01;
        double eChi_g01=eg01;
        double eChi_tao=etao;
        
        if(chi_red > 1)
            {
                eChi_a01=eChi_a01*sqrt(chi_red);
                eChi_b01=eChi_b01*sqrt(chi_red);
                eChi_c01=eChi_c01*sqrt(chi_red);
                eChi_e01=eChi_e01*sqrt(chi_red);
                eChi_f01=eChi_f01*sqrt(chi_red);
                eChi_g01=eChi_g01*sqrt(chi_red);
                eChi_tao=eChi_tao*sqrt(chi_red);
                        
            }
        
        
        Results.precision(6);
        
        Results << "\n chi_red  :\t " << chi_red  << "\t NDF  : \t"  << FIT_REACT->GetNDF() <<  endl;
        

        Results << "\n kis  :\t " << a01 << "\t +/- \t"  << ea01 <<  "\t eChi \t"  << eChi_a01 <<   endl;
        Results << "\n ksi  :\t " << b01 << "\t +/- \t"  << eb01 <<  "\t eChi \t"  << eChi_b01 <<   endl;
        Results << "\n kpi  :\t " << c01 << "\t +/- \t"  << ec01 <<  "\t eChi \t"  << eChi_c01 <<   endl;
        Results << "\n kip  :\t " << e01 << "\t +/- \t"  << ee01 <<  "\t eChi \t"  << eChi_e01 <<   endl;
        Results << "\n ksp  :\t " << f01 << "\t +/- \t"  << ef01 <<  "\t eChi \t"  << eChi_f01 <<   endl;
        Results << "\n kps  :\t " << g01 << "\t +/- \t"  << eg01 <<  "\t eChi \t"  << eChi_g01 <<   endl;
        Results << "\n t0   :\t " << tao << "\t +/- \t"  << etao <<  "\t eChi \t"  << eChi_tao <<   endl;

        Results << "\n \t\t\t\t -------------------- \t \t " << endl;    


    
        timer.Stop();
    
       double cputime = timer.CpuTime();
       int hours=cputime/3600;
       int minutes=(cputime-3600*hours)/60;
       double seconds=cputime-3600*hours-60*minutes;
    cout << "\n\nCPU time (h:m:s)\t\t" << hours << ":" << minutes << ":" << seconds << "\n\n" << endl;
    
}


