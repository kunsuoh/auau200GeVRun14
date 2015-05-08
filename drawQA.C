#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"

const int nbin = 101;
double ptbin[nbin+1] = {0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.05, 1.1,  1.15,  1.2,  1.25,  1.3,  1.35,  1.4,  1.45,  1.5,  1.55,  1.6,  1.65,  1.7,  1.75,  1.8,  1.85,  1.9,  1.95, 2,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5};
double ptbin2[nbin+1] = {0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 1.0, 1.5, 2.5, 3.5, 5, 7, 10};
//double ptbin2[nbin+1] = {1.0, 1.5, 2.5, 3.5, 5, 7, 10};
int centbin[10] = {11,23,44,77,125,190,275,381,443,650}; // Guannan Xie's definition at HF PWG meeting Apr. 23
double zdcxbin[7] = {15000, 20000, 25000, 30000, 35000, 40000, 45000};

void drawQA(int selCent = 1,unsigned short opt = 0x2000){
    
    if (opt>>0 & 1) {
        TNtuple * ntuple = (TNtuple*)infile->Get("nt");
        TCanvas *c1 = new TCanvas("c1","c1",400,0,1000,680);
        c1->Divide(5,3);
        
        //drawing options
        TString vars[15] = {"pt", "nsige", "eoverp", "nphi", "neta",
            "phiDist", "zDist", "etaTowDist", "phiTowDist", "beta",
            "pairCharge", "pairMass", "pairDca", "positionX", "positionY"};
        TString xname[15] = {"pt", "nsige", "eoverp", "nphi", "neta",
            "phiDist", "zDist", "etaTowDist", "phiTowDist", "beta",
            "pairCharge", "pairMass", "pairDca", "positionX", "positionY"};
        int nbin[15]={100,289,100,10,10,
            100,100,100,100,100,
            3,100,100,100,100};
        float low[15]={0,-13,0,0,0,
            -0.2,-30,-2,-2,0.8,
            -3,0,0,-10,-10};
        float high[15]={12,13,3,10,10,
            0.2,30,2,2,1.2,
            3,1,3,10,10};
        
        TString cut_pt = "pt > 0.2";
        TString cut_mass = cut_pt + " && pairMass < .3";
        TString cut_dca =  cut_mass + " && pairDca < 0.2";
        TString cut_dca2=  cut_mass + " && pairDca < 0.1";
        TString cut_dca3=  cut_mass + " && pairDca < 0.01";
        TString cut_dca4=  cut_mass + " && pairDca < 0.001";
        TString cut_dca5=  cut_mass + " && pairDca < 0.0001";
        
        TString cut_tof = cut_mass + "&& abs(1-beta) < 0.2";
        TString cut_emc = cut_mass + "&& abs(phiTowDist) < 0.05 && abs(etaTowDist) < 0.05 && eoverp > 0.5 && eoverp < 2 && neta > 1 && nphi > 1 && abs(zDist) < 5 && abs(phiDist) < 0.01 ";
        
        TH1F * h[15];
        //draw histograms
        for (int i=0; i<15; i++) {
            h[i] = new TH1F(Form("h%d",i),Form("h%d",i),nbin[i],low[i],high[i]);
            h[i]->GetXaxis()->SetTitle(xname[i]);
            h[i]->SetMinimum(0.7);
            h[i]->SetFillColor(45);
            if (i==1) h[i]->GetXaxis()->SetRangeUser(-3,3);
            
            c1->cd(i+1);
            c1->cd(i+1)->SetLogy();
            
            ntuple->Draw(vars[i]+Form(">>h%d",i),cut_pt);
            
            c1->Update();
            cout << "drawing " << vars[i] << "-mass cut" << endl;
            ntuple->SetFillColor(38);
            ntuple->Draw(vars[i],cut_mass,"same");
            
            c1->Update();
            cout << "drawing " << vars[i] << "-dca cut" << endl;
            ntuple->SetFillColor(36);
            ntuple->Draw(vars[i],cut_dca,"same");
            
            c1->Update();
            cout << "drawing " << vars[i] << "-dca2 cut" << endl;
            ntuple->SetFillColor(34);
            ntuple->Draw(vars[i],cut_dca2,"same");
            
            c1->Update();
            cout << "drawing " << vars[i] << "-dca3 cut" << endl;
            ntuple->SetFillColor(30);
            ntuple->Draw(vars[i],cut_dca3,"same");
            
            c1->Update();
            cout << "drawing " << vars[i] << "-dca4 cut" << endl;
            ntuple->SetFillColor(24);
            ntuple->Draw(vars[i],cut_dca4,"same");
            
            c1->Update();
            cout << "drawing " << vars[i] << "-dca5 cut" << endl;
            ntuple->SetFillColor(12);
            ntuple->Draw(vars[i],cut_dca5,"same");
            
            
            //setting histograms
            c1->cd();
            c1->Update();
        }
    }
    if (opt>>1 & 1) {
        
        TTree * ntuple = (TNtuple*)infile->Get("t1");
        //drawing options
        TString vars[16] = {"pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY"};
        TString xname[16] = {"Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]"};
        TString yname[16] = {"Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]"};
        int xnbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float xlow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float xhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        int ynbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float ylow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float yhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        TString baseCut = "pairCharge==0 && ";//nsige < 2 && nsige > -1 && pairCharge==0 && abs(1-beta) < 0.025 &&";
        TString baseCut2 = "pairCharge!=0 && ";//nsige < 2 && nsige > -1 && pairCharge!=0 && abs(1-beta) < 0.025 &&";
        TString cut[16] = {"pairDca < 0.1", "pairDca < 0.01", "pairDca < 0.001", "pairDca < 0.0001",
            "pairDca < 0.1 && pt < 0.5", "pairDca < 0.01 && pt < 0.5", "pairDca < 0.001 && pt < 0.5", "pairDca < 0.0001 && pt < 0.5",
            "pairDca < 0.1 && pt > 0.5 && pt < 1 ", "pairDca < 0.01 && pt > 0.5 && pt < 1 ", "pairDca < 0.001 && pt > 0.5 && pt < 1 ", "pairDca < 0.0001 && pt > 0.5 && pt < 1 ",
            "pairDca < 0.1 && pt > 1", "pairDca < 0.01 && pt > 1", "pairDca < 0.001 && pt > 1", "pairDca < 0.0001 && pt > 1"
        };
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=16; i<101; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;pair mass;Entries [#]",i),200,0,0.5);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),200,0,0.5);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),200,0,0.5);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("pairMass>>h3%d",i),baseCut + cut[1] + Form("&& pt > %f && pt < %f", ptbin[i], ptbin[i+1]),"p");        // unlike sign
            ntuple->Draw(Form("pairMass>>h4%d",i),baseCut2 + cut[1] + Form("&& pt > %f && pt < %f", ptbin[i], ptbin[i+1]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/phe_%.2f_%.2f.pdf",ptbin[i],ptbin[i+1]));
            
            cout << i << " " << h5[i]->GetEntries() << " " << h5[i]->Integral()/h3[i]->Integral() << endl;
        }
        
        c3->cd();
        c3->SaveAs("~/Desktop/out_PtbyMassDca.pdf");
        
    }
    if (opt>>2 & 1) {
        TNtuple * ntuple = (TNtuple*)infile->Get("nt_pico");
        TCanvas *c4 = new TCanvas("c4","c4",200,0,800,800);
        c4->Divide(10,10);
        
        //drawing options
        const int nbin = 101;
        double ptbin[nbin+1] = {0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.05, 1.1,  1.15,  1.2,  1.25,  1.3,  1.35,  1.4,  1.45,  1.5,  1.55,  1.6,  1.65,  1.7,  1.75,  1.8,  1.85,  1.9,  1.95, 2,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5};
        TString vars[nbin] = {"nsige:(1-1/beta*TMath::Sqrt(1-0.000511*0.000511/pt/pt/TMath::CosH(eta)/TMath::CosH(eta)))"};
        TString xname[nbin] = {"dbeta"};
        TString yname[nbin] = {"nsige"};
        int ynbin[nbin]={289};
        float ylow[nbin]={-13};
        float yhigh[nbin]={13};
        int xnbin[nbin]={800};
        float xlow[nbin]={-0.2};
        float xhigh[nbin]={0.6};
        TString baseCut = "";
        TString cut[nbin] = {""};
        
        TH2F * h6[nbin];
        TPaveText * pt3[nbin];
        TFile * fout = new TFile("out.root","recreate");
        fout->cd();
        //draw histograms
        for (int j=0; j<5; j++) {
            for (int i=0; i<nbin; i++) {
                if (i == 100) continue;
                h6[i] = new TH2F(Form("h6%d",i),Form("h6%d%d",j,i),xnbin[0],xlow[0],xhigh[0],ynbin[0],ylow[0],yhigh[0]);
                h6[i]->GetXaxis()->SetTitle(xname[0]);
                h6[i]->GetYaxis()->SetTitle(yname[0]);
                //   h6[i]->GetYaxis()->SetRangeUser(-3,3);
                c4->cd(i+1);
                c4->cd(i+1)->SetLogz();
                
                TString cutForElectron = " && ((1-1/beta*TMath::Sqrt(1-0.000511*0.000511/pt/pt/TMath::CosH(eta)/TMath::CosH(eta)))**2/0.02**2 + nsige**2/9 < 1)";
                ntuple->Draw(vars[0]+Form(">>h6%d",i),Form("pt > %f && pt < %f && abs(eta) > %f  && abs(eta) < %f ", ptbin[i], ptbin[i+1], (float)j*0.1, (float)j*0.1 + 0.1 ) ,"col2");
                pt3[i] = new TPaveText(-0.09,-10,0.09,-8);
                pt3[i]->AddText(Form("%.2f < pT < %.2f",ptbin[i], ptbin[i+1] ));
                pt3[i]->SetFillColor(0);
                pt3[i]->Draw("same");
                
                c4->Update();
                c4->cd();
                h6[i]->Write();
                cout << i << endl;
            }
        }
        fout->Close();
    }
    if (opt>>3 & 1) {
        
        TTree * ntuple = (TTree*)infile->Get("t1");
        TCanvas *c2 = new TCanvas("c2","c2",200,0,1000,1000);
        c2->Divide(4,4);
        
        //drawing options
        TString vars[16] = {"pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY"};
        TString xname[16] = {"Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]"};
        TString yname[16] = {"Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]"};
        int xnbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float xlow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float xhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        int ynbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float ylow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float yhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        TString baseCut = "nsige < 2 && nsige > -1 && pairCharge==0 && pairAngle3d > 0.015 && pt > 0.5 && ";
        TString cut[16] = {
            "pairMass < 0.3  && pairDca < 0.1 && pairDca > 0.01", "pairMass < 0.3  && pairDca < 0.01 && pairDca > 0.001", "pairMass < 0.3  && pairDca < 0.001 && pairDca > 0.0001", "pairMass < 0.3  && pairDca < 0.0001",
            "pairMass < 0.1  && pairDca < 0.1 && pairDca > 0.01", "pairMass < 0.1  && pairDca < 0.01 && pairDca > 0.001", "pairMass < 0.1  && pairDca < 0.001 && pairDca > 0.0001", "pairMass < 0.1  && pairDca < 0.0001",
            "pairMass < 0.05 && pairDca < 0.1 && pairDca > 0.01", "pairMass < 0.05 && pairDca < 0.01 && pairDca > 0.001", "pairMass < 0.05 && pairDca < 0.001 && pairDca > 0.0001", "pairMass < 0.05 && pairDca < 0.0001",
            "pairMass < 0.01 && pairDca < 0.1 && pairDca > 0.01", "pairMass < 0.01 && pairDca < 0.01 && pairDca > 0.001", "pairMass < 0.01 && pairDca < 0.001 && pairDca > 0.0001", "pairMass < 0.01 && pairDca < 0.0001"
        };
        
        TH2F * h2[16];
        TH1F * h3[16];
        TH1F * h4[16];
        TH1F * h5[16];
        TPaveText * pt[16];
        TPaveText * pt2[16];
        
        //draw histograms
        for (int i=0; i<16; i++) {
            h2[i] = new TH2F(Form("h2%d",i),Form("h2%d",i),xnbin[i],xlow[i],xhigh[i],ynbin[i],ylow[i],yhigh[i]);
            h2[i]->GetXaxis()->SetTitle(xname[i]);
            h2[i]->GetYaxis()->SetTitle(yname[i]);
            
            c2->cd(i+1);
            c2->cd(i+1)->SetLogz();
            
            ntuple->Draw(vars[i]+Form(">>h2%d",i),baseCut + cut[i],"col2");
            pt[i] = new TPaveText(-10,-13,1,-11);
            pt[i]->AddText(Form("# of pairs: %.0f",h2[i]->GetEntries()));
            pt[i]->SetFillColor(0);
            pt[i]->Draw("same");
            
            c2->Update();
            
            
            cout << i << " " << h2[i]->GetEntries() << endl;
        }
        
        c2->cd();
        
        c2->SaveAs("~/Desktop/out_PositionXYbyMassDca.pdf");
        
    }
    if (opt>>4 & 1) {
        TString infilename = "Ana_12X_8.root";
        TFile * infile = new TFile(infilename);

        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString vars[16] = {"pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY",
            "pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY","pairPositionX:pairPositionY"};
        TString xname[16] = {"Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]",
            "Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]","Conversion X [cm]"};
        TString yname[16] = {"Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]",
            "Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]","Conversion Y [cm]"};
        int xnbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float xlow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float xhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        int ynbin[16]={1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,
            1000,1000,1000,1000,1000,1000};
        float ylow[16]={-10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,
            -10,-10,-10,-10,-10,-10};
        float yhigh[16]={10,10,10,10,10,
            10,10,10,10,10,
            10,10,10,10,10,10};
        TString baseCut = "pairCharge==0 && ";//nsige < 2 && nsige > -1 && pairCharge==0 && abs(1-beta) < 0.025 &&";
        TString baseCut2 = "pairCharge!=0 && ";//nsige < 2 && nsige > -1 && pairCharge!=0 && abs(1-beta) < 0.025 &&";
        TString cut[16] = {"pairDca < 0.1", "pairDca < 0.01", "pairDca < 0.001", "pairDca < 0.0001",
            "pairDca < 0.1 && pt < 0.5", "pairDca < 0.01 && pt < 0.5", "pairDca < 0.001 && pt < 0.5", "pairDca < 0.0001 && pt < 0.5",
            "pairDca < 0.1 && pt > 0.5 && pt < 1 ", "pairDca < 0.01 && pt > 0.5 && pt < 1 ", "pairDca < 0.001 && pt > 0.5 && pt < 1 ", "pairDca < 0.0001 && pt > 0.5 && pt < 1 ",
            "pairDca < 0.1 && pt > 1", "pairDca < 0.01 && pt > 1", "pairDca < 0.001 && pt > 1", "pairDca < 0.0001 && pt > 1"
        };
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=0; i<20; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;pair mass;Entries [#]",i),40,0,0.2);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),40,0,0.2);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),40,0,0.2);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("pairMass>>h3%d",i),baseCut +  Form(" pairDca < %f", dcacut[i]),"p");        // unlike sign
            ntuple->Draw(Form("pairMass>>h4%d",i),baseCut2 + Form(" pairDca < %f", dcacut[i]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/MassbyDca/phe_dca%f.pdf",dcacut[i]));
            
            cout << dcacut[i] << " " << h5[i]->Integral(1,10) << " " << h5[i]->Integral(1,10)/h3[i]->Integral(1,10) << endl;
        }
    }
    
    // partner electorn nSigE distribution by pairDca cut : 32
    if (opt>>5 & 1) {
        TString infilename = "Ana_12X_8.root";
        TFile * infile = new TFile(infilename);
        
        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString baseCut = "pairCharge==0 ";
        TString baseCut2 = "pairCharge!=0 ";
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        
        int nbin = 289;
        float xmin = -13;
        float xmax = 13;
        
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=0; i<20; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;Partner electron nSigE;Entries [#]",i),nbin,xmin,xmax);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),nbin,xmin,xmax);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),nbin,xmin,xmax);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMaximum(2e5);
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("partner_nsige>>h3%d",i),baseCut +  Form("&& pairDca < %f", dcacut[i]),"p");        // unlike sign
            ntuple->Draw(Form("partner_nsige>>h4%d",i),baseCut2 + Form("&& pairDca < %f", dcacut[i]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            TF1 * funFit = new TF1("funFit","gaus(0)",-3,3);
            funFit->SetParameters( 4.11648e+03,
                                  -2.57344e-01,
                                  9.30512e-01);
            h5[i]->Fit(funFit,"NOR");
            funFit->SetLineColor(2);
            funFit->SetLineStyle(2);
            funFit->SetRange(-13,13);
            funFit->Draw("same");
            
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/PartnerNSigE/pairDca_%f.pdf",dcacut[i]));
            
        }
    }
    
    // draw mass distribution with w/ and w/o pairAngle3d cut by pT : 64
    if (opt>>6 & 1) {
        TString infilename = "Ana_12X_8.root";
        TFile * infile = new TFile(infilename);
        
        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString baseCut = "pairCharge==0 ";
        TString baseCut2 = "pairCharge!=0 ";
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        
        int nbin = 100;
        float xmin = 0;
        float xmax = 0.5;
        
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TH1D * h3_2[101];
        TH1D * h4_2[101];
        TH1D * h5_2[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=0; i<11; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;pair mass;Entries [#]",i),nbin,xmin,xmax);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),nbin,xmin,xmax);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),nbin,xmin,xmax);
            h3_2[i] = new TH1D(Form("h3_2%d",i),Form("h3_2%d;Partner electron nSigE;Entries [#]",i),nbin,xmin,xmax);
            h4_2[i] = new TH1D(Form("h4_2%d",i),Form("h4_2%d",i),nbin,xmin,xmax);
            h5_2[i] = new TH1D(Form("h5_2%d",i),Form("h5_2%d",i),nbin,xmin,xmax);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(24);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            h3_2[i]->SetMarkerStyle(20);
            h3_2[i]->SetMarkerColor(4);
            h3_2[i]->SetMarkerSize(0.5);
            h4_2[i]->SetMarkerStyle(20);
            h4_2[i]->SetMarkerColor(2);
            h4_2[i]->SetMarkerSize(0.5);
            h5_2[i]->SetMarkerStyle(20);
            h5_2[i]->SetMarkerColor(1);
            h5_2[i]->SetMarkerSize(0.5);
            h3_2[i]->Sumw2();
            h4_2[i]->Sumw2();
            h5_2[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
         //   h3[i]->SetMaximum(2e5);
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("pairMass>>h3%d",i),baseCut +  Form("&& pt > %f&& pt < %f", ptbin2[i], ptbin2[i+1]),"p");        // unlike sign
            ntuple->Draw(Form("pairMass>>h4%d",i),baseCut2 + Form("&& pt > %f&& pt < %f", ptbin2[i], ptbin2[i+1]),"psame");   // like sign background
            
            ntuple->Draw(Form("pairMass>>h3_2%d",i),baseCut +  Form("&& pairAngle3d < 0.2 && pt > %f&& pt < %f", ptbin2[i], ptbin2[i+1]),"psame");        // unlike sign
            ntuple->Draw(Form("pairMass>>h4_2%d",i),baseCut2 + Form("&& pairAngle3d < 0.2 && pt > %f&& pt < %f", ptbin2[i], ptbin2[i+1]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            h5_2[i]->Add(h4_2[i],h3_2[i],-1.,1.);
            h5_2[i]->Draw("psame");
            
            cout << h5[i]->Integral() << " " << h5_2[i]->Integral() <<  h5_2[i]->Integral()/h5[i]->Integral() << endl;
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/MassByPtAngleCut/Pt_%.2f_%.2f.pdf",ptbin2[i],ptbin2[i+1]));
            
        }
    
    }
    
    // partner electorn nSigE distribution by partner pT : 128
    if (opt>>7 & 1) {
        TString infilename = "Ana_12X_8.root";
        TFile * infile = new TFile(infilename);
        
        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString baseCut = "pairCharge==0 ";
        TString baseCut2 = "pairCharge!=0 ";
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        
        int nbin = 289;
        float xmin = -13;
        float xmax = 13;
        
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=0; i<11; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;Partner electron pT;Entries [#]",i),nbin,xmin,xmax);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),nbin,xmin,xmax);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),nbin,xmin,xmax);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMaximum(2e5);
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("partner_nsige>>h3%d",i),baseCut +  Form("&& partner_pt > %f && partner_pt < %f", ptbin2[i], ptbin2[i+1]),"p");        // unlike sign
            ntuple->Draw(Form("partner_nsige>>h4%d",i),baseCut2 + Form("&& partner_pt > %f && partner_pt < %f", ptbin2[i], ptbin2[i+1]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            TF1 * funFit = new TF1("funFit","gaus(0)",-3,3);
            funFit->SetParameters( 4.11648e+03,
                                  -2.57344e-01,
                                  9.30512e-01);
            h5[i]->Fit(funFit,"NOR");
            funFit->SetLineColor(2);
            funFit->SetLineStyle(2);
            funFit->SetRange(-13,13);
            funFit->Draw("same");
            
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/PartnerNSigE/partner_pt_%f_%f.pdf",ptbin2[i],ptbin2[i+1]));
            
        }
    }
    
    // partner electorn nSigE distribution by tagged pT : 256
    if (opt>>8 & 1) {
        TString infilename = "Ana_12X_8.root";
        TFile * infile = new TFile(infilename);
        
        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString baseCut = "pairCharge==0 ";
        TString baseCut2 = "pairCharge!=0 ";
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        
        int nbin = 289;
        float xmin = -13;
        float xmax = 13;
        
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        TCanvas *c3[101];
        
        //draw histograms
        for (int i=0; i<11; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;Tagged electron pT;Entries [#]",i),nbin,xmin,xmax);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),nbin,xmin,xmax);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),nbin,xmin,xmax);
            
            h3[i]->SetMarkerStyle(24);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(0.5);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMaximum(2e5);
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("partner_nsige>>h3%d",i),baseCut +  Form("&& pt > %f && pt < %f", ptbin2[i], ptbin2[i+1]),"p");        // unlike sign
            ntuple->Draw(Form("partner_nsige>>h4%d",i),baseCut2 + Form("&& pt > %f && pt < %f", ptbin2[i], ptbin2[i+1]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("psame");
            TF1 * funFit = new TF1("funFit","gaus(0)",-3,3);
            funFit->SetParameters( 4.11648e+03,
                                  -2.57344e-01,
                                  9.30512e-01);
            h5[i]->Fit(funFit,"NOR");
            funFit->SetLineColor(2);
            funFit->SetLineStyle(2);
            funFit->SetRange(-13,13);
            funFit->Draw("same");
            
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/PartnerNSigE/pt_%f_%f.pdf",ptbin2[i],ptbin2[i+1]));
            
        }
    }
    
    // pairMass distribution by partner nSigE cut : 512
    if (opt>>9 & 1) {
        TString infilename = "Ana_12X_12.root";
        TFile * infile = new TFile(infilename);
        
        TTree * ntuple = (TNtuple*)infile->Get("tPhE");
        //drawing options
        TString baseCut = "abs(1-beta) < 0.025 && pairDca < 1";
        TString chargeCut = "&& pairCharge==0";
        TString chargeCut2 = "&& pairCharge!=0";
        float dcacut[20] = {1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3,0.2 ,0.1, 0.05, 0.01, 0.005, 0.001, 0.0001};
        float partnerNSigECut[20] = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 13};
        
        const int nbin = 19;
        const int nMassBin = 400;
        float xmin = 0.;
        float xmax = 0.4;
        
        TH2F * h2[101];
        TH1D * h3[101];
        TH1D * h4[101];
        TH1D * h5[101];
        
        TCanvas *c3[nbin];
        float mass[nbin];

        //draw histograms
        for (int i=0; i<nbin; i++) {
            c3[i] = new TCanvas(Form("c3%d",i),Form("c3%d",i),0,0,500,400);
            
            h3[i] = new TH1D(Form("h3%d",i),Form("h3%d;electron pair mass;Entries [#]",i),nMassBin,xmin,xmax);
            h4[i] = new TH1D(Form("h4%d",i),Form("h4%d",i),nMassBin,xmin,xmax);
            h5[i] = new TH1D(Form("h5%d",i),Form("h5%d",i),nMassBin,xmin,xmax);
            
            h3[i]->SetMarkerStyle(20);
            h3[i]->SetMarkerColor(4);
            h3[i]->SetMarkerSize(0.5);
            h4[i]->SetMarkerStyle(24);
            h4[i]->SetMarkerColor(2);
            h4[i]->SetMarkerSize(0.5);
            h5[i]->SetMarkerStyle(20);
            h5[i]->SetMarkerColor(1);
            h5[i]->SetMarkerSize(1);
            h3[i]->Sumw2();
            h4[i]->Sumw2();
            h5[i]->Sumw2();
            
            c3[i]->cd();
            c3[i]->SetLogy();
            
            h3[i]->SetMaximum(1e4);
            h3[i]->SetMinimum(0.7);
            
            ntuple->Draw(Form("pairMass>>h3%d",i),baseCut + chargeCut +  Form("&& abs(partner_nsige) < %f", partnerNSigECut[i]),"p");        // unlike sign
            ntuple->Draw(Form("pairMass>>h4%d",i),baseCut + chargeCut2 + Form("&& abs(partner_nsige) < %f", partnerNSigECut[i]),"psame");   // like sign background
            
            h5[i]->Add(h4[i],h3[i],-1.,1.);
            h5[i]->Draw("same");
            c3[i]->Update();
            c3[i]->SaveAs(Form("~/Desktop/pairMassByPartnerNSigECut_%.1f.pdf",partnerNSigECut[i]));
            
            cout << h5[i]->Integral(1,200) << endl;
            mass[i] = h3[i]->Integral(1,10)-h4[i]->Integral(1,10)*h3[i]->Integral(201,400)/h4[i]->Integral(201,400); // 30 MeV/c
        }
        TCanvas * cc=  new TCanvas("cc","cc",500,400);
        TGraph * gra = new TGraph(nbin, partnerNSigECut, mass);
        gra->SetMarkerSize(1);
        gra->SetMarkerColor(1);
        gra->SetMarkerStyle(20);
        gra->Draw("");
        cc->SaveAs("~/Desktop/ElectronSignalYieldBynSigECut.pdf");
    }
    
    // hadrons dca distribution by pT
    if (opt>>10 & 1) {
        cout << "hadrons dca distribution by pT. " << endl; 
        TString infilename = "Ana_7.root";
//        infilename="out_166_15166046.root";
        TFile * infile = new TFile(infilename);
        TTree * tPion = (TTree*)infile->Get("tIncPion");
        TTree * tE = (TTree*)infile->Get("tInc");
        
        const int nbin      = 13;
        const int xnbin     = 98;      // pt
        const float xmin    = 0.2;
        const float xmax    = 10;
        const int ynbin     = 101;      // dca
        const float ymin    = -0.1;
        const float ymax    =  0.1;
        
        TCanvas *cPion[nbin];
        TCanvas *cE[nbin];
        TH2F * h2dPion[nbin];
        TH2F * h2dE[nbin];
        TH1F * h1dPion[nbin];
        TH1F * h1dE[nbin];

        
        double mean_pion[nbin];
        double mean_pionerr[nbin];
        double rms_pion[nbin];
        double rms_pionerr[nbin];
        double mean_electron[nbin];
        double mean_electronerr[nbin];
        double rms_electron[nbin];
        double rms_electronerr[nbin];
        
        double pt2[nbin];
        double pt2err[nbin];
        TString baseCut = "pt > 0.2";
        
        
        
        
        for (int i=0; i<nbin; i++) {
            TString cutE = baseCut + " && nsige > 0 && ((pt < 1. && abs(1-beta) < 0.025) || pt > 1.)";
            TString cutPion = baseCut;

            cPion[i] = new TCanvas(Form("cPion%d",i),Form("cPion%d",i),0,0,500,400);
            if (i==0) cE[i] = new TCanvas(Form("cE%d",i),Form("cE%d",i),500,0,500,400);
            cPion[i]->SetLogz();
            if (i==0) cE[i]->SetLogz();

            if (i==0) {
                h2dPion[i] = new TH2F(Form("h2dPion%d",i),Form("h2dPion%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                h2dE[i] = new TH2F(Form("h2dE%d",i),Form("h2dE%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                
                cPion[i]->cd();
                tPion->Draw(Form("dca:pt>>h2dPion%d",i),cutPion,"col2");
                
                cE[i]->cd();
                tE->Draw(Form("dca:pt>>h2dE%d",i),cutE,"col2");
            }
            
            else {
                TString ptcut = Form("&& pt > %f && pt < %f",ptbin2[i-1],ptbin2[i]);
                
                cPion[i]->SetLogy();

                h1dPion[i] = new TH1F(Form("h1dPion%d",i),Form("h1dPion%d;DCA_{XY} [cm];Entries [#]",i),ynbin,ymin,ymax);
                h1dE[i] = new TH1F(Form("h1dE%d",i),Form("h1dE%d",i),ynbin,ymin,ymax);

                cPion[i]->cd();
                
                TString centcut = Form("&& mRefMult > %d && mRefMult <= %d", centbin[selCent-1], centbin[9]);
                tPion->Draw(Form("dca>>h1dPion%d",i),cutPion + ptcut + centcut,"");
                tE->Draw(Form("dca>>h1dE%d",i),cutE + ptcut,"");
                h1dPion[i]->Sumw2();
                h1dE[i]->Sumw2();
                
                h1dPion[i]->Scale(1./h1dPion[i]->GetMaximum());
                h1dE[i]->Scale(1./h1dE[i]->GetMaximum());
                
                h1dE[i]->SetMarkerStyle(25);
                h1dE[i]->SetMarkerColor(1);
                h1dE[i]->SetMarkerSize(0.3);
                
                h1dPion[i]->SetMarkerStyle(24);
                h1dPion[i]->SetMarkerColor(2);
                h1dPion[i]->SetMarkerSize(0.3);
                
                h1dE[i]->SetLineColor(1);
                h1dPion[i]->SetLineColor(2);
                
                h1dPion[i]->SetMaximum(3);
                h1dPion[i]->SetMinimum(1e-5);

                h1dPion[i]->Draw("");
                h1dE[i]->Draw("same");
                h1dPion[i]->Draw("same");
                
                
                
                pt3 = new TPaveText(-0.05,1e-2,0.05,1e-4);
                pt3->AddText(Form("%.2f < pT < %.2f",ptbin2[i-1], ptbin2[i] ));
                pt3->SetFillStyle(0);
                pt3->SetLineColorAlpha(kBlue, 0);
                pt3->Draw("same");

                mean_pion[i-1] =      h1dPion[i]->GetMean();
                mean_pionerr[i-1] =   h1dPion[i]->GetMeanError();
                rms_pion[i-1] =       h1dPion[i]->GetRMS();
                rms_pionerr[i-1] =    h1dPion[i]->GetRMSError();

                mean_electron[i-1] =  h1dE[i]->GetMean();
                mean_electronerr[i-1] = h1dE[i]->GetMeanError();
                rms_electron[i-1] = h1dE[i]->GetRMS();
                rms_electronerr[i-1] = h1dE[i]->GetRMSError();

                pt2[i-1] = (ptbin2[i]+ptbin2[i-1])/2;
                pt2err[i-1] = ptbin2[i]-pt2[i-1];
                
                
            }
            if (i==0) {cE[i]->Update();cE[i]->SaveAs(Form("~/Desktop/electronDcaByPt_Cent%d.pdf",selCent));}
            cPion[i]->Update();
            if (i==0) cPion[i]->SaveAs(Form("~/Desktop/pionDcaByPt_Cent%d.pdf",selCent));
            else cPion[i]->SaveAs(Form("~/Desktop/pionAndElectronDca_Pt%.1fPt%.1f_Cent%d.pdf",ptbin2[i-1],ptbin2[i],selCent));
        }
        
        TCanvas * cc=  new TCanvas("cc","cc",500,400);
        TH1F * dum = new TH1F("dum","dum;p_{T} [GeV/c];Meam [cm]",100,0.2,10);
        dum->SetMaximum(0.002);
        dum->SetMinimum(-0.002);
        dum->SetLineStyle(2);
        dum->SetLineColor(2);
        dum->Draw();
        TGraphErrors * gra2 = new TGraphErrors(nbin-1, pt2, mean_pion,pt2err,mean_pionerr);
        gra2->SetMarkerSize(1);
        gra2->SetMarkerColor(2);
        gra2->SetMarkerStyle(20);
        gra2->Draw("plsame");
        TGraphErrors * gra = new TGraphErrors(nbin-1, pt2, mean_electron,pt2err,mean_electronerr);
        gra->SetMarkerSize(1);
        gra->SetMarkerColor(1);
        gra->SetMarkerStyle(20);
        gra->Draw("plsame");
        cc->SaveAs(Form("~/Desktop/mean_Cent%d.pdf",selCent));

        
        
        TCanvas * cc2=  new TCanvas("cc2","cc2",500,400);
        TH1F * dum2 = new TH1F("dum2","dum2;p_{T} [GeV/c];RMS [cm]",100,0.2,10);
        dum2->SetMaximum(0.04);
        dum2->SetMinimum(0);
        dum2->Draw();
        TGraphErrors * gra3 = new TGraphErrors(nbin-1, pt2, rms_pion,pt2err,rms_pionerr);
        gra3->SetMarkerSize(1);
        gra3->SetMarkerColor(2);
        gra3->SetMarkerStyle(20);
        gra3->Draw("plsame");
        TGraphErrors * gra4 = new TGraphErrors(nbin-1, pt2, rms_electron,pt2err,rms_electronerr);
        gra4->SetMarkerSize(1);
        gra4->SetMarkerColor(1);
        gra4->SetMarkerStyle(20);
        gra4->Draw("plsame");
        cc2->SaveAs(Form("~/Desktop/rms_Cent%d.pdf",selCent));

    }

    // hadrons dca distribution by Centrality
    if (opt>>11 & 1) {
        TString infilename = "Ana_6.root";
        //        infilename="out_166_15166046.root";
        TFile * infile = new TFile(infilename);
        TTree * tPion = (TTree*)infile->Get("tIncPion");
        TTree * tE = (TTree*)infile->Get("tInc");
        
        const int nbin      = 8;
        const int xnbin     = 98;      // pt
        const float xmin    = 0.2;
        const float xmax    = 10;
        const int ynbin     = 101;      // dca
        const float ymin    = -0.1;
        const float ymax    =  0.1;
        
        TCanvas *cPion[nbin];
        TCanvas *cE[nbin];
        TH2F * h2dPion[nbin];
        TH2F * h2dE[nbin];
        TH1F * h1dPion[nbin];
        TH1F * h1dE[nbin];
        
        
        double mean_pion[nbin];
        double mean_pionerr[nbin];
        double rms_pion[nbin];
        double rms_pionerr[nbin];
        double mean_electron[nbin];
        double mean_electronerr[nbin];
        double rms_electron[nbin];
        double rms_electronerr[nbin];
        
        double pt2[nbin];
        double pt2err[nbin];
        
        
        TString baseCut = "pt < 1";
        
        
        
        
        for (int i=0; i<nbin; i++) {
            TString cutE = baseCut + " && abs(nsige) < 2";
            TString cutPion = baseCut;
            
            cPion[i] = new TCanvas(Form("cPion%d",i),Form("cPion%d",i),0,0,500,400);
            if (i==0) cE[i] = new TCanvas(Form("cE%d",i),Form("cE%d",i),500,0,500,400);
            cPion[i]->SetLogz();
            if (i==0) cE[i]->SetLogz();
            
            if (i==0) {
                h2dPion[i] = new TH2F(Form("h2dPion%d",i),Form("h2dPion%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                h2dE[i] = new TH2F(Form("h2dE%d",i),Form("h2dE%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                
                cPion[i]->cd();
                tPion->Draw(Form("dca:pt>>h2dPion%d",i),cutPion,"col2");
                
                cE[i]->cd();
                tE->Draw(Form("dca:pt>>h2dE%d",i),cutE,"col2");
            }
            
            else {
                cPion[i]->SetLogy();
                
                h1dPion[i] = new TH1F(Form("h1dPion%d",i),Form("h1dPion%d;DCA_{XY} [cm];Entries [#]",i),ynbin,ymin,ymax);
                h1dE[i] = new TH1F(Form("h1dE%d",i),Form("h1dE%d",i),ynbin,ymin,ymax);
                
                cPion[i]->cd();
                
                TString centcut = Form("&& mRefMult > %d && mRefMult <= %d", centbin[i-1], centbin[i]);
                tPion->Draw(Form("dca>>h1dPion%d",i),cutPion + centcut,"");
                tE->Draw(Form("dca>>h1dE%d",i),cutE + centcut,"");
                h1dPion[i]->Sumw2();
                h1dE[i]->Sumw2();
                
                h1dPion[i]->Scale(1./h1dPion[i]->GetMaximum());
                h1dE[i]->Scale(1./h1dE[i]->GetMaximum());
                
                h1dE[i]->SetMarkerStyle(25);
                h1dE[i]->SetMarkerColor(1);
                h1dE[i]->SetMarkerSize(0.3);
                
                h1dPion[i]->SetMarkerStyle(24);
                h1dPion[i]->SetMarkerColor(2);
                h1dPion[i]->SetMarkerSize(0.3);
                
                h1dE[i]->SetLineColor(1);
                h1dPion[i]->SetLineColor(2);
                
                h1dPion[i]->SetMaximum(3);
                h1dPion[i]->SetMinimum(1e-5);
                
                h1dPion[i]->Draw("");
                h1dE[i]->Draw("same");
                h1dPion[i]->Draw("same");
                
                
                
                pt3 = new TPaveText(-0.05,1e-2,0.05,1e-4);
                pt3->AddText(Form("%d < mRefMult #leq %d",centbin[i-1], centbin[i] ));
                pt3->SetFillStyle(0);
                pt3->SetLineColorAlpha(kBlue, 0);
                pt3->Draw("same");
                
                mean_pion[i-1] =      h1dPion[i]->GetMean();
                mean_pionerr[i-1] =   h1dPion[i]->GetMeanError();
                rms_pion[i-1] =       h1dPion[i]->GetRMS();
                rms_pionerr[i-1] =    h1dPion[i]->GetRMSError();
                
                mean_electron[i-1] =  h1dE[i]->GetMean();
                mean_electronerr[i-1] = h1dE[i]->GetMeanError();
                rms_electron[i-1] = h1dE[i]->GetRMS();
                rms_electronerr[i-1] = h1dE[i]->GetRMSError();
                
                pt2[i-1] = (centbin[i]+centbin[i-1])/2;
                pt2err[i-1] = centbin[i]-pt2[i-1];
                
                
            }
            if (i==0) {cE[i]->Update();cE[i]->SaveAs(Form("~/Desktop/electronDcaByCent.pdf"));}
            cPion[i]->Update();
            if (i==0) cPion[i]->SaveAs(Form("~/Desktop/pionDcaByCent.pdf"));
            else cPion[i]->SaveAs(Form("~/Desktop/pionAndElectronDca_mRefMult%d_%d.pdf",centbin[i-1],centbin[i]));
        }

        TCanvas * cc=  new TCanvas("cc","cc",500,400);
        TH1F * dum = new TH1F("dum","dum;mRefMult;Meam [cm]",100,0,443);
        dum->SetMaximum(0.002);
        dum->SetMinimum(-0.002);
        dum->SetLineStyle(2);
        dum->SetLineColor(2);
        dum->Draw();
        TGraphErrors * gra2 = new TGraphErrors(nbin-1, pt2, mean_pion,pt2err,mean_pionerr);
        gra2->SetMarkerSize(1);
        gra2->SetMarkerColor(2);
        gra2->SetMarkerStyle(20);
        gra2->Draw("plsame");
        TGraphErrors * gra = new TGraphErrors(nbin-1, pt2, mean_electron,pt2err,mean_electronerr);
        gra->SetMarkerSize(1);
        gra->SetMarkerColor(1);
        gra->SetMarkerStyle(20);
        gra->Draw("plsame");
        cc->SaveAs(Form("~/Desktop/mean_Cent%d.pdf",selCent));
        
        
        
        TCanvas * cc2=  new TCanvas("cc2","cc2",500,400);
        TH1F * dum2 = new TH1F("dum2","dum2;mRefMult;RMS [cm]",100,0,443);
        dum2->SetMaximum(0.04);
        dum2->SetMinimum(0);
        dum2->Draw();
        TGraphErrors * gra3 = new TGraphErrors(nbin-1, pt2, rms_pion,pt2err,rms_pionerr);
        gra3->SetMarkerSize(1);
        gra3->SetMarkerColor(2);
        gra3->SetMarkerStyle(20);
        gra3->Draw("plsame");
        TGraphErrors * gra4 = new TGraphErrors(nbin-1, pt2, rms_electron,pt2err,rms_electronerr);
        gra4->SetMarkerSize(1);
        gra4->SetMarkerColor(1);
        gra4->SetMarkerStyle(20);
        gra4->Draw("plsame");
        cc2->SaveAs(Form("~/Desktop/rms_Cent%d.pdf",selCent));
        
    }

    // hadrons dca distribution by ZDCx
    if (opt>>12 & 1) {
        TString infilename = "Ana_7.root";
        //        infilename="out_166_15166046.root";
        TFile * infile = new TFile(infilename);
        TTree * tPion = (TTree*)infile->Get("tIncPion");
        TTree * tE = (TTree*)infile->Get("tInc");
        
        const int nbin      = 7;
        const int xnbin     = 98;      // pt
        const float xmin    = 0.2;
        const float xmax    = 10;
        const int ynbin     = 101;      // dca
        const float ymin    = -0.1;
        const float ymax    =  0.1;
        
        TCanvas *cPion[nbin];
        TCanvas *cE[nbin];
        TH2F * h2dPion[nbin];
        TH2F * h2dE[nbin];
        TH1F * h1dPion[nbin];
        TH1F * h1dE[nbin];
        
        
        double mean_pion[nbin];
        double mean_pionerr[nbin];
        double rms_pion[nbin];
        double rms_pionerr[nbin];
        double mean_electron[nbin];
        double mean_electronerr[nbin];
        double rms_electron[nbin];
        double rms_electronerr[nbin];
        
        double pt2[nbin];
        double pt2err[nbin];
        
        
        TString baseCut = "pt < 1 && pt > 0.2 ";
        
        for (int i=0; i<nbin; i++) {
            TString cutE = baseCut + " && abs(nsige) < 2";
            TString cutPion = baseCut;
            
            cPion[i] = new TCanvas(Form("cPion%d",i),Form("cPion%d",i),0,0,500,400);
            if (i==0) cE[i] = new TCanvas(Form("cE%d",i),Form("cE%d",i),500,0,500,400);
            cPion[i]->SetLogz();
            if (i==0) cE[i]->SetLogz();
            
            if (i==0) {
                h2dPion[i] = new TH2F(Form("h2dPion%d",i),Form("h2dPion%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                h2dE[i] = new TH2F(Form("h2dE%d",i),Form("h2dE%d;p_{T} [GeV/c];DCA_{XY} [cm]",i),xnbin,xmin,xmax,ynbin,ymin,ymax);
                
                cPion[i]->cd();
                tPion->Draw(Form("dca:pt>>h2dPion%d",i),cutPion,"col2");
                
                cE[i]->cd();
                tE->Draw(Form("dca:pt>>h2dE%d",i),cutE,"col2");
            }
            
            else {
                cPion[i]->SetLogy();
                
                h1dPion[i] = new TH1F(Form("h1dPion%d",i),Form("h1dPion%d;DCA_{XY} [cm];Entries [#]",i),ynbin,ymin,ymax);
                h1dE[i] = new TH1F(Form("h1dE%d",i),Form("h1dE%d",i),ynbin,ymin,ymax);
                
                cPion[i]->cd();
                
                TString centcut = Form("&& mZDCx > %f && mZDCx <= %f", zdcxbin[i-1], zdcxbin[i]);
                tPion->Draw(Form("dca>>h1dPion%d",i),cutPion + centcut,"");
                tE->Draw(Form("dca>>h1dE%d",i),cutE + centcut,"");
                h1dPion[i]->Sumw2();
                h1dE[i]->Sumw2();
                
                h1dPion[i]->Scale(1./h1dPion[i]->GetMaximum());
                h1dE[i]->Scale(1./h1dE[i]->GetMaximum());
                
                h1dE[i]->SetMarkerStyle(25);
                h1dE[i]->SetMarkerColor(1);
                h1dE[i]->SetMarkerSize(0.3);
                
                h1dPion[i]->SetMarkerStyle(24);
                h1dPion[i]->SetMarkerColor(2);
                h1dPion[i]->SetMarkerSize(0.3);
                
                h1dE[i]->SetLineColor(1);
                h1dPion[i]->SetLineColor(2);
                
                h1dPion[i]->SetMaximum(3);
                h1dPion[i]->SetMinimum(1e-5);
                
                h1dPion[i]->Draw("");
                h1dE[i]->Draw("same");
                h1dPion[i]->Draw("same");
                
                
                
                pt3 = new TPaveText(-0.05,1e-2,0.05,1e-4);
                pt3->AddText(Form("%.0fk < mZDCx #leq %.0fk",zdcxbin[i-1]/1000, zdcxbin[i]/1000 ));
                pt3->SetFillStyle(0);
                pt3->SetLineColorAlpha(kBlue, 0);
                pt3->Draw("same");
                
                mean_pion[i-1] =      h1dPion[i]->GetMean();
                mean_pionerr[i-1] =   h1dPion[i]->GetMeanError();
                rms_pion[i-1] =       h1dPion[i]->GetRMS();
                rms_pionerr[i-1] =    h1dPion[i]->GetRMSError();
                
                mean_electron[i-1] =  h1dE[i]->GetMean();
                mean_electronerr[i-1] = h1dE[i]->GetMeanError();
                rms_electron[i-1] = h1dE[i]->GetRMS();
                rms_electronerr[i-1] = h1dE[i]->GetRMSError();
                
                pt2[i-1] = (zdcxbin[i]+zdcxbin[i-1])/2;
                pt2err[i-1] = zdcxbin[i]-pt2[i-1];
                
                
            }
            if (i==0) {cE[i]->Update();cE[i]->SaveAs(Form("~/Desktop/electronDcaByZDCx_Pt2.pdf"));}
            cPion[i]->Update();
            if (i==0) cPion[i]->SaveAs(Form("~/Desktop/pionDcaByZDCx_Pt2.pdf"));
            else cPion[i]->SaveAs(Form("~/Desktop/pionAndElectronDca_mZDCx%.0fk_%.0fk_Pt2.pdf",zdcxbin[i-1]/1000,zdcxbin[i]/1000));
        }
        
        TCanvas * cc=  new TCanvas("cc","cc",500,400);
        TH1F * dum = new TH1F("dum","dum;mZDCx;Meam [cm]",100,zdcxbin[0],zdcxbin[6]);
        dum->SetMaximum(0.002);
        dum->SetMinimum(-0.002);
        dum->SetLineStyle(2);
        dum->SetLineColor(2);
        dum->Draw();
        TGraphErrors * gra2 = new TGraphErrors(nbin-1, pt2, mean_pion,pt2err,mean_pionerr);
        gra2->SetMarkerSize(1);
        gra2->SetMarkerColor(2);
        gra2->SetMarkerStyle(20);
        gra2->Draw("plsame");
        TGraphErrors * gra = new TGraphErrors(nbin-1, pt2, mean_electron,pt2err,mean_electronerr);
        gra->SetMarkerSize(1);
        gra->SetMarkerColor(1);
        gra->SetMarkerStyle(20);
        gra->Draw("plsame");
        cc->SaveAs(Form("~/Desktop/mean_zdc%d_Pt2.pdf",selCent));
        
        
        
        TCanvas * cc2=  new TCanvas("cc2","cc2",500,400);
        TH1F * dum2 = new TH1F("dum2","dum2;mZDCx;RMS [cm]",100,zdcxbin[0],zdcxbin[6]);
        dum2->SetMaximum(0.04);
        dum2->SetMinimum(0);
        dum2->Draw();
        TGraphErrors * gra3 = new TGraphErrors(nbin-1, pt2, rms_pion,pt2err,rms_pionerr);
        gra3->SetMarkerSize(1);
        gra3->SetMarkerColor(2);
        gra3->SetMarkerStyle(20);
        gra3->Draw("plsame");
        TGraphErrors * gra4 = new TGraphErrors(nbin-1, pt2, rms_electron,pt2err,rms_electronerr);
        gra4->SetMarkerSize(1);
        gra4->SetMarkerColor(1);
        gra4->SetMarkerStyle(20);
        gra4->Draw("plsame");
        cc2->SaveAs(Form("~/Desktop/rms_zdc%d_Pt2.pdf",selCent));
        
    }

    // electron TOF pid : pairMass, nsige, 
    if (opt>>13 & 1) {
        cout << "electron TOF PID study" << endl;
        const int n = 12;
        
        TString base = "abs(eta) < 0.7 && nsige > -1";
        TString basecutUL = base + "&& pairCharge==0";
        TString basecutLS = base + "&& pairCharge!=0";
        TString ptcut = "&& pt > 1.5";
        TString drawingObj[n] = {"nsige","e0/pt/TMath::CosH(eta)","e/pt/TMath::CosH(eta)","zDist","phiDist","etaTowDist","phiTowDist","neta","nphi","pairMass","partner_nsige","dca"};
        
        TString infilename = "Ana_2.root";
        TFile * infile = new TFile(infilename);
        TTree * tPureE = (TTree*)infile->Get("tPureE");


        TCanvas * cc = new TCanvas("cc","cc",500,500);
        cc->SetLogz();
        TH2D * hPosition = new TH2D("hPosition","hPosition",500,-10,10,500,-10,10);
        tPureE->Draw("pairPositionX:pairPositionY>>hPosition",basecutUL + ptcut,"col2");
        cc->Update();
        TCanvas * cc2 = new TCanvas("cc2","cc2",500,500);
        cc2->SetLogz();
        TH2D * hPosition2 = new TH2D("hPosition2","hPosition2",500,-10,10,500,-10,10);
        tPureE->Draw("pairPositionX:pairPositionY>>hPosition2",basecutLS + ptcut,"col2");
        cc2->Update();
        
        
        int nbin[n] = {289, 100, 100, 100, 100, 100, 100, 10, 10, 100, 289, 100};
        double xmin[n] = {-13, 0, 0, -20, -0.1, -0.1, -0.1, 0, 0, 0, -13, -0.1};
        double xmax[n] = {13, 3, 3, 20, 0.1, 0.1, 0.1, 10, 10, 0.01, 13, 0.1};
        
        TCanvas * c[n];
        TH1D * h0;
        TH1D * h1;
        TH1D * h2;
        
        
        for (int i = 0; i<n; i++) {
            
            h0 = new TH1D("h0",drawingObj[i],nbin[i],xmin[i],xmax[i]);
            h1 = new TH1D("h1",drawingObj[i],nbin[i],xmin[i],xmax[i]);
            h2 = new TH1D("h2",drawingObj[i],nbin[i],xmin[i],xmax[i]);

            
            
            c[i] = new TCanvas(Form("c%d",i),Form("c%d",i),500,400);
         //   c[i]->SetLogy();
            
            tPureE->Draw(drawingObj[i] + ">>h1",basecutUL + ptcut,"");
            tPureE->Draw(drawingObj[i] + ">>h2",basecutLS + ptcut,"same");
            h1->Sumw2();
            h2->Sumw2();
            h0->Add(h1,h2,1,-1);
            h0->Draw("same");
            
            
            h1->SetLineColor(1);
            h2->SetLineColor(4);
            h0->SetLineColor(2);
            if (i==0){
                TF1 * fitfun = new TF1("fitfun","gaus(0)",-2,3);
                fitfun->SetParameters(  6.04752e+03, -4.30945e-01, 9.59738e-01);
                h0->Fit(fitfun,"LNOR");
                
                fitfun->SetLineColor(2);
                fitfun->SetLineStyle(2);
                fitfun->SetRange(-13,13);
                fitfun->Draw("same");
                
                cout << fitfun->GetChisquare() << "/" << fitfun->GetNDF() << endl;
            }
            if (i==1){
                TF1 * fitfun = new TF1("fitfun","gaus(0)",0.5,1.5);
                fitfun->SetParameters(100, 1., 1.);
                h0->Fit(fitfun,"LNOR");
                
                fitfun->SetLineColor(2);
                fitfun->SetLineStyle(2);
                fitfun->SetRange(-13,13);
                fitfun->Draw("same");
                
                cout << fitfun->GetChisquare() << "/" << fitfun->GetNDF() << endl;
            }        
            if (i==2){
                TF1 * fitfun = new TF1("fitfun","gaus(0)",0.5,1.5);
                fitfun->SetParameters(100, 1., 1.);
                h0->Fit(fitfun,"LNOR");
                
                fitfun->SetLineColor(2);
                fitfun->SetLineStyle(2);
                fitfun->SetRange(-13,13);
                fitfun->Draw("same");
                
                cout << fitfun->GetChisquare() << "/" << fitfun->GetNDF() << endl;
            }            
            if (i==10){
                TF1 * fitfun = new TF1("fitfun","gaus(0)",-3,3);
                fitfun->SetParameters(100, 0., 1.);
                h0->Fit(fitfun,"LNOR");
                
                fitfun->SetLineColor(2);
                fitfun->SetLineStyle(2);
                fitfun->SetRange(-13,13);
                fitfun->Draw("same");
                
                cout << fitfun->GetChisquare() << "/" << fitfun->GetNDF() << endl;
            }
            c[i]->SaveAs(Form("~/Desktop/c%d.pdf",i));
            cout << h0->Integral() << " " << h0->Integral(1,nbin[i]) << " " << h0->GetEntries() << endl;
        }
        

    }
        // hadron pid 
    if (opt>>14 & 1) {
        cout << "hadron PID study" << endl;
        const int n = 9;
        
        TString base = "abs(1-beta) < 0.025 && nsige > 0 && abs(eta) < 0.7";
        TString basecutUL = base + "&& pairCharge==0";
        TString basecutLS = base + "&& pairCharge!=0";
        TString ptcut = "&& pt > 1. && pt < 2.";
        TString drawingObj[n] = {"nsige","e0/pt/TMath::CosH(eta)","e/pt/TMath::CosH(eta)","zDist","phiDist","etaTowDist","phiTowDist","neta","nphi"};
        
        TString infilename = "Ana_2.root";
        TFile * infile = new TFile(infilename);
        TTree * tPureE = (TTree*)infile->Get("tPureE");

        int nbin[n] = {289, 100, 100, 100, 100, 100, 100, 10, 10};
        double xmin[n] = {-13, 0, 0, -20, -0.1, -0.1, -0.1, 0, 0};
        double xmax[n] = {13, 3, 3, 20, 0.1, 0.1, 0.1, 10, 10};
        
        TCanvas * c[n];
        TH1D * h0;
        TH1D * h1;
        TH1D * h2;
        
        for (int i = 0; i<n; i++) {
            
            h0 = new TH1D("h0","h0",nbin[i],xmin[i],xmax[i]);
            h1 = new TH1D("h1","h1",nbin[i],xmin[i],xmax[i]);
            h2 = new TH1D("h2","h2",nbin[i],xmin[i],xmax[i]);

            
            
            c[i] = new TCanvas(Form("c%d",i),Form("c%d",i),500,400);
         //   c[i]->SetLogy();
            
            tPureE->Draw(drawingObj[i] + ">>h1",basecutUL + ptcut,"");
            tPureE->Draw(drawingObj[i] + ">>h2",basecutLS + ptcut,"same");
            h1->Sumw2();
            h2->Sumw2();
            h0->Add(h1,h2,1,-1);
            h0->Draw("same");
            
            
            h1->SetLineColor(1);
            h2->SetLineColor(4);
            h0->SetLineColor(2);
            if (i==0){
                TF1 * fitfun = new TF1("fitfun","gaus(0)",-3,3);
                fitfun->SetParameters(  6.04752e+03, -4.30945e-01, 9.59738e-01);
                h0->Fit(fitfun,"LNOR");
                
                fitfun->SetLineColor(2);
                fitfun->SetLineStyle(2);
                fitfun->SetRange(-13,13);
                fitfun->Draw("same");
                
                cout << fitfun->GetChisquare() << "/" << fitfun->GetNDF() << endl;
            }
        }
        

    }

}
