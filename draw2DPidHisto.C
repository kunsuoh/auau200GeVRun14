const int nbin = 101;
double ptbin[nbin+1] = {0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.05, 1.1,  1.15,  1.2,  1.25,  1.3,  1.35,  1.4,  1.45,  1.5,  1.55,  1.6,  1.65,  1.7,  1.75,  1.8,  1.85,  1.9,  1.95, 2,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5};


float xMeanE[5][101];
float xSigmaE[5][101];
float yMeanE[5][101];
float ySigmaE[5][101];


float xMeanPi[5][101];
float xSigmaPi[5][101];
float yMeanPi[5][101];
float ySigmaPi[5][101];

void draw2DPidHisto(TString infilename = "Ana/out_165_15165044.root"){
    TFile * infile = new TFile(infilename);
    THnSparse * hs = infile->Get("hsTrackPid");
    loadEIDcut();

    TH2D * histo2d[5][nbin];
    TH2D * histo2dCut[5][nbin];
    for (unsigned long i = 0; i<5; i++) {
        for (int j = 0 ; j<nbin; j++) {
            histo2d[i][j] = new TH2D(Form("histo_%d_%d",i,j),Form("histo_%d_%d",i,j),800,-0.2,0.6,289,-13,13);
            histo2dCut[i][j] = new TH2D(Form("histoCut_%d_%d",i,j),Form("histoCut_%d_%d",i,j),800,-0.2,0.6,289,-13,13);
        }
    }
    
    int bins[4] =   {20, 1000,   289, 800};
    double mins[4] = {-1, 0,      -13, -0.2};
    double maxs[4] = {1, 10,     13,  0.6};
    
    
    TH1D * h[4];
    double var[4];
    for (int j=0;j<4;j++) {
        h[j] = new TH1D(Form("h%d",j),Form("h%d",j),bins[j],mins[j],maxs[j]);
    }
    
    
    unsigned long nentries = hs->GetNbins();
  //  nentries= 1000000;
    cout << "Total: " << nentries << " tracks" << endl;
    for (unsigned long i=0;i<nentries;i++){
        Int_t coord[4]={0};
        hs->GetBinContent(i, coord);
        
        for (int j=0;j<4;j++) {
            var[j] = ((double)coord[j]- 0.5)/bins[j]*(maxs[j]-mins[j]) + mins[j];
            h[j]->Fill(var[j]);
        }
        if (i%100000 == 0) cout << i << " " << var[0]  << " " << var[1]  << " " << var[2]  << " " << var[3] <<endl;
        
        // set iPt and iEta
        int iPt=0;
        int iEta=0;
        double Eta = var[0];
        double Pt = var[1];
        
        if (Pt > 5.) continue;
        if (Eta > 0.5 || Eta < -0.5) continue;
        
        for (int k=0;k<nbin+1;k++) {
            if (Pt < ptbin[k]) {
                iPt=k-1;
                break;
            }
        }
        for (int k=0;k<5;k++) if (TMath::Abs(Eta) < k*0.1 + 0.1) {iEta=k;break;}
        //       cout << iPt << " " << iEta << " " << var[0] << " " << var[1] << " " << var[2] << " " << var[3] << endl;
        histo2d[iEta][iPt]->Fill(var[3],var[2]);
        
        float x = var[3];
        float y = var[2];
        
        float x0 = xMeanE[iEta][iPt];
        float y0 = yMeanE[iEta][iPt];
        float a0 = xSigmaE[iEta][iPt] * 3;
        float b0 = ySigmaE[iEta][iPt] * 3;
        
        float x1 = xMeanPi[iEta][iPt];
        float y1 = yMeanPi[iEta][iPt];
        float a1 = xSigmaPi[iEta][iPt] * 3.5;
        float b1 = ySigmaPi[iEta][iPt] * 3.5;

        if(((x-x0)/a0)**2 + ((y-y0)/b0)**2 < 1 && ((x-x1)/a1)**2 + ((y-y1)/b1)**2 > 1 ) {
            histo2dCut[iEta][iPt]->Fill(var[3],var[2]);
        }
    }
    cout << "Save histograms..." << endl;
    TFile * fout = new TFile("out_2dHisto.root","recreate");
    fout->cd();
    for (unsigned long i = 0; i<5; i++) {
        for (int j = 0 ; j<nbin; j++) {
            histo2d[i][j]->Write();
            histo2dCut[i][j]->Write();
        }
    }
    fout->Close();
}

Int_t loadEIDcut(){
    TString p0;
    float xMean, xSigma, yMean, ySigma;
    
    for (int j=0; j<5; j++) {
        TString inparaname = Form("inpara/parameter_ver13_Cent4_Eta%d.txt",j);
        ifstream inpara(inparaname);
        if (!inpara) {
            cout << "Failed to load parameter... " << endl;
            continue;
        } else cout << inparaname << endl;
        float dum;
        for (int k=0; k<101; k++){
            for (int l=0; l<8; l++) {
                inpara >> p0;
                inpara >> xMean;
                inpara >> xSigma;
                inpara >> yMean;
                inpara >> ySigma;
                inpara >> dum;
                inpara >> dum;
                inpara >> dum;
                inpara >> dum;
                inpara >> dum;
                //         cout << i << " " << j << " " << k << " " << l << " " << p0 << " " << xMean << " " << xSigma << " " << yMean << " " << ySigma << " /// dum=" << dum << endl;
                
                if (l==2){
                    xMeanPi[j][k]=xMean;
                    xSigmaPi[j][k]=xSigma;
                    yMeanPi[j][k]=yMean;
                    ySigmaPi[j][k]=ySigma;
                }
                if (l==0){
                    xMeanE[j][k]=xMean;
                    xSigmaE[j][k]=xSigma;
                    yMeanE[j][k]=yMean;
                    ySigmaE[j][k]=ySigma;
                }
            }
        }
    }
    
    return 1;
}