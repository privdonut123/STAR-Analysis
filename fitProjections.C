//#include <iostream>
#include <vector>
using namespace std;
#include <cmath>
// #include <TSystem.h>
// #include <TSystemDirectory.h>
// #include <TSystemFile.h>
// #include <TList.h>

class FitAndEnergy: public TObject
    {
        public:
            FitAndEnergy() {}
            FitAndEnergy(double amplitude_fit, double mean_fit, double mean_fit_error, double std_dev_fit, double ecal_energy = -1, double nEntries = -1, string fit_type = "gaussian")
                {
                    this->ecal_energy = ecal_energy;
                    this->amplitude_fit = amplitude_fit;
                    this->mean_fit = mean_fit;
                    this->mean_fit_error = mean_fit_error;
                    this->std_dev_fit = std_dev_fit;
                    this->nEntries = nEntries;
                }
            FitAndEnergy(double amplitude_landau_fit, double mpv_landau_fit, double scale_landau_fit, double ecal_energy = -1, double nEntries = -1, string fit_type = "landau")
                {
                    this->ecal_energy = ecal_energy;
                    this->amplitude_landau_fit = amplitude_landau_fit;
                    this->mpv_landau_fit = mpv_landau_fit;
                    this->scale_landau_fit = scale_landau_fit;
                    this->nEntries = nEntries;
                }
            
            // Accessors
            double GetEcalEnergy() const { return ecal_energy; }
            double GetAmplitudeFit() const { return amplitude_fit; }
            double GetSkewnessFit() const { return skewness_fit; }
            double GetMeanFit() const { return mean_fit; }
            double GetMeanFitError() const { return mean_fit_error; }
            double GetStdDevFit() const { return std_dev_fit; }
            double GetAmplitudeLandauFit() const { return amplitude_landau_fit; }
            double GetMpvLandauFit() const { return mpv_landau_fit; }
            double GetScaleLandauFit() const { return scale_landau_fit; }
            double GetNEntries() const { return nEntries; }

            // Mutators
            void SetEcalEnergy(double energy) { this->ecal_energy = energy; }
            void SetAmplitudeFit(double amplitude) { this->amplitude_fit = amplitude; }
            void SetSkewnessFit(double skewness) { this->skewness_fit = skewness; }
            void SetMeanFit(double mean) { this->mean_fit = mean; }
            void SetMeanFitError(double mean_error) { this->mean_fit_error = mean_error; }
            void SetStdDevFit(double std_dev) { this->std_dev_fit = std_dev; }
            void SetAmplitudeLandauFit(double amplitude_landau) { this->amplitude_landau_fit = amplitude_landau; }
            void SetMpvLandauFit(double mpv_landau) { this->mpv_landau_fit = mpv_landau; }
            void SetScaleLandauFit(double scale_landau) { this->scale_landau_fit = scale_landau; }
            void SetNEntries(double entries) { this->nEntries = entries; }

            void skew_gauss_dummy_to_full_fit(double amplitude_guess, double mean_guess, double std_dev_guess, TH1F *h1);
            void landau_dummy_to_full_fit(double amplitude_guess, double mpv_guess, double scale_guess, TH1F *h1);

        private:
            double ecal_energy;
            
            //Fit parameters for skewed gaussian fit
            double skewness_fit;
            double amplitude_fit;
            double mean_fit;
            double mean_fit_error;
            double std_dev_fit;

            //Fit parameters for landau fit
            double amplitude_landau_fit;
            double mpv_landau_fit;
            double scale_landau_fit;

            double nEntries;
            //ClassDef(FitAndEnergy, 1);
            

    };

void FitAndEnergy::skew_gauss_dummy_to_full_fit(double amplitude_guess, double mean_guess, double std_dev_guess, TH1F *h1)
    {
        TF1 *fit_dummy = new TF1("fit_dummy", "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", mean_guess - std_dev_guess, mean_guess + std_dev_guess);
        fit_dummy->SetParameters(amplitude_guess, mean_guess, std_dev_guess, 0);
        fit_dummy->SetParLimits(0, 0, amplitude_guess*5);
        fit_dummy->SetParLimits(1, mean_guess - std_dev_guess, mean_guess + std_dev_guess);
        fit_dummy->SetParLimits(2, 0, std_dev_guess*2);
        fit_dummy->SetParLimits(3, -10, 10); // Adjust the limits for the skewness parameter
        fit_dummy->FixParameter(3, 0);
        h1->Fit(fit_dummy, "RNQ");
        this->amplitude_fit = fit_dummy->GetParameter(0);
        this->mean_fit = fit_dummy->GetParameter(1);
        this->mean_fit_error = fit_dummy->GetParError(1);
        this->std_dev_fit = fit_dummy->GetParameter(2);
        this->skewness_fit = fit_dummy->GetParameter(3); // Add the skewness parameter
        this->nEntries = h1->GetEntries();
        delete fit_dummy;
    }

void FitAndEnergy::landau_dummy_to_full_fit(double amplitude_guess, double mpv_guess, double scale_guess, TH1F *h1)
    {
        TF1 *fit_dummy = new TF1("fit_dummy", "[0]*TMath::Landau(x, [1], [2])", mpv_guess - 1, mpv_guess + 1); // No reason for choosing 3, just a guess
        fit_dummy->SetParameters(amplitude_guess, mpv_guess, scale_guess);
        fit_dummy->SetParLimits(0, 0, amplitude_guess*5);
        fit_dummy->SetParLimits(1, mpv_guess, mpv_guess + 1);
        fit_dummy->SetParLimits(2, 0, 2.5);
        h1->Fit(fit_dummy, "RNQ");
        this->amplitude_landau_fit= fit_dummy->GetParameter(0);
        this->mpv_landau_fit = fit_dummy->GetParameter(1);
        this->scale_landau_fit = fit_dummy->GetParameter(2);
        this->nEntries = h1->GetEntries();
        delete fit_dummy;
    }

//ClassImp(FitAndEnergy);

void fitProjections() 
    {
        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
        loadSharedLibraries();
        gSystem->Load("StFcsDbMaker");
        gSystem->Load("libStFcsWaveformFitMaker.so");
        //PeakAna *pA = new PeakAna();
        
        TSystemDirectory dir("dir", "input/fit_proj_input/");
        TList *files = dir.GetListOfFiles();

        if (files) 
            {
                TSystemFile *file;
                TString fname;
                TIter next(files);
                TFile *out_file = new TFile("./output/output_projections.root", "RECREATE"); // Create the output file outside the loop
                while ((file=(TSystemFile*)next())) 
                {
                    fname = file->GetName();
                    if (!file->IsDirectory() && fname.EndsWith(".root")) 
                    {
                        TString filePath = "input/fit_proj_input/" + fname;
                        TFile *in_file = TFile::Open(filePath); // Open the file

                        TDirectory *subdir = out_file->mkdir(fname); // Create a subdirectory for this file
                        subdir->cd(); // Change to the subdirectory

                        TH2F *h2_ecal_hcal_cluster_energy = (TH2F*)in_file->Get("h2_ecal_hcal_cluster_energy"); // replace with your histogram name
                        TH2F *h2_ecal_hcal_hit_energy = (TH2F*)in_file->Get("h2_ecal_hcal_hit_energy"); // replace with your histogram name
                        cout << "Current subdirectory in out_file: " << subdir->GetName() << endl;
                        projections_and_fits(subdir, h2_ecal_hcal_cluster_energy, "cluster");
                        subdir->GetMotherDir()->cd(); // Change back to the parent directory
                        cout << "Current directory in out_file: " << subdir->GetMotherDir()->GetName() << endl;
                        projections_and_fits(subdir, h2_ecal_hcal_hit_energy, "hit");
                        
                        in_file->Close();
                        
                    }
                    
                }
                out_file->Close(); // Close the output file after the loop
            }
        // TFile *file = TFile::Open("output/output_pi-e60.root"); // replace with your file name


        // TFile *outfile = new TFile("./output/output_projections.root", "RECREATE");
        // TH2F *h2_ecal_hcal_cluster_energy = (TH2F*)file->Get("h2_ecal_hcal_cluster_energy"); // replace with your histogram name
        // TH2F *h2_ecal_hcal_hit_energy = (TH2F*)file->Get("h2_ecal_hcal_hit_energy"); // replace with your histogram name
        // projections_and_fits(outfile, h2_ecal_hcal_cluster_energy, "cluster");
        // projections_and_fits(outfile, h2_ecal_hcal_hit_energy, "hit");
        // outfile->Close();
        // file->Close();
    }

void projections_and_fits(TDirectory *dir, TH2F *h2, TString cluster_or_hit)
    { 
        PeakAna *pA = new PeakAna();
        TDirectory *subdir;
        if (cluster_or_hit == "cluster")
        {
            subdir = dir->mkdir("Clusters");
        }
        else if (cluster_or_hit == "hit")
        {
            subdir = dir->mkdir("Hits");
        }
        subdir->cd();

        TDirectory *subdir1 = dir->mkdir("Preliminary_Summed_Projections_Ecal_on_Hcal"); // Change type to TDirectory*
        TDirectory *subdir2 = dir->mkdir("Peak_Finder_Summed_Projections_Ecal_on_Hcal"); // Change type to TDirectory*

        
        int nbinsX = h2->GetXaxis()->GetNbins();
        TClonesArray *fits_and_energy = new TClonesArray("FitAndEnergy");
        /*
        // Create a dynamic 2D array
        double** ecal_energies_and_fits = new double*[3];
        for(int i = 0; i < 3; ++i)
            {
                ecal_energies_and_fits[i] = new double[nbinsX];
            }
        */
        // Loop over the x-axis bins
        //vector<FitAndEnergy> fits_and_energy;

        for (int i = 1; i <= nbinsX; ++i) 
            {
                subdir1->cd();
                // Create a 1D histogram of the y-axis projection for this bin
                TH1D *h1 = h2->ProjectionY(Form("h1_%d", i), i, i);
                if (h1->GetEntries() == 0)
                    {
                        continue;
                    }

                double amplitude_guess = h1->GetMaximum();
                double mean_guess = h1->GetBinCenter(h1->GetMaximumBin());
                double std_dev_guess = h1->GetRMS(); // Get the RMS as an initial guess for the standard deviation
                double skewness_guess = 0;
                /*
                for (int i = 1; i <= h1->GetNbinsX(); ++i) 
                    {
                        double standardized_value = (i - mean_guess) / std_dev_guess;
                        double cubed_value = pow(standardized_value, 3);
                        double contribution = cubed_value * h1->GetBinContent(i);
                        skewness_guess += contribution;
                    }
                */
                // cout << "Estimated amplitude: " << amplitude_guess << endl;
                // cout << "Estimated mean: " << mean_guess << endl;
                // cout << "Estimated standard deviation: " << std_dev_guess << endl;
                // cout << "Estimated skewness: " << skewness_guess << endl;
                // Fith1_clean
                TF1 *fit_dummy = new TF1(Form("fit1_%d", i), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", mean_guess - std_dev_guess, mean_guess + std_dev_guess);
                fit_dummy->SetParameters(amplitude_guess, mean_guess, std_dev_guess, skewness_guess);
                fit_dummy->SetParLimits(0, 0, amplitude_guess*5);
                fit_dummy->SetParLimits(1, h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                fit_dummy->SetParLimits(2, 0, std_dev_guess*5);
                fit_dummy->SetParLimits(3, skewness_guess-10, 0);
                fit_dummy->FixParameter(3, 0);
                h1->Fit(fit_dummy, "RNQ");
                FitAndEnergy *fit_and_energy = (FitAndEnergy*)fits_and_energy->ConstructedAt(i-1);
                fit_and_energy->SetEcalEnergy(h2->GetXaxis()->GetBinCenter(i));
                fit_and_energy->SetAmplitudeFit(fit_dummy->GetParameter(0));
                fit_and_energy->SetMeanFit(fit_dummy->GetParameter(1));
                fit_and_energy->SetMeanFitError(fit_dummy->GetParError(1));
                fit_and_energy->SetStdDevFit(fit_dummy->GetParameter(2));
                fit_and_energy->SetNEntries(h1->GetEntries());
                double xmin = fit_and_energy->GetMeanFit() - fit_and_energy->GetStdDevFit();
                double xmax = fit_and_energy->GetMeanFit() + fit_and_energy->GetStdDevFit();
                
                TF1 *fit1 = new TF1(Form("fit1_%d", i), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                fit1->SetParameters(fit_and_energy->GetAmplitudeFit(), fit_and_energy->GetMeanFit(), fit_and_energy->GetStdDevFit(), fit_and_energy->GetSkewnessFit());
                fit1->FixParameter(0, fit_and_energy->GetAmplitudeFit());
                fit1->FixParameter(1, fit_and_energy->GetMeanFit());
                fit1->FixParameter(2, fit_and_energy->GetStdDevFit());
                fit1->FixParameter(3, fit_and_energy->GetSkewnessFit());
                fit1->SetNpx(1000);

                // Create a canvas to draw the histogram and fit
                double ecal_energy = h2->GetXaxis()->GetBinCenter(i);
                //TCanvas *c1 = new TCanvas("c1","c1", 800, 600);
                h1->SetTitle(Form("Projected ECAL Energy on HCal: %.2f (Gev)", ecal_energy));
                TCanvas *c1 = new TCanvas(Form("ECAL_Energy_%.2f", ecal_energy), Form("ECAL_Energy:_%.2f", ecal_energy), 800, 600);
                //gStyle->SetOptStat(0); 
                c1->SetTitle(Form("ECAL Energy: %.2f", ecal_energy));
                // Draw the histogram
                // Create a legend
                TH1F *h1_new = (TH1F*)h1->Clone(Form("h1_new_%d", i));
                h1_new->GetYaxis()->SetTitle("Counts");
                h1_new->GetXaxis()->SetTitle("Hcal Summed Energy (GeV)");
                h1_new->Fit(fit1,"RQ");
                //h1_new->Fit(fit2,"R+"); // Add the "+" option to add the second fit to the same histogram
                h1_new->SetStats(0);
                h1_new->Draw();
                //fit1->Draw("SAME");
                //fit2->Draw("SAME");
                // new TH1F(Form("h1_new_%d", i), "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                
                //h1_new->GetListOfFunctions()->Add(fit1);
                double chi2_ndf_whole_range = fit1->GetNDF() != 0 ? fit1->GetChisquare()/fit1->GetNDF() : -1;
                double chi2_ndf_fit_range = fit_dummy->GetNDF() != 0 ? fit_dummy->GetChisquare()/fit_dummy->GetNDF() : -1;
                //fit1->SetLineColor(kRed);
                
                TLegend *leg = new TLegend(0.5, 0.5, 0.9, 0.9); // coordinates are in the pad's fraction
                leg->SetTextSize(0.10); // Set text size to 4% of the pad height
                leg->AddEntry("Skewed gaussian fit parameters");
                leg->AddEntry((TObject*)0, Form("Amplitude (A): %.3f", fit1->GetParameter(0)), "");
                leg->AddEntry((TObject*)0, Form("Mean (mu): %.3f", fit1->GetParameter(1)), "");
                leg->AddEntry((TObject*)0, Form("Mean fit: %.3f", fit1->Mean(fit1->GetXmin(), fit1->GetXmax())), "");
                leg->AddEntry((TObject*)0, Form("Standard Deviation (sigma): %.3f", fit1->GetParameter(2)), "");
                leg->AddEntry((TObject*)0, Form("Skewness (alpha): %.3f", fit1->GetParameter(3)), "");
                //leg->AddEntry((TObject*)0, Form("Chi^2: %.3f", fit1->GetChisquare()), "");
                //leg->AddEntry((TObject*)0, Form("NDF: %.3f", fit1->GetNDF()), "");
                leg->AddEntry((TObject*)0, Form("Chi^2/NDF whole range: %.3f", chi2_ndf_whole_range), "");
                leg->AddEntry((TObject*)0, Form("Chi^2/NDF fit range: %.3f", chi2_ndf_fit_range), "");
                leg->AddEntry((TObject*)0, Form("Histogram Mean: %.3f", mean_guess), "");
                leg->AddEntry((TObject*)0, Form("Histogram RMS: %.3f", std_dev_guess), "");
                leg->Draw();
                
                if(i==5)
                    {
                        c1->SaveAs("test.png");
                    }
                // Write the canvas to the ROOT file
                c1->Write();
                delete c1;
                dir->cd();
                
                subdir2->cd();

                TCanvas *c1_alt = new TCanvas(Form("Prelim_ECAL_Energy_%i", (int)ecal_energy), Form("Prelim_ECAL_Energy:_%i", (int)ecal_energy), 800, 600);
                TH1F *h1_clean = new TH1F(Form("h1_clean_%d", i), "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                for (int bin = 1; bin <= h1->GetNbinsX(); ++bin)
                    {
                        h1_clean->SetBinContent(bin, h1->GetBinContent(bin));
                        //h1_new->SetBinError(bin, h1->GetBinError(bin));
                    }
                //gStyle->SetOptStat(0); 
                //c1_alt->SetTitle(Form("ECAL Energy: %.2f", ecal_energy));
                //h1_clean->Draw();
                pA->SetData(h1_clean);
                pA->SetBaseline(0,.75);
                pA->SetRange(0, 0, 100, 1000000);
                pA->SetFilter(2,5);
                pA->SetTunnelThreshold(0.1);
                pA->SetTunnelScale(0.1);
                pA->SetTunnelSigma(10);
                pA->AnalyzeForPeak();
                //pA->Print();
                for(int i_peak = 0; i_peak < pA->NPeaks(); i_peak++)
                    {
                        FitAndEnergy *peak_fit_param = new FitAndEnergy();
                        if (i_peak == 0 && i_peak < pA->NPeaks() - 1)
                            {
                                cout << "MIP found at " << pA->GetPeak(i_peak).mPeakX << " with amplitude " << pA->GetPeak(i_peak).mPeakY << endl;
                                cout << "Attempting to fit MIP with Landau distribution" << endl;
                                peak_fit_param->landau_dummy_to_full_fit(h1_clean->GetBinContent(1), 1, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                                cout << "Amplitude: " << peak_fit_param->GetAmplitudeLandauFit() << endl;
                                TF1 *peak_fit = new TF1(Form("peak_fit_%d", i_peak), "[0]*TMath::Landau(x, [1], [2])", h1_clean->GetXaxis()->GetXmin(), h1_clean->GetXaxis()->GetXmax());
                                peak_fit->FixParameter(0, peak_fit_param->GetAmplitudeLandauFit());
                                peak_fit->FixParameter(1, peak_fit_param->GetMpvLandauFit());
                                peak_fit->FixParameter(2, peak_fit_param->GetScaleLandauFit());
                            }
                        else
                            {
                                cout << "Attempting to fit peak " << i_peak << " with skewed gaussian distribution" << endl;
                                peak_fit_param->skew_gauss_dummy_to_full_fit(pA->GetPeak(i_peak).mPeakY, pA->GetPeak(i_peak).mPeakX, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                                TF1 *peak_fit = new TF1(Form("peak_fit_%d", i_peak), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1_clean->GetXaxis()->GetXmin(), h1_clean->GetXaxis()->GetXmax());
                                peak_fit->FixParameter(0, peak_fit_param->GetAmplitudeFit());
                                peak_fit->FixParameter(1, peak_fit_param->GetMeanFit());
                                peak_fit->FixParameter(2, peak_fit_param->GetStdDevFit());
                                peak_fit->FixParameter(3, 0);
                            }
                        //peak_fit_param->dummy_fit_to_full_fit(pA->GetPeak(i_peak).mPeakY, pA->GetPeak(i_peak).mPeakX, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                        // peak_fit_param->ecal_energy = ecal_energy;
                        // peak_fit_param->amplitude_fit = pA->GetPeak(i_peak).mPeakY;
                        // peak_fit_param->mean_fit = pA->GetPeak(i_peak).mPeakX;
                        // peak_fit_param->std_dev_fit = (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 9; //I have no reason why I chose 9, just that visually the fits seem to approximate the peaks well
                        // peak_fit_param->nEntries = h1_clean->GetEntries();

                        
                        peak_fit->SetLineColor(kBlue);
                        // peak_fit->FixParameter(0, amplitude_fit);
                        // peak_fit->FixParameter(1, mean_fit);
                        // peak_fit->FixParameter(2, std_dev_fit);
                        // peak_fit->FixParameter(3, skewness_fit);
                        peak_fit->SetNpx(1000);
                        if (pA->PeakTunnel(pA->GetPeak(i_peak)) == false)
                            {
                                cout << "Not-A-Peak Probability: " << pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) << endl;
                                // cout << "Peak " << i_peak << " is not a good peak" << endl;
                                peak_fit->SetLineColor(kRed);
                                //continue;
                            }
                        // cout << "Peak Window: " << pA->GetPeak(i_peak).mStartX << " to " << pA->GetPeak(i_peak).mEndX << endl;
                        // cout << "Fitting peak " << i_peak << " with parameters: " << pA->GetPeak(i_peak).mPeakY << ", " << pA->GetPeak(i_peak).mPeakX << ", " << pA->GetPeak(i_peak).mEndX-pA->GetPeak(i_peak).mStartX << ", " << skewness_fit << endl;
                        h1_clean->Fit(peak_fit, "RQ+");
                        //peak_fit->Draw("LSAME");
                        // cout << "Peak " << i_peak << " at (" << pA->GetPeak(i_peak).mPeakX << ", " << pA->GetPeak(i_peak).mPeakY << ")" << endl;
                        //delete peak_fit_param;
                        delete peak_fit_param;
                    }
                h1_clean->SetStats(0);
                h1_clean->SetTitle(Form("ECAL Energy: %.2f (Gev)", ecal_energy));
                h1_clean->GetXaxis()->SetTitle("Hcal Summed Energy (GeV)");
                h1_clean->GetYaxis()->SetTitle("Counts");
                h1_clean->Draw();
                // pA->Draw("PL;A");
                TPaveText *pt = new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC");
                pA->AddPeakStats(pt, "ad");
                pt->Draw();

                if (i == 4)
                    {
                        c1_alt->SaveAs("test2.png");
                    }
                
                
                //pA->Get_mPeaks()
                // pA->Print();

                // pA->Draw("PL;A");
                // if(i == 4)
                //     {
                //         pA->Write();
                //     }
                
                // cout << "Saving Peak finding histogram" << endl;
                //c1_alt->Update();
                c1_alt->Write();
                //c1_alt->Print("test.png")
                dir->cd();
                delete c1_alt;
            }
        //outfile->cd();
        TF1 *f2 = new TF1("fit_2", "[0] + [1]*x", 0, 100);
        dir->cd();
        //h2_ecal_hcal_cluster_energy->Fit(f2);
        if (cluster_or_hit == "cluster")
            {
                TCanvas *c2 = new TCanvas("h2_ecal_hcal_cluster_energy","h2_ecal_hcal_cluster_energy", 800, 600);
            }
        else if (cluster_or_hit == "hit")
            {
                TCanvas *c2 = new TCanvas("h2_ecal_hcal_hit_energy","h2_ecal_hcal_hit_energy", 800, 600);
            }
        else
            {
                cout << "Error: cluster_or_hit must be either \"cluster\" or \"hit\"" << endl;
            }
        
        
        h2->SetStats(0);
        h2->Draw("COLZ");

        // Create a TGraphErrors object
        TGraphErrors *gr = new TGraphErrors();

        // Loop over the bins
        
        for (int j = 1; j <= nbinsX; ++j) 
            {
                FitAndEnergy *fit_and_energy = (FitAndEnergy*)fits_and_energy->At(j);
                // Set the point and its error in the TGraphErrors object
                if (fit_and_energy == 0 ||(fit_and_energy != 0 && fit_and_energy->GetMeanFit() == 0))
                    {
                        continue;
                    }

                gr->SetPoint(gr->GetN(), fit_and_energy->GetEcalEnergy(), fit_and_energy->GetMeanFit());
                double weight = 1.0 / sqrt(fit_and_energy->GetNEntries()); // calculate weight
                gr->SetPointError(gr->GetN()-1, 0, fit_and_energy->GetMeanFitError());
            }

        // Draw the TGraphErrors object on the canvas
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kViolet); // Change marker color to purple
        gr->SetLineColor(kViolet); // Change error bars color to purple
        gr->Draw("P same"); // "P" option for points, "same" to draw on the same canvas
        gr->Fit(f2,"Q"); // "Q" option for quiet mode
        //f2->Draw("same");

        // Add legend for f2 parameters
        
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(f2, "f2 Parameters", "l");
        legend->AddEntry((TObject*)0, Form("y-intercept: %.2f", f2->GetParameter(0)), "");
        legend->AddEntry((TObject*)0, Form("slope: %.2f", f2->GetParameter(1)), "");
        legend->AddEntry((TObject*)0, Form("chi^2: %f", f2->GetChisquare()), "");
        legend->AddEntry((TObject*)0, Form("NDF: %i", f2->GetNDF()), "");
        legend->AddEntry((TObject*)0, Form("chi^2/NDF: %f", f2->GetChisquare()/f2->GetNDF()), "");
        legend->Draw();
        
        //c2->cd();
        c2->Write();
        // Close the ROOT files
        dir->cd();
            
    }