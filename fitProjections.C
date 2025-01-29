#include <ostream>
#include <vector>
using namespace std;
#include <cmath>
// #include <TSystem.h>
// #include <TSystemDirectory.h>
// #include <TSystemFile.h>
// #include <TList.h>
// #include <iostream>
const string red("\033[0;31m");
const string green("\033[1;32m");
const string yellow("\033[1;33m");
const string cyan("\033[0;36m");
const string magenta("\033[0;35m");
const string reset("\033[0m");

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
            double GetChi2NdfFitRange() const { return chi2_ndf_fit_range; }

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

            void skew_gauss_dummy_to_full_fit(double amplitude_guess, double mean_guess, double std_dev_guess, TH1 *h1);
            void landau_dummy_to_full_fit(double amplitude_guess, double mpv_guess, double scale_guess, TH1 *h1);

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
            double chi2_ndf_fit_range;
            //ClassDef(FitAndEnergy, 1);
            

    };

void FitAndEnergy::skew_gauss_dummy_to_full_fit(double amplitude_guess, double mean_guess, double std_dev_guess, TH1 *h1)
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
        this->chi2_ndf_fit_range = fit_dummy->GetNDF() != 0 ? fit_dummy->GetChisquare()/fit_dummy->GetNDF() : -1;
        delete fit_dummy;
    }

void FitAndEnergy::landau_dummy_to_full_fit(double amplitude_guess, double mpv_guess, double scale_guess, TH1 *h1)
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
        
        TSystemDirectory dir("dir", "input/track_qa/eta_3_0/");
        TList *files = dir.GetListOfFiles();

        if (files) 
            {
                TSystemFile *file;
                TString file_name;
                TIter next(files);
                TFile *out_file = new TFile("./output/output_projections.root", "RECREATE"); // Create the output file outside the loop
                while ((file=(TSystemFile*)next())) 
                {
                    file_name = file->GetName();
                    if (!file->IsDirectory() && file_name.EndsWith(".root")) 
                    {
                        TString filePath = "input/track_qa/eta_3_0/" + file_name;
                        TFile *in_file = TFile::Open(filePath); // Open the file
                        
                        TDirectory *subdir_pid = out_file->mkdir(file_name.ReplaceAll(".root", "")); // Create a subdirectory for this file
                        subdir_pid->cd(); // Change to the subdirectory

                        TH2F *h2_ecal_hcal_cluster_energy = (TH2F*)in_file->Get("h2_ecal_hcal_cluster_energy"); // replace with your histogram name
                        TH2F *h2_ecal_hcal_hit_energy = (TH2F*)in_file->Get("h2_ecal_hcal_hit_energy"); // replace with your histogram name
                        TH1D *h1_fwd_pt_res = (TH1D*)in_file->Get("h1_fwd_pt_res");
                        cout << "Current subdirectory in out_file: " << subdir_pid->GetName() << endl;
                        projections_and_fits(subdir_pid, h2_ecal_hcal_cluster_energy, h1_fwd_pt_res, "cluster");
                        subdir_pid->GetMotherDir()->cd(); // Change back to the parent directory
                        cout << "Current directory in out_file: " << subdir_pid->GetMotherDir()->GetName() << endl;
                        projections_and_fits(subdir_pid, h2_ecal_hcal_hit_energy, h1_fwd_pt_res, "hit");
                        subdir_pid->GetMotherDir()->cd(); // Change back to the parent directory
                        cout << "Current directory in out_file: " << subdir_pid->GetMotherDir()->GetName() << endl;
                        cout << "Current subdirectory in out_file: " << subdir_pid->GetName() << endl;
                        projections_and_fits(subdir_pid, h2_ecal_hcal_cluster_energy, h1_fwd_pt_res, "track_pt");
                        
                        in_file->Close();
                        
                    }
                    
                }
                out_file->Close(); // Close the output file after the loop
            }
    }

void projections_and_fits(TDirectory *dir, TH2F *h2, TH1D *h1_fwd_track, TString cluster_or_hit)
    { 
        PeakAna *pA = new PeakAna();
        //pA->SetDebug(2);
        TDirectory *subdir;
        if (cluster_or_hit == "cluster")
            {
                subdir = dir->mkdir("Clusters");
            }
        else if (cluster_or_hit == "hit")
            {
                subdir = dir->mkdir("Hits");
            }
        else if (cluster_or_hit == "track_pt")
            {
                subdir = dir->mkdir("Fwd_track_pt");
                subdir->cd();
                PeakAna *pA_pt = new PeakAna();

                pA_pt->SetData(h1_fwd_track);
                // pA_pt->SetBaseline(0,.75);
                pA_pt->SetRange(-1, 0, 1, 1000000);
                pA_pt->SetFilter(2,5);
                
                //pA->SetTunnelThreshold(0.1);
                pA_pt->SetTunnelScale(0.1);
                pA_pt->SetTunnelSigma(1);
                pA_pt->SetBaseline(10000,5);
                pA_pt->AnalyzeForPeak();
                pA_pt->Print();
                if (pA_pt->NPeaks() == 0)
                    {
                        cout << "No peaks found in the forward track pT histogram" << endl;
                        return;
                    }
                else
                    {
                        cout << "Peaks found in the forward track pT histogram" << endl;
                        
                        for(int i_peak = 0; i_peak < pA_pt->NPeaks(); i_peak++)
                            {
                                if(pA_pt->PeakProb(pA_pt->GetPeak(i_peak), pA_pt->TunnelScale(), pA_pt->TunnelSigma()) < 1e-10)
                                    {
                                        cout << "Preparing for fit of peak "  << 1e-10 << endl;
                                        FitAndEnergy *pt_fit_param = new FitAndEnergy();
                                        pt_fit_param->skew_gauss_dummy_to_full_fit(pA_pt->GetPeak(i_peak).mPeakY, pA_pt->GetPeak(i_peak).mPeakX, (pA_pt->GetPeak(i_peak).mEndX - pA_pt->GetPeak(i_peak).mStartX) / 2, h1_fwd_track);
                                        // pt_fit_param->skew_gauss_dummy_to_full_fit(h1_fwd_track->GetBinContent(h1_fwd_track->GetMaximumBin()), h1_fwd_track->GetBinCenter(h1_fwd_track->GetMaximumBin()), .15, h1_fwd_track);
                                        cout << "Peak guess: " << h1_fwd_track->GetBinCenter(h1_fwd_track->GetMaximumBin()) << " with amplitude: " << h1_fwd_track->GetBinContent(h1_fwd_track->GetMaximumBin()) << endl;
                                        TF1 *peak_fit = new TF1(Form("peak_fit_%d", i_peak), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1_fwd_track->GetXaxis()->GetXmin(), h1_fwd_track->GetXaxis()->GetXmax());
                                        peak_fit->FixParameter(0, pt_fit_param->GetAmplitudeFit());
                                        peak_fit->FixParameter(1, pt_fit_param->GetMeanFit());
                                        peak_fit->FixParameter(2, pt_fit_param->GetStdDevFit());
                                        peak_fit->FixParameter(3, pt_fit_param->GetSkewnessFit());
                                        peak_fit->SetLineColor(kBlue);
                                        peak_fit->SetNpx(1000);
                                        h1_fwd_track->Fit(peak_fit, "RQ+");
                                        cout << "Peak " << i_peak << " with window: " << pA_pt->GetPeak(i_peak).mStartX << " to " << pA_pt->GetPeak(i_peak).mEndX << endl;
                                        cout << "Peak " << i_peak << " with amplitude: " << pA_pt->GetPeak(i_peak).mPeakY << endl;
                                    }
                            }
                    }
                TCanvas *c_fwd_track = new TCanvas("Fwd_Track_pt", "Fwd_Track_pt", 800, 600);
                h1_fwd_track->SetTitle("Forward Track pT");
                h1_fwd_track->GetXaxis()->SetTitle("pT (GeV)");
                h1_fwd_track->GetYaxis()->SetTitle("Counts");
                h1_fwd_track->SetStats(0);
                h1_fwd_track->Draw();

                TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9); // coordinates are in the pad's fraction
                leg->SetTextSize(0.05); // Set text size to 4% of the pad height
                leg->AddEntry("Skewed gaussian fit parameters");
                leg->AddEntry((TObject*)0, Form("Amplitude (A): %.3f", peak_fit->GetParameter(0)), "");
                leg->AddEntry((TObject*)0, Form("Mean (mu): %.3f", peak_fit->GetParameter(1)), "");
                leg->AddEntry((TObject*)0, Form("Mean fit: %.3f", peak_fit->Mean(peak_fit->GetXmin(), peak_fit->GetXmax())), "");
                leg->AddEntry((TObject*)0, Form("Standard Deviation (sigma): %.3f", peak_fit->GetParameter(2)), "");
                leg->AddEntry((TObject*)0, Form("Skewness (alpha): %.3f", peak_fit->GetParameter(3)), "");
                leg->AddEntry((TObject*)0, Form("Chi^2/NDF fit range: %.3f", pt_fit_param->GetChi2NdfFitRange()), "");
                leg->Draw();

                c_fwd_track->Write();
                delete c_fwd_track;
                delete pA_pt;
                delete leg;
                delete pt_fit_param;
                dir->cd();
                return;
            }
        subdir->cd();
        cout << "Current subdirectory in out_file: " << subdir->GetName() << endl;
        // TDirectory *subdir_pelim_cal = subdir->mkdir("Preliminary_Summed_Projections_Ecal_on_Hcal"); // Change type to TDirectory*
        TDirectory *subdir_peak_find_cal = subdir->mkdir("Peak_Finder_Summed_Projections_Ecal_on_Hcal"); // Change type to TDirectory*

        // Calorimeter response
        int nbinsX = h2->GetXaxis()->GetNbins();
        TClonesArray *fits_and_energy = new TClonesArray("FitAndEnergy");

        for (int i = 1; i <= nbinsX; ++i) 
            {
                //subdir_pelim_cal->cd();
                // Create a 1D histogram of the y-axis projection for this bin
                TH1D *h1 = h2->ProjectionY(Form("h1_%d", i), i, i);
                if (h1->GetEntries() == 0)
                    {
                        continue;
                    }

                double ecal_energy = h2->GetXaxis()->GetBinCenter(i);
                
                subdir_peak_find_cal->cd();

                TCanvas *c1_alt = new TCanvas(Form("Peak_Find_ECAL_Energy_%i", (int)ecal_energy), Form("Peak_Find_ECAL_Energy:_%i", (int)ecal_energy), 800, 600);
                TH1D *h1_clean = new TH1D(Form("h1_clean_%d", i), "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                for (int bin = 1; bin <= h1->GetNbinsX(); ++bin)
                    {
                        h1_clean->SetBinContent(bin, h1->GetBinContent(bin));
                    }
                pA->SetData(h1_clean);
                pA->SetBaseline(0,.75);
                pA->SetRange(0, 0, 150, 1000000);
                pA->SetFilter(2,5);
                //pA->SetTunnelThreshold(0.1);
                pA->SetTunnelScale(0.1);
                pA->SetTunnelSigma(10);
                pA->AnalyzeForPeak();
                double min_not_a_peak_prob = 1;
                
                
                for (int i_peak = 0; i_peak < pA->NPeaks(); i_peak++)
                    {
                        if (pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) < min_not_a_peak_prob)
                            {
                                // cout << cyan << "New min not-a-peak probability: " << pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) << reset << endl;
                                min_not_a_peak_prob = pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma());
                            }
                    }
                //cout << "Max Not-A-Peak Probability: " << max_not_a_peak_prob << endl;
                // cout << "Fitting peaks for bin " << i << " with ecal energy " << ecal_energy << " GeV" << endl;
                TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9); // coordinates are in the pad's fraction
                leg->SetTextSize(0.05); // Set text size to 4% of the pad height
                TGraph *mipMarker = new TGraph();
                mipMarker->SetMarkerStyle(20); // Set marker style
                mipMarker->SetMarkerColor(kRed); // Set marker color

                TGraph *peakMarker = new TGraph();
                peakMarker->SetMarkerStyle(21); // Set marker style
                peakMarker->SetMarkerColor(kGreen); // Set marker color

                for(int i_peak = 0; i_peak < pA->NPeaks() && h1_clean->GetEntries() > 100; i_peak++) 
                    {
                        // leg->AddEntry("Skewed gaussian fit parameters");
                        // cout << "Fitting peaks for energy: " << ecal_energy << " GeV for peak " << i_peak << endl;
                        // cout << red << "Peak " << i_peak << " with window: " << pA->GetPeak(i_peak).mStartX << " to " << pA->GetPeak(i_peak).mEndX << reset << endl;
                        if (i_peak == 0 && (i_peak < (pA->NPeaks() - 1)))
                            {
                                FitAndEnergy *mip_fit_param = new FitAndEnergy();
                                // if(ecal_energy == 0.5)
                                //     {
                                //         cout << red << "MIP found. Min not-a-peak probability: " << min_not_a_peak_prob << reset << endl;
                                //         // cout << pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) << endl;
                                //     }
                                // cout << "MIP found at " << pA->GetPeak(i_peak).mPeakX << " with amplitude " << pA->GetPeak(i_peak).mPeakY << endl;
                                mip_fit_param->skew_gauss_dummy_to_full_fit(h1_clean->GetBinContent(1), 1, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                                // cout << red << "Error in dummy fit: " << reset << endl;
                                TF1 *peak_fit = new TF1(Form("mip_fit_%d", i_peak), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1_clean->GetXaxis()->GetXmin(), h1_clean->GetXaxis()->GetXmax());
                                peak_fit->FixParameter(0, mip_fit_param->GetAmplitudeFit());
                                peak_fit->FixParameter(1, mip_fit_param->GetMeanFit());
                                peak_fit->FixParameter(2, mip_fit_param->GetStdDevFit());
                                peak_fit->FixParameter(3, mip_fit_param->GetSkewnessFit());
                                peak_fit->SetLineColor(kRed);
                                // TLegendEntry *entry = leg->AddEntry((TObject*)0, Form("MIP Skew Gaussian fit parameters"), "");
                                // entry->SetMarkerStyle(20);
                                // entry->SetMakerColor(kRed);
                                leg->AddEntry(mipMarker, "MIP Skew Gaussian fit parameters", "p");
                                leg->AddEntry((TObject*)0, Form("Amplitude (A): %.3f", peak_fit->GetParameter(0)), "");
                                leg->AddEntry((TObject*)0, Form("Mean: %.3f", peak_fit->GetParameter(1)), "");
                                leg->AddEntry((TObject*)0, Form("StdDev: %.3f", peak_fit->GetParameter(2)), "");
                                leg->AddEntry((TObject*)0, Form("Skewness: %.3f", peak_fit->GetParameter(3)), "");
                                delete mip_fit_param;
                            }
                        else if (i_peak > 0 && (i_peak < (pA->NPeaks() - 1)) && pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) != min_not_a_peak_prob)
                            {
                                FitAndEnergy *peak_fit_param = new FitAndEnergy();
                                // cout << "Secondary peak found. Attempting to fit peak " << i_peak << " with skewed gaussian distribution" << endl;
                                peak_fit_param->skew_gauss_dummy_to_full_fit(pA->GetPeak(i_peak).mPeakY, pA->GetPeak(i_peak).mPeakX, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                                TF1 *peak_fit = new TF1(Form("secondary_peak_%d", i_peak), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1_clean->GetXaxis()->GetXmin(), h1_clean->GetXaxis()->GetXmax());
                                peak_fit->FixParameter(0, peak_fit_param->GetAmplitudeFit());
                                peak_fit->FixParameter(1, peak_fit_param->GetMeanFit());
                                peak_fit->FixParameter(2, peak_fit_param->GetStdDevFit());
                                peak_fit->FixParameter(3, peak_fit_param->GetSkewnessFit());
                                peak_fit->SetLineColor(kBlue);
                                delete peak_fit_param;
                            }
                        else
                            {
                                FitAndEnergy *peak_fit_param = (FitAndEnergy*)fits_and_energy->ConstructedAt(i);
                                // cout << "Attempting to fit peak " << i_peak << " with skewed gaussian distribution" << endl;
                                peak_fit_param->skew_gauss_dummy_to_full_fit(pA->GetPeak(i_peak).mPeakY, pA->GetPeak(i_peak).mPeakX, (pA->GetPeak(i_peak).mEndX - pA->GetPeak(i_peak).mStartX) / 8, h1_clean);
                                peak_fit_param->SetEcalEnergy(ecal_energy);
                                TF1 *peak_fit = new TF1(Form("peak_fit_%d", i_peak), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1_clean->GetXaxis()->GetXmin(), h1_clean->GetXaxis()->GetXmax());
                                peak_fit->FixParameter(0, peak_fit_param->GetAmplitudeFit());
                                peak_fit->FixParameter(1, peak_fit_param->GetMeanFit());
                                peak_fit->FixParameter(2, peak_fit_param->GetStdDevFit());
                                peak_fit->FixParameter(3, peak_fit_param->GetSkewnessFit());
                                peak_fit->SetLineColor(kGreen);
                                leg->AddEntry(peakMarker, Form("Peak %d Skewed gaussian fit parameters", i_peak), "p");
                                leg->AddEntry((TObject*)0, Form("Amplitude (A): %.3f", peak_fit->GetParameter(0)), "");
                                leg->AddEntry((TObject*)0, Form("Mean (mu): %.3f", peak_fit->GetParameter(1)), "");
                                leg->AddEntry((TObject*)0, Form("Mean fit: %.3f", peak_fit->Mean(peak_fit->GetXmin(), peak_fit->GetXmax())), "");
                                leg->AddEntry((TObject*)0, Form("Standard Deviation (sigma): %.3f", peak_fit->GetParameter(2)), "");
                                leg->AddEntry((TObject*)0, Form("Skewness (alpha): %.3f", peak_fit->GetParameter(3)), "");
                                leg->AddEntry((TObject*)0, Form("Chi^2/NDF fit range: %.3f", peak_fit_param->GetChi2NdfFitRange()), "");
                                leg->AddEntry((TObject*)0, Form("Not-A-Peak Probability: %.3e", pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma())), "");
                            }
                        //(FitAndEnergy*)fits_and_energy->ConstructedAt(i-1) = peak_fit_param;
                        
                        peak_fit->SetNpx(1000);
                        pA->PeakTunnel(pA->GetPeak(i_peak));
                        // if (pA->PeakTunnel(pA->GetPeak(i_peak)) == false)
                        //     {
                        //         //pA->Print();
                        //         //cout << "Not-A-Peak Probability: " << pA->PeakProb(pA->GetPeak(i_peak), pA->TunnelScale(), pA->TunnelSigma()) << endl;
                        //         peak_fit->SetLineColor(kRed);
                        //     }
                        h1_clean->Fit(peak_fit, "RQ+");
                        //delete peak_fit_param;
                    }
                
                h1_clean->SetStats(0);
                h1_clean->SetTitle(Form("ECAL Energy: %.2f (Gev)", ecal_energy));
                h1_clean->GetXaxis()->SetTitle("Hcal Summed Energy (GeV)");
                h1_clean->GetYaxis()->SetTitle("Counts");
                h1_clean->Draw();
                leg->Draw();
                // TPaveText *pt = new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC");
                // pA->AddPeakStats(pt, "ad");
                // pt->Draw();
                c1_alt->Write();
                dir->cd();
                delete mipMarker;
                delete peakMarker;
                delete c1_alt;
            }
        TF1 *f2 = new TF1("fit_2", "[0] + [1]*x", 0, 100);
        dir->cd();
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
        
        for (int j = 0; j <= nbinsX; ++j) 
            {
                FitAndEnergy *fit_and_energy = (FitAndEnergy*)fits_and_energy->At(j);
                // Set the point and its error in the TGraphErrors object
                if (fit_and_energy == 0 ||(fit_and_energy != 0 && fit_and_energy->GetMeanFit() == 0))
                    {
                        continue;
                    }
                // cout << "Setting point " << gr->GetN() << " with energy " << fit_and_energy->GetEcalEnergy() << " and mean fit " << fit_and_energy->GetMeanFit() << endl;
                gr->SetPoint(gr->GetN(), fit_and_energy->GetEcalEnergy(), fit_and_energy->GetMeanFit());
                double weight = 1.0 / sqrt(fit_and_energy->GetNEntries()); // calculate weight
                gr->SetPointError(gr->GetN()-1, 0, fit_and_energy->GetMeanFitError());
            }

        // Draw the TGraphErrors object on the canvas
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kViolet); // Change marker color to purple
        gr->SetLineColor(kViolet); // Change error bars color to purple
        gr->Draw("P same"); // "P" option for points, "same" to draw on the same canvas
        //cout << "Fitting the TGraphErrors object" << endl;
        gr->Fit(f2,"Q"); // "Q" option for quiet mode
        //cout << "Fitting complete" << endl;

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
        delete c2;
        // Close the ROOT files
        dir->cd();
            
    }