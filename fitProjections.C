//#include <iostream>
#include <vector>
using namespace std;
#include <cmath>

class FitAndEnergy: public TObject 
    {
    public:
        FitAndEnergy() {}
        double ecal_energy;
        double mean_fit;
        double mean_fit_error;
        double std_dev_fit;
        double nEntries;

    };


void fitProjections() 
    {
        TFile *file = TFile::Open("output/output_pi-e60.root"); // replace with your file name
        TFile *outfile = new TFile("./output/output_projections.root", "RECREATE");
        TH2F *h2_ecal_hcal_cluster_energy = (TH2F*)file->Get("h2_ecal_hcal_cluster_energy"); // replace with your histogram name
        TH2F *h2_ecal_hcal_hit_energy = (TH2F*)file->Get("h2_ecal_hcal_hit_energy"); // replace with your histogram name
        projections_and_fits(outfile, h2_ecal_hcal_cluster_energy, "cluster");
        projections_and_fits(outfile, h2_ecal_hcal_hit_energy, "hit");
        outfile->Close();
        file->Close();
    }
void projections_and_fits(TFile *outfile, TH2F *h2, TString cluster_or_hit)
    {
        TDirectory *dir; // Change type to TDirectory*
        
        // Load libraries
        // Load the ROOT file and get the histogram
        
        if (cluster_or_hit == "cluster")
            {
                dir = outfile->mkdir("Clusters");
                dir->cd();
            }
            else if (cluster_or_hit == "hit")
            {
                dir = outfile->mkdir("Hits");
                dir->cd();
            }
            
        TDirectory *subdir = dir->mkdir("Summed_Projections_Ecal_on_Hcal"); // Change type to TDirectory*
        
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
                subdir->cd();
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
                cout << "Estimated amplitude: " << amplitude_guess << endl;
                cout << "Estimated mean: " << mean_guess << endl;
                cout << "Estimated standard deviation: " << std_dev_guess << endl;
                cout << "Estimated skewness: " << skewness_guess << endl;
                // Fit the histogram with a Gaussian function
                TF1 *fit_dummy = new TF1(Form("fit1_%d", i), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", mean_guess - std_dev_guess, mean_guess + std_dev_guess);
                fit_dummy->SetParameters(amplitude_guess, mean_guess, std_dev_guess, skewness_guess);
                fit_dummy->SetParLimits(0, 0, amplitude_guess*5);
                fit_dummy->SetParLimits(1, h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                fit_dummy->SetParLimits(2, 0, std_dev_guess*5);
                fit_dummy->SetParLimits(3, skewness_guess-10, 0);
                fit_dummy->FixParameter(3, 0);
                h1->Fit(fit_dummy, "R");
                double amplitude_fit = fit_dummy->GetParameter(0);
                double mean_fit = fit_dummy->GetParameter(1);
                double mean_fit_error = fit_dummy->GetParError(1);
                double std_dev_fit = fit_dummy->GetParameter(2);
                double skewness_fit = fit_dummy->GetParameter(3);
                double ecal_energy = h2->GetXaxis()->GetBinCenter(i);
                double xmin = mean_fit - std_dev_fit;
                double xmax = mean_fit + std_dev_fit;
                FitAndEnergy *fit_and_energy = fits_and_energy->ConstructedAt(i-1);
                fit_and_energy->ecal_energy = ecal_energy;
                fit_and_energy->mean_fit = mean_fit;
                fit_and_energy->mean_fit_error = mean_fit_error;
                fit_and_energy->std_dev_fit = std_dev_fit;
                fit_and_energy->nEntries = h1->GetEntries();
                
                TF1 *fit1 = new TF1(Form("fit1_%d", i), "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]))", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                fit1->SetParameters(amplitude_fit, mean_fit, std_dev_fit, skewness_fit);
                fit1->FixParameter(0, amplitude_fit);
                fit1->FixParameter(1, mean_fit);
                fit1->FixParameter(2, std_dev_fit);
                fit1->FixParameter(3, skewness_fit);
                fit1->SetNpx(10000);
                
                // Create a canvas to draw the histogram and fit
                double ecal_energy = h2->GetXaxis()->GetBinCenter(i);
                //TCanvas *c1 = new TCanvas("c1","c1", 800, 600);
                h1->SetTitle(Form("Projected ECAL Energy on HCal: %.2f (Gev)", ecal_energy));
                TCanvas *c1 = new TCanvas(Form("ECAL_Energy:_%.2f", ecal_energy), Form("ECAL_Energy:_%.2f", ecal_energy), 800, 600);
                //gStyle->SetOptStat(0); 
                c1->SetTitle(Form("ECAL Energy: %.2f", ecal_energy));
                // Draw the histogram
                // Create a legend
                TH1F *h1_new = new TH1F(Form("h1_new_%d", i), "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
                for (int bin = 1; bin <= h1->GetNbinsX(); ++bin)
                    {
                        h1_new->SetBinContent(bin, h1->GetBinContent(bin));
                        //h1_new->SetBinError(bin, h1->GetBinError(bin));
                    }
                h1_new->SetLineColor(h1->GetLineColor());
                h1_new->SetLineStyle(h1->GetLineStyle());
                h1_new->SetLineWidth(h1->GetLineWidth());
                h1_new->SetFillColor(h1->GetFillColor());
                h1_new->SetFillStyle(h1->GetFillStyle());
                h1_new->Fit(fit1, "R"); //This just wastes computation time since the parameters are already fixed. Its only here for the graphical display
                h1_new->GetYaxis()->SetTitle("Counts");
                h1_new->GetXaxis()->SetTitle("Hcal Summed Energy (GeV)");
                h1_new->SetTitle(Form("Projected ECAL Energy on HCal: %.2f (Gev)", ecal_energy));
                //h1_new->GetListOfFunctions()->Add(fit1);
                double chi2_ndf_whole_range = fit1->GetNDF() != 0 ? fit1->GetChisquare()/fit1->GetNDF() : -1;
                double chi2_ndf_fit_range = fit_dummy->GetNDF() != 0 ? fit_dummy->GetChisquare()/fit_dummy->GetNDF() : -1;
                h1_new->SetStats(0);
                h1_new->Draw();
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
                //dir->cd();
                // Write the canvas to the ROOT file
                c1->Write();
                dir->cd();
                cout << i << endl;
                //delete c1;
            }
        //outfile->cd();
        TF1 *f2 = new TF1("fit2", "[0] + [1]*x", 0, 100);
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
                if (fit_and_energy == 0 ||(fit_and_energy != 0 && fit_and_energy->mean_fit == 0))
                    {
                        continue;
                    }

                gr->SetPoint(gr->GetN(), fit_and_energy->ecal_energy, fit_and_energy->mean_fit);
                double weight = 1.0 / sqrt(fit_and_energy->nEntries); // calculate weight
                gr->SetPointError(gr->GetN()-1, 0, fit_and_energy->mean_fit_error);
            }

        // Draw the TGraphErrors object on the canvas
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kViolet); // Change marker color to purple
        gr->SetLineColor(kViolet); // Change error bars color to purple
        gr->Draw("P same"); // "P" option for points, "same" to draw on the same canvas
        gr->Fit(f2); // "Q" option for quiet mode
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
        outfile->cd();
            
    }