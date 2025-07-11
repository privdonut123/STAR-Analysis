#include <iostream>
#include <vector>
using namespace std;
#include <cmath>
#include <string>
#include <cctype>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFile.h>
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


void sortData(std::vector<double>& energy, std::vector<double>& mean, std::vector<double>& std_dev, vector<double>& cal_hit_mean, vector<double>& cal_hit_std_dev)
    {
        vector<int> indices(energy.size());
        for (int i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        
        for (int i = 0; i < energy.size(); ++i) {
            for (int j = i + 1; j < energy.size(); ++j) {
                if (energy[i] > energy[j]) {
                    swap(energy[i], energy[j]);
                    swap(mean[i], mean[j]);
                    swap(std_dev[i], std_dev[j]);
                    swap(cal_hit_mean[i], cal_hit_mean[j]);
                    swap(cal_hit_std_dev[i], cal_hit_std_dev[j]);
                    swap(indices[i], indices[j]);
                }
            }
        }
    }

void pt_qa(TString rapidity = "3.0")
    {
        cout << cyan << "Opening ./output/output_projections.root file" << reset << endl;
        TFile *in_file = TFile::Open("./output/output_projections.root"); // Open the file
        if (!in_file)
            {
                cout << red << "Error: Could not open ./output/output_projections.root file" << reset << endl;
                return;
            }
        
        //TH1F* pt_res = (TH1F*)in_file->Get("mu-e15/Fwd_track_pt/Fwd_track_pt");
        // Get the list of keys (subdirectories) in the ROOT file
        TIter next(in_file->GetListOfKeys());
        TKey *key;

        vector<double> pi_energy;
        vector<double> pi_pt_mean;
        vector<double> pi_pt_std_dev;
        vector<double> pi_cal_mean;
        vector<double> pi_cal_std_dev;

        vector<double> mu_energy;
        vector<double> mu_pt_mean;
        vector<double> mu_pt_std_dev;
        vector<double> mu_cal_mean;
        vector<double> mu_cal_std_dev;

        // Loop over each key (subdirectory)
        while ((key = (TKey*)next())) {
        TDirectory *subdir = (TDirectory*)key->ReadObj();
        if (subdir) {
            cout << "Subdirectory " << key->GetName() << " found" << endl;

            // Extract the number from the subdirectory name
            string subdir_name = key->GetName();
            string number_str;
            for (size_t i = 0; i < subdir_name.length(); ++i) {
                if (isdigit(subdir_name[i]) || subdir_name[i] == '.') {
                    number_str += subdir_name[i];
                }
            }

            if (!number_str.empty()) {
                stringstream ss(number_str);
                double subdir_number;
                ss >> subdir_number;
                // cout << red << "Subdirectory number: " << subdir_number << reset << endl;
                // Check the first three characters of the subdirectory name
                if (subdir_name.substr(0, 3) == "pi-") {
                    pi_energy.push_back(subdir_number);
                } else if (subdir_name.substr(0, 3) == "mu-") {
                    mu_energy.push_back(subdir_number);
                } else {
                    cout << "Error: Subdirectory name does not start with 'pi-' or 'mu-': " << subdir_name << endl;
                    continue;
                }
            } else {
                cout << "Error: Could not find a number in subdirectory name " << subdir_name << endl;
                continue;
            }

            Double_t pt_res_mean = 0;
            Double_t pt_res_std_dev = 0;

            TDirectoryFile *subdir_pt = subdir->GetDirectory("Resolution_qa");
            
            TCanvas *pt_canvas = (TCanvas*)subdir_pt->Get("Fwd_Track_pt");
            if (!pt_canvas) 
                {
                    cout << "Error: Could not find canvas in subdirectory " << key->GetName() << endl;
                    if (subdir_name.substr(0, 3) == "pi-") 
                        {
                            pi_pt_mean.push_back(pt_res_mean);
                            pi_pt_std_dev.push_back(pt_res_std_dev);
                        } 
                    else if (subdir_name.substr(0, 3) == "mu-") 
                        {
                            mu_pt_mean.push_back(pt_res_mean);
                            mu_pt_std_dev.push_back(pt_res_std_dev);
                        }
                    continue;
                }

            TH1F *hist = (TH1F*)pt_canvas->GetListOfPrimitives()->At(0);
            if (!hist) 
                {
                    cout << "Error: Could not find histogram in subdirectory " << key->GetName() << endl;
                    continue;
                }

            TF1 *fitFunc = (TF1*)hist->GetListOfFunctions()->At(0);
            if (!fitFunc) 
                {
                    cout << "Error: Could not find fit function in subdirectory " << key->GetName() << endl;
                    continue;
                }

            pt_res_mean = fitFunc->GetParameter(1);
            pt_res_std_dev = fitFunc->GetParameter(2);

            if (subdir_name.substr(0, 3) == "pi-") 
                {
                    pi_pt_mean.push_back(pt_res_mean);
                    pi_pt_std_dev.push_back(pt_res_std_dev);
                } 
            else if (subdir_name.substr(0, 3) == "mu-") 
                {
                    mu_pt_mean.push_back(pt_res_mean);
                    mu_pt_std_dev.push_back(pt_res_std_dev);
                }

            Double_t cal_res_mean = 0;
            Double_t cal_res_std_dev = 0;

            TCanvas *cal_hit_canvas = (TCanvas*)subdir_pt->Get("Fwd_Cal_hit");
            if (!cal_hit_canvas) 
                {
                    cout << "Error: Could not find canvas in subdirectory " << key->GetName() << endl;
                    if (subdir_name.substr(0, 3) == "pi-") 
                        {
                            pi_cal_mean.push_back(cal_res_mean);
                            pi_cal_std_dev.push_back(cal_res_std_dev);
                        } 
                    else if (subdir_name.substr(0, 3) == "mu-") 
                        {
                            mu_cal_mean.push_back(cal_res_mean);
                            mu_cal_std_dev.push_back(cal_res_std_dev);
                        }
                    continue;
                }

            TH1F *hist = (TH1F*)cal_hit_canvas->GetListOfPrimitives()->At(0);
            if (!hist) {
                cout << "Error: Could not find histogram in subdirectory " << key->GetName() << endl;
                continue;
            }

            TF1 *fitFunc = (TF1*)hist->GetListOfFunctions()->At(0);
            if (!fitFunc) {
                cout << "Error: Could not find fit function in subdirectory " << key->GetName() << endl;
                continue;
            }

            Double_t cal_res_mean = fitFunc->GetParameter(1);
            Double_t cal_res_std_dev = fitFunc->GetParameter(2);

            // Sort the mean and std_dev values into the appropriate vectors
            if (subdir_name.substr(0, 3) == "pi-") {
                pi_cal_mean.push_back(cal_res_mean);
                pi_cal_std_dev.push_back(cal_res_std_dev);
            } else if (subdir_name.substr(0, 3) == "mu-") {
                mu_cal_mean.push_back(cal_res_mean);
                mu_cal_std_dev.push_back(cal_res_std_dev);
            }
        } else {
            cout << "Subdirectory " << key->GetName() << " not found" << endl;
        }

        delete subdir;
    }
    
    sortData(pi_energy, pi_pt_mean, pi_pt_std_dev, pi_cal_mean, pi_cal_std_dev);
    sortData(mu_energy, mu_pt_mean, mu_pt_std_dev, mu_cal_mean, mu_cal_std_dev);

    TGraphErrors *pi_pt_graph = new TGraphErrors(pi_energy.size(), &pi_energy[0], &pi_pt_mean[0], 0, &pi_pt_std_dev[0]);
    TGraphErrors *mu_pt_graph = new TGraphErrors(mu_energy.size(), &mu_energy[0], &mu_pt_mean[0], 0, &mu_pt_std_dev[0]);
    TGraphErrors *pi_cal_graph = new TGraphErrors(pi_energy.size(), &pi_energy[0], &pi_cal_mean[0], 0, &pi_cal_std_dev[0]);
    TGraphErrors *mu_cal_graph = new TGraphErrors(mu_energy.size(), &mu_energy[0], &mu_cal_mean[0], 0, &mu_cal_std_dev[0]);

    TGraphErrors *pi_pt_std_dev_graph = new TGraphErrors(pi_energy.size(), &pi_energy[0], &pi_pt_std_dev[0]);
    TGraphErrors *mu_pt_std_dev_graph = new TGraphErrors(mu_energy.size(), &mu_energy[0], &mu_pt_std_dev[0]);
    TGraphErrors *pi_cal_std_dev_graph = new TGraphErrors(pi_energy.size(), &pi_energy[0], &pi_cal_std_dev[0]);
    TGraphErrors *mu_cal_std_dev_graph = new TGraphErrors(mu_energy.size(), &mu_energy[0], &mu_cal_std_dev[0]);

    TCanvas *c1 = new TCanvas("c1", "Graphs", 800, 600);

    // Draw the pi_pt_graph
    pi_pt_graph->SetMarkerStyle(20);
    pi_pt_graph->SetMarkerColor(kRed);
    pi_pt_graph->SetLineColor(kRed);
    pi_pt_graph->SetTitle("FTS p_{T} and FCAL Energy resolution, #eta = " + rapidity + ";Particle Energy; Resolution Mean");
    pi_pt_graph->GetYaxis()->SetRangeUser(-1.2, 1.2);
    pi_pt_graph->Draw("AP");

    // Draw the mu_pt_graph on the same canvas
    mu_pt_graph->SetMarkerStyle(21);
    mu_pt_graph->SetMarkerColor(kBlue);
    mu_pt_graph->SetLineColor(kBlue);
    mu_pt_graph->Draw("P SAME");

    // Draw the pi_cal_graph on the same canvas
    pi_cal_graph->SetMarkerStyle(22);
    pi_cal_graph->SetMarkerColor(kGreen);
    pi_cal_graph->SetLineColor(kGreen);
    pi_cal_graph->Draw("P SAME");

    // Draw the mu_cal_graph on the same canvas
    mu_cal_graph->SetMarkerStyle(23);
    mu_cal_graph->SetMarkerColor(kBlack);
    mu_cal_graph->SetLineColor(kBlack);
    mu_cal_graph->Draw("P SAME");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(pi_pt_graph, "#pi^{-} FTS p_{T}", "p");
    legend->AddEntry(mu_pt_graph, "#mu^{-} FTS p_{T}", "p");
    legend->AddEntry(pi_cal_graph, "#pi^{-} FCAL Energy", "p");
    legend->AddEntry(mu_cal_graph, "#mu^{-} FCAL Energy", "p");
    legend->Draw();

    // Create a ROOT file and save the canvas
    TFile *file = new TFile("output/pt_qa.root", "RECREATE");
    c1->Write();

    TCanvas *c2 = new TCanvas("c2", "Graphs", 800, 600);
    // Draw the pi_pt_std_dev_graph
    pi_pt_std_dev_graph->SetMarkerStyle(20);
    pi_pt_std_dev_graph->SetMarkerColor(kRed);
    pi_pt_std_dev_graph->SetLineColor(kRed);
    pi_pt_std_dev_graph->SetTitle("FTS p_{T} and FCAL Energy resolution, #eta = " + rapidity + ";Particle Energy; Resolution Std Dev");
    pi_pt_std_dev_graph->GetYaxis()->SetRangeUser(-10, 10);
    pi_pt_std_dev_graph->Draw("AP");
    // Draw the mu_pt_std_dev_graph on the same canvas
    mu_pt_std_dev_graph->SetMarkerStyle(21);
    mu_pt_std_dev_graph->SetMarkerColor(kBlue);
    mu_pt_std_dev_graph->SetLineColor(kBlue);
    mu_pt_std_dev_graph->Draw("P SAME");
    // Draw the pi_cal_std_dev_graph on the same canvas
    pi_cal_std_dev_graph->SetMarkerStyle(22);
    pi_cal_std_dev_graph->SetMarkerColor(kGreen);
    pi_cal_std_dev_graph->SetLineColor(kGreen);
    pi_cal_std_dev_graph->Draw("P SAME");
    // Draw the mu_cal_std_dev_graph on the same canvas
    mu_cal_std_dev_graph->SetMarkerStyle(23);
    mu_cal_std_dev_graph->SetMarkerColor(kBlack);
    mu_cal_std_dev_graph->SetLineColor(kBlack);
    mu_cal_std_dev_graph->Draw("P SAME");
    TLegend *legend2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend2->AddEntry(pi_pt_std_dev_graph, "#pi^{-} FTS p_{T}", "p");
    legend2->AddEntry(mu_pt_std_dev_graph, "#mu^{-} FTS p_{T}", "p");
    legend2->AddEntry(pi_cal_std_dev_graph, "#pi^{-} FCAL Energy", "p");
    legend2->AddEntry(mu_cal_std_dev_graph, "#mu^{-} FCAL Energy", "p");
    legend2->Draw();

    c2->Write();

    file->Close();

    // Clean up
    delete c1;
    delete c2;
    delete pi_pt_graph;
    delete mu_pt_graph;
    delete pi_cal_graph;
    delete mu_cal_graph;
    delete pi_pt_std_dev_graph;
    delete mu_pt_std_dev_graph;
    delete pi_cal_std_dev_graph;
    delete mu_cal_std_dev_graph;
    delete file;
    

    // cout << "Mean values for pi_energy: " << endl;
    // for (size_t j = 0; j < pi_energy.size(); j++) {
    //     cout << pi_energy[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for pi_pt_mean: " << endl;
    // for (size_t j = 0; j < pi_pt_mean.size(); j++) {
    //     cout << pi_pt_mean[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for pi_pt_std_dev: " << endl;
    // for (size_t j = 0; j < pi_pt_std_dev.size(); j++) {
    //     cout << pi_pt_std_dev[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for pi_cal_mean: " << endl;
    // for (size_t j = 0; j < pi_cal_mean.size(); j++) {
    //     cout << pi_cal_mean[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for pi_cal_std_dev: " << endl;
    // for (size_t j = 0; j < pi_cal_std_dev.size(); j++) {
    //     cout << pi_cal_std_dev[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for mu_energy: " << endl;
    // for (size_t j = 0; j < mu_energy.size(); j++) {
    //     cout << mu_energy[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for mu_pt_mean: " << endl;
    // for (size_t j = 0; j < mu_pt_mean.size(); j++) {
    //     cout << mu_pt_mean[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for mu_pt_std_dev: " << endl;
    // for (size_t j = 0; j < mu_pt_std_dev.size(); j++) {
    //     cout << mu_pt_std_dev[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for mu_cal_mean: " << endl;
    // for (size_t j = 0; j < mu_cal_mean.size(); j++) {
    //     cout << mu_cal_mean[j] << " ";
    // }
    // cout << endl;

    // cout << "Mean values for mu_cal_std_dev: " << endl;
    // for (size_t j = 0; j < mu_cal_std_dev.size(); j++) {
    //     cout << mu_cal_std_dev[j] << " ";
    // }   
    // cout << endl;

    }