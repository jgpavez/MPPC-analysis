#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TROOT.h" 
#include "TGraph.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <iostream>
#include <fstream>
#include <map>
#include <sys/stat.h>
std::map<int, int> data_map;

Long_t Directory = 20120301;
Int_t Run = 1;
char MainDir[] = "/lustre/atlas/orlando";
Int_t temp = 5;
Int_t Box = -1;
const Int_t N_points = 16;
Int_t MPPC_map[32];
Double_t *cell_list = new Double_t[N_points];
Double_t *volt_list = new Double_t[N_points];
Double_t *gain_list = new Double_t[N_points];
Double_t *gain_chan_list = new Double_t[N_points];
Double_t *cross_talk_list = new Double_t[N_points];
Double_t *dark_rate_list = new Double_t[N_points];
Double_t *mean_phe_list = new Double_t[N_points];

using namespace std;

void plot( const Int_t position) {
  for(Int_t n=0; n<N_points; n++) {
    cell_list[n] = 0.;
    volt_list[n] = 0.;
    gain_list[n] = 0.;
    gain_chan_list[n] = 0.;
    cross_talk_list[n] = 0.;
    dark_rate_list[n] = 0.;
    mean_phe_list[n] = 0.;
  }

  Int_t mppc_cell;
  Double_t volt, gain, gain_chan, cross_talk, dark_rate, mean_phe;
  ifstream in;
  in.open(Form("%s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory, Run, position, MPPC_map[position]));
  if(!in.is_open()) {cout << "position " << position << " is not included!" << endl; return;}
  char first_line[100];
  in.getline(first_line,100);
  in.getline(first_line,100);
  Int_t i=0;
  while(in >> mppc_cell >> volt >> gain >> gain_chan >> cross_talk >> mean_phe >> dark_rate) {
    cell_list[i] = i+1;
    volt_list[i] = volt-70.;
    gain_list[i] = gain;
    gain_chan_list[i] = gain_chan/10.;
    cross_talk_list[i] = 10*cross_talk;
    dark_rate_list[i] = 10*dark_rate;
    mean_phe_list[i] = mean_phe;
    i++;
  }

  Double_t mean_VOP = 0.;
  Double_t mean_gain = 0.;
  Double_t mean_CT = 0.;
  Double_t mean_DR  = 0.;
  Double_t mean_NPh  = 0.;
  Double_t err_VOP = 0.;
  Double_t err_gain = 0.;
  Double_t err_CT = 0.;
  Double_t err_DR  = 0.;
  Double_t err_NPh  = 0.;
  for(Int_t k = 0; k < N_points; k++) {
    mean_VOP = mean_VOP + volt_list[k] + 70.;
    mean_gain = mean_gain + gain_list[k];
    mean_CT = mean_CT + cross_talk_list[k]/10.;
    mean_DR = mean_DR + dark_rate_list[k]/10;
    mean_NPh = mean_NPh + mean_phe_list[k];
  }
  mean_VOP = mean_VOP/N_points;
  mean_gain = mean_gain/N_points;
  mean_CT = mean_CT/N_points;
  mean_DR = mean_DR/N_points;
  mean_NPh = mean_NPh/N_points;
  
  for(Int_t k = 0; k < N_points; k++) {
    err_VOP = err_VOP + TMath::Power((volt_list[k] + 70. - mean_VOP),2);
    err_gain = err_gain + TMath::Power((gain_list[k] - mean_gain),2);
    err_CT = err_CT + TMath::Power((cross_talk_list[k]/10. - mean_CT),2);
    err_DR = err_DR + TMath::Power((dark_rate_list[k]/10 - mean_DR),2);
    err_NPh = err_NPh + TMath::Power((mean_phe_list[k] - mean_NPh),2);
  }
  err_VOP = TMath::Sqrt(err_VOP/N_points);
  err_gain = TMath::Sqrt(err_gain/N_points);
  err_CT = TMath::Sqrt(err_CT/N_points);
  err_DR = TMath::Sqrt(err_DR/N_points);
  err_NPh = TMath::Sqrt(err_NPh/N_points);
  
  cout << "MPPC ID = " << MPPC_map[position] << " at position = " << position << endl;
  cout << "Average operational volatage = " << mean_VOP << "    " << err_VOP << endl;
  cout << "Average gain = " << mean_gain << "    " << err_gain << endl;
  cout << "Average cross talk = " << mean_CT << "    " << err_CT << endl;
  cout << "Average dark rate = " << mean_DR << "    " << err_DR << endl;
  cout << "Average NPhe = " << mean_NPh << "    " << err_NPh << endl;

  ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/report_T%02d.txt",MainDir,Box,temp,Directory,  Run, temp),ios::app);
  out << MPPC_map[position] << "\t" << position << "\t" << mean_VOP << "\t" << err_VOP << "\t" << mean_gain << "\t" << err_gain << "\t" << mean_CT << "\t" << err_CT << "\t" << mean_DR << "\t" << err_DR << "\t" << mean_NPh << "\t" << err_NPh << endl;;
  out.close();

  TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
  c1->SetFillColor(42);
  c1->SetGrid();

  TGraph *gr1 = new TGraph(N_points, cell_list, volt_list);
  gr1->SetMinimum(0);
  gr1->SetMaximum(7);
  gr1->SetLineColor(2);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(1.1);
  gr1->Draw("ACP");

  TGraph *gr2 = new TGraph(N_points, cell_list, gain_list);
  gr2->SetLineColor(5);
  gr2->SetMarkerColor(5);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.1);
  gr2->Draw("CP:same");

  TGraph *gr3 = new TGraph(N_points, cell_list, gain_chan_list);
  gr3->SetLineColor(4);
  gr3->SetMarkerColor(4);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(1.1);
  gr3->Draw("CP:same");

  TPaveText *pt = new TPaveText(0.00862069,0.9533898,0.03735632,0.9957627,"blNDC");
  pt->SetName("title");
  pt->SetBorderSize(2);
  TText *text = pt->AddText("Graph");
  pt->Draw();
  TLatex *tex = new TLatex(1.70977,7.193527,"Vop-70");
  tex->SetTextColor(2);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(5.324174,7.193527,"Gain(chan/10)");
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(11.17008,7.226619,"Gain(compxE5)");
  tex->SetTextColor(5);
  tex->SetLineWidth(2);
  tex->Draw();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  c1->ToggleToolBar();
  
  system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/gain_voltage_%d_%d.gif",MainDir, Box, temp , Directory, Run, position, MPPC_map[position]));
  c1->SaveAs(Form("%s/USM%03d/T%02d_%d_run%03d/gain_voltage_%d_%d.gif",MainDir, Box, temp , Directory, Run, position, MPPC_map[position]));
  
  TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",200,10,700,500);
  c2->SetFillColor(42);
  c2->SetGrid();

  TGraph *gr4 = new TGraph(N_points, cell_list, cross_talk_list);
  gr4->SetMinimum(0);
  gr4->SetMaximum(4);
  gr4->SetLineColor(3);
  gr4->SetMarkerColor(3);
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerSize(1.1);
  gr4->Draw("ACP");

  TGraph *gr5 = new TGraph(N_points, cell_list, dark_rate_list);
  gr5->SetLineColor(1);
  gr5->SetMarkerColor(1);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(1.1);
  gr5->Draw("CP:same");

  TGraph *gr6 = new TGraph(N_points, cell_list, mean_phe_list);
  gr6->SetLineColor(6);
  gr6->SetMarkerColor(6);
  gr6->SetMarkerStyle(21);
  gr6->SetMarkerSize(1.1);
  gr6->Draw("CP:same");

  pt = new TPaveText(0.01,0.9390678,0.121954,0.995,"blNDC");
  pt->SetName("title");
  pt->SetBorderSize(2);
  text = pt->AddText("Graph");
  pt->Draw();
  tex = new TLatex(1.67834,4.161017,"NPE");
  tex->SetTextColor(6);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(5.261315,4.161017,"(CrossTalk/10)%");
  tex->SetTextColor(3);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(11.54723,4.139831,"(DarkRate/10)%");
  tex->SetLineWidth(2);
  tex->Draw();
  c2->Modified();
  c2->cd();
  c2->SetSelected(c2);
  c2->ToggleToolBar();
  
  system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/phe_dr_ct_%d_%d.gif",MainDir, Box, temp , Directory,Run, position, MPPC_map[position]));
  c2->SaveAs(Form("%s/USM%03d/T%02d_%d_run%03d/phe_dr_ct_%d_%d.gif",MainDir, Box, temp , Directory,Run, position, MPPC_map[position]));
  
  c1->Close();
  c2->Close();
  return;
}

bool checkDirectory()
{
    struct stat st,st2;
    if (stat(Form("%s/USM%03d",MainDir,Box),&st) == 0){
        if (stat(Form("%s/USM%03d/T%02d_%d_run%03d",MainDir,Box,temp,Directory,Run),&st2) == 0){                   
		cout << "Directory es already created" <<endl;
                return false;
        }else{
                //mkdir(Form("/lustre/atlas/jgpavez/%d/T%02d_run%03d",Directory,temp),0777);
                cout<<"You should run read_data.C first"<<endl;
	}
    }else{
        //mkdir(Form("/lustre/atlas/jgpavez/%d",Directory),0777);
        //mkdir(Form("/lustre/atlas/jgpavez/%d/T%02d_run%03d",Directory,temp),0777);
        cout<<"You should run read_data.C first"<<endl;
        return true;
    }


}

bool readParameters()
{
		char header[50];
    ifstream parameters;
    parameters.open("../parameters.txt");
    if ( parameters.is_open()){
				parameters>>header>>Box;
        parameters>>header>>Directory;
        parameters>>header>>Run;
        parameters>>header>>temp;
	parameters>>header>>MainDir;
	parameters.close();
        return true;
    }else{
        cout<<"You need to set the parameters in parameters.txt"<<endl;
    }
    parameters.close();
    return false;
}

void read_mppc_id() {
    char first_line[100];
    Int_t slot, channel, position, mppc;
    Double_t biasVoltage;
    ifstream in;
    in.open(Form("/lustre/atlas/data/%d/etc/biasvoltage.run%03d.map",Directory,Run));
    in.getline(first_line,100);
    while(in >> slot >> channel >> position >> mppc >> biasVoltage) {
      if(position==0) if(mppc==0) continue;
      MPPC_map[position] = mppc;
      cout<< position << " " << mppc << endl;
    }
}


int main(int argc, char *argv[]) {
  if(!readParameters()) return 0;
  if(checkDirectory()) return 0;
  read_mppc_id();
  Int_t position;
  int min_position = 0, max_position = 31;
  system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/report_T%02d.txt",MainDir,Box,temp,Directory,  Run, temp));
  ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/report_T%02d.txt",MainDir,Box,temp,Directory,  Run, temp),ios::app);
  out << "MPPC_ID\tPOS\tVOP\tVOP_err\tgain\tgain_err\tCT\tCT_err\tDR\tDR_err\tNPhe_\tNPhe_err" <<endl;;
  out.close();
  for ( int i = min_position; i <= max_position; i++)
    plot(position);
  return 0;
}
