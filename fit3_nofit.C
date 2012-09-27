#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <map>
#include <sys/stat.h>
#include <TDatime.h>

char MainDir[] = "/lustre/atlas/orlando";
Int_t Box = -1;
Long_t Directory = 20120313;
Int_t Run = 1;
Int_t temp = 5;
Double_t delta_oper_volt = 0.9;
Int_t min_volt = 0, max_volt = 12;
Long_t BAD = 0;
Int_t MIN_POS = 0;
Int_t MAX_POS = 31;
Double_t s_fit = 8.;


Int_t MPPC_map[32];
std::map<double, double> VOLT_map;
std::map<double, double> AMPF_map;

const Int_t N_points = max_volt-min_volt+1;
Int_t N_points_act = N_points;

Double_t *gain_list = new Double_t[N_points];
Double_t *lambda_list = new Double_t[N_points];
Double_t *volt_list = new Double_t[N_points];
Double_t *cross_talk_list = new Double_t[N_points];
Double_t *dark_rate_list = new Double_t[N_points];
Double_t *gain_list_err = new Double_t[N_points];
Double_t *volt_list_err = new Double_t[N_points];

Int_t npeaks = 30;
Int_t npeaks_ped = 30;

using namespace std;

void gaussianFitPed(int pos, int pin, int volt, TFile *file, Double_t DR_gain, Double_t DR_crosstalk)
{
    Double_t par[3000];
    TFile *ftemp = new TFile(Form("%s/USM%03d/T%02d_%d_run%03d/gaussian_fit.root",MainDir,Box,temp,Directory,Run),"READ");
    TH1F *h = (TH1F*)gROOT->FindObject(Form("hist_Pedest_pos%d_pin%d_volt%d",pos,pin,volt));

    TSpectrum *s = new TSpectrum(2*npeaks_ped);
    Int_t nfound = s->Search(h,5.,"",0.005);
    printf("Found %d candidate peaks to fit\n",nfound);
    // Estimate background using TSpectrum::Background
    TH1 *hb = s->Background(h,15,"goff");
    //estimate linear background using a fitting method
    TF1 *fline = new TF1("fline","pol1",0,1000);
    h->Fit("fline","qn");
    //Loop on all found peaks. Eliminate peaks at the background level
    par[0] = fline->GetParameter(0);
    par[1] = fline->GetParameter(1);
    npeaks_ped = 0;
    Float_t *xpeaks = s->GetPositionX();
    Double_t Mean[5];
    for (Int_t p=0; p<nfound; p++) {
        Float_t xp = xpeaks[p];
        Int_t bin = h->GetXaxis()->FindBin(xp);
        Float_t yp = h->GetBinContent(bin);
        if ( xp < 100. ) continue;
        if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
        Mean[npeaks_ped] = xp;
        npeaks_ped++;
    }
    printf("Found %d useful peaks to fit\n",npeaks_ped);

    Int_t N_mean = npeaks_ped;
    //Double_t Mean[5];
    Double_t Minimum;
    Int_t index[5];
    TMath::Sort(N_mean,Mean,index,0);

    Double_t N0 = h->Integral(1,Int_t(Mean[0]+DR_gain)-1);
    Double_t NR = h->Integral(Int_t(Mean[0]+DR_gain),4096);
    Double_t corr_factor = NR/(h->GetEntries());

    Double_t VAR = TMath::Power((h->GetRMS()/DR_gain),2);
    Double_t MEAN = ( (h->GetMean() - Mean[0])/DR_gain );
    Double_t corr_factor1 = MEAN / (1 + DR_crosstalk);
    Double_t corr_factor2 = MEAN * TMath::Exp(DR_crosstalk) / (1 + DR_crosstalk);
    Double_t corr_factor3 = VAR / (TMath::Power((1. + DR_crosstalk),2) + DR_crosstalk);
    Double_t corr_factor4 = VAR * TMath::Exp(DR_crosstalk) / (TMath::Power((1. + DR_crosstalk),2) + DR_crosstalk);
    Double_t lambda = (-(VAR-MEAN) + sqrt( (VAR-MEAN)*(VAR-MEAN) + 4*MEAN*MEAN) )/2;
    Double_t cross_talk = MEAN/lambda - 1;
    Double_t corr_factor5 = lambda;
    Double_t corr_factor6 = lambda * TMath::Exp(cross_talk);

    dark_rate_list[volt-min_volt] = corr_factor;

    ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    //if(DR_gain>0) out << "correction_factor1=" << corr_factor1 << " cf2=" << corr_factor2 <<  " cf3=" << corr_factor3 << " cf4=" << corr_factor4 << " cf5=" << corr_factor5 << " cf6=" << corr_factor6 << "\n";
    if(DR_gain>0) out << corr_factor << "\n";
    out.close();

    file->cd();
    h->Write(Form("hist_ped_pos%d_pin%d_volt%d",pos,pin,volt));
    delete ftemp;
    delete s;
    delete fline;
}

void gaussianFitLED(int pos, int pin, int volt, TFile *file, Double_t& DR_gain, Double_t& DR_crosstalk)
{
    Double_t par[3000];
    TFile *ftemp = new TFile(Form("%s/USM%03d/T%02d_%d_run%03d/gaussian_fit.root",MainDir,Box,temp,Directory,Run),"READ");
    TH1F *h = (TH1F*)gROOT->FindObject(Form("hist_wLight_pos%d_pin%d_volt%d",pos,pin,volt));

    //TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    //h->Draw();
    TSpectrum *s = new TSpectrum(2*npeaks);
    Int_t nfound = s->Search(h,s_fit,"",0.005);
    printf("Found %d candidate peaks to fit\n",nfound);
    // Estimate background using TSpectrum::Background
    TH1 *hb = s->Background(h,15,"goff");
    Float_t *xpeaks = s->GetPositionX();
    Int_t N_mean = nfound;
    Int_t index[50];
    TMath::Sort(N_mean,xpeaks,index,0);
    if(N_mean>=4) N_mean=4;
    printf("Found %d useful peaks to fit\n",npeaks);
    printf("Now fitting: Be patient\n");

    Double_t Mean[4],Minimum[4];
    for(Int_t k=0; k<N_mean; k++) Mean[k] = xpeaks[index[k]];
    Double_t gain=0;
    Double_t cross_talk = 0;
    if(N_mean==4 || N_mean==3) gain = Mean[2] - Mean[1];
    else if (N_mean==2) gain = Mean[1] - Mean[0];
    Double_t VAR = TMath::Power((h->GetRMS()/gain),2);
    Double_t MEAN = ( (h->GetMean() - Mean[0])/gain );
    Double_t lambda = (-(VAR-MEAN) + sqrt( (VAR-MEAN)*(VAR-MEAN) + 4*MEAN*MEAN) )/2;
    cross_talk = MEAN/lambda - 1;
    DR_gain = gain;
    DR_crosstalk = cross_talk;
    cout << "gain=" << gain << endl;
    cout << "lambda=" << lambda << endl;
    cout << "cross_talk=" << cross_talk << endl;

    Double_t real_volt = VOLT_map[100*pos + volt];

    ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    if(gain>0) out << pin << "\t" << real_volt << "\t" << gain << "\t" << cross_talk << "\t" << lambda << "\t";
    out.close();

    volt_list[volt-min_volt] = real_volt;
    gain_list[volt-min_volt] = gain;
    lambda_list[volt-min_volt] = lambda;
    cross_talk_list[volt-min_volt] = cross_talk;
    volt_list_err[volt-min_volt] = 0.;
    gain_list_err[volt-min_volt] = 0.;
    //gain_list_err[volt-12] = TMath::Sqrt(TMath::Power(fit->GetParameter(3+3*index[2]),2)+TMath::Power(fit->GetParameter(3+3*index[1]),2));

    h->SetTitle(Form("gain=%.2f cross_talk=%.2f",gain,cross_talk));
    file->cd();
    h->Write(Form("hist_wLight_pos%d_pin%d_volt%d",pos,pin,volt));

    ftemp->Close();
    delete ftemp;
    ftemp = 0L;
    delete s;
}

void draw_gain_volt(int pos, int pin, Double_t& BreakDownVoltage) {

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();

    Int_t low_limit_ind; 
    for(Int_t i=0; i<N_points_act; i++) {
      if(gain_list[i]<1) break;
      low_limit_ind=i;
    }

    TGraphErrors *gr1 = new TGraphErrors(N_points_act, volt_list, gain_list, volt_list_err, gain_list_err);
    gr1->SetLineColor(38);
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(1.1);
    gr1->Fit("pol1","","",volt_list[low_limit_ind-1],volt_list[1]);
    TF1 *myfunc = gr1->GetFunction("pol1");
    if(myfunc->GetParameter(1)>0) BreakDownVoltage = (-1) * myfunc->GetParameter(0) / myfunc->GetParameter(1);
    else BreakDownVoltage=0;
    cout<<"Break Down Voltage = " << BreakDownVoltage << endl;
    gr1->SetTitle(Form("Break Down Voltage = %.2f",BreakDownVoltage));
    //gr1->SetMinimum(0);
    //gr1->SetMaximum(130);
    gr1->Draw("ACP");

    c1->SaveAs(Form("%s/USM%03d/T%02d_%d_run%03d/gain_vs_volt_pos%d_pin%d_%d.gif",MainDir,Box,temp,Directory,Run,pos,pin,MPPC_map[pos]));

    ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    Double_t gain_oper_chan = (myfunc->GetParameter(1)) * (BreakDownVoltage + delta_oper_volt) + (myfunc->GetParameter(0));
    Double_t gain_oper = gain_oper_chan * 10. / 1.6 / AMPF_map[100*pos+pin];
    Int_t oper_volt_index = 0;
    while (BreakDownVoltage + delta_oper_volt < volt_list[oper_volt_index]) oper_volt_index++;
    Double_t cross_talk_oper = (cross_talk_list[oper_volt_index] + cross_talk_list[oper_volt_index-1]) / 2.;
    Double_t dark_rate_oper = (dark_rate_list[oper_volt_index] + dark_rate_list[oper_volt_index-1]) / 2.;
    out << pin << "\t" << BreakDownVoltage + delta_oper_volt << "\t" << gain_oper << "\t" << gain_oper_chan << "\t" << cross_talk_oper << "\t";// << dark_rate_oper <<"\t";
    out.close();

    ofstream out1(Form("%s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    out1 << pin << "\t" << BreakDownVoltage + delta_oper_volt << "\t" << gain_oper << "\t" << gain_oper_chan << "\t" << cross_talk_oper << "\t";// << dark_rate_oper <<"\t";
    out1.close();

    delete gr1;
    delete c1;
}

void draw_lambda_volt(int pos, int pin, Double_t BreakDownVoltage) {

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();

    Int_t low_limit_ind; 
    for(Int_t i=0; i<N_points_act; i++) {
      if(gain_list[i]<1) break;
      low_limit_ind=i;
    }

    TGraphErrors *gr1 = new TGraphErrors(N_points_act, volt_list, lambda_list, volt_list_err, gain_list_err);
    gr1->SetLineColor(38);
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(1.1);
    //gr1->Fit("pol2","","",volt_list[0],volt_list[N_points_act]);
    gr1->Fit("pol1","","",volt_list[low_limit_ind-1], volt_list[1]);
    //TF1 *myfunc = gr1->GetFunction("pol2");
    TF1 *myfunc = gr1->GetFunction("pol1");
    //Double_t Num_Phe =  myfunc->GetParameter(2) * (BreakDownVoltage + delta_oper_volt) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(1) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(0);
    Double_t Num_Phe = myfunc->GetParameter(1) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(0);
    Double_t ChiSqr = myfunc->GetChisquare() / myfunc->GetNDF();
    cout<<"Number of photoelectrons = " << Num_Phe << endl;
    gr1->SetTitle(Form("#chi^{2}=%.2f Num_Phe = %.2f",ChiSqr,Num_Phe));
    gr1->SetMinimum(0);
    gr1->SetMaximum(3.8);
    gr1->Draw("ACP");
    c1->SaveAs(Form("%s/USM%03d/T%02d_%d_run%03d/lambda_vs_volt_pos%d_pin%d_%d.gif",MainDir,Box,temp,Directory,Run,pos,pin,MPPC_map[pos]));

    ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    out << Num_Phe << "\t";
    out.close();

    ofstream out1(Form("%s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    out1 << Num_Phe << "\t";
    out1.close();

    delete gr1;
    delete c1;
}

void draw_DR_volt(int pos, int pin, Double_t BreakDownVoltage) {

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();

    Int_t low_limit_ind; 
    for(Int_t i=0; i<N_points_act; i++) {
      if(gain_list[i]<1) break;
      low_limit_ind=i;
    }

    TGraphErrors *gr1 = new TGraphErrors(N_points_act, volt_list, dark_rate_list, volt_list_err, gain_list_err);
    gr1->SetLineColor(38);
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(1.1);
    //gr1->Fit("pol2","","",volt_list[0],volt_list[N_points_act]);
    gr1->Fit("pol1","","",volt_list[low_limit_ind-3], volt_list[1]);
    //TF1 *myfunc = gr1->GetFunction("pol2");
    TF1 *myfunc = gr1->GetFunction("pol1");
    //Double_t Num_Phe =  myfunc->GetParameter(2) * (BreakDownVoltage + delta_oper_volt) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(1) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(0);
    Double_t Dark_Rate = myfunc->GetParameter(1) * (BreakDownVoltage + delta_oper_volt) + myfunc->GetParameter(0);
    Double_t ChiSqr = myfunc->GetChisquare() / myfunc->GetNDF();
    cout<<"Dark Rate = " << Dark_Rate << endl;
    gr1->SetTitle(Form("#chi^{2}=%.2f Dark_Rate = %.2f",ChiSqr,Dark_Rate));
    gr1->SetMinimum(0);
    gr1->SetMaximum(0.1);
    gr1->Draw("ACP");
    c1->SaveAs(Form("%s/USM%03d/T%02d_%d_run%03d/DR_vs_volt_pos%d_pin%d_%d.gif",MainDir,Box,temp,Directory,Run,pos,pin,MPPC_map[pos]));

    ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    out << Dark_Rate << "\n";
    out.close();

    ofstream out1(Form("%s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,pos,MPPC_map[pos]),ios::app);
    out1 << Dark_Rate << "\n";
    out1.close();

    delete gr1;
    delete c1;
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

void read_volt() { // 100*mppc_pos + volt_index -> ket_volt
    //char first_line[100];
    Int_t temp, volt_index, mppc_pos;
    Double_t set_volt, caen_volt, caen_curr, ket_volt, ket_rms;
    ifstream in;
    //in.open(Form("/lustre/atlas/data/%d/QDC/run%03d_OperationalData.dat",Directory,Run));
    in.open(Form("/lustre/atlas/data/%d/etc/run%03d_OperationalData.dat",Directory,Run));
    //in.getline(first_line,100);
    while(in >> temp >> volt_index >> mppc_pos >> set_volt >> caen_volt >> caen_curr >> ket_volt >> ket_rms) {
      VOLT_map[100*mppc_pos + volt_index] = ket_volt;
     }
    cout<<"volt map is created"<<endl;
}

void read_amplifier() { // 100*mppc_pos + mppc_cell -> amp_gain
    char first_line[100];
    Int_t mppc_cell, mppc_pos, qdc_channel, addr;
    Double_t pedes, count, output_Q, amp_gain;
    ifstream in;
    in.open(Form("/lustre/atlas/data/%d/etc/AmplifiersGain.csv",Directory));
    in.getline(first_line,100);
    while(in >> addr >> mppc_cell >> mppc_pos >> qdc_channel >> pedes >> count >> output_Q >> amp_gain) {
        AMPF_map[100*mppc_pos + mppc_cell] = amp_gain;
    }
    cout<<"amplifier map is created"<<endl;
}

bool checkDirectory()
{
   struct stat st,st2;
  if (stat(Form("%s/USM%03d",MainDir,Box),&st) == 0){
		if (stat(Form("%s/USM%03d/T%02d_%d_run%03d",MainDir,Box,temp,Directory,Run),&st2) == 0){
			cout << "Directory is already created" <<endl;
			return false;
	}
	else{
		//mkdir(Form("/lustre/atlas/jgpavez/%d/T%02d_run%03d",temp,Directory),0777);
		cout<<Form("%s/USM%03d/T%02d_%d_run%03d",MainDir,Box,temp,Directory,Run)<<endl;
		cout<<"You should run read_data.C first"<<endl;
  }
  }
	else{
	//mkdir(Form("/lustre/atlas/jgpavez/%d",Directory),0777);
	//mkdir(Form("/lustre/atlas/jgpavez/%d/T%02d_run%03d",temp,Directory),0777);
	cout<<Form("%s/%03d",MainDir,Box)<<endl;
	cout<<"You should run read_data.C first"<<endl;
	return true;
  }


}

bool readParameters()
{
char	header[50];
int bs;
    ifstream parameters;
    parameters.open("parameters.txt");
    if ( parameters.is_open()){
	parameters>>header>>Box;
	parameters>>header>>Directory;
	parameters>>header>>Run;
	parameters>>header>>temp;
	parameters>>header>>MainDir;
	parameters>>header>>delta_oper_volt;
	parameters>>header>>min_volt;
	parameters>>header>>max_volt;
	parameters>>header>>bs;
	while (bs>-1){
		BAD |= ((BAD>>bs)|1)<<bs;
		parameters>>bs;	
	}
	parameters>>header>>MIN_POS;
	parameters>>header>>MAX_POS;
	parameters>>header>>s_fit;

	parameters.close();
	N_points_act = max_volt-min_volt + 1; 
    	return true;
    }else{
	cout<<"You need to set the parameters in parameters.txt"<<endl;
    }
    parameters.close();
    return false;
}

int main( int argc, char *argv[]) {
    TDatime start;
    if (!readParameters()) return 0;
    if (checkDirectory()) return 0;  
    //system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/*.gif",MainDir,Box,temp,Directory,Run));
    read_mppc_id();
    read_volt();
    read_amplifier();
    const int min_pos = MIN_POS, max_pos = MAX_POS , min_pin = 0, max_pin = 15;
    //for ( int i_pos = min_pos; i_pos <=  max_pos; i_pos++){
    for ( int i_pos = min_pos; i_pos <=  max_pos; i_pos++){ 
      if( (BAD>>i_pos)&1 ) continue;
      TFile * file = new TFile(Form("%s/USM%03d/T%02d_%d_run%03d/gaussian_fit_MPPC_%d_%d.root",MainDir,Box,temp,Directory,Run,i_pos,MPPC_map[i_pos]),"RECREATE");
      system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,i_pos,MPPC_map[i_pos]));
      ofstream out(Form("%s/USM%03d/T%02d_%d_run%03d/BreakDownVoltage_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,i_pos,MPPC_map[i_pos]),ios::app);
      out << "MPPC_ID: " << MPPC_map[i_pos] << " at position: " << i_pos << endl;
      out << "PIN" << "\t" << "VOP" << "\t" << "\t" << "gain_ch" << "\t" << "CT" << "\t" << "NPhe" <<"\t" << "DR" << endl;;
      out.close();
      system(Form("rm -f %s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,i_pos,MPPC_map[i_pos]));
      ofstream out1(Form("%s/USM%03d/T%02d_%d_run%03d/results_MPPC_%d_%d.txt",MainDir,Box,temp,Directory,Run,i_pos,MPPC_map[i_pos]),ios::app);
      out1 << "MPPC_ID: " << MPPC_map[i_pos] << " at position: " << i_pos << endl;
      out1 << "PIN" << "\t" << "VOP" << "\t" << "gain" << "\t" << "gain_ch" << "\t" << "CT" << "\t" << "NPhe" <<"\t" << "DR" << endl;;
      out1.close();
      for ( int i_pin = min_pin; i_pin <= max_pin ; i_pin++) {
	for ( int volt = min_volt; volt <= max_volt; volt++ ) {
	  cout<<"mppc_pos: "<<i_pos<<" mppc_pin: "<<i_pin<<" volt_index: "<<volt<<endl;
	  Double_t DR_gain, DR_crosstalk;
	  gaussianFitLED(i_pos, i_pin, volt, file, DR_gain, DR_crosstalk);
	  gaussianFitPed(i_pos, i_pin, volt, file, DR_gain, DR_crosstalk);
	}
	Double_t BreakDownVoltage;
	draw_gain_volt(i_pos, i_pin, BreakDownVoltage);
	draw_lambda_volt(i_pos, i_pin, BreakDownVoltage);
	draw_DR_volt(i_pos, i_pin, BreakDownVoltage);
      }  
      file->Close();
      delete file;
      file = 0L;
    }
    TDatime end;
  	double time_use=(end.Convert()-start.Convert());
		cout<<"elapsed time: "<<time_use<<endl;
    return 0;
}
