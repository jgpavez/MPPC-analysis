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

using namespace std;

Int_t Directory = 20120312;
char MainDir[] = "/lustre/atlas/orlando";
Int_t Run = 8;
Int_t temp = 5;

std::map<int, int> QDC_map;
std::map<int, int> reverse_map;
std::map<int, int> HIST_map;
std::map<int, int> reverse_HIST_map;

const Int_t N_MPPC = 32;
const Int_t N_pins = 16;
const Int_t min_volt = 0;
const Int_t max_volt = 12;
Long_t Box;
TH1F *h[N_MPPC*N_pins];

void read_qdc() { // 10*qdc_ch + addr -> 100*mppc_pos + mppc_pin
    char first_line[100];
    int addr, qdc_ch, mppc_pos, mppc_pin, slot, BVch;
    ifstream in;
    in.open("/lustre/atlas/etc/qdcconfigNEW.map");
    in.getline(first_line,100);
    while(in >> addr >> qdc_ch >> mppc_pos >> mppc_pin >> slot >> BVch) {
        QDC_map[10*qdc_ch + addr] = 100*mppc_pos + mppc_pin;
        reverse_map[100*mppc_pos + mppc_pin] = 10*qdc_ch + addr;
    }
    cout<<"qdc map is created"<<endl;
    cout<<"reverse map is created"<<endl;
}

void hist_maping(TH1F *hist[N_MPPC*N_pins]) {//100*mppc + pin -> index
  int index=0;
  for(int mppc=0; mppc<N_MPPC; mppc++){
    for(int pin=0; pin<N_pins; pin++){
      HIST_map[100*mppc + pin] = index;
      reverse_HIST_map[index] = 100*mppc + pin;
      hist[index++]=new TH1F(Form("hist_mppc%d_cell%d",mppc,pin),"",4096,0,4096);
    }
  }
  cout<<"hist map is created"<<endl;
  cout<<"hist reverse map is created"<<endl;
}

void readFiles(int Run,int temp, bool light, int volt, TFile *file){
  
  TH1F *hist[N_MPPC*N_pins];
  hist_maping(hist);
  
  for (Int_t addr=0; addr<8; addr++){
    std::ifstream infile(Form("/lustre/atlas/data/%d/QDC/run%03d_T%02d_index%02d_%s-addr%d.dat",Directory,Run,temp,volt,(light==1)?"wLight":"Pedest",addr));
    Int_t event,channel,count;	
    while(infile>>event>>channel>>count){
      Int_t mpos_mcell = QDC_map[10*channel+addr];
      Int_t index = HIST_map[100*(mpos_mcell/100) + mpos_mcell%100];
      hist[index]->Fill(count);
    }
    infile.close();
  }
  
  for(Int_t index=0; index<N_MPPC*N_pins; index++){
    Int_t mpos_mcell = reverse_HIST_map[index];
    hist[index] -> Write(Form("hist_%s_pos%d_pin%d_volt%d",(light==1)?"wLight":"Pedest",mpos_mcell/100,mpos_mcell%100,volt));
    delete hist[index];
  }
  
}

bool checkDirectory()
{
    struct stat st,st2;
    if (stat(Form("%s/USM%03d",MainDir,Box),&st) == 0){
        if (stat(Form("%s/USM%03d/T%02d",MainDir,Box,temp,Directory),&st2) == 0){
                cout << "Directory is already created" <<endl;
                return false;
        }else{
	  mkdir(Form("%s/USM%03d/T%02d_%d_run%03d",MainDir,Box,temp,Directory,Run),0777);
        }
    }else{
        mkdir(Form("%s/USM%03d",MainDir,Box),0777);
        mkdir(Form("%s/USM%03d/T%02d_%d_run%03d",MainDir,Box,temp,Directory,Run),0777);
        return true;
    }


}


bool readParameters()
{
		char header[50];
    ifstream parameters;
    parameters.open("../parameters.txt");
    if ( parameters.is_open()) {
				parameters>>header>>Box;        
				parameters>>header>>Directory;
        parameters>>header>>Run;
        parameters>>header>>temp;
				parameters>>header>>MainDir;
        parameters.close();
        return true;
    } else {
        cout<<"You need to set the parameters in parameters.txt"<<endl;
    }
    parameters.close();
    return false;
}


int main() {
  TDatime start;
  readParameters();
  checkDirectory();  
  
  TFile * file = new TFile(Form("%s/USM%03d/T%02d_%d_run%03d/gaussian_fit.root",MainDir,Box,temp,Directory,Run),"RECREATE");
  read_qdc();
  
  for ( int volt = min_volt; volt <= max_volt; volt++ ) {
    try {
      readFiles(Run,temp,0,volt,file);//for the Pedestal run
    } catch (int e) {
      cout<<"Error de lectura... Se continua"<<endl;
    }
    try {
      readFiles(Run,temp,1,volt,file);//for the wLight run
    } catch(int e) {
      cout<<"Error de lectura... Se continua"<<endl;
    }
  }
  
  file->Close();
  delete file;
  file = 0L;
  TDatime end;
 	double time_use=(end.Convert()-start.Convert());
	cout<<"elapsed time: "<<time_use<<"[s]"<<endl;
  return 0;
  
}
