// 
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: 2015.04.21
// 

#include <TApplication.h> 


// TCanvas* drawPt(vector<TString> inputFiles,
// 		TString histType){

// }


void print_usage(){
  printf("NAME\n\tdrawPt - draw transverse momentum spectrum\n");
  printf("\nSYNOPSIS\n\tdrawPt label input.root\n "); 
  // printf("\nSYNOPSIS\n\tdrawPt [-t hist-type ] [-opt draw-option]\n "); 
  // printf("\t[-h hist-name ] [-vmax max-value] [-npad num-pad] [-b] input1 input2 ...\n");
  printf("\nOPTIONS\n");
  printf("\t%-5s  %-40s\n", "label", "options: pho");
  // printf("\t%-5s  %-40s\n", "-t", "hist type [TH1D, TH2D]");
  // printf("\n\t%-5s  %-40s\n", "-opt", "draw option for histgram [colz, surf2]");
  // printf("\n\t%-5s  %-40s\n", "-h", "histogram name");
  // printf("\n\t%-5s  %-40s\n", "-vmax", "limit the histogram within vmax-value");
  // printf("\n\t%-5s  %-40s\n", "-vmin", "limit the histogram within vmin-value");
  // printf("\n\t%-5s  %-40s\n", "-p", "print pixel values to stdout");
  printf("\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n");
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

  TApplication theApp("App", 0, 0);
  theApp.SetReturnFromRun(true);
  //drawPt(inputFiles);
  theApp.Run();
}

