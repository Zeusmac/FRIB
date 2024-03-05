/*



Mac Wheeler

Aug 28, 2023
*/

#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <cstring>

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TChain.h"
#include "TCutG.h"
#include "TF1.h"
#include "TString.h"

using namespace std;

// defining max size of arrays
const int Maxsize = 200;

// defining tree branch vars
Int_t multi;            // multiplitiy of arrays
Int_t detID[Maxsize];   // Detector ID
Double_t e[Maxsize];    // energy kev
ULong64_t e_t[Maxsize]; // time stamp(4 ns)
Int_t cfd[Maxsize];     // constant fraction discrimnation

struct YSOHIT
{

  int anode;
  int energy;
  int time;
};

struct YSO
{

  int nhits;
  YSOHIT hits[20];
  int anode;

  void reset()
  {
    nhits = 0;
  }
  void add_hits(int anode, int energy, int time)
  {
    hits[nhits].anode = anode;
    hits[nhits].energy = energy;
    hits[nhits].time = time;
    nhits += 1;
  }
};

struct dega
{

  YSO hsides[8];
  YSO lsides[8];
  YSO vetof;
  YSO vetor;
  int nhside;
  int nlside;

  int YSOID[2][6];

  void Pull_YSOID()
  {

    ifstream YSO("YSOID.txt");

    string IDs[2];
    int h_l = 0;

    if (YSO.is_open())
    {

      string myline;
      while (getline(YSO, myline))
      {

        string delim = " ";
        size_t pos = 0;

        

        while ((pos = myline.find(delim)) != string::npos)
        {

          IDs[0] = myline.substr(0, pos);
          myline.erase(0, pos + delim.length());
          IDs[1] = myline;
        }
        YSOID[h_l][stoi(IDs[0])] = stoi(IDs[1]);
        if (stoi(IDs[0]) == 4)
        {
          h_l++;
        }
      }
    }
    else
    {
      cerr << "could not open file\n";
      exit(1);
    }

    YSO.close();
  }

  int IDYSO(int ID)
  {
    for (int h = 0; h < 10; h++)
      for (int k = 0; k < 2; k++)
      {
        {
          if (YSOID[k][h] == ID)
          {
            return h;
          }
        }
      }
    return 0;
  }
};

// start of func main
int main(int argc, char *argv[])
{

  int indx[4][4];

  indx[0][0] = -1;
  indx[1][1] = -1;
  indx[3][3] = -1;
  indx[4][4] = -1;

  indx[0][1] = 1;
  indx[1][0] = 1;
  indx[2][0] = 2;
  indx[0][2] = 2;
  indx[3][0] = 3;
  indx[0][3] = 3;

  indx[1][2] = 4;
  indx[2][1] = 4;
  indx[1][3] = 5;
  indx[3][1] = 5;

  indx[2][3] = 6;
  indx[3][2] = 6;

  dega array; // defining data type dega

  array.Pull_YSOID();

  // starting tree chain
  TChain *ch = new TChain("tree");

  // opening file with the data
  ifstream myfile("250rundata.txt");

  string myline;

  if (myfile.is_open())
  {

    while ( getline(myfile, myline))
    {

      TString line(myline);

      ch->Add(line); // adding files to chain
    }
  }
  else
  {
    cout << "Couldn't open file 250\n";
  }

  myfile.close();

  // Making histos
  TH2D *rawEn = new TH2D("rawEn", "2D histogram of energy of Scints", 4096, 0.0, 65536.0, 416, 0, 416.0); // Raw energy Histo

  TH2D *YSOhigh = new TH2D("YSOhigh", "3D histogram of the high gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0, 1.0);
  TH2D *YSOlow = new TH2D("YSOlow", "3D histogram of the low gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0.0, 1.0);
  TH2D *YSOhEn = new TH2D("YSOhEn", "2D histo of the high gain anodes and dynode vs energy", 4096, 0.0, 65536.0, 5, 0, 5);
  TH2D *YSOlEn = new TH2D("YSOlEn", "2D histo of the low gain anodes and dynode vs energy", 4096, 0.0, 65536.0, 5, 0, 5);

  TH2D *tdiffh = new TH2D("tdiffh", "2D histogram of high gain time differences between anodes", 100, -500, 500, 6, 1, 7);
  TH2D *tdiffl = new TH2D("tdiffl", "2D histogram of low gain time differences between anodes", 100, -500, 500, 6, 1, 7);

  TH2D *vetoEn = new TH2D("vetoEn","Histo of front and rear veto det vs energy", 4096, 0.0 , 65536.0, 2, 0, 2);
  TH2D *vetordiff = new TH2D("vetordiff","Time diffs of rear veto between high gain and low gain", 100, -500 , 500, 8, 0, 8);
  TH2D *vetofdiff = new TH2D("vetofdiff","Time diffs of front veto between high gain and low gain", 100, -500 , 500 , 8, 0, 8);
  // making a new file to write histograms onto
  TFile *f2 = new TFile(argv[1], "RECREATE");

  // settingaddress of branches to vars
  ch->SetBranchAddress("multi", &multi);
  ch->SetBranchAddress("detID", detID);
  ch->SetBranchAddress("e", e);
  ch->SetBranchAddress("e_t", e_t);
  ch->SetBranchAddress("cfd", cfd);

  // getting the number of enteries
  Int_t nentries = (Int_t)ch->GetEntries();

  // looping through all the entries
  for (int i = 0; i < nentries; i++)
  {

    /*if ((i % 75) == 0)
    {
      cout << i << " " << nentries << '\n';
      double f = nentries;
      double g = i * 100.0;

      cout << (g / f) << "%" << '\n';
    }*/

    ch->GetEntry(i); // getting entries

    array.nhside = 0;
    array.nlside = 0;

    //setting nhits to zero for veto dets
    array.vetof.reset();
    array.vetor.reset();

    if (multi < 4)
    {
      continue;
    } // skiping entries with multi=0
    if (multi > Maxsize)
    {

      cerr << multi << " Error:increase maxsize!\n";
      exit(1);
    }


    // looping through all events
    for (int k = 0; k < multi; k++)
    {
      if (e[k] > 60000)
      {
        continue;
      }

      rawEn->Fill(e[k], detID[k]);

      int anode;
      int found = 0;
      if ((detID[k] <= 255 && detID[k] >= 250)&& e[k] > 4000)
      {
        // setting to high gain anodes
        anode = array.IDYSO(detID[k]);
        // cout << anode << endl;

        // find out if this clover already exits in this event
        for (int j = 0; j < array.nhside; ++j)
        {
          if (array.hsides[j].anode == anode)
          {
            found = 1;
            array.hsides[j].add_hits(anode, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
          }
        }

        if (found == 0)
        { // if cant find det ID
          array.hsides[array.nhside].reset();
          array.hsides[array.nhside].add_hits(anode, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
          array.hsides[array.nhside].anode = anode;
          array.nhside += 1;
        }
      }
      if ((detID[k] <= 205 && detID[k] >= 200) && e[k] > 207)
      {
        // low gain annodes
        anode = array.IDYSO(detID[k]);

        for (int j = 0; j < array.nlside; ++j)
        {
          if (array.lsides[j].anode == anode)
          {
            found = 1;
            array.lsides[j].add_hits(anode, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
          }
        }

        if (found == 0)
        { // if cant find det ID
          array.lsides[array.nlside].reset();
          array.lsides[array.nlside].add_hits(anode, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
          array.lsides[array.nlside].anode = anode;
          array.nlside += 1;
        }
      }
      if((detID[k] == 290) && (e[k] > 100))//front veto scint
      { 
        array.vetof.add_hits(0,e[k],(e_t[k] + cfd[k] / 16384.0) * 4.0);
      } 
      if(detID[k] == 291 && (e[k] > 100))//rear veto scint
      {
        array.vetor.add_hits(1,e[k],(e_t[k] + cfd[k] / 16384.0) * 4.0);
      }
      

      float vEn[2] = {0};
      float vtime[2];

      for(int v = 0; v < array.vetof.nhits;v++)
      {
        if(vEn[0] < array.vetof.hits[v].energy)
        {
          vEn[0] = array.vetof.hits[v].energy;
          vtime[0] = array.vetof.hits[v].time;
        }
      }
      for(int o = 0; o < array.vetor.nhits; o++)
      {
        if(vEn[1] < array.vetor.hits[o].energy){
          vEn[1] = array.vetor.hits[o].energy;
          vtime[1] = array.vetof.hits[o].time;
        }
      }
      if( vEn[0] > 0){
      vetoEn -> Fill(vEn[0],0);
      }
      if(vEn[1] > 0){
      vetoEn -> Fill(vEn[1],1);
      }


      double maxtdh = 0;//max time diff between anodes for high gain
      double maxtdl = 0;//max time diff between anodes for low gainß
      float vfdiff; //time diff front veto - high and low gain
      float vrdiff;// time diff rear veto - high and low gain
      
      if ((array.nhside == 5)&& !(array.vetor.nhits > 0 || array.vetof.nhits > 0))
      { // make sure all four anodes of high gain hit. Also making sure both rear and front veto dont hit

        float ah[5] = {0}; // array for energy of each anode
        float timeh[5];
        float tdh;

        for (int g = 0; g < array.nhside; g++)
        {

          int el = array.hsides[g].anode;

          for (int p = 0; p < array.hsides[g].nhits; p++)
          {

            if (ah[el] < array.hsides[g].hits[p].energy)
            {

              ah[el] = array.hsides[g].hits[p].energy;
              timeh[el] = array.hsides[g].hits[p].time;
            }
          }
          YSOhEn->Fill(ah[el], el);
        }
        for (int f = 0; f < 4; f++)
        {
         
          for (int u = f + 1; u < 4; u++)
          {
            tdh = timeh[u] - timeh[f];
            if (abs(tdh) > maxtdh)
            {
              maxtdh = abs(tdh);
            }
            if ((u > f))
            {
              tdiffh->Fill(tdh, indx[u][f]);
            }
            else if (f > u)
            {
              tdiffh->Fill(-tdh, indx[u][f]);
            }
            
          }
          vfdiff = vtime[0] - timeh[f];
          vrdiff = vtime[1] - timeh[f];
          vetofdiff -> Fill(vfdiff,f+4);
          vetordiff -> Fill(vfdiff,f+4);
        }

        float sumh = ah[0] + ah[1] + ah[2] + ah[3];
        float yh = (ah[3] + ah[2]) / sumh;
        float xh = (ah[1] + ah[2]) / sumh;
        float avg = sumh / 4.0;

        if ((xh > 1) || (yh > 1))
        {
          cout << "x: " << xh << "\n";
          cout << "y: " << yh << "\n";
          cout << "something went wrong high gain\n";
        }
        if (maxtdh < 10)
        {
          YSOhigh->Fill(xh, yh);
        }
      }
      if (((array.nlside == 5)&& !(array.vetor.nhits > 0)) && (array.vetof.nhits > 0))
      { // make sure all anodes four of low gain hit and there no rear veto hit

        float al[5] = {0}; // array for energy of each anode
        float timel[5];
        float tdl;

        for (int g = 0; g < array.nlside; g++)
        {

          int ell = array.lsides[g].anode;

          for (int p = 0; p < array.lsides[g].nhits; p++)
          {

            if (al[ell] < array.lsides[g].hits[p].energy)
            {

              al[ell] = array.lsides[g].hits[p].energy;
              timel[ell] = array.lsides[g].hits[p].time;
            }
          }
          YSOlEn->Fill(al[ell], ell);
        }
        for (int f = 0; f < 4; f++)
        {
          for (int u = f + 1; u < 4; u++)
          {
            if (abs(tdl) > maxtdl)
            {
              maxtdl = abs(tdl);
            }
            tdl = timel[u] - timel[f];
            if ((u > f))
            {
              tdiffl->Fill(tdl, indx[u][f]);
            }
            else if (f > u)
            {
              tdiffl->Fill(-tdl, indx[u][f]);
            }
          }
          vfdiff = vtime[0] - timel[f];
          vrdiff = vtime[1] - timel[f];
          vetofdiff -> Fill(vfdiff,f);
          vetordiff -> Fill(vfdiff,f);
        }

        float suml = al[0] + al[1] + al[2] + al[3];
        float yl = (al[3] + al[2]) / suml;
        float xl = (al[1] + al[2]) / suml;
        float avg = suml / 4.0;

        if ((xl > 1) || (yl > 1))
        {
          cout << "x: " << xl << "\n";
          cout << "y: " << yl << "\n";
          cout << "something went wrong low gain\n";
        }
        if (maxtdl < 10)
        {
          YSOlow->Fill(xl, yl);
        }
      }
    }
  } 
  rawEn->Write();

  vetoEn ->Write();

  YSOhEn->Write();
  YSOlEn->Write();
  YSOhigh->Write();
  YSOlow->Write();
  
  tdiffh->Write();
  tdiffl->Write();
  vetordiff -> Write();
  vetofdiff -> Write();


  f2->Close();

  return 0; // end of function main
}