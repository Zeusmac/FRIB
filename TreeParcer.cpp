/* TreeParcer makes a root file with histograms of energy, time of flight and energy vs TOF

USE
./TreeParcer data.txt fileName.root

##OUTPUT##

makes fileName.root


####

Mac Wheeler
Jun 8, 2023
*/

#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TString.h"
#include "TCutG.h"

using namespace std;

// defining max size of arrays
const int Maxsize = 200;

// defining tree branch vars
Int_t multi;            // multiplitiy of arrays
Int_t detID[Maxsize];   // Detector ID
Double_t e[Maxsize];    // energy kev
ULong64_t e_t[Maxsize]; // time stamp(4 ns)
Int_t cfd[Maxsize];     // constant fraction discrimnation

struct hit
{
  int ID;
  int energy;
  double time;
};

struct PIDhits
{
  int nhits;
  hit hits[20];
  int ID;

  void reset()
  {
    nhits = 0;
  }

  void add_hit(int ID, int energy, double time)
  {
    hits[nhits].ID = ID;
    hits[nhits].energy = energy;
    hits[nhits].time = time;
    hits[nhits];
    nhits++;
  }
};

struct det
{

  PIDhits msx[2];
  PIDhits cross_scint;
  PIDhits image_scint[4];
  int ncross;
  int nimage;
  int nms;

  TCutG *cuts[100];
  TH2D *hcuts[100];
  TString cut_names[100];
  int cutID[100];
  int ncuts;

  void read_cuts(string cut_txt)
  {
    /*cut_txt formatting:
    #comments
    IDnumber    npoints
    X   y
    X   y
    X   y
    X   y
    */
    fstream cutfile(cut_txt);

    if (cutfile.is_open())
    {

      string line;
      stringstream ss;

      while (getline(cutfile, line))
      {
        char f = line[0];
        if (line.size() == 0)
        {
          continue;
        }
        if (f == '#')
        {
          continue;
        }
        if (f == ' ')
        {
          continue;
        }

        ss.clear();
        ss.str(line);

        int ID, npoints;
        string name;
        ss >> ID;
        ss >> npoints;
        ss >> name;

        cout << ID << " " << npoints << endl;

        cuts[ncuts] = new TCutG();
        cut_names[ncuts] = name;
        cutID[ncuts] = ID;

        cout << "loading cuts" << ID << endl;
        for (int h = 0; h < npoints; h++)
        {
          getline(cutfile, line);
          if (line.size() == 0)
          {
            continue;
          }
          if (f == '#')
          {
            continue;
          }
          if (f == ' ')
          {
            continue;
          }

          ss.clear();
          ss.str(line);

          double x, y;
          ss >> x >> y;
          cout << x << " " << y << endl;

          cuts[ncuts]->AddPoint(x, y);
        }
        hcuts[ncuts] = new TH2D(cut_names[], "2D histogram cut", 500, -150.0, -50.0, 500, 0.0, 6000);
        ncuts++;
      }
      
    }
    else
    {
      cerr << " could not open cuts sfile\n";
      exit(1);
    }
    cutfile.close();
  }
  void Fill_cuts(double x, double y)
  {
    for (int h = 0; h < ncuts; h++)
    {
      if (cuts[h]->IsInside(x, y))
      {
        hcuts[h]->Fill(x, y);
      }
    }
  }
  void Write_cuts()
  {
    for (int y = 0; y < ncuts; y++)
    {
      hcuts[y]->Write();
    }
  }
};

int main(int argc, char *argv[])
{

  // making a new file to write histograms onto
  TFile *f2 = new TFile(argv[2], "RECREATE");

  det array;

  // starting a Tchaing for seg 250
  TChain *ch = new TChain("tree");

  ifstream myfile(argv[1]);

  string myline;

  if (myfile.is_open())
  {

    while (getline(myfile, myline))
    {
      // cout << myline << endl;
      TString line(myline);

      ch->Add(line); // adding files to chain
    }
  }
  else
  {
    cout << "Couldn't open run data file\n";
  }

  myfile.close();

  // Setting up histograms
  TH2D *PID100 = new TH2D("PID100", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000);
  TH2D *PID40 = new TH2D("PID40", "2D histogram of msx40 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000); // energy vs TOF
  TH2D *rawEn = new TH2D("rawEn", "2D histogram of energy of Scints", 4096, 0.0, 65536.0, 416, 0, 416.0);     // Raw energy Histo

  array.read_cuts("PIDcuts.txt");

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

    array.nms = 0;
    array.ncross = 0;
    array.nimage = 0;

    array.cross_scint.reset();

    ch->GetEntry(i); // getting entries

    /* if ((i % 75) == 0)
     {
       cout << i << " " << nentries << '\n';
       double f = nentries;
       double g = i * 100.0;

       cout << (g / f) << "%" << '\n';
     }*/

    if (multi < 3)
    {
      continue;
    }

    if (multi > Maxsize)
    {

      cerr << multi << " Error:increase maxsize!\n";
      exit(1);
    }

    // looping through all events
    for (int k = 0; k < multi; k++)
    {

      int found;
      int ID = detID[k];
      rawEn->Fill(e[k], ID); // filling histo

      if (((detID[k] == 300 || detID[k] == 301) && e[k] > 5000) && e[k] < 25000)
      { // R or L scint
        found = 0;
        for (int m = 0; m < array.nimage; m++)
        {
          if (array.image_scint[m].ID = ID)
          {
            found = 1;
            array.image_scint[m].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }
        if (found == 0)
        {
          array.image_scint[array.nimage].reset();
          array.image_scint[array.nimage].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          array.image_scint[array.nimage].ID = ID;
          array.nimage++;
        }
      }
      if (detID[k] == 307)
      { // B2 cross scint

        array.cross_scint.add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
        array.ncross = 1;
      }
      if (detID[k] == 308 || detID[k] == 309)
      { // msx40
        found = 0;
        ID = detID[k];
        for (int m = 0; m < array.nimage; m++)
        {
          if (array.msx[m].ID = ID)
          {
            found = 1;
            array.msx[m].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }
        if (found == 0)
        {
          array.msx[array.nms].reset();
          array.msx[array.nms].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          array.msx[array.nms].ID = ID;
          array.nms++;
        }
      }
    }
    if ((array.nms >= 1 && array.nimage >= 1) && array.ncross >= 1)
    {

      int energy[4] = {0};
      double stime;
      double etime;
      double time[2];

      for (int k = 0; k < array.nimage; k++)
      {
        for (int l = 0; l < array.image_scint[k].nhits; l++)
        {
          if (energy[0] < array.image_scint[k].hits[l].energy)
          {
            energy[0] = array.image_scint[k].hits[l].energy;
            stime = array.image_scint[k].hits[l].time;
          }
        }
      }
      for (int l = 0; l < array.cross_scint.nhits; l++)
      {
        if (energy[1] < array.cross_scint.hits[l].energy)
        {
          energy[1] = array.cross_scint.hits[l].energy;
          etime = array.cross_scint.hits[l].time;
        }
      }
      for (int k = 0; k < array.nms; k++)
      {
        for (int l = 0; l < array.msx[k].nhits; l++)
        {
          if (array.msx[k].hits[l].ID == 308)
          {
            if (energy[2] < array.msx[k].hits[l].energy)
            {
              energy[2] = array.msx[k].hits[l].energy;
              time[0] = array.msx[k].hits[l].time;
            }
          }
          if (array.msx[k].hits[l].ID == 309)
          {
            if (energy[3] < array.msx[k].hits[l].energy)
            {
              energy[3] = array.msx[k].hits[l].energy;
              time[1] = array.msx[k].hits[l].time;
            }
          }
        }
      }

      double TOF = (etime - stime); // to ns
      PID40->Fill(TOF, energy[2]);
      PID100->Fill(TOF, energy[3]);

      array.Fill_cuts(TOF, energy[2]);
    }
  }

  array.Write_cuts();

  rawEn->Write();
  PID100->Write();
  PID40->Write();

  f2->Close();
  return 0; // end of function main
}
