/* TreeParcer makes a root file with histograms of energy, time of flight and energy vs TOF

USE
./Cloverdata filename.root data.txt calibrationfile.txt

##OUTPUT##

makes fileName.root
in current directory

####

Mac Wheeler
Jul 26, 2023
*/

#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <cstring>
#include <sstream>

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH2D.h"
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
Double_t e[Maxsize];    // energy
ULong64_t e_t[Maxsize]; // time stamp(4 ns)
Int_t cfd[Maxsize];     // constant fraction discrimnation

struct cloverHIT
{

  int cloverID;
  int crystalID;
  int crystalside;
  unsigned int energy;
  double time;
};

struct clover
{

  int nhits;
  cloverHIT hits[20];
  int cloverID;
  int naddback = 0;
  int ncrystalhits = 0;

  void reset()
  {
    nhits = 0;
  }
  void add_hits(int cloverID, int crystalID, int energy, int side, double time)
  {
    if (energy == 0)
    {
      cout << "zero hit added\n";
    }
    hits[nhits].cloverID = cloverID;
    hits[nhits].crystalID = crystalID;
    hits[nhits].energy = energy;
    hits[nhits].crystalside = side;
    hits[nhits].time = time;
    nhits += 1;
  }
};

struct dega
{

  clover detectors[16];
  int ndet;

  double offset[16][4];
  double gain[16][4];
  int clovers[16][4];

  void PullcloverIDs()
  {

    ifstream clover("CloverID.txt");

    string ID[2];
    string dud;

    if (clover.is_open())
    {

      getline(clover, dud);
      string myline;

      while (getline(clover, myline))
      {

        string delim = "\t\t\t";
        size_t pos = 0;
        int i = 0;

        while ((pos = myline.find(delim)) != string::npos)
        {

          ID[i] = myline.substr(0, pos);
          myline.erase(0, pos + delim.length());
          i++;
        }
        cout << "myline:" << myline << endl;
        clovers[stoi(ID[0]) - 1][stoi(myline)] = stoi(ID[1]);
      }
    }
    else
    {
      cerr << "could not open file CloverID\n";
      exit(1);
    }

    clover.close();
  }

  int IDclover(int crystalID, int &cloverID, int &crystalside)
  {
    for (int j = 0; j < 13; j++)
    {
      for (int y = 0; y < 4; y++)
      {
        if (clovers[j][y] == crystalID)
        {
          cloverID = j + 1;
          crystalside = y;
          return 0;
        }
      }
    }
    return 0;
  }

  void Pullcalibration(string CloverEnCal)
  {

    fstream CalEn(CloverEnCal);

    if (CalEn.is_open())
    {
      std::string En;
      std::stringstream ss;

      while (getline(CalEn, En))
      {
        if (En.size() == 0)
        {
          continue;
        }

        ss.clear();
        ss.str(En);

        int ID;
        double intercept, slope;

        ss >> ID >> intercept >> slope;

        cout << ID << " " << intercept << " " << slope << endl;

        for (int g = 0; g < 13; g++)
        {
          for (int y = 0; y < 4; y++)
          {
            if (clovers[g][y] == (ID - 1))
            {
              offset[g][y] = intercept;
              gain[g][y] = slope;
            }
          }
        }
      }
    }
    else
    {
      cerr << "CalEn couldnt open cloverEncal250\n";
      exit(1);
    }
    CalEn.close();
  }

  double calibrate(double energy, int crystalside, int clover)
  {
    // cout << " energy " << energy << " crystalside " << crystalside << " clover " << clover << " ";
    energy = offset[clover - 1][crystalside] +
             (gain[clover - 1][crystalside] * energy);
    // cout << " offset " << offset[clover][crystalside] << " gain " << gain[clover][crystalside] << endl;
    return energy;
  }
};

// start of func main
int main(int argc, char *argv[])
{

  dega array; // defining data type dega

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

  // starting tree chain
  TChain *ch = new TChain("tree");

  // opening file with the data
  ifstream myfile(argv[2]);

  string myline;

  if (myfile.is_open())
  {

    while (myfile.good())
    {

      getline(myfile, myline);

      TString line(myline);

      ch->Add(line); // adding files to chain
    }
  }
  else
  {
    cerr << "Couldn't open rundata file\n";
    exit(1);
  }

  myfile.close();

  array.PullcloverIDs();
  cout << "did it make it here?\n";
  array.Pullcalibration(argv[3]);
  for (int g = 0; g < 13; g++)
  {
    for (int y = 0; y < 4; y++)
    {

      cout << "clover:" << array.clovers[g][y] << endl;
      cout << "offset:" << array.offset[g][y] << endl;
      cout << "gain:" << array.gain[g][y] << endl;
    }
  }

  // looping through and adding the trees of files to chain

  TH2D *rawEn = new TH2D("rawEn", "2D histogram of raw data from DAQ", 4096, 0.0, 65536.0, 416, 0, 416.0); // Raw energy Histo
  TH2D *gamEn = new TH2D("gamEn", "2D histogram of crystal vs energy", 8192, 0.0, 8192, 52, 0, 52);        // histo of crystal vs energy
  TH2D *gamaddback = new TH2D("gamaddback", "2D histogram of clover vs addback", 8192, 0.0, 8192, 13, 1, 14);
  TH2D *hittdiff = new TH2D("hittdiff", "2D Histo of time difference between crystals", 1000, -1000, 1000, 6, 0, 6);

  // making a new file to write histograms onto
  TFile *f2 = new TFile(argv[1], "RECREATE"); // making a new file

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

    if ((i % 75) == 0)
    {
      cout << i << " " << nentries << '\n';
      double f = nentries;
      double g = i *100.0;

      cout << (g/f) << "%" <<'\n';
    }

    ch->GetEntry(i); // getting entries

    array.ndet = 0;

    if ((multi == 0))
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

      if ((e[k] > 62500 || e[k] == 0))
      {
        continue;
      }
      if ((int)e[k] == 0)
      {
        cout << e[k] << "entery:" << k << " " << i << " why are there zeros?\n";
      }

      rawEn->Fill(e[k], detID[k]);

      if ((detID[k] < 52 && detID[k] >= 0))
      {
        // find which clover and which crystal the ID corrisponds to
        int cloverID, crystalside;

        array.IDclover(detID[k], cloverID, crystalside);

        int found = 0;
        // find out if this clover already exits in this event
        for (int j = 0; j < array.ndet; ++j)
        {

          if (array.detectors[j].cloverID == cloverID)
          {

            // add this entry to the detector

            found = 1;
            if (e[k] == 0)
            {
              cout << "why zeros?\n";
            }
            array.detectors[j].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }

        if (found == 0)
        { // if cant find det ID
          array.detectors[array.ndet].reset();
          if ((int)e[k] == 0)
          {
            cout << "why not zeros?\n";
          }
          array.detectors[array.ndet].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          array.detectors[array.ndet].cloverID = cloverID;
          array.ndet += 1;
        }
      }
    }
    // loop through all clovers
    for (int g = 0; g < array.ndet; g++)
    {
      int crystal;
      int crystalside;
      double energyr;
      double energyp;
      double timediff;
      double gammas[10];
      int addback[10] = {0};
      int clover = array.detectors[g].cloverID;
      int ngammas = 0;

      if ((clover == 8) || (clover == 4))
      {
        continue;
      }

      // to see if there is any gamma rays that hit mutliple crystals in one clover
      for (int p = 0; p < array.detectors[g].nhits; p++)
      {

        int side1 = array.detectors[g].hits[p].crystalside;
        int crystal = array.detectors[g].hits[p].crystalID;
        energyp = array.calibrate(array.detectors[g].hits[p].energy, side1, clover);
        gamEn->Fill(energyp, crystal);
        for (int r = 0; r < array.detectors[g].nhits; r++)
        {

          int side2 = array.detectors[g].hits[r].crystalside;
          energyr = array.calibrate(array.detectors[g].hits[r].energy, side2, clover);

          if (r == p)
          {
            continue;
          }

          if ((addback[r] == 1) || (addback[p] == 1))
          {
            continue;
          }
          if ((energyp < 30) || (energyr < 30))
          {
            continue;
          }

          timediff = array.detectors[g].hits[p].time - array.detectors[g].hits[r].time;

          if ((side1 > side2))
          {
            hittdiff->Fill(timediff, indx[side1 - 1][side2 - 1]);
          }
          else if (side2 > side1)
          {
            hittdiff->Fill(-timediff, indx[side1 - 1][side2 - 1]);
          }
          if ((sqrt(pow(timediff, 2)) < 200))
          {
            addback[r] = 1;
            addback[p] = 1;
            gammas[ngammas] = array.calibrate(array.detectors[g].hits[r].energy, side2, clover) + array.calibrate(array.detectors[g].hits[p].energy, side1, clover);
            ngammas++;
          }
        }
        if (addback[p] == 0)
        {
          gammas[ngammas] = array.calibrate(array.detectors[g].hits[p].energy, side1, clover);
          ngammas++;
        }
      }

      for (int y = 0; y < ngammas; y++)
      {
        gamaddback->Fill(gammas[y], clover);
      }
    }
  }

  hittdiff->Write();
  gamaddback->Write();
  rawEn->Write();
  gamEn->Write();

  f2->Close();

  return 0; // end of function main
}