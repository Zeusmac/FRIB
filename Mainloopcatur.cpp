/*
Same as Mainloop.cpp but made for Caturs event builder instead of ryans

hello
####
./Mainloop histfile.root Datafile.txt clovercalibrationfile.txt

Mac Wheeler
Feb 23, 2024
*/
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TString.h"
#include "TCutG.h"

using namespace std;

int gatedcut40 = 5; // var for cut num of the gate for msx40

// defining max size of arrays
const int Maxsize = 400;

// defining tree branch vars
Int_t multi;            // multiplitiy of arrays
Int_t detID[Maxsize];   // Detector ID
Int_t e[Maxsize];    // energy kev
Long64_t e_t[Maxsize]; // time stamp(4 ns)
Int_t cfd[Maxsize];     // constant fraction discrimnation

fstream()


struct hit
{
  int ID;
  int energy;
  double time;
};

struct PIDhits
{
  int nhits;
  hit hits[200];
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
  int ncuts = 0;
  bool yncut;

  void read_cuts(std::string cut_txt, bool cut = 0)
  {
    /*cut_txt formatting:
    #comments
    IDnumber    npoints
    X   y
    X   y
    X   y
    X   y
    */
    yncut = cut;
    if (yncut == 0)
    {
      return;
    }
    std::fstream cutfile(cut_txt);

    if (cutfile.is_open())
    {

      std::string line;
      std::stringstream ss;

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
        std::string name;

        ss >> ID;
        ss >> npoints;
        ss >> name;

        cuts[ID] = new TCutG();
        cut_names[ID] = name;

        std::cout << "loading cut: " << ID << std::endl;
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

          cuts[ID]->AddPoint(x, y);
        }
        hcuts[ID] = new TH2D(cut_names[ID], "2D histogram cut", 500, -150.0, -50.0, 500, 0.0, 6000);
        ncuts = ID + 1;
      }
    }
    else
    {
      std::cerr << " could not open cuts sfile\n";
      exit(1);
    }
    cutfile.close();
  }
  void Fill_cuts(double x, double y)
  {
    if (yncut == 0)
    {
      return;
    }
    else
    {
      for (int h = 0; h < ncuts; h++)
      {
        if (cuts[h]->IsInside(x, y) == 1)
        {
          hcuts[h]->Fill(x, y);
        }
      }
    }
  }
  int PID_gate(int cutnum, double x, double y)
  {
    return (cuts[cutnum]->IsInside(x, y));
  }
  void Write_cuts()
  {
    if (yncut == 0)
    {
      return;
    }
    for (int y = 0; y < ncuts; y++)
    {
      hcuts[y]->Write();
    }
  }
};

struct YSOHIT
{

  int anode;
  int energy;
  double time;
};

struct YSOHITS
{

  int nhits;
  YSOHIT hits[200];
  int anode;

  void reset()
  {
    nhits = 0;
  }
  void add_hits(int anode, int energy, double time)
  {
    hits[nhits].anode = anode;
    hits[nhits].energy = energy;
    hits[nhits].time = time;
    nhits += 1;
  }
};

struct YSOdet
{
  YSOHITS hsides[8];
  YSOHITS lsides[8];
  YSOHITS vetof;
  YSOHITS vetor;
  int nhside;
  int nlside;

  int YSOID[2][6];

  void Pull_YSOID()
  {

    std::ifstream YSO("YSOIDcopy.txt");

    std::string IDs[2];
    int h_l = 0;

    if (YSO.is_open())
    {

      std::string myline;
      while (getline(YSO, myline))
      {

        std::string delim = " ";
        size_t pos = 0;

        while ((pos = myline.find(delim)) != std::string::npos)
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
      std::cerr << "could not open file\n";
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

struct cloverHIT
{

  int cloverID;
  int crystalID;
  int crystalside;
  int energy;
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
    // cout << "reset" << endl;
    nhits = 0;
  }
  void add_hits(int cloverID, int crystalID, int energy, int side, double time)
  {
    // cout << "add hit" << endl;
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

  double offset[60];
  double gain[60];
  int clovers[13][4];

  void PullcloverIDs()
  {

    std::ifstream clover("CloverID2.txt");

    std::string dud;
    std::stringstream ss;

    if (clover.is_open())
    {

      getline(clover, dud);
      std::string myline;

      while (getline(clover, myline))
      {

        if (myline.size() == 0)
        {
          continue;
        }

        ss.clear();
        ss.str(myline);

        int clID, crystal, side;

        ss >> clID >> crystal >> side;

        clovers[clID - 1][side] = crystal;
        //cout << (clID - 1) << " " << side << " " << crystal << std::endl;
      }
    }
    else
    {
      std::cerr << "could not open file CloverID\n";
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

  void Pullcalibration(std::string CloverEnCal)
  {

    std::fstream CalEn(CloverEnCal);

    std::string cal[3];

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
        offset[(ID-1)] = intercept;
        gain[(ID -1)] = slope;
        }
      }
    else
    {
      std::cerr << "CalEn couldnt open cloverEncal250\n";
      exit(1);
    }
    CalEn.close();
  }

  double calibrate(double energy, int cloverID)
  {
    //cout << " energy " << energy << " crystalside " << crystalside << " clover " << clover << " ";
    energy = offset[cloverID-1] +
             (gain[cloverID] * energy);
    //cout << "cloverID" << cloverID+1 << " offset " << offset[cloverID] << " gain " << gain[cloverID] << endl;
    return energy;
  }
};

struct store
{
  vector<double> Imptime;
  vector<double> Impx;
  vector<double> Impy;

  vector<double> betatime;
  vector<double> betax;
  vector<double> betay;

  int found = 0;

  void add_imp(double time, double x, double y)
  {
    Imptime.push_back(time);
    Impx.push_back(x);
    Impy.push_back(y);
  }
  void add_beta(double time, double x, double y)
  {
    betatime.push_back(time);
    betax.push_back(x);
    betay.push_back(y);
  }
};

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

  // making a new file to write histograms onto
  TFile *f2 = new TFile(argv[1], "RECREATE");

  det PID;
  YSOdet YSO;
  dega gamma;
  store store40;
  store store100;

  // starting a Tchaing for seg 250
  TChain *ch = new TChain("Pixie16");

  ifstream myfile(argv[2]);

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
  TH2D *rawEn = new TH2D("rawEn", "2D histogram of energy of Scints", 4096, 0.0, 65536.0, 416, 0, 416.0); // Raw energy Histo

  TH2D *PID100 = new TH2D("PID100", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000);
  TH2D *PID40 = new TH2D("PID40", "2D histogram of msx40 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000); // energy vs TOF

  TH2D *YSOhigh = new TH2D("YSOhigh", "3D histogram of the high gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0, 1.0);
  TH2D *YSOlow = new TH2D("YSOlow", "3D histogram of the low gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0.0, 1.0);
  TH2D *GatedYSO100 = new TH2D("GatedYSO100", "3D histogram of the low gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0, 1.0);
  TH2D *GatedYSO40 = new TH2D("GatedYSO40", "3D histogram of the low gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0.0, 1.0);
  TH2D *YSOhEn = new TH2D("YSOhEn", "2D histo of the high gain anodes and dynode vs energy", 4096, 0.0, 65536.0, 5, 0, 5);
  TH2D *YSOlEn = new TH2D("YSOlEn", "2D histo of the low gain anodes and dynode vs energy", 4096, 0.0, 65536.0, 5, 0, 5);
  TH2D *YSOhtfilt = new TH2D("YSOhfilt", "2D histogram of the high gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0, 1.0);
  TH2D *YSOhtspacefilt = new TH2D("YSOhtspacefilt", "2D histogram of the high gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0, 1.0);

  TH2D *tdiffh = new TH2D("tdiffh", "2D histogram of high gain time differences between anodes", 100, -500, 500, 6, 1, 7);
  TH2D *tdiffl = new TH2D("tdiffl", "2D histogram of low gain time differences between anodes", 100, -500, 500, 6, 1, 7);
  TH1D *tdifflh40 = new TH1D("tdifflh40", "1d histo of the difference between high gain and low gain points", 800, -100e6, 100e6);
  TH1D *tdifflh100 = new TH1D("tdifflh100", "1d histo of the difference between high gain and low gain points", 1000, -100e6, 100e6);

  TH2D *vetoEn = new TH2D("vetoEn", "Histo of front and rear veto det vs energy", 4096, 0.0, 65536.0, 2, 0, 2);
  TH1D *vetofgated = new TH1D("vetofgated", "Histo of front veto vs energy", 4096, 0.0, 65536.0);
  TH2D *vetordiff = new TH2D("vetordiff", "Time diffs of rear veto between high gain and low gain", 100, -500, 500, 8, 0, 8);
  TH2D *vetofdiff = new TH2D("vetofdiff", "Time diffs of front veto between high gain and low gain", 100, -500, 500, 8, 0, 8);

  TH2D *gamEn = new TH2D("gamEn", "2D histogram of crystal vs energy", 8192, 0.0, 8192, 52, 0, 52); // histo of crystal vs energy
  TH2D *gamaddback = new TH2D("gamaddback", "2D histogram of clover vs addback", 8192, 0.0, 8192, 13, 1, 14);
  TH2D *hittdiff = new TH2D("hittdiff", "2D Histo of time difference between crystals", 1000, -1000, 1000, 6, 1, 7);

  gamma.PullcloverIDs();
  gamma.Pullcalibration(argv[3]);
  PID.read_cuts("PIDcuts.txt", true);
  YSO.Pull_YSOID();

  // settingaddress of branches to vars
  ch->SetBranchAddress("mult", &multi);
  ch->SetBranchAddress("id", detID);
  ch->SetBranchAddress("energy", e);
  ch->SetBranchAddress("pxitime", e_t);
  ch->SetBranchAddress("cfdtime", cfd);

  // getting the number of enteries
  Int_t nentries = (Int_t)ch->GetEntries();

  // looping through all the entries
  for (int i = 0; i < nentries; i++)
  {

    PID.nms = 0;
    PID.ncross = 0;
    PID.nimage = 0;
    YSO.nhside = 0;
    YSO.nlside = 0;
    gamma.ndet = 0;

    // setting nhits to zero
    PID.cross_scint.reset();

    YSO.vetof.reset();
    YSO.vetor.reset();

    ch->GetEntry(i); // getting entries

    // if ((i % 75) == 0)
    //  {
    //    cout << i << " " << nentries << '\n';
    //    double f = nentries;
    //    double g = i * 100.0;

    //    cout << (g / f) << "%" << '\n';
    //  }

    if (multi < 3)
    {
      continue;
    }

    if (multi > Maxsize)
    {

      cerr << multi << " Error:increase maxsize!\n";
      exit(1);
    }

    // looping through all hits in the event
    for (int k = 0; k < multi; k++)
    {

      int foundp;
      int ID = detID[k];
      rawEn->Fill(e[k], ID); // filling histo

      if (((detID[k] == 234 || detID[k] == 235) && e[k] > 5000) && e[k] < 25000)
      { // R or L scint for TOF
        foundp = 0;
        for (int m = 0; m < PID.nimage; m++)
        {
          if (PID.image_scint[m].ID = ID)
          {
            foundp = 1;
            PID.image_scint[m].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }
        if (foundp == 0)
        {
          PID.image_scint[PID.nimage].reset();
          PID.image_scint[PID.nimage].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          PID.image_scint[PID.nimage].ID = ID;
          PID.nimage++;
        }
      }
      if (detID[k] == 248)
      { // B2 cross scint

        PID.cross_scint.add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
        PID.ncross = 1;
      }
      if (detID[k] == 240 || detID[k] == 241)
      { // msx40 or msx100 for energy
        foundp = 0;
        ID = detID[k];
        for (int m = 0; m < PID.nimage; m++)
        {
          if (PID.msx[m].ID = ID)
          {
            foundp = 1;
            PID.msx[m].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }
        if (foundp == 0)
        {
          PID.msx[PID.nms].reset();
          PID.msx[PID.nms].add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          PID.msx[PID.nms].ID = ID;
          PID.nms++;
        }
      }

      int anode;
      int foundy = 0;
      if (((detID[k] <= 215 && detID[k] >= 212) || (detID[k] = 208))) //&& (e[k] > 4000 && e[k] < 60000))
      {
        //cout << "high\n";
        // setting to high gain anodes for YSO
        anode = YSO.IDYSO(detID[k]);
        // cout << anode << endl;

        // find out if this clover already exits in this event
        for (int j = 0; j < YSO.nhside; ++j)
        {
          if (YSO.hsides[j].anode == anode)
          {
            foundy = 1;
            YSO.hsides[j].add_hits(anode, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }

        if (foundy == 0)
        { // if cant find det ID
          YSO.hsides[YSO.nhside].reset();
          YSO.hsides[YSO.nhside].add_hits(anode, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          YSO.hsides[YSO.nhside].anode = anode;
          YSO.nhside += 1;
        }
      }
      if (((detID[k] <= 219 && detID[k] >= 216) || detID[k] == 209)) //&& (e[k] > 207 && e[k] < 60000))
      {
        cout << "low:"<< e[k] << "\n";
        // low gain annodes for YSO
        anode = YSO.IDYSO(detID[k]);

        for (int j = 0; j < YSO.nlside; ++j)
        {
          if (YSO.lsides[j].anode == anode)
          {
            foundy = 1;
            YSO.lsides[j].add_hits(anode, e[k], ((e_t[k] + cfd[k] / 16384.0) * 4.0));
          }
        }

        if (foundy == 0)
        { // if cant find det ID
          YSO.lsides[YSO.nlside].reset();
          YSO.lsides[YSO.nlside].add_hits(anode, e[k], ((e_t[k] + cfd[k] / 16384.0) * 4.0));
          YSO.lsides[YSO.nlside].anode = anode;
          YSO.nlside += 1;
        }
      }
      if ((detID[k] == 336) && (e[k] > 100 && e[k] < 60000)) // front veto scint
      {
        YSO.vetof.add_hits(0, e[k], (e_t[k])); //+ cfd[k] / 16384.0) * 4.0);
      }
      if (detID[k] == 337 && (e[k] > 100 && e[k] < 60000)) // rear veto scint
      {
        YSO.vetor.add_hits(1, e[k], (e_t[k])); //+ cfd[k] / 16384.0) * 4.0);
      }
      if ((detID[k] <= 341 && detID[k] >= 256)) //For filling the clovers
      {
        if ((e[k] > 62500 || e[k] == 0))
        {
          continue;
        }
        // find which clover and which crystal the ID corrisponds to
        int cloverID, crystalside;
        int foundg = 0;
        gamma.IDclover(detID[k], cloverID, crystalside);

        // cout << "cloverID:" << cloverID << " side:" << crystalside << endl;

        // find out if this clover already exits in this event
        for (int j = 0; j < gamma.ndet; ++j)
        {
          if (gamma.detectors[j].cloverID == cloverID)
          {
            // add this entry to the detector

            foundg = 1;
            gamma.detectors[j].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k]+ (cfd[k] / 16384.0)) * 4.0);
          }
        }

        if (foundg == 0)
        { // if cant find det ID
          gamma.detectors[gamma.ndet].reset();
          gamma.detectors[gamma.ndet].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k]+ (cfd[k] / 16384.0)) * 4.0);
          gamma.detectors[gamma.ndet].cloverID = cloverID;
          gamma.ndet += 1;
        }
      }
    }

    double TOF;
    int energy[4] = {0};
    if ((PID.nms >= 1 && PID.nimage >= 1) && PID.ncross >= 1)
    {
      double stime;
      double etime;
      double time[2];

      for (int k = 0; k < PID.nimage; k++)
      {
        for (int l = 0; l < PID.image_scint[k].nhits; l++)
        {
          if (energy[0] < PID.image_scint[k].hits[l].energy)
          {
            energy[0] = PID.image_scint[k].hits[l].energy;
            stime = PID.image_scint[k].hits[l].time;
          }
        }
      }
      for (int l = 0; l < PID.cross_scint.nhits; l++)
      {
        if (energy[1] < PID.cross_scint.hits[l].energy)
        {
          energy[1] = PID.cross_scint.hits[l].energy;
          etime = PID.cross_scint.hits[l].time;
        }
      }
      for (int k = 0; k < PID.nms; k++)
      {
        for (int l = 0; l < PID.msx[k].nhits; l++)
        {
          if (PID.msx[k].hits[l].ID == 240) // msx40 energy
          {
            if (energy[2] < PID.msx[k].hits[l].energy)
            {
              energy[2] = PID.msx[k].hits[l].energy;
              time[0] = PID.msx[k].hits[l].time;
            }
          }
          if (PID.msx[k].hits[l].ID == 241) // msx100 energy
          {
            if (energy[3] < PID.msx[k].hits[l].energy)
            {
              energy[3] = PID.msx[k].hits[l].energy;
              time[1] = PID.msx[k].hits[l].time;
            }
          }
        }
      }

      TOF = (etime - stime);       // to ns
      PID40->Fill(TOF, energy[2]); // Filling MX40  PID

      // cout << "TOF: " << TOF << " Energy" << energy[3] << endl;
      PID100->Fill(TOF, energy[3]); // Filling MX100 PID

      PID.Fill_cuts(TOF, energy[2]); // Filling PDI cuts
    }

    float vEn[2] = {0};
    float vtime[2];

    for (int v = 0; v < YSO.vetof.nhits; v++)
    {
      if (vEn[0] < YSO.vetof.hits[v].energy)
      {
        vEn[0] = YSO.vetof.hits[v].energy;
        vtime[0] = YSO.vetof.hits[v].time;
      }
    }
    for (int o = 0; o < YSO.vetor.nhits; o++)
    {
      if (vEn[1] < YSO.vetor.hits[o].energy)
      {
        vEn[1] = YSO.vetor.hits[o].energy;
        vtime[1] = YSO.vetof.hits[o].time;
      }
    }
    if (vEn[0] > 0)
    {
      if(PID.PID_gate(gatedcut40, TOF, energy[2]) == 1){
        vetofgated -> Fill(vEn[0]);
      }
      vetoEn->Fill(vEn[0], 0);
    }
    if (vEn[1] > 0)
    {
      vetoEn->Fill(vEn[1], 1);
    }

    double maxtdh = 0; // max time diff between anodes for high gain
    double maxtdl = 0; // max time diff between anodes for low gainß
    float vfdiff;      // time diff front veto - high and low gain
    float vrdiff;      // time diff rear veto - high and low gain

    float al[5] = {0};    // array for energy of each anode
    float timel[5] = {0}; // array for time of each low energy hit

    if (((YSO.nlside == 4))) //&& !(YSO.vetor.nhits > 0)) && (YSO.vetof.nhits > 0))
    { // make sure all anodes four of low gain hit and there no rear veto hit

      float tdl;

      for (int g = 0; g < YSO.nlside; g++)
      {

        int ell = YSO.lsides[g].anode;

        for (int p = 0; p < YSO.lsides[g].nhits; p++)
        {

          if (al[ell] < YSO.lsides[g].hits[p].energy)
          {

            al[ell] = YSO.lsides[g].hits[p].energy;
            timel[ell] = YSO.lsides[g].hits[p].time;
            // cout << timel[ell] << endl;
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
        vetofdiff->Fill(vfdiff, f);
        vetordiff->Fill(vfdiff, f);
      }

      float suml = al[0] + al[1] + al[2] + al[3];
      float yl = (al[3] + al[2]) / suml;
      float xl = (al[1] + al[2]) / suml;

      if ((xl > 1) || (yl > 1))
      {
        cout << "x: " << xl << "\n";
        cout << "y: " << yl << "\n";
        cout << "something went wrong low gain\n";
      }
      if (maxtdl < 10000)
      {
        if (PID.PID_gate(4, TOF, energy[3]) == 1) // gate for the msx100
        {
          GatedYSO100->Fill(xl, yl);
          // cout << "add100" << endl;
          store100.add_imp(timel[4], xl, yl);
        }
        if (PID.PID_gate(gatedcut40, TOF, energy[2]) == 1) // gate for the msx40
        {
          GatedYSO40->Fill(xl, yl);
          store40.add_imp(timel[4], xl, yl);
          //cout << "add40" << endl;
        }
        YSOlow->Fill(xl, yl);
      }
    }

    float ah[5] = {0}; // array for energy of each anode
    float timeh[5] = {0};

    if ((YSO.nhside == 5)) //&& !(YSO.vetor.nhits > 0 || YSO.vetof.nhits > 0))
    { // make sure all four anodes of high gain hit. Also making sure both rear and front veto dont hits

      float tdh;

      for (int g = 0; g < YSO.nhside; g++)
      {

        int el = YSO.hsides[g].anode;

        for (int p = 0; p < YSO.hsides[g].nhits; p++) // looping through the hits and picking the highest energy hit
        {

          if (ah[el] < YSO.hsides[g].hits[p].energy)
          {
            ah[el] = YSO.hsides[g].hits[p].energy;
            timeh[el] = YSO.hsides[g].hits[p].time;
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
        vetofdiff->Fill(vfdiff, f + 4);
        vetordiff->Fill(vfdiff, f + 4);
      }

      float sumh = ah[0] + ah[1] + ah[2] + ah[3]; // equations for finding YSO x and y
      float yh = (ah[3] + ah[2]) / sumh;
      float xh = (ah[1] + ah[2]) / sumh;

      if ((xh > 1) || (yh > 1))
      {
        cout << "x: " << xh << "\n";
        cout << "y: " << yh << "\n";
        cout << "something went wrong high gain\n";
      }

      if (maxtdh < 10000)
      {
        YSOhigh->Fill(xh, yh);
        store40.add_beta(timeh[4], xh, yh);
        store100.add_beta(timeh[4], xh, yh);
      }
    }
    for (int g = 0; g < gamma.ndet; g++)
    {
      int crystal;
      int crystalside;
      double energyr;
      double energyp;
      double timediff;
      double gammas[10];
      int addback[10] = {0};
      int clover = gamma.detectors[g].cloverID;
      int ngammas = 0;

      if ((clover == 8) || (clover == 4))// there are no clover 4 and 8
      {
        continue;
      }

      // to see if there is any gamma rays that hit mutliple crystals in one clover
      for (int p = 0; p < gamma.detectors[g].nhits; p++)
      {
        int side1 = gamma.detectors[g].hits[p].crystalside;
        int crystal = gamma.detectors[g].hits[p].crystalID;
        // cout << " energy " << gamma.detectors[g].hits[p].energy << " crystalside " << side1 << " clover " << clover << " " << endl;
        energyp = gamma.calibrate(gamma.detectors[g].hits[p].energy, (clover-1)*4 + side1 );
        //cout <<  "energyp" << energyp << endl;
        gamEn->Fill(energyp, (clover-1)*4+side1);
        for (int r = 0; r < gamma.detectors[g].nhits; r++)
        {

          int side2 = gamma.detectors[g].hits[r].crystalside;
          energyr = gamma.calibrate(gamma.detectors[g].hits[r].energy, (clover-1)*4 + side2);

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

          timediff = gamma.detectors[g].hits[p].time - gamma.detectors[g].hits[r].time;

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
            gammas[ngammas] = gamma.calibrate(gamma.detectors[g].hits[r].energy, (clover-1)*4 + side2) + gamma.calibrate(gamma.detectors[g].hits[p].energy, (clover-1)*4 + side1);
            ngammas++;
          }
        }
        if (addback[p] == 0)
        {
          gammas[ngammas] = gamma.calibrate(gamma.detectors[g].hits[p].energy,(clover-1)*4 + side1);
          ngammas++;
        }
      }

      for (int y = 0; y < ngammas; y++)
      {
        gamaddback->Fill(gammas[y], clover);
      }
    }
  }
  for (int l = 0; l < store40.betatime.size();l++){
    for(int k = 0; k < store40.Imptime.size();k++) {
     if ((abs(store40.betatime[l]-store40.Imptime[k]) < 100e6)) { //
      double r = r = sqrt(pow((store40.Impx[k] - store40.betax[l]), 2) + pow((store40.Impy[k] - store40.betay[l]), 2));
      YSOhtspacefilt -> Fill(store40.betax[l],store40.betay[l]);
      if (r < .25){
        // cout << "made it " << store40.betax[l]<< endl;
        // cout << "diff " << store40.betatime[l]-store40.Imptime[k] << endl;
        
        tdifflh40 -> Fill(store40.betatime[l]-store40.Imptime[k]);
      }
     }
    }
  }
  for (int l = 0; l < store100.betatime.size();l++){
    for(int k = 0; k < store100.Imptime.size();k++) {
     if ((abs(store100.betatime[l]-store100.Imptime[k]) < 100e6)) {
      double r = r = sqrt(pow((store100.Impx[k] - store100.betax[l]), 2) + pow((store100.Impy[k] - store100.betay[l]), 2));
      if (r < .25){
        tdifflh100 -> Fill(store100.betatime[l]-store100.Imptime[k]);
      }
     }
    }
  }

  PID.Write_cuts();

  rawEn->Write();

  GatedYSO100->Write();
  GatedYSO40->GetXaxis()->SetTitle("x and y");
  GatedYSO40->GetYaxis()->SetTitle("x and y");
  GatedYSO40->Write();
  YSOhtfilt->Write();
  YSOhtspacefilt->Write();

  YSOhEn->Write();
  YSOlEn->Write();
  YSOhigh->Write();
  YSOlow->Write();

  tdiffh->Write();
  tdiffl->Write();
  tdifflh40->Write();
  tdifflh100->Write();
  cout << "tdifflh40:" << (tdifflh40->Integral()) << endl;
  cout << "tdifflh100:" << (tdifflh100->Integral()) << endl;

  vetoEn->Write();
  vetofgated -> Write();
  vetordiff->Write();
  vetofdiff->Write();

  PID100->Write();
  PID40->Write();

  hittdiff->Write();
  gamaddback->Write();
  gamEn->Write();

  f2->Close();
  return 0; // end of function main
}