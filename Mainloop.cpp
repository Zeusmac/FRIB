/*

Na32-> 885

####
./Mainloop histfile.root Datafile.txt clovercalibrationfile.txt

Mac Wheeler
Nov 8, 2023
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

// defining max size of arrays
const int Maxsize = 200;

// defining tree branch vars
Int_t multi;            // multiplitiy of arrays
Int_t detID[Maxsize];   // Detector ID
Double_t e[Maxsize];    // energy kev
ULong64_t e_t[Maxsize]; // time stamp(4 ns)
Int_t cfd[Maxsize];     // constant fraction discrimnationx

double decaytime = 400e6; // ns time window for decay curve

struct hit
{
  int ID;
  int energy;
  double time;
};

struct PIDhits
{
  int nhits;
  hit hits[100];
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

  TCutG *cuts[30];
  TH2D *hcuts[30];
  TH2D *YSOlcut[30];
  TH1D *tdifflh[30];
  TString cut_names[30];
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
        YSOlcut[ID] = new TH2D(cut_names[ID] + "YSOhigh", "3D histogram of the low gain YSO Surface and energy", 4096, 0.0, 1.0, 4096, 0.0, 1.0);
        tdifflh[ID] = new TH1D(cut_names[ID] + "tdifflh", "1d histo of the difference between high gain and low gain points", 200, -(decaytime), decaytime);

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
  void Fill_cut(double x, double y)
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
      YSOlcut[y]->Write();
      tdifflh[y]->Write();
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
  YSOHIT hits[100];
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
    // cout << "anode " << anode << " energy " << energy << " time " << endl;
  }
};

struct YSOdet
{
  YSOHITS hsides[5];
  YSOHITS lsides[5];
  YSOHITS vetof;
  YSOHITS vetor;
  int nhside;
  int nlside;

  int YSOlowID[5];
  int YSOhighID[5];

  int IDYSO(int ID)
  {
    for (int h = 0; h < 5; h++)
    {
      if (YSOlowID[h] == ID)
      {
        return h;
      }
      if (YSOhighID[h] == ID)
      {
        return h;
      }
    }
    cerr << ID << "YSOID problem" << endl;
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
  cloverHIT hits[100];
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

  vector<double> offset;
  vector<double> gain;
  int clovers[13][4];

  int IDclover(int crystalID, int &cloverID, int &crystalside)
  {
    for (int j = 0; j < 13; j++)
    {
      for (int u = 0; u < 4; u++)
      {
        if (clovers[j][u] == crystalID)
        {
          cloverID = j + 1;
          crystalside = u;
          // cout <<  u << endl;
          // cout << "ID " << j+1 << " side " << y << " crystal " << crystalID << endl;
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
        offset.push_back(intercept);
        gain.push_back(slope);
      }
    }
    else
    {
      std::cerr << "CalEn couldnt open cloverEncal250\n";
      exit(1);
    }
    CalEn.close();
  }

  double calibrate(double energy, int cl)
  {
    // cout << " energy " << energy << " crystalside " << crystalside << " clover " << clover << " ";
    energy = offset[cl] +
             (gain[cl] * energy);
    // cout << cl << " offset " << offset[cl-1] << " gain " << gain[cl-1] << endl;
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

  det PID;
  YSOdet YSO;
  dega gamma;
  store store40;
  store store100[10];
  vector<int> gatecuts;
  int cloverminID; // min and max clover IDs for the if statement
  int clovermaxID;
  int TOFID[3]; // stime1, stime2 then etime IDS
  int msx40ID;
  int msx100ID;
  int vetofID; // front veto ID
  int vetorID; // rear veto ID

  ifstream Map("IDMap.txt");

  if (Map.is_open()) // unpacking ID map
  {
    string mapline;
    stringstream ssmap;
    string dud;

    getline(Map, dud);
    getline(Map, dud);

    ssmap.clear();
    ssmap.str(mapline);

    string junk;
    int clovernum;

    ssmap >> junk >> clovernum;

    getline(Map, dud);
    for (int i = 0; i < 52; i++)
    {
      getline(Map, mapline);
      ssmap.clear();
      ssmap.str(mapline);
      int clnum;
      int side;
      int clID;

      ssmap >> clnum >> clID >> side;

      gamma.clovers[clnum - 1][side] = clID;
      // cout << gamma.clovers[clnum-1][side] << endl;
    }

    int k = 0;
    while (getline(Map, mapline))
    {
      if (k == 17)
      {
        break;
      }
      if (((mapline.size() == 0) || (mapline[0] == '#')) || (mapline[0] == ' '))
      {
        continue;
      }
      ssmap.clear();
      ssmap.str(mapline);
      int ID;
      ssmap >> junk >> ID;
      if (k <= 2)
      {
        TOFID[k] = ID;
        k++;
        continue;
      }
      if ((k == 3))
      {
        msx40ID = ID;
        k++;
        continue;
      }
      if ((k == 4))
      {
        msx100ID = ID;
        k++;
        continue;
      }
      if ((k == 5))
      {
        vetofID = ID;
        k++;
        continue;
      }
      if ((k == 6))
      {
        vetorID = ID;
        k++;
        continue;
      }
      if ((k == 7))
      {
        YSO.YSOhighID[4] = ID;
        k++;
        continue;
      }
      if ((k == 8))
      {
        YSO.YSOlowID[4] = ID;
        k++;
        continue;
      }
      if ((k >= 9) && (k <= 12))
      {
        YSO.YSOhighID[k - 9] = ID;
        k++;
        continue;
      }
      if ((k >= 13) && (k <= 16))
      {
        YSO.YSOlowID[k - 13] = ID;
        k++;
        continue;
      }
    }
    while (getline(Map, mapline))
    {
      if (((mapline.size() == 0) || (mapline[0] == '#')) || (mapline[0] == ' '))
      {
        continue;
      }
      ssmap.clear();
      ssmap.str(mapline);

      int gatecut;
      ssmap >> gatecut;
      gatecuts.push_back(gatecut);
    }
    Map.close();
  }
  else
  {
    cerr << "couldnt open map\n";
  }

  cloverminID = gamma.clovers[0][0];
  clovermaxID = gamma.clovers[0][0];
  for (int l = 0; l < 13; l++)
  {
    for (int o = 0; o < 4; o++)
    {
      if (cloverminID > gamma.clovers[l][o])
      {
        cloverminID = gamma.clovers[l][o];
      }
      if (clovermaxID < gamma.clovers[l][o])
      {
        clovermaxID = gamma.clovers[l][o];
      }
    }
  }
  vector<double> avg_TOF; // avg TOF per run
  vector<int> entries;    // number of enteries per run

  ifstream avgTOF("avgTOFfile.txt"); // pulling avg TOF info for each run

  if (avgTOF.is_open())
  {
    
    string TOFline;
    stringstream ss;
    while (getline(avgTOF, TOFline))
    {
      ss.clear();
      ss.str(TOFline);
      double avgTOF;
      ss >> avgTOF;
      avg_TOF.push_back(avgTOF); // time of flight corrections
    }
  }
  else
  {
    cerr << "couldnt open avgTOFfile.txt" << endl;
  }
  avgTOF.close();

  // making a new file to write histograms onto
  TFile *f2 = new TFile(argv[1], "RECREATE");

  // starting a Tchaing for seg 250
  TChain *ch = new TChain("tree");

  ifstream myfile(argv[2], ios::in);

  string myline;

  if (myfile.is_open())
  {
    int p = 0;
    int runent;
    while (getline(myfile, myline))
    {
      if ((myline[19] == '0') && (myline[18] == '0'))
      {
        cout << myline << endl;
        runent = (Int_t)ch->GetEntries();
        if (runent > 0)
        {
          entries.push_back(runent);
        }
        cout << runent << endl;
      }
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

  int indx[4][4]; // index for the clover crystal time diff histo

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

  TH1D *TOFHist[50]; // for finding the avg TOF perun
  TString name;

  for (int y = 0; y < (entries.size() + 1); y++)
  {
    name = "TOF" + to_string(y);
    TOFHist[y] = new TH1D(name, "TOF", 500, -150.0, -50.0); // making the a his for each run
  }

  // Setting up histograms
  TH2D *rawEn = new TH2D("rawEn", "2D histogram of energy of Scints", 4096, 0.0, 65536.0, 416, 0, 416.0); // Raw energy Histo

  TH2D *PID100Imp = new TH2D("PID100Imp", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000);
  TH2D *PID100ImpBox = new TH2D("PID100ImpBox", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 1000.0, 6000);
  TH2D *PID40ImpBox = new TH2D("PID40ImpBox", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 1000.0, 6000);
  TH2D *PID100 = new TH2D("PID100", "Histo of msx100 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000);
  TH2D *PID40 = new TH2D("PID40", "2D histogram of msx40 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000);       // energy vs TOF
  TH2D *PID40Imp = new TH2D("PID40Imp", "2D histogram of msx40 energy vs TOF", 500, -150.0, -50.0, 500, 0.0, 6000); // energy vs TOF

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
  TH1D *tdifflh40 = new TH1D("tdifflh40", "1d histo of the difference between high gain and low gain points", 200, -(decaytime), decaytime);

  TH2D *vetoEn = new TH2D("vetoEn", "Histo of front and rear veto det vs energy", 4096, 0.0, 65536.0, 2, 0, 2);
  TH1D *vetofgated = new TH1D("vetofgated", "Histo of front veto vs energy", 4096, 0.0, 65536.0);
  TH2D *vetordiff = new TH2D("vetordiff", "Time diffs of rear veto between high gain and low gain", 100, -500, 500, 8, 0, 8);
  TH2D *vetofdiff = new TH2D("vetofdiff", "Time diffs of front veto between high gain and low gain", 100, -500, 500, 8, 0, 8);

  TH2D *gamEn = new TH2D("gamEn", "2D histogram of crystal vs energy", 8192, 0.0, 8192, 52, 0, 52); // histo of crystal vs energy
  TH2D *gamaddback = new TH2D("gamaddback", "2D histogram of clover vs addback", 8192, 0.0, 8192, 13, 1, 14);
  TH2D *hittdiff = new TH2D("hittdiff", "2D Histo of time difference between crystals", 1000, -1000, 1000, 6, 1, 7);

  gamma.Pullcalibration(argv[3]);
  PID.read_cuts("PIDcuts2.txt", true);

  // settingaddress of branches to vars
  ch->SetBranchAddress("multi", &multi);
  ch->SetBranchAddress("detID", detID);
  ch->SetBranchAddress("e", e);
  ch->SetBranchAddress("e_t", e_t);
  ch->SetBranchAddress("cfd", cfd);

  // getting the number of enteries
  Int_t nentries = (Int_t)ch->GetEntries();

  entries.push_back(nentries);

  int this_run = 0; // var to determine the run number

  // looping through all the entries
  for (int i = 0; i < nentries; i++)
  {
    ch->GetEntry(i); // getting entries

    if (multi > Maxsize)
    {

      cerr << multi << " Error:increase maxsize!\n";
      exit(1);
    }
    if (entries[this_run] <= i) // determining run number
    {
      this_run++;
    }
    if (multi < 3)
    {
      continue;
    }

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

    if ((i % 175) == 0)
    {
      // cout << i << " " << nentries << '\n';
      double f = nentries;
      double g = i * 100.0;
      printf("\r%.2f%%", (g / f));
    }

    // cout << "fill" << endl;
    //  looping through all hits in the event
    for (int k = 0; k < multi; k++)
    {

      int foundp;
      int ID = detID[k];
      rawEn->Fill(e[k], ID); // filling histo

      if (((detID[k] == TOFID[0] || detID[k] == TOFID[1]) && e[k] > 5000) && e[k] < 25000)
      { // R or L scint
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
      if (detID[k] == TOFID[2])
      { // B2 cross scint

        PID.cross_scint.add_hit(ID, e[k], (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
        PID.ncross = 1;
      }
      if (detID[k] == msx40ID || detID[k] == msx100ID)
      { // msx40 and msx100
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
      if ((detID[k] <= YSO.YSOhighID[3] && detID[k] >= YSO.YSOhighID[0]) || (detID[k] == YSO.YSOhighID[4]))
      {
        if ((e[k] < 840 || e[k] > 60000))
        { // gateing on the dynode
          continue;
        }
        // setting to high gain anodes
        anode = YSO.IDYSO(detID[k]);
        // cout << detID[k] << endl;

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
      if ((detID[k] <= YSO.YSOlowID[3] && detID[k] >= YSO.YSOlowID[0]) || (detID[k] == YSO.YSOlowID[4]))
      {
        if ((e[k] < 230 || e[k] > 60000))
        { // gateing on the dynode
          continue;
        }
        // low gain annodes
        anode = YSO.IDYSO(detID[k]);
        // cout << "lanode " << anode << "detID " << detID[k] << endl;
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
      if ((detID[k] == vetofID) && (e[k] > 100 && e[k] < 60000)) // front veto scint
      {
        YSO.vetof.add_hits(0, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
      }
      if (detID[k] == vetorID && (e[k] > 100 && e[k] < 60000)) // rear veto scint
      {
        YSO.vetor.add_hits(1, e[k], (e_t[k] + cfd[k] / 16384.0) * 4.0);
      }
      if ((detID[k] < clovermaxID && detID[k] >= cloverminID))
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
            gamma.detectors[j].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          }
        }

        if (foundg == 0)
        { // if cant find det ID
          gamma.detectors[gamma.ndet].reset();
          gamma.detectors[gamma.ndet].add_hits(cloverID, detID[k], e[k], crystalside, (e_t[k] + (cfd[k] / 16384.0)) * 4.0);
          gamma.detectors[gamma.ndet].cloverID = cloverID;
          gamma.ndet += 1;
        }
      }
    }

    double TOF;    // time of flight
    double TOFcor; // corrected time of flight (by run)
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
          if (PID.msx[k].hits[l].ID == msx40ID) // msx40 energy
          {
            if (energy[2] < PID.msx[k].hits[l].energy)
            {
              energy[2] = PID.msx[k].hits[l].energy;
              time[0] = PID.msx[k].hits[l].time;
            }
          }
          if (PID.msx[k].hits[l].ID == msx100ID) // msx100 energy
          {
            if (energy[3] < PID.msx[k].hits[l].energy)
            {
              energy[3] = PID.msx[k].hits[l].energy;
              time[1] = PID.msx[k].hits[l].time;
            }
          }
        }

      }

      TOF = (etime - stime); // to ns
      if(avg_TOF.size() > 0){
      TOFcor = TOF + (avg_TOF[0] - avg_TOF[this_run]);
      }
      else{
        TOFcor = TOF;
      }

      PID40->Fill(TOF, energy[2]); // Filling MX40  PID

      // cout << "TOF: " << TOF << " Energy" << energy[3] << endl;
      PID100->Fill(TOF, energy[3]); // Filling MX100 PID
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
      // if (PID.PID_gate(gatecuts[0], TOF, energy[2]) == 1)
      // {
      //   vetofgated->Fill(vEn[0]);
      // }
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
    // cout << "ysol" << endl;

    if (((YSO.nlside == 5) && !(YSO.vetor.nhits > 0)) && ((YSO.vetof.nhits > 0) && (PID.nms >= 1 && PID.ncross >= 1)) && (PID.nimage >= 1))
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
      if (maxtdl < 10)
      {
        PID100Imp->Fill(TOFcor, energy[3]);
        if ((abs(TOF) > 105 && abs(TOF) < 135) && energy[3] > 1000)
        {
          PID100ImpBox->Fill(TOFcor, energy[3]);
          TOFHist[this_run]->Fill(TOF); // FIlling TOF histo of this run
        }
        for (int y = 0; y < PID.ncuts; y++)
        {
          if (PID.PID_gate(y, TOFcor, energy[3]) == 1) // gate for the msx100
          {
            store100[y].add_imp(timel[4], xl, yl);
            PID.YSOlcut[y]->Fill(xl, yl);
          }
        }
        PID40Imp->Fill(TOFcor, energy[2]);
        if ((TOF < -105 && TOF > -130) && energy[3] > 1000)
        {
          PID40ImpBox->Fill(TOFcor, energy[3]);
        }
        if (PID.PID_gate(gatecuts[0], TOFcor, energy[2]) == 1) // gate for the msx40
        {
          GatedYSO40->Fill(xl, yl);
          store40.add_imp(timel[4], xl, yl);
          // cout << "add40" << endl;
        }
        YSOlow->Fill(xl, yl);
      }
    }

    float ah[5] = {0}; // array for energy of each anode
    float timeh[5] = {0};
    // cout << "ysoh" << endl;
    if (((YSO.nhside == 5) && ((PID.nms == 0) && PID.nimage == 0)) && (PID.ncross == 0 && !(YSO.vetor.nhits > 0)))
    { // make sure all four anodes of high gain hit. Also making sure both rear and front veto dont hits
      // cout << "helo"<< endl;
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

      if (maxtdh < 10)
      {
        YSOhigh->Fill(xh, yh);
        store40.add_beta(timeh[4], xh, yh);
        store100[0].add_beta(timeh[4], xh, yh);
      }
    }
    // cout << "gamma" << endl;
    for (int g = 0; g < gamma.ndet; g++)
    {
      int crystal;
      int crystalside;
      double energyr;
      double energyp;
      double timediff;
      double gammas[10];
      int addback[10] = {0};
      int cl = gamma.detectors[g].cloverID;
      int ngammas = 0;

      if ((cl == 8) || (cl == 4))
      {
        continue;
      }

      // to see if there is any gamma rays that hit mutliple crystals in one clover
      for (int p = 0; p < gamma.detectors[g].nhits; p++)
      {
        int side1 = gamma.detectors[g].hits[p].crystalside;
        // cout << "side1" << side1 << endl;
        // cout << "clover" << cl << endl;
        //  cout << " energy " << gamma.detectors[g].hits[p].energy << " crystalside " << side1 << " clover " << clover << " " << endl;
        energyp = gamma.calibrate(gamma.detectors[g].hits[p].energy, ((cl - 1) * 4 + side1));
        // cout << "energyp" << energyp << endl;
        gamEn->Fill(energyp, (cl - 1) * 4 + side1 + 1);
        for (int r = 0; r < gamma.detectors[g].nhits; r++)
        {

          int side2 = gamma.detectors[g].hits[r].crystalside;
          if (side2 > 3)
          {
            // cout << "what?? " << cl << " " << side2 << endl;
          }
          // cout << "side2" << side2 << endl;
          // cout << "clover" << cl << endl;
          energyr = gamma.calibrate(gamma.detectors[g].hits[r].energy, ((cl - 1) * 4 + side2));

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
            gammas[ngammas] = gamma.calibrate(gamma.detectors[g].hits[r].energy, ((cl - 1) * 4 + side2)) + gamma.calibrate(gamma.detectors[g].hits[p].energy, ((cl - 1) * 4 + side1));
            ngammas++;
          }
        }
        if (addback[p] == 0)
        {
          gammas[ngammas] = gamma.calibrate(gamma.detectors[g].hits[p].energy, ((cl - 1) * 4 + side1));
          ngammas++;
        }
      }

      for (int y = 0; y < ngammas; y++)
      {
        gamaddback->Fill(gammas[y], cl);
      }
    }
  }

  printf("\rEntries are done looping\n");
  // if(avg_TOF.empty() == true)
  // {
  //   ofstream TOFfile("avgTOFfile.txt", ios::trunc); 
  //   cout << avg_TOF.size() << endl;
  //   // storing avg TOFs in file
  //   if (TOFfile.is_open())
  //   {
  //     cout << "is it open?" << endl;
  //     for (int t = 0; t < (entries.size()); t++)
  //     {
  //       TOFfile << TOFHist[t]->GetMean(1) << endl;
  //       printf("\rstoring TOF %i\n", t);
  //     }
  //     cout << "closing" << endl;
  //      TOFfile.close();
  //   }
  //   else
  //   {
  //     cerr << "cant open TOFfile\n";
  //   }
  // }

  cout << "\rNow doing the beta implant corrilations\n";
  for (int l = 0; l < store40.betatime.size(); l++)
  {
    for (int k = 0; k < store40.Imptime.size(); k++)
    {
      if ((abs(store40.betatime[l] - store40.Imptime[k]) < decaytime))
      {
        double r = r = sqrt(pow((store40.Impx[k] - store40.betax[l]), 2) + pow((store40.Impy[k] - store40.betay[l]), 2));
        YSOhtspacefilt->Fill(store40.betax[l], store40.betay[l]);
        if (r < .25)
        {
          // cout << "made it " << store40.betax[l]<< endl;
          // cout << "diff " << store40.betatime[l]-store40.Imptime[k] << endl;
          tdifflh40->Fill(store40.betatime[l] - store40.Imptime[k]);
        }
      }
    }
  }
  for (int f = 0; f < PID.ncuts; f++)
  {
    for (int l = 0; l < store100[f].betatime.size(); l++)
    {
      for (int k = 0; k < store100[f].Imptime.size(); k++)
      {
        if ((abs(store100[0].betatime[l] - store100[f].Imptime[k]) < decaytime))
        {
          double r = r = sqrt(pow((store100[f].Impx[k] - store100[0].betax[l]), 2) + pow((store100[f].Impy[k] - store100[0].betay[l]), 2));
          if (r < .25)
          {
            PID.tdifflh[f]->Fill(store100[0].betatime[l] - store100[f].Imptime[k]);
          }
        }
      }
    }
  }
  printf("r/beta corrilations are done\n");

  PID.Write_cuts();

  rawEn->Write();

  for (int t = 0; t < (entries.size()); t++)
  {
    TOFHist[t]->Write();
    printf("\rstoring TOF %i\n", t);
  }
  // TCanvas * canvas = new TCanvas();
  // canvas->divid(3, 3);

  // canvas->cd(1);
  // canvas->cd(1)->setLogz();
  // //hist1->Draw();
  // canvas->cd(2)
  // //hist2->Draw();

  // canvas->Write();

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
  cout << "tdifflh40:" << (tdifflh40->Integral()) << endl;

  vetoEn->Write();
  vetofgated->Write();
  vetordiff->Write();
  vetofdiff->Write();

  PID100->Write();
  PID100Imp->Write();
  PID100ImpBox->Write();
  cout << "Mean of TOF,PID100: " << PID100ImpBox->GetMean(1) << endl;
  cout << "Standard error of TOF,PID100: " << PID100ImpBox->GetMean(13) << endl;
  PID40->Write();
  PID40Imp->Write();
  PID40ImpBox->Write();

  hittdiff->Write();
  gamaddback->Write();
  gamEn->Write();

  f2->Close();
  printf("Everything is done\n");
  return 0; // end of function main
}