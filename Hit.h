#include <vector>


struct HIT
{
vector<int> ID;
vector<int> energy;
vector<double> time;

  
  int ID;

  void reset()
  {
    ID.clear();
    energy.clear();
    time.clear();
  }

  void add_hit(int ID, int energy, double time)
  {
    ID.push_back(ID);
    energy.push_back(energy);
    time.push_back(time);  
    int nhits = ID.size();
  }
  

};


