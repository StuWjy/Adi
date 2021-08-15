#include "adi-satellite-container.h"
#include "adi-station-container.h"
#include "adi-link-helper.h"
#include "adi-constant.h"
#include "adi-util.h"
#include <ctime>

#define ISL

using namespace adi;

int main ()
{
  clock_t cpuStart, cpuStop;
  cpuStart = clock ();
  PointingBound bound = PointingBound {
    Bound {adi::DegreesToRadians (-70.0), adi::DegreesToRadians (70.0)},
    Bound {adi::DegreesToRadians (-70.0), adi::DegreesToRadians (70.0)}
  };
  DateTime epoch = DateTime (2022, 1, 1, 0, 0, 0);
  DateTime start = DateTime (2022, 1, 1, 0, 0, 0);
  DateTime stop = DateTime (2023, 1, 1, 0, 0, 0);
  Interval interval;
  double sma = K_RE + 500.0;
  double ecc = 0.0;
  double inc = 97.406 * M_PI / 180.0;
  double raan = 85.0 * M_PI / 180.0;
  double aop = 0.0;
  double ma = 35.0 * M_PI / 180.0;
  size_t N_PLANE = 3;
  size_t N_SAT = 40;
  std::vector<SatelliteContainer> Constellation;
  for (size_t plane = 0;plane < N_PLANE;++plane)
  {
    Constellation.push_back (SatelliteContainer::DistributeEvenly (Satellite::Element {sma, ecc, inc, raan, aop, ma}, N_SAT));
    Constellation[plane].Install (Bottom | Left | Right | Front | Back, bound);
    raan += adi::DegreesToRadians (15);
  }
  #ifdef S2G
  size_t N_STA = 40;
  StationContainer nodes = StationContainer::Create (N_STA);
  nodes[0]->SetElement (-35.3009,  149.1234,  0.5587,  "Canberra");
  nodes[1]->SetElement (-25.7445,   28.2287,  1.3801,  "Pretoria");
  nodes[2]->SetElement (-15.8366,  -47.9065,  1.0300,  "Brasilia");
  nodes[3]->SetElement (  1.3034,  103.8208,  0.0288,  "Singapore");
  nodes[4]->SetElement (  3.1596,  101.7278,  0.0471,  "KualaLumpur");
  nodes[5]->SetElement ( 13.7627,  100.5671,  0.0032,  "Bangkok");
  nodes[6]->SetElement ( 14.5422,  121.0297,  0.0225,  "Manila");
  nodes[7]->SetElement ( 16.7888,   96.1346,  0.0242,  "Yangon");
  nodes[8]->SetElement ( 18.3130,  109.3110,  0.0180,  "Sanya");
  nodes[9]->SetElement ( 21.0324,  105.8380,  0.0184,  "HaNoi");
  nodes[10]->SetElement ( 19.3370,  -99.1978,  2.3109,  "MexicoCity");
  nodes[11]->SetElement ( 23.1167,  113.2500,  0.0434,  "Guangzhou");
  nodes[12]->SetElement ( 24.6735,   46.6263,  0.6595,  "Riyadh");
  nodes[13]->SetElement ( 24.4409,   54.3712,  0.0068,  "AbuDhabi");
  nodes[14]->SetElement ( 28.5986,   77.1897,  0.2264,  "NewDelhi");
  nodes[15]->SetElement ( 30.0682,   31.2188,  0.0315,  "Cairo");
  nodes[16]->SetElement ( 32.0917,   34.7747,  0.0171,  "Haifa");
  nodes[17]->SetElement ( 30.6667,  104.0667,  0.4850,  "Chengdu");
  nodes[18]->SetElement ( 31.1094,  121.3681,  0.0040,  "Shanghai");
  nodes[19]->SetElement ( 33.2976,   44.3103,  0.0402,  "Baghdad");
  nodes[20]->SetElement ( 33.7344,   73.1358,  0.5610,  "Islamabad");
  nodes[21]->SetElement ( 35.8045,   51.4775,  1.5842,  "Teheran");
  nodes[22]->SetElement ( 37.5822,  126.9710,  0.0462,  "Seoul");
  nodes[23]->SetElement ( 35.6563,  139.7272,  0.0307,  "Tokyo");
  nodes[24]->SetElement ( 39.8992,   32.8727,  0.9860,  "Ankara");
  nodes[25]->SetElement ( 40.1172,  116.2280,  0.0380,  "Beijing");
  nodes[26]->SetElement ( 39.0544,  125.7526,  0.3397,  "Pyongyang");
  nodes[27]->SetElement ( 38.9428,  -77.0664,  0.0982,  "Washington");
  nodes[28]->SetElement ( 41.9229,   12.4996,  0.0592,  "Rome");
  nodes[29]->SetElement ( 40.4483,   -3.6417,  0.6922,  "Madrid");
  nodes[30]->SetElement ( 45.4364,  -75.6847,  0.0587,  "Ottawa");
  nodes[31]->SetElement ( 46.9395,    7.4654,  0.5541,  "Berne");
  nodes[32]->SetElement ( 47.5093,   19.0757,  0.1061,  "Budapest");
  nodes[33]->SetElement ( 48.8502,    2.3158,  0.0362,  "Paris");
  nodes[34]->SetElement ( 50.4442,   30.5407,  0.1962,  "Kyiv");
  nodes[35]->SetElement ( 52.5131,   13.4172,  0.0387,  "Berlin");
  nodes[36]->SetElement ( 52.2507,   21.0024,  0.1057,  "Warsaw");
  nodes[37]->SetElement ( 51.5211,   -0.1457,  0.0292,  "London");
  nodes[38]->SetElement ( 55.7109,   37.5164,  0.1789,  "Moscow");
  nodes[39]->SetElement ( 64.1500,  -21.9500,  0.0390,  "Reykjavik");

  nodes.Install (Top, bound);
  LinkHelper linkHelper;
  for (uint32_t plane = 0;plane < N_PLANE;++plane)
  {
    for (uint32_t sat = 0;sat < N_SAT;++sat)
    {
      for (uint32_t sta = 0;sta < N_STA;++sta)
      {
        linkHelper.AddLink (
          Link (
            Constellation[plane][sat]->GetTurntable (Bottom),
            nodes[sta]->GetTurntable (Top),
            SRC2DST | DST2SRC,
            SRC_DAY | DST_DAY | BEYOND_DISTANCE
          )
        );
      }
    }
  }
  for (uint32_t day = 0;day < 365;++day)
  {
    interval.SetStart (start);
    interval.SetStop (start + Day);
    linkHelper.SetInterval (interval);
    LinkInfoList& linkDatas = linkHelper.CalcLink ();
    for (size_t i = 0;i < linkDatas.size ();i++)
    {
      DateTime start = linkDatas[i].linkDatas.front ().time;
      DateTime stop  = linkDatas[i].linkDatas.back ().time;
      std::cout << linkDatas[i].src->GetObject ()->GetName () << " to " << linkDatas[i].dst->GetObject ()->GetName () << " " << Interval (start, stop) << std::endl;
    }
    start = start + Day;
  }
  #endif
  #ifdef ISL
  LinkHelper linkHelper;
  linkHelper.AddLink (
    Link (
      Constellation[0][0]->GetTurntable (Front),
      Constellation[0][1]->GetTurntable (Back),
      SRC2DST | DST2SRC,
      SRC_DAY | DST_DAY | BEYOND_DISTANCE
    )
  );
  interval.SetStart (start);
  interval.SetStop (start + Day);
  linkHelper.SetInterval (interval);
  LinkInfoList& linkDatas = linkHelper.CalcLink ();
  for (const LinkInfo& link : linkDatas)
  {
    std::cout << link;
  }
  #endif
  cpuStop = clock ();
  std::cout << (double)(cpuStop - cpuStart) / CLOCKS_PER_SEC << "sec";
  return 0;
}