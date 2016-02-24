#ifndef __DGSEvent__
#define __DGSEvent__
#include "TObject.h"

enum Channels {
  GESignal,
  GEReset,
  GRFPSignal,
  BGOSignal,
  TestSignal,
  POBSignal
};


class DGSEvent : public TObject {
 public:
   DGSEvent();
  ~DGSEvent();
  float                   ehi, DetAngle;            // Calculated Energy (Fang add angle for doppler correction)
  float                   ECal;
  short int               id;
  unsigned short int      tpe, tid; // Ge/BGO/Si etc  
  unsigned short int      flag;           // clean/dirty flag = 0/1 
  
  unsigned short int      chan_id;  
  unsigned short int      board_id;
  unsigned short int      geo_addr;
  unsigned short int      packet_length;
  
  unsigned short int      header_type;
  unsigned short int      event_type;
  unsigned short int      header_length;

  //34
  
  unsigned long long int  event_timestamp;
  unsigned long long int  last_disc_timestamp;
  unsigned long long int  peak_timestamp;

  //58
  
  unsigned short int      timestamp_match_flag;
  unsigned short int      external_disc_flag;
  unsigned short int      cfd_valid_flag;
  unsigned short int      pileup_only_flag;  
  unsigned short int      offset_flag;
  unsigned short int      sync_error_flag; 
  unsigned short int      general_error_flag; 

  unsigned short int      peak_valid_flag;

  unsigned short int      pileup_flag ;
  
  int                     sampled_baseline;
  int                     cfd_sample_0;
  int                     cfd_sample_1;
  int                     cfd_sample_2;
  int                     pre_rise_energy;
  int                     post_rise_energy;

  //100
  
  unsigned short int      m1_begin_sample;
  unsigned short int      m1_end_sample;
  unsigned short int      m2_begin_sample;
  unsigned short int      m2_end_sample;
  unsigned short int      peak_sample;
  unsigned short int      base_sample;
  
  int                     baseline;
  
  unsigned short int      traceLen;
  //short int               trace[MAXTRACELEN];

  //120

  void Clear();
  
  ClassDef(DGSEvent, 1);
};


#endif //__DGSEvent__
