#include "DGSEvent.h"

DGSEvent::DGSEvent(){
  this->Clear();
}
DGSEvent::~DGSEvent(){;}

void DGSEvent::Clear() {
  ehi = 0.0;
  DetAngle = 0.0;
  ECal = 0.0;
  id = 0;
  tpe=0, tid=0; 
  flag=0;       

  chan_id=0;  
  board_id=0;
  geo_addr=0;
  packet_length=0;

  header_type=0;
  event_type=0;
  header_length=0;

  event_timestamp=0;
  last_disc_timestamp=0;
  peak_timestamp=0;

  timestamp_match_flag=0;
  external_disc_flag=0;
  cfd_valid_flag=0;
  pileup_only_flag=0;  
  offset_flag=0;
  sync_error_flag=0; 
  general_error_flag=0; 

  peak_valid_flag=0;

  pileup_flag =0;

  sampled_baseline=0;
  cfd_sample_0=0;
  cfd_sample_1=0;
  cfd_sample_2=0;
  pre_rise_energy=0;
  post_rise_energy=0;

  m1_begin_sample=0;
  m1_end_sample=0;
  m2_begin_sample=0;
  m2_end_sample=0;
  peak_sample=0;
  base_sample=0;

  baseline=0;

  traceLen=0;
}
