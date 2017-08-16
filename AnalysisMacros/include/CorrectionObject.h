#ifndef CORRECTIONOBJECT_H
#define CORRECTIONOBJECT_H

#pragma once

#include <cmath>
#include <iostream>
#include <TString.h>
#include <TFile.h>

#include <string>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

  class CorrectionObject {



  public:

    // Constructors, destructor
    CorrectionObject(const TString & runnr, const TString & generator, const TString & collection,const TString & input_path, const TString & input_path_MC, const TString & weight_path, const bool & closuretest = false,const bool & trigger_fwd = false,const bool & trigger_central = false);
    CorrectionObject(const CorrectionObject &) = default;
    CorrectionObject & operator = (const CorrectionObject &) = default;
    ~CorrectionObject() = default;
    inline void CloseFiles(){_MCFile->Close(); _DATAFile->Close();};
    
    // Setter and getter functions
    inline TString runnr(){ return _runnr;}
    inline TString collection(){ return _collection;}
    inline TString input_path(){return _input_path;}
    inline TString input_path_MC(){return _input_path_MC;}    
    inline TString weight_path(){return _weight_path;}
    inline TString generator(){ return _generator;}
    inline TString jettag(){ return _jettag;}
    inline bool closuretest(){return _closuretest;}
    inline bool trigger_fwd(){return _trigger_fwd;}
    inline bool trigger_central(){return _trigger_central;}
    inline TString MCPath(){return _MCpath;}
    inline TString DATAPath(){return _DATApath;}
    inline TString OutPath(){return _outpath;}
    inline TString lumitag(){return _lumitag;}
    inline const TString runnr() const{return _runnr;}
    inline const TString collection() const{return _collection;}
    inline const TString input_path()const{return _input_path;}
    inline const TString input_path_MC()const{return _input_path_MC;}
    inline const TString weight_path()const{return _weight_path;}
    inline const TString generator() const{return _generator;}
    inline const TString jettag() const{return _jettag;}
    inline const bool closuretest() const {return _closuretest;}
    inline const bool trigger_fwd() const {return _trigger_fwd;}
    inline const bool trigger_central() const {return _trigger_central;}
    inline const TString MCPath() const {return _MCpath;}
    inline const TString DATAPath() const {return _DATApath;}
    inline const TString OutPath() const {return _outpath;}
    inline const TString lumitag() const {return _lumitag;}
    inline void set_runnr(TString x){_runnr = x;}
    inline void set_collection(TString x){_collection = x;}
    inline void set_input_path(TString x){_input_path = x;}
    inline void set_weight_path(TString x){_weight_path = x;}
    inline void set_generator(TString x){_generator = x;}
    inline void set_jettag(TString x){_jettag = x;}
    inline void set_closuretest(bool x){_closuretest = x;}
    inline void set_trigger_fwd(bool x){_trigger_fwd = x;}
    inline void set_trigger_central(bool x){_trigger_central = x;}
    inline void set_MCPath(TString x){_MCpath = x; _MCFile->Close(); _MCFile = new TFile(_MCpath,"READ");}
    inline void set_DATAPath(TString x){_DATApath = x; _DATAFile->Close(); _DATAFile = new TFile(_DATApath,"READ");}
    inline void set_outpath(TString x){_outpath = x;}
    inline void set_lumitag(TString x){_lumitag = x;}

    //Main functions for calculating L2 residuals, defined in CorrectionObject.cc
    void ControlPlots();
    void kFSR();
    void kFSR_CorrectFormulae();
    void kFSR_CorrectFormulae_eta();
    void Pt_Extrapolation(bool mpfMethod = true);
    void Pt_Extrapolation_Alternative(bool mpfMethod = true);
    void Pt_Extrapolation_Alternative_CorrectFormulae(bool mpfMethod = true);
    void Pt_Extrapolation_Alternative_CorrectFormulae_eta(bool mpfMethod = true);
    void L2ResOutput();
    void L2ResOutput_eta();
    void L2ResAllRuns();
    void L2Res_JEC();
    void L2ResOverlay(bool is_MPF);
    void InputForGlobalFit();
    void InputForGlobalFit_eta_0_13();
    void FinalControlPlots();
    void FinalControlPlots_CorrectFormulae();
    void FinalControlPlots_CorrectFormulae_eta();
    void CalculateMCWeights();
    void CalculateMCWeights_TriggerThresholds(bool CentralTriggers);
    void FullCycle_CorrectFormulae();
    void FullCycle_CorrectFormulae_eta();
    void MatchingPlots();
    
static bool make_path(std::string path_name){
    int status = 1;

    std::string path_base = "";
    std::string dir_name ="";
    std::string path_remain = path_name;

    size_t last_pos = path_remain.find("/",1);
    while(last_pos != std::string::npos){
      dir_name = path_remain.substr(0,last_pos);
      path_remain = path_remain.substr(last_pos);

      path_base += dir_name;
      status *=  mkdir(path_base.c_str() , S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

      last_pos = path_remain.find("/",1);
    }

    return status!=0;
  }

static bool make_path(TString path_name){
  return CorrectionObject::make_path(std::string(path_name.Data()));
  }

  private:
    TString _runnr;
    TString _collection;
    TString _generator, _generator_tag;
    TString _input_path;
    TString _input_path_MC;   
    TString _weight_path;
    TString _jettag;
    TString _lumitag;
    TString _MCpath, _MCpath_ForWeights, _MCpath_ForWeights_FLAT, _MCpath_ForWeights_FWD, _DATApath, _DATApath_ForWeights, _DATApath_ForWeights_FLAT, _DATApath_ForWeights_FWD;
    TString _outpath;
    TString _weightpath_FLAT, _weightpath_FWD;
    TFile*  _MCFile;
    TFile*  _DATAFile;
    TString _weightpath;
    bool    _closuretest;
    bool    _trigger_fwd;
    bool    _trigger_central;
   

  }; // end of class CorrectionObject



// Equality operators

inline bool operator == (const CorrectionObject & a, const CorrectionObject & b) {
  return a.runnr() == b.runnr() && a.collection() == b.collection() && a.generator() == b.generator() && a.closuretest() == b.closuretest();
}


inline bool operator != (const CorrectionObject & a, const CorrectionObject & b) {
  return a.runnr() != b.runnr() || a.collection() != b.collection() || a.generator() != b.generator();
}

// I/O operations

inline ostream & operator << (ostream & os, const CorrectionObject & q) {
  return os << "(" << q.runnr() << "," << endl 
                  << q.generator() << ","  << endl
                  << q.collection() << "," << endl 
                  << q.closuretest() << "," << endl
                  << q.MCPath() << ","  << endl 
	          << q.DATAPath() << "," << endl
	          << q.OutPath() << ")";
}


#endif
