#include "../include/CorrectionObject.h"
#include <TString.h>
#include <TFile.h>
// #include "../../include/helper_fkt.h"

CorrectionObject::CorrectionObject(const TString & runnr, const TString & generator, const TString & collection, const TString & input_path, const TString & input_path_MC, const TString & weight_path, const bool & closuretest,const bool & trigger_fwd,const bool & trigger_central):
  _runnr(runnr), _collection(collection), _generator(generator),_input_path(input_path),_input_path_MC(input_path_MC),_weight_path(weight_path), _closuretest(closuretest), _trigger_fwd(trigger_fwd), _trigger_central(trigger_central)
    {
      TString s_tmp = _collection;
      s_tmp.ReplaceAll("CHS", "PFchs");
      s_tmp.ReplaceAll("Puppi", "PFpuppi");
      _jettag = s_tmp;
      TString inputPath;
      TString inputPath_MC;
      
      //Declare path for input, weights, output
      inputPath  = _input_path;//+"/";
      inputPath_MC  = _input_path_MC;
      _weightpath_FLAT = _weight_path;//+"/CENTRAL/";
      _weightpath_FWD  = _weight_path;//+"/FWD/";
     _weightpath      = _weight_path;

     TString input_base = inputPath;
     input_base.Resize(input_base.Last('/')+1);
     _outpath    = input_base +"Run" + _runnr + "/";

          
// "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS.root
     //FIXME hardcoded 2016 MC paths
     //For QCD pT binned samples 
      if(_generator == "pythia"){
	if(_trigger_central && !_trigger_fwd)     { _MCpath = inputPath_MC;}
	else if(!_trigger_central && _trigger_fwd){ _MCpath = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS.root";}
       	//else if(_trigger_central && _trigger_fwd) {_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_" + _collection  +".root";}
       	else if(_trigger_central && _trigger_fwd)  {_MCpath = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS.root";}
	else throw runtime_error("In Correction Object: No valid Trigger-Flag (main.C) was set.");
	
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_" + _collection + ".root";//  "_Flat.root";
	_MCpath_ForWeights_FWD  = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_" + _collection  + ".root";// "_Fwd.root";
	_MCpath_ForWeights  = _weightpath + "uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_" + _collection  + ".root";
	_generator_tag = "pythia8";
      }
     
      //For flat MC samples:
      /*
      if(_generator == "pythia"){
	if(_trigger_central && !_trigger_fwd)     { _MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Flat.root";}
	else if(!_trigger_central && _trigger_fwd){ _MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Fwd.root";}
	//else if(_trigger_central && _trigger_fwd)  {_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  +".root";}
	else if(_trigger_central && _trigger_fwd)  {_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Full_Run" + _runnr  + ".root";}
	else throw runtime_error("In Correction Object: No valid Trigger-Flag (main.C) was set.");
	
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD  = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Fwd.root";
	_MCpath_ForWeights  = _weight_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + ".root";
	_generator_tag = "pythia8";
      }
      */

      else if(_generator == "herwig"){
	_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_"+ _collection  +".root";
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_" + _collection  + "_Fwd.root";
	_generator_tag = "herwigpp";
      }
      else if(_generator == "madgraph"){
	_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_"+ _collection  +".root";
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_" + _collection  + "_Fwd.root";
	_generator_tag = "madgraphMLM";
      }

      //DATA 
      _DATApath = input_path // + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root"
	;
      _DATApath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      _DATApath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      _DATApath_ForWeights = _weight_path + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      
      //Check if files are in place:
      cout << "Opening MC file:   " << _MCpath << endl;
      cout << "Opening DATA file: " << _DATApath << endl << endl;
      cout << "Output path is: " <<_outpath << endl;
      if(!  CorrectionObject::make_path(std::string(_outpath.Data()))){
      	cout << "No new Directory was created" << endl;
      }
      cout << endl;
      
      _MCFile = new TFile(_MCpath,"READ");
      _DATAFile = new TFile(_DATApath,"READ");

      if(_MCFile->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _MCpath+" does not exist!");
      if(_DATAFile->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _DATApath+" does not exist!");


      //lumitags
      if(_runnr == "BCD")    _lumitag      = "RunBCD  12.9 fb^{-1}";
      else if(_runnr == "G") _lumitag      = "RunG  7.6 fb^{-1}";
      else if(_runnr == "B") _lumitag      = "RunB  TODO fb^{-1}";
      //FIXME differentiate between 2016 (5.8 fb^{-1}) and 2017
      else if(_runnr == "C") _lumitag      = "RunC  2.6 fb^{-1}";
      else if(_runnr == "D") _lumitag      = "RunD  4.3 fb^{-1}";
      else if(_runnr == "E") _lumitag      = "RunE  4.1 fb^{-1}";
      else if(_runnr == "F") _lumitag      = "RunF  3.2 fb^{-1}";
      else if(_runnr == "G") _lumitag      = "RunG  7.6 fb^{-1}";
      else if(_runnr == "Fearly") _lumitag = "RunFearly  2.7 fb^{-1}";
      else if(_runnr == "EFearly") _lumitag = "RunEFearly  6.8 fb^{-1}";
      else if(_runnr == "Flate") _lumitag = "RunFlate  0.4 fb^{-1}";
      else if(_runnr == "FlateG") _lumitag = "RunFlateG  8.0 fb^{-1}";
      else if(_runnr == "FlateGH") _lumitag = "RunFlateGH  X.X fb^{-1}";
      else if(_runnr == "H") _lumitag = "RunH  8.5 fb^{-1}";
      else if(_runnr == "BCDEFearly") _lumitag = "RunBCDEFearly  19.7 fb^{-1}";
      else if(_runnr == "BCDEFGH") _lumitag = "RunBCDEFGH  36.8 fb^{-1}";
      else throw runtime_error("In constructor: Invalid RunNr. specified.");
    }

void CorrectionObject::FullCycle_CorrectFormulae(){
  std::cout<<"Doing the Full Cycle"<<std::endl;
  std::cout<<"\nStarting ControlPlots()\n"<<std::endl;
  CorrectionObject::ControlPlots();  
  std::cout<<"\nStarting kFSR_CorrectFormulae()\n"<<std::endl;
  CorrectionObject::kFSR_CorrectFormulae();
  std::cout<<"\nStarting Pt_Extrapolation_Alternative_CorrectFormulae(true)\n"<<std::endl;
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae(true);
  std::cout<<"\nStarting Pt_Extrapolation_Alternative_CorrectFormulae(false)\n"<<std::endl;
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae(false);
  std::cout<<"\nStarting L2ResOutput()\n"<<std::endl;
  CorrectionObject::L2ResOutput();
  std::cout<<"\nStarting FinalControlPlots_CorrectFormulae()\n"<<std::endl;
  CorrectionObject::FinalControlPlots_CorrectFormulae();
}

//Full cycle to calculate L2Res with extended eta range. negative Values
void CorrectionObject::FullCycle_CorrectFormulae_eta(){
  CorrectionObject::ControlPlots();
  CorrectionObject::kFSR_CorrectFormulae_eta();
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae_eta(true);
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae_eta(false);
  CorrectionObject::L2ResOutput_eta();
  CorrectionObject::FinalControlPlots_CorrectFormulae_eta();
}

