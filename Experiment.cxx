void Experiment (TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
   gSystem->CompileMacro("Detector.cxx",opt.Data());
  gSystem->CompileMacro("Point.cxx",opt.Data());
  gSystem->CompileMacro("Vertex.cxx",opt.Data());
  gSystem->CompileMacro("Track.cxx",opt.Data());
  gSystem->CompileMacro("DoExperiment.cxx",opt.Data());
  gSystem->CompileMacro("ReconstructVertex.cxx",opt.Data());
  gSystem->CompileMacro("ResultsStudy.cxx",opt.Data());
}