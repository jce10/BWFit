#include "SpectrumCalibration.h"
#include "LineShapes.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <TCanvas.h>
#include <TDecompChol.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVectorD.h>

namespace bwfit {
namespace {
using Parameter = CalibrationParameter;
std::map<std::string, CalibrationShape>& Shapes() {
  static std::map<std::string, CalibrationShape> shapes{
    {"gaussian", [](double x, double m, double w) { return Gaussian(x,m,w,1); }},
    {"bw", [](double x, double m, double w) { return BreitWigner(x,m,w,1); }}
  };
  return shapes;
}
void Require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error("Calibration: " + message);
}
bool Finite(double x) { return std::isfinite(x); }
void ValidateParameter(const Parameter& p, const std::string& name) {
  Require(Finite(p.value) && Finite(p.min) && Finite(p.max) && p.min <= p.value &&
          p.value <= p.max, "invalid value/bounds for " + name);
}
void Validate(const CalibrationSpec& s) {
  Require(Finite(s.xmin_mm) && Finite(s.xmax_mm) && s.xmin_mm < s.xmax_mm &&
          s.bins > 0, "invalid raw histogram range/bins");
  Require(s.order == 1 || s.order == 2, "ORDER must be 1 or 2");
  Require(s.background_order >= -1 && s.background_order <= 2,
          "BACKGROUND_ORDER must be -1 (none), 0, 1 or 2");
  Require(s.background.size() == static_cast<size_t>(s.background_order+1),
          "supply one BG row per background coefficient");
  for (const auto& p : s.background) {
    ValidateParameter(p, "background");
    Require(p.min >= 0, "background density must be nonnegative");
  }
  int anchors = 0;
  std::set<std::string> names;
  for (const auto& p : s.peaks) {
    Require(!p.name.empty() && names.insert(p.name).second, "peak names must be unique");
    Require(p.name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")==std::string::npos,
            "peak names must use letters, digits, underscore or hyphen");
    Require(Shapes().count(p.shape), "unknown shape: " + p.shape);
    ValidateParameter(p.yield, p.name+" yield");
    ValidateParameter(p.centroid, p.name+" centroid");
    ValidateParameter(p.width, p.name+" width");
    Require(p.yield.min >= 0 && p.width.min > 0, "yield >= 0 and width > 0 required");
    Require(p.centroid.min >= s.xmin_mm && p.centroid.max <= s.xmax_mm,
            "centroid bounds outside fit range: " + p.name);
    Require(Finite(p.energy_MeV) && Finite(p.energy_error_MeV) &&
            p.energy_error_MeV >= 0, "invalid reference energy/error");
    if (p.anchor) {
      ++anchors;
      Require(p.centroid.min < p.centroid.max, "anchor centroid must float: " + p.name);
    }
  }
  Require(anchors >= s.order+2, "need at least ORDER+2 anchors for residual checks");
  const int i = s.interference_first, j = s.interference_second;
  if (i != -1 || j != -1) {
    Require(i >= 0 && j >= 0 && i != j && i < static_cast<int>(s.peaks.size()) &&
            j < static_cast<int>(s.peaks.size()), "invalid interfering peak indices");
    Require(s.peaks[i].shape == "bw" && s.peaks[j].shape == "bw",
            "interference requires two nonrelativistic bw components");
    ValidateParameter(s.phase, "phase");
  }
}
int BackgroundIndex(const CalibrationSpec& s) { return 3*static_cast<int>(s.peaks.size()); }
int PhaseIndex(const CalibrationSpec& s) { return BackgroundIndex(s)+s.background_order+1; }
int ParameterCount(const CalibrationSpec& s) {
  return PhaseIndex(s)+(s.interference_first >= 0 ? 1 : 0);
}
double Background(double x, const double* p, const CalibrationSpec& s) {
  const double t = (x-s.xmin_mm)/(s.xmax_mm-s.xmin_mm);
  const int b = BackgroundIndex(s);
  switch (s.background_order) {
    case 0: return p[b];
    case 1: return p[b]*(1-t)+p[b+1]*t;
    case 2: return p[b]*(1-t)*(1-t)+2*p[b+1]*t*(1-t)+p[b+2]*t*t;
    default: return 0;
  }
}
// TF1 gives counts per bin; with ROOT option I the bin-averaged function
// equals the integral of the density over each equal-width bin.
class RawModel {
  CalibrationSpec spec_;
  std::vector<CalibrationShape> shapes_;
  double bin_width_;
public:
  RawModel(const CalibrationSpec& s, double w): spec_(s), bin_width_(w) {
    for (const auto& peak : s.peaks) shapes_.push_back(Shapes().at(peak.shape));
  }
  double operator()(double* x, double* p) const {
    double y = Background(x[0], p, spec_);
    for (size_t k=0; k<shapes_.size(); ++k)
      y += p[3*k]*shapes_[k](x[0],p[3*k+1],p[3*k+2]);
    if (spec_.interference_first >= 0) {
      const int i=3*spec_.interference_first, j=3*spec_.interference_second;
      const double X=(x[0]-p[i+1])*(x[0]-p[j+1])+p[i+2]*p[j+2]/4;
      const double Y=p[i+2]/2*(x[0]-p[j+1])-p[j+2]/2*(x[0]-p[i+1]);
      const double phase=p[PhaseIndex(spec_)];
      y += std::sqrt(p[i]*p[j]*p[i+2]*p[j+2])/std::acos(-1.0)*
           (X*std::cos(phase)+Y*std::sin(phase))/(X*X+Y*Y);
    }
    return y*bin_width_;
  }
};
void SetParameter(TF1& f, int i, const Parameter& p, const std::string& name) {
  f.SetParName(i,name.c_str()); f.SetParameter(i,p.value);
  if (p.min == p.max) f.FixParameter(i,p.value);
  else f.SetParLimits(i,p.min,p.max);
}
TMatrixDSym Inverse(const TMatrixDSym& a) {
  TDecompChol factor(a);
  Bool_t ok = false;
  TMatrixDSym inv = factor.Invert(ok);
  Require(ok, "singular covariance/design matrix; check anchors and their span");
  return inv;
}
void FitPolynomial(CalibrationResult& r, const CalibrationSpec& s) {
  std::vector<int> index;
  for (size_t k=0;k<r.peaks.size();++k) if(r.peaks[k].input.anchor) index.push_back(k);
  const int n=index.size(), m=s.order+1;
  double lo=r.peaks[index[0]].centroid_mm, hi=lo;
  for (int k:index) { lo=std::min(lo,r.peaks[k].centroid_mm); hi=std::max(hi,r.peaks[k].centroid_mm); }
  Require(hi>lo,"anchors have no position span");
  const double x0=(lo+hi)/2, scale=(hi-lo)/2;
  TMatrixD design(n,m);
  TVectorD energy(n);
  TMatrixDSym cx(n);
  for(int i=0;i<n;++i) {
    const double z=(r.peaks[index[i]].centroid_mm-x0)/scale;
    design(i,0)=1; design(i,1)=z; if(m==3) design(i,2)=z*z;
    energy[i]=r.peaks[index[i]].input.energy_MeV;
    for(int j=0;j<n;++j) cx(i,j)=r.spectrum_covariance(3*index[i]+1,3*index[j]+1);
  }
  TMatrixD transpose(TMatrixD::kTransposed,design);
  TVectorD beta(m); beta.Zero();
  TMatrixDSym weights(n), covariance(m);
  bool converged=false;
  for(int iteration=0;iteration<100;++iteration) {
    TMatrixDSym v(n);
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) {
      const double di=(beta[1]+(m==3 ? 2*beta[2]*design(i,1) : 0))/scale;
      const double dj=(beta[1]+(m==3 ? 2*beta[2]*design(j,1) : 0))/scale;
      v(i,j)=iteration==0 ? (i==j ? 1.0 : 0.0) : di*cx(i,j)*dj;
      if(iteration>0 && i==j) v(i,j)+=std::pow(r.peaks[index[i]].input.energy_error_MeV,2);
    }
    weights=Inverse(v);
    TMatrixD normal=transpose*TMatrixD(weights)*design;
    TMatrixDSym ns(m);
    for(int i=0;i<m;++i) for(int j=0;j<m;++j) ns(i,j)=normal(i,j);
    covariance=Inverse(ns);
    TVectorD next=TMatrixD(covariance)*(transpose*(TMatrixD(weights)*energy));
    double change=0;
    for(int i=0;i<m;++i) change=std::max(change,std::abs(next[i]-beta[i]));
    beta=next;
    if(iteration>0 && change<1e-11) { converged=true; break; }
  }
  Require(converged,"effective-variance calibration did not converge");
  TMatrixD transform(3,m); transform.Zero();
  transform(0,0)=1; transform(0,1)=-x0/scale; transform(1,1)=1/scale;
  if(m==3) {
    transform(0,2)=x0*x0/(scale*scale); transform(1,2)=-2*x0/(scale*scale);
    transform(2,2)=1/(scale*scale);
  }
  TVectorD a=transform*beta;
  TMatrixD tt(TMatrixD::kTransposed,transform);
  TMatrixD ca=transform*TMatrixD(covariance)*tt;
  for(int i=0;i<3;++i) {
    r.a[i]=a[i];
    for(int j=0;j<3;++j) r.coefficient_covariance(i,j)=ca(i,j);
  }
  const double dmin=r.a[1]+2*r.a[2]*s.xmin_mm;
  const double dmax=r.a[1]+2*r.a[2]*s.xmax_mm;
  Require(dmin*dmax>0,"calibration is not monotonic over the raw fit range");
  const TVectorD residual=energy-design*beta;
  r.calibration_chi2=residual*(TMatrixD(weights)*residual);
  r.calibration_ndf=n-m; r.order=s.order;
  for(auto& p:r.peaks) p.residual_keV=1000*(p.input.energy_MeV-r.Energy(p.centroid_mm));
}
std::string Quote(const std::string& value) {
  std::ostringstream o; o << '"';
  for(unsigned char c:value) {
    if(c=='"' || c=='\\') o << '\\' << c;
    else if(c<32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else o << c;
  }
  o << '"'; return o.str();
}
void MatrixJSON(std::ostream& out, const TMatrixDSym& a) {
  out << '[';
  for(int i=0;i<a.GetNrows();++i) {
    if(i) out << ',';
    out << '[';
    for(int j=0;j<a.GetNcols();++j) { if(j) out << ','; out << a(i,j); }
    out << ']';
  }
  out << ']';
}
} // namespace

double CalibrationResult::Energy(double x) const { return a[0]+x*(a[1]+x*a[2]); }
double CalibrationResult::EnergyError(double x) const {
  const double v[]={1,x,x*x}; double variance=0;
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) variance+=v[i]*coefficient_covariance(i,j)*v[j];
  return std::sqrt(std::max(0.0,variance));
}
void RegisterCalibrationShape(const std::string& name, CalibrationShape density) {
  Require(!name.empty() && bool(density),"empty custom shape/name");
  Require(!Shapes().count(name),"shape name already registered: "+name);
  Shapes().emplace(name,std::move(density));
}
CalibrationSpec ReadCalibrationSpec(const std::string& path) {
  std::ifstream input(path); Require(bool(input),"cannot open card: "+path);
  CalibrationSpec s;
  std::set<std::string> seen;
  std::string line; int number=0;
  while(std::getline(input,line)) {
    ++number; line=line.substr(0,line.find('#'));
    std::istringstream row(line); std::string key;
    if(!(row>>key)) continue;
    if(key!="PEAK" && key!="BG") Require(seen.insert(key).second,"duplicate key: "+key);
    if(key=="ROOTFILE") row>>s.rootfile;
    else if(key=="TREE") row>>s.tree;
    else if(key=="BRANCH") row>>s.branch;
    else if(key=="HISTOGRAM") row>>s.histogram;
    else if(key=="OUTPUT_PREFIX") row>>s.output_prefix;
    else if(key=="RANGE_MM") row>>s.xmin_mm>>s.xmax_mm;
    else if(key=="BINS") row>>s.bins;
    else if(key=="ORDER") row>>s.order;
    else if(key=="BACKGROUND_ORDER") row>>s.background_order;
    else if(key=="STATISTIC") {
      std::string stat; row>>stat;
      Require(stat=="poisson" || stat=="chi2","STATISTIC must be poisson or chi2");
      s.poisson=stat=="poisson";
    } else if(key=="BG") {
      Parameter p; row>>p.value>>p.min>>p.max; s.background.push_back(p);
    } else if(key=="INTERFERENCE") {
      row>>s.interference_first>>s.interference_second>>s.phase.value>>s.phase.min>>s.phase.max;
    } else if(key=="PEAK") {
      CalibrationPeak p; int anchor=0;
      row>>p.name>>p.shape>>anchor>>p.energy_MeV>>p.energy_error_MeV
         >>p.yield.value>>p.yield.min>>p.yield.max
         >>p.centroid.value>>p.centroid.min>>p.centroid.max
         >>p.width.value>>p.width.min>>p.width.max;
      Require(anchor==0 || anchor==1,"anchor flag must be 0 or 1");
      p.anchor=anchor==1; s.peaks.push_back(p);
    } else Require(false,"unknown card key: "+key);
    Require(!row.fail(),"malformed card row "+std::to_string(number));
    std::string extra; Require(!(row>>extra),"extra tokens on card row "+std::to_string(number));
  }
  Require(!s.rootfile.empty(),"ROOTFILE is required");
  Validate(s); return s;
}
std::unique_ptr<TH1> LoadRawSpectrum(const CalibrationSpec& s) {
  Validate(s);
  TFile f(s.rootfile.c_str(),"READ"); Require(!f.IsZombie(),"cannot open ROOT file: "+s.rootfile);
  if(!s.histogram.empty()) {
    auto* h=dynamic_cast<TH1*>(f.Get(s.histogram.c_str()));
    Require(h && h->GetDimension()==1,"missing one-dimensional histogram: "+s.histogram);
    auto copy=std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone("raw_calibration")));
    copy->SetDirectory(nullptr); return copy;
  }
  auto* tree=dynamic_cast<TTree*>(f.Get(s.tree.c_str())); Require(tree,"missing tree: "+s.tree);
  auto* leaf=tree->GetLeaf(s.branch.c_str());
  Require(leaf && leaf->GetLenStatic()==1 && !leaf->GetLeafCount(),"BRANCH must be a scalar numeric leaf");
  const std::string type=leaf->GetTypeName();
  Require(type=="Double_t" || type=="Float_t" || type=="Int_t" || type=="UInt_t" ||
          type=="Long64_t" || type=="ULong64_t" || type=="Short_t" || type=="UShort_t",
          "unsupported scalar leaf type: "+type);
  auto h=std::make_unique<TH1D>("raw_calibration","",s.bins,s.xmin_mm,s.xmax_mm);
  h->SetDirectory(nullptr);
  for(Long64_t i=0;i<tree->GetEntries();++i) {
    Require(tree->GetEntry(i)>0,"failed reading tree entry");
    const double x=leaf->GetValue(); Require(Finite(x),"nonfinite raw position"); h->Fill(x);
  }
  return h;
}
CalibrationResult FitRawCalibrationSpectrum(TH1& h, const CalibrationSpec& s) {
  Validate(s); Require(h.GetDimension()==1,"expected 1D histogram");
  const double w=h.GetBinWidth(1);
  Require(std::abs(h.GetXaxis()->GetXmin()-s.xmin_mm)<1e-8 &&
          std::abs(h.GetXaxis()->GetXmax()-s.xmax_mm)<1e-8,
          "histogram axis must match RANGE_MM (all regular bins are fitted)");
  h.GetXaxis()->SetRange(0,0);
  for(int b=1;b<=h.GetNbinsX();++b) {
    Require(std::abs(h.GetBinWidth(b)-w)<1e-8*w,"uniform bin widths required");
    const double y=h.GetBinContent(b), e=h.GetBinError(b);
    Require(Finite(y) && Finite(e) && e>=0,"nonfinite bin content/error");
    if(s.poisson) {
      Require(y>=0 && std::abs(y-std::round(y))<1e-7,
              "poisson requires unsubtracted, unweighted integer counts");
      Require(std::abs(e*e-y)<1e-5*std::max(1.0,y),"poisson bin errors indicate weighted data");
    }
  }
  Require(h.Integral()>0,"empty histogram");
  RawModel model(s,w);
  TF1 fit("raw_position_fit",model,s.xmin_mm,s.xmax_mm,ParameterCount(s));
  fit.SetNpx(std::max(2000,h.GetNbinsX()*2));
  for(size_t k=0;k<s.peaks.size();++k) {
    SetParameter(fit,3*k,s.peaks[k].yield,s.peaks[k].name+"_yield");
    SetParameter(fit,3*k+1,s.peaks[k].centroid,s.peaks[k].name+"_centroid_mm");
    SetParameter(fit,3*k+2,s.peaks[k].width,s.peaks[k].name+"_FWHM_mm");
  }
  for(size_t k=0;k<s.background.size();++k)
    SetParameter(fit,BackgroundIndex(s)+k,s.background[k],"bg_"+std::to_string(k));
  if(s.interference_first>=0) SetParameter(fit,PhaseIndex(s),s.phase,"phase_rad");
  TFitResultPtr result=h.Fit(&fit,s.poisson ? "RSNILQ" : "RSNIQ");
  Require(result.Get() && int(result)==0 && result->IsValid(),"raw fit failed; adjust seeds/bounds/model");
  Require(result->CovMatrixStatus()==3,"raw covariance not accurate; do not calibrate from this fit");
  CalibrationResult r;
  r.spectrum_status=int(result); r.covariance_status=result->CovMatrixStatus();
  r.spectrum_ndf=fit.GetNDF(); r.spectrum_chi2=fit.GetChisquare();
  Require(r.spectrum_ndf>0,"raw fit has no residual degrees of freedom");
  r.spectrum_covariance.ResizeTo(fit.GetNpar(),fit.GetNpar());
  r.spectrum_covariance=result->GetCovarianceMatrix();
  for(int k=0;k<fit.GetNpar();++k) {
    r.spectrum_parameters.push_back(fit.GetParameter(k)); r.spectrum_errors.push_back(fit.GetParError(k));
    Require(Finite(fit.GetParameter(k)) && Finite(fit.GetParError(k)),"nonfinite fitted parameter/error");
    double lo=0,hi=0; fit.GetParLimits(k,lo,hi);
    if(hi>lo && std::min(fit.GetParameter(k)-lo,hi-fit.GetParameter(k))<1e-4*(hi-lo))
      r.warnings.push_back(std::string("Parameter near bound: ")+fit.GetParName(k));
  }
  for(size_t k=0;k<s.peaks.size();++k) {
    CalibrationPeakResult p; p.input=s.peaks[k];
    p.yield=fit.GetParameter(3*k); p.yield_error=fit.GetParError(3*k);
    p.centroid_mm=fit.GetParameter(3*k+1); p.centroid_error_mm=fit.GetParError(3*k+1);
    p.fwhm_mm=fit.GetParameter(3*k+2); p.fwhm_error_mm=fit.GetParError(3*k+2);
    if(p.input.anchor) {
      Require(p.centroid_error_mm>0,"anchor centroid error is zero");
      const double edge=std::min(p.centroid_mm-p.input.centroid.min,p.input.centroid.max-p.centroid_mm);
      Require(edge>1e-4*(p.input.centroid.max-p.input.centroid.min),"anchor centroid at bound: "+p.input.name);
    }
    r.peaks.push_back(p);
  }
  // Explicit Poisson deviance evaluated using the same bin integrals as the fit.
  r.spectrum_deviance=0;
  for(int b=1;b<=h.GetNbinsX();++b) {
    const double mu=fit.Integral(h.GetXaxis()->GetBinLowEdge(b),h.GetXaxis()->GetBinUpEdge(b))/w;
    const double y=h.GetBinContent(b); Require(mu>0 && Finite(mu),"nonpositive model bin expectation");
    if(y>=0) r.spectrum_deviance+=2*(mu-y+(y>0 ? y*std::log(y/mu) : 0));
  }
  FitPolynomial(r,s); return r;
}

void WriteCalibrationReports(TH1& h, const CalibrationSpec& s, const CalibrationResult& r) {
  const std::filesystem::path prefix(s.output_prefix);
  if(prefix.has_parent_path()) std::filesystem::create_directories(prefix.parent_path());
  const std::string stem=prefix.string();
  std::ofstream json(stem+".json"), csv(stem+"_peaks.csv"), card(stem+"_coefficients.txt"), input(stem+"_input.txt");
  Require(bool(json)&&bool(csv)&&bool(card)&&bool(input),"cannot create output reports");
  json << std::setprecision(17); csv << std::setprecision(17); card << std::setprecision(17);
  input << std::setprecision(17) << "ROOTFILE " << s.rootfile << "\nTREE " << s.tree
        << "\nBRANCH " << s.branch << "\nOUTPUT_PREFIX " << s.output_prefix
        << "\nRANGE_MM " << s.xmin_mm << ' ' << s.xmax_mm << "\nBINS " << h.GetNbinsX()
        << "\nORDER " << s.order << "\nBACKGROUND_ORDER " << s.background_order
        << "\nSTATISTIC " << (s.poisson ? "poisson" : "chi2") << '\n';
  if(!s.histogram.empty()) input << "HISTOGRAM " << s.histogram << '\n';
  for(const auto& p:s.background) input << "BG " << p.value << ' ' << p.min << ' ' << p.max << '\n';
  for(const auto& p:s.peaks) input << "PEAK " << p.name << ' ' << p.shape << ' ' << p.anchor << ' '
    << p.energy_MeV << ' ' << p.energy_error_MeV << ' ' << p.yield.value << ' ' << p.yield.min << ' ' << p.yield.max << ' '
    << p.centroid.value << ' ' << p.centroid.min << ' ' << p.centroid.max << ' '
    << p.width.value << ' ' << p.width.min << ' ' << p.width.max << '\n';
  if(s.interference_first>=0) input << "INTERFERENCE " << s.interference_first << ' ' << s.interference_second << ' '
    << s.phase.value << ' ' << s.phase.min << ' ' << s.phase.max << '\n';
  json << "{\n  \"schema_version\":1,\n  \"rootfile\":" << Quote(s.rootfile)
       << ",\n  \"tree\":" << Quote(s.tree) << ",\n  \"branch\":" << Quote(s.branch)
       << ",\n  \"histogram\":" << Quote(s.histogram)
       << ",\n  \"position_unit\":\"mm\",\n  \"energy_unit\":\"MeV\",\n  \"width_convention\":\"FWHM\",\n"
       << "  \"statistic\":" << Quote(s.poisson ? "poisson" : "chi2")
       << ",\n  \"range_mm\":[" << s.xmin_mm << ',' << s.xmax_mm << "],\n"
       << "  \"spectrum_status\":" << r.spectrum_status << ",\n  \"covariance_status\":" << r.covariance_status
       << ",\n  \"spectrum_chi2_ROOT\":" << r.spectrum_chi2 << ",\n  \"spectrum_deviance\":" << r.spectrum_deviance
       << ",\n  \"spectrum_ndf\":" << r.spectrum_ndf << ",\n  \"spectrum_parameters\":[";
  for(size_t k=0;k<r.spectrum_parameters.size();++k) { if(k) json<<','; json<<r.spectrum_parameters[k]; }
  json << "],\n  \"spectrum_errors\":[";
  for(size_t k=0;k<r.spectrum_errors.size();++k) { if(k) json<<','; json<<r.spectrum_errors[k]; }
  json << "],\n  \"spectrum_covariance\":"; MatrixJSON(json,r.spectrum_covariance);
  json << ",\n  \"calibration_method\":\"iterative effective-variance GLS, first-order covariance\",\n"
       << "  \"order\":" << r.order << ",\n  \"a\":[" << r.a[0] << ',' << r.a[1] << ',' << r.a[2]
       << "],\n  \"coefficient_covariance\":"; MatrixJSON(json,r.coefficient_covariance);
  json << ",\n  \"calibration_chi2\":" << r.calibration_chi2 << ",\n  \"calibration_ndf\":" << r.calibration_ndf
       << ",\n  \"interference_indices\":[" << s.interference_first << ',' << s.interference_second << "],\n  \"peaks\":[\n";
  csv << "name,shape,anchor,energy_MeV,energy_error_MeV,yield,yield_error,centroid_mm,centroid_error_mm,FWHM_mm,FWHM_error_mm,known_minus_calibrated_keV\n";
  for(size_t k=0;k<r.peaks.size();++k) {
    const auto& p=r.peaks[k]; if(k) json << ",\n";
    json << "    {\"name\":" << Quote(p.input.name) << ",\"shape\":" << Quote(p.input.shape)
         << ",\"anchor\":" << (p.input.anchor ? "true":"false") << ",\"energy_MeV\":" << p.input.energy_MeV
         << ",\"energy_error_MeV\":" << p.input.energy_error_MeV << ",\"yield\":" << p.yield
         << ",\"yield_error\":" << p.yield_error << ",\"centroid_mm\":" << p.centroid_mm
         << ",\"centroid_error_mm\":" << p.centroid_error_mm << ",\"FWHM_mm\":" << p.fwhm_mm
         << ",\"FWHM_error_mm\":" << p.fwhm_error_mm << ",\"known_minus_calibrated_keV\":" << p.residual_keV << '}';
    csv << p.input.name << ',' << p.input.shape << ',' << p.input.anchor << ',' << p.input.energy_MeV << ','
        << p.input.energy_error_MeV << ',' << p.yield << ',' << p.yield_error << ',' << p.centroid_mm << ','
        << p.centroid_error_mm << ',' << p.fwhm_mm << ',' << p.fwhm_error_mm << ',' << p.residual_keV << '\n';
  }
  json << "\n  ],\n  \"warnings\":[";
  for(size_t k=0;k<r.warnings.size();++k) { if(k) json<<','; json<<Quote(r.warnings[k]); }
  json << "]\n}\n";
  card << "# E(MeV) = intercept + slope*x(mm) + quadratic*x(mm)^2\n"
       << "CALIB_SLOPE: " << r.a[1] << "\nCALIB_INTERCEPT: " << r.a[0]
       << "\nCALIB_QUADRATIC: " << r.a[2] << '\n';
  json.flush(); csv.flush(); card.flush(); input.flush();
  Require(bool(json)&&bool(csv)&&bool(card)&&bool(input),"failed writing output reports");

  RawModel model(s,h.GetBinWidth(1));
  TF1 f("raw_model_report",model,s.xmin_mm,s.xmax_mm,ParameterCount(s));
  f.SetParameters(r.spectrum_parameters.data());
  TGraph total, background, pull;
  const int n=std::max(2000,h.GetNbinsX()*2);
  std::vector<TGraph> components(s.peaks.size());
  TGraph cross;
  for(int k=0;k<n;++k) {
    const double x=s.xmin_mm+(s.xmax_mm-s.xmin_mm)*k/(n-1), w=h.GetBinWidth(1);
    total.SetPoint(k,x,f.Eval(x));
    const double bg=w*Background(x,r.spectrum_parameters.data(),s);
    background.SetPoint(k,x,bg);
    double incoherent=bg;
    for(size_t j=0;j<s.peaks.size();++j) {
      const double y=w*r.peaks[j].yield*Shapes().at(s.peaks[j].shape)(x,r.peaks[j].centroid_mm,r.peaks[j].fwhm_mm);
      components[j].SetPoint(k,x,y); incoherent+=y;
    }
    cross.SetPoint(k,x,f.Eval(x)-incoherent);
  }
  int ip=0;
  for(int b=1;b<=h.GetNbinsX();++b) {
    const double mu=f.Integral(h.GetXaxis()->GetBinLowEdge(b),h.GetXaxis()->GetBinUpEdge(b))/h.GetBinWidth(b);
    const double sigma=s.poisson ? std::sqrt(mu) : h.GetBinError(b);
    if(sigma>0) pull.SetPoint(ip++,h.GetBinCenter(b),(h.GetBinContent(b)-mu)/sigma);
  }
  TFile out((stem+".root").c_str(),"RECREATE"); Require(!out.IsZombie(),"cannot create ROOT report");
  h.Write("raw_histogram"); total.Write("total_fit"); background.Write("background"); pull.Write("pulls");
  for(size_t j=0;j<components.size();++j) components[j].Write(("component_"+s.peaks[j].name).c_str());
  if(s.interference_first>=0) cross.Write("interference_cross_term");
  r.coefficient_covariance.Write("coefficient_covariance"); r.spectrum_covariance.Write("spectrum_covariance");
  TCanvas raw("raw_calibration_canvas","Raw position fit",1100,850); raw.Divide(1,2);
  raw.cd(1); h.SetStats(false); h.SetTitle("");
  h.GetXaxis()->SetTitle("Focal-plane position (mm)"); h.GetYaxis()->SetTitle("Counts / bin"); h.Draw("hist");
  total.SetLineColor(2); total.SetLineWidth(2); total.Draw("L SAME");
  background.SetLineStyle(2); background.Draw("L SAME");
  TLegend legend(.12,.64,.38,.88); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.AddEntry(&total,"Total fit","l"); legend.AddEntry(&background,"Background","l");
  for(size_t j=0;j<components.size();++j) {
    components[j].SetLineColor(3+j); components[j].Draw("L SAME");
    legend.AddEntry(&components[j],s.peaks[j].name.c_str(),"l");
  }
  legend.Draw();
  raw.cd(2); pull.SetTitle(";Focal-plane position (mm);(Data-model) / uncertainty"); pull.SetMarkerStyle(7); pull.Draw("AP");
  TLine zero(s.xmin_mm,0,s.xmax_mm,0); zero.SetLineStyle(2); zero.Draw();
  raw.Write(); raw.SaveAs((stem+"_spectrum.pdf").c_str());

  TGraphErrors anchors, residuals;
  int ia=0;
  for(const auto& p:r.peaks) if(p.input.anchor) {
    anchors.SetPoint(ia,p.centroid_mm,p.input.energy_MeV);
    anchors.SetPointError(ia,p.centroid_error_mm,p.input.energy_error_MeV);
    residuals.SetPoint(ia,p.centroid_mm,p.residual_keV);
    // These are input errors, NOT uncertainties of fitted residuals.
    const double d=r.a[1]+2*r.a[2]*p.centroid_mm;
    residuals.SetPointError(ia,0,1000*std::hypot(p.input.energy_error_MeV,d*p.centroid_error_mm));
    ++ia;
  }
  TF1 polynomial("energy_calibration","[0]+[1]*x+[2]*x*x",s.xmin_mm,s.xmax_mm);
  polynomial.SetParameters(r.a.data());
  TCanvas cal("energy_calibration_canvas","Energy calibration",1100,850); cal.Divide(1,2);
  cal.cd(1); anchors.SetTitle(";Focal-plane position (mm);Known excitation energy (MeV)");
  anchors.SetMarkerStyle(20); anchors.Draw("AP"); polynomial.Draw("SAME");
  anchors.GetXaxis()->SetLimits(s.xmin_mm,s.xmax_mm);
  cal.cd(2); residuals.SetTitle("Input-error bars;Focal-plane position (mm);Known - calibrated (keV)");
  residuals.SetMarkerStyle(20); residuals.Draw("AP");
  residuals.GetXaxis()->SetLimits(s.xmin_mm,s.xmax_mm); zero.Draw();
  anchors.Write("calibration_anchors"); residuals.Write("calibration_residuals"); polynomial.Write();
  cal.Write(); cal.SaveAs((stem+"_calibration.pdf").c_str());
  out.Flush();
  Require(!out.TestBit(TFile::kWriteError),"failed writing ROOT report");
  out.Close();
}
CalibrationResult RunCalibration(const CalibrationSpec& s) {
  auto h=LoadRawSpectrum(s); auto r=FitRawCalibrationSpectrum(*h,s);
  WriteCalibrationReports(*h,s,r);
  std::cout << std::setprecision(12) << "\nRaw spectrum fit in mm:\n";
  for(const auto& p:r.peaks)
    std::cout << "  " << p.input.name << " (" << p.input.shape << (p.input.anchor ? ", anchor" : ", model only")
              << "): x=" << p.centroid_mm << " +/- " << p.centroid_error_mm
              << " mm; FWHM=" << p.fwhm_mm << " +/- " << p.fwhm_error_mm << " mm\n";
  std::cout << "Calibration E(x)=a0+a1*x+a2*x^2:\n";
  for(int k=0;k<3;++k) std::cout << "  a" << k << " = " << r.a[k] << " +/- "
    << std::sqrt(std::max(0.0,r.coefficient_covariance(k,k))) << '\n';
  std::cout << "Calibration chi2/NDF = " << r.calibration_chi2 << '/' << r.calibration_ndf
            << "\nReports: " << s.output_prefix << "*\n";
  for(const auto& w:r.warnings) std::cerr << "Calibration warning: " << w << '\n';
  return r;
}
} // namespace bwfit

// Compile this file with -DBWFIT_CALIBRATION_MAIN to get a separate command.
// Without the macro it has no main(), so it can be linked into bwsingf.
#ifdef BWFIT_CALIBRATION_MAIN
int main(int argc,char** argv) {
  if(argc!=2) { std::cerr << "Usage: bwcalibrate config/calibration.txt\n"; return 1; }
  gROOT->SetBatch(true);
  try { bwfit::RunCalibration(bwfit::ReadCalibrationSpec(argv[1])); return 0; }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
#endif
