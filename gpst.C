/**
 * This is the Gamma-ray Polarimetry Simulation Toolkit GPST
 * =========================================================
 *
 * Based on the X-Calibur Simulation Toolkit (XST)
 * by Henric Krawczynski and Fabian Kislat
 *
 * A detailed manual can be found at 
 *    https://sites.physics.wustl.edu/xcalibur/wiki/index.php/X-Calibur_observation_simulation_code
 *
 * Change Log
 * ----------
 * 06/05/17 - Version 14  - Initial adaptation from XST version 13 - FK
 */


#define GPST_VERSION 14
#define GPST_VERSION_STRING "14.0"
#define GPST_RELEASE_DATE "06/05/2017"


#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TInterpreter.h>

#include <iostream>
#include <string>


struct Result {
  Result();
  Result(const std::string &inputfile, TCanvas *canvas_,
         TGraphAsymmErrors *flux_, TGraphAsymmErrors *fraction_,
         TGraphAsymmErrors *angle_);

  std::string file_prefix;
  TCanvas *canvas;
  TGraphAsymmErrors *flux;
  TGraphAsymmErrors *fraction;
  TGraphAsymmErrors *angle;

  void info();
  void save();
  void saveAs(const std::string &filename);
};


void gpst(const char *name, bool plot_clean, bool large, bool debug);
void gpst_version();
void save();
void saveAs(const std::string &filename);

Result result;

// Let's not bother CINT with the details of this code
#ifndef __CINT__

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>


using namespace std;

/* Mathematical Constants */
#define LOG10   2.302585093
#define PI      M_PI
#define MaxBins 10000

static int random_seed = 0;

enum FluxUnit {
  CGS,
  CGSe9
};

enum FractionUnit {
  Frac,
  Percent
};

enum Panel {
  Flux,
  Fraction,
  Angle
};


double Crab(double energy);
double getChi(double q, double u);
void getValuesAndErrors(double Nsrc, double Nbg, double mu, double aT,
                        double phiT, double &fluxM, double &aM, double &phiM,
                        double &fluxE, double &aELow, double &aEHigh,
                        double &phiE);
double prob(double ratio,double aT,double mu,double aM,double deltaPhi);
void loadSensitivity(double *ns, double *nbg, double *mu, int numbins,
                     double modulation_factor=0.52, double start = 2.,
                     double deltaE = 0.2, string sensitivityFile = "",
                     string backgroundFile = "", string modulationFile = "",
                     double scale = 1);
void interpolate_stokes(double *xs, double *Is, double *Qs, double *Us,
                        int numbins, double x, double &I, double &Q, double &U);
double get_table(double *t_x,double *t_v,int numbins,double x);
void linearIntegralQUI(double xMin, double xMax, int nPoints, double *x,
                       double *flux, double *pi, double *chi,
                       double &qsum, double &usum, double &isum);
void style();
bool read_version_info(const std::string &line, int &version, int linecount);
bool interpretBoolean(const std::string &value);
Style_t interpretStyle(const std::string &value);
Color_t interpretColor(const std::string &value);
FluxUnit interpretFluxUnit(const std::string &value);
FractionUnit interpretFractionUnit(const std::string &value);
Panel interpretPanel(const std::string &value);
string make_filename(const std::string &base, std::string input);


static std::set<std::string> SPECIAL_VALUES;

#define GPST_ADD_SPECIAL_VALUE(map, key, value) \
  map[key] = value;                            \
  SPECIAL_VALUES.insert(key)

inline bool is_special_string(const std::string &value) {
  return SPECIAL_VALUES.find(value) != SPECIAL_VALUES.end();
}

std::string prepare_special_value(const std::string &value);

bool is_special_value(const std::string &value) {
  if (is_special_string(value)) return true;

  try {
    interpretStyle(value);
    return true;
  } catch (...) {}

  try {
    interpretColor(value);
    return true;
  } catch (...) {}

  try {
    interpretBoolean(value);
    return true;
  } catch (...) {}

  try {
    interpretFractionUnit(value);
    return true;
  } catch (...) {}

  try {
    interpretFluxUnit(value);
    return true;
  } catch (...) {}

  try {
    interpretPanel(value);
    return true;
  } catch (...) {}

  return false;
}


static bool do_debug = false;


void __attribute__((__format__(printf, 1, 3)))
log_debug_impl(const char *format, int line, ...)
{
  if (!do_debug) return;

  printf("GPST (%d) DEBUG: ", line);

  string format_str = string(format) + "\n";

  va_list va;
  va_start(va, line);
  vprintf(format_str.c_str(), va);
  va_end(va);
}

#define log_debug(format, ...)                          \
  log_debug_impl(format, __LINE__, ##__VA_ARGS__)

template <typename T>
void debug_print_vector(const std::string &label, const vector<T> &v) {
  if (!do_debug) return;

  stringstream ss;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0) ss << ", ";
    ss << v[i];
  }
  log_debug("%s: [ %s ]", label.c_str(), ss.str().c_str());
}

void __attribute__((__format__(printf, 1, 3)))
log_info_impl(const char *format, int line, ...)
{
  printf("\x1b[1mGPST (%d) INFO: \x1b[0m", line);
  
  string format_str = string(format) + "\n";

  va_list va;
  va_start(va, line);
  vprintf(format_str.c_str(), va);
  va_end(va);
}

#define log_info(format, ...)                           \
  log_info_impl(format, __LINE__, ##__VA_ARGS__)

void __attribute__((__format__(printf, 1, 3)))
log_warn_impl(const char *format, int line, ...)
{
  printf("\x1b[33mGPST (%d) WARNING: \x1b[0m", line);
  
  string format_str = string(format) + "\n";

  va_list va;
  va_start(va, line);
  vprintf(format_str.c_str(), va);
  va_end(va); 
}

#define log_warn(format, ...)                           \
  log_warn_impl(format, __LINE__, ##__VA_ARGS__)

void __attribute__((__format__(printf, 1, 3)))
log_error_impl(const char *format, int line, ...)
{
  fprintf(stderr, "\x1b[1m\x1b[31mGPST (%d) ERROR: \x1b[0m", line);
  
  string format_str = string(format) + "\n";

  va_list va;
  va_start(va, line);
  vfprintf(stderr, format_str.c_str(), va);
  va_end(va);
}

#define log_error(format, ...)                          \
  log_error_impl(format, __LINE__, ##__VA_ARGS__)

void __attribute__((__format__(printf, 1, 3))) __attribute__((__noreturn__))
log_fatal_impl(const char *format, int line, ...)
{
  fprintf(stderr, "\x1b[1m\x1b[31mGPST (%d), FATAL: ", line);
  
  string format_str = string(format) + "\x1b[0m\n";

  va_list va;
  va_start(va, line);
  vfprintf(stderr, format_str.c_str(), va);
  va_end(va);

  abort();
}

#define log_fatal(format, ...)                          \
  log_fatal_impl(format, __LINE__, ##__VA_ARGS__)


inline TRandom3& getRNG() {
  static TRandom3 rng(random_seed);
  return rng;
}


enum ParserState {
  ExpectVersionInfo,
  ExpectAssignment,
  ExpectArray
};

enum ParserFlags {
  NoParserFlags = 0,
  ExpectComponentEnergies = 1,
  ExpectComponentFlux = 1<<1,
  ExpectComponentPolarization = 1<<2,
  ExpectComponentAngle = 1<<3,
  HaveComponentPoints = 1<<4,
  PreparePhase = 1<<5,
  PrepareExtraModel = 1<<6,
  ExpectModelEnergiesOrPhases = 1<<7
};


namespace strings {

  void strip(string &s);
  void strip(string &s, char c);
  void split(const string &str, char delim, vector<string> &v);
  bool split_assignment(const string &str, string &key, string &value);

  template <typename T>
  T to(const std::string &from)
  {
    stringstream ss(from);
    ss.exceptions(ios_base::badbit | ios_base::failbit);
    T v;
    ss >> v;
    return v;
  }

}


namespace math {
  
  template <typename T>
  inline T sqr(const T &x) {
    return x*x;
  }

  inline bool isnan(double x) {
    return x != x;
  }
  
}


struct Dimensions {
  int frameHeight;
  int frameWidth;
  int leftMargin;
  int rightMargin;
  int bottomMargin;
  int topMargin;
  int textHeight;
  int width;
  int height;
  int yoffset;
};


struct PanelConfig {
  int height;
  double fontsize;
  double bottom;
  double top;
  double bottomMargin;
  double topMargin;
  double yoffset;
};


struct ComponentInfo {
  ComponentInfo()
    : line(-1), name(""), flags(0), color(-1), style(-1), width(2)
  {}

  int line;
  string name;
  int flags;
  int color;
  int style;
  Width_t width;
  int legend;
};


struct ExtraModel {
  ExtraModel() : line(-1), name(""), color(-1), style(-1), width(2) {}
  
  int line;
  string name;
  int color;
  int style;
  Width_t width;
  vector<double> energy;
  vector<double> phase;
  vector<double> flux;
  vector<double> fraction;
  vector<double> angle;
  int legend;

  void clear() {
    line = -1;
    name = "";
    color = -1;
    style = -1;
    width = 2;
    energy.clear();
    phase.clear();
    flux.clear();
    fraction.clear();
    angle.clear();
    legend = -1;
  }
};


struct LegendConfig {
  double x1;
  double y1;
  double x2;
  double y2;
  Panel panel;
  unsigned panel_index;
  double legendfont;
  int columns;
  TLegend *legend;
};


PanelConfig makePanelConfig(unsigned index, unsigned nPads, double &top,
                            const Dimensions &d, bool large);

double calculateTitleOffset(int pixels, double xmin, double xmax,
                            double fontsize, bool logx, const Dimensions &d);


Result::Result()
  : file_prefix(""), canvas(NULL), flux(NULL), fraction(NULL), angle(NULL)
{}


Result::Result(const std::string &inputfile, TCanvas *canvas_,
               TGraphAsymmErrors *flux_, TGraphAsymmErrors *fraction_,
               TGraphAsymmErrors *angle_)
  : file_prefix(""), canvas(canvas_), flux(flux_), fraction(fraction_),
    angle(angle_)
{
  size_t slash = inputfile.rfind('/');
  string filename = slash != string::npos ? inputfile.substr(slash+1) : inputfile;
  string dirname = slash != string::npos ? inputfile.substr(0, slash+1) : "";

  size_t dot = filename.rfind('.');
  if (dot != string::npos && dot > 0)
    filename = filename.substr(0, dot);

  if (filename.length() > 4 && filename.substr(0, 4) == "data")
    filename.replace(0, 4, "canvas");
  else
    filename = "canvas" + filename;

  file_prefix = dirname + filename;
}


void Result::info()
{
  if (!canvas) {
    log_error("No canvas. Result invalid.");
  } else {
    log_info("Calling save() will save the canvas as");
    log_info("    %s.root", file_prefix.c_str());
    log_info("and");
    log_info("    %s.eps", file_prefix.c_str());
    log_info("You can also use saveAs(filename) to save under a custom file");
    log_info("name and format. Supported formats are:");
    log_info("    root, C, png, eps, pdf");
    log_info("All other file extensions will store the simulated data in a");
    log_info("text file.");
  }
}


void Result::save()
{
  if (canvas) {
    canvas->SaveAs((file_prefix + ".root").c_str());
    canvas->Print((file_prefix + ".eps").c_str());
  } else {
    log_error("No canvas. Result cannot be saved.");
  }
}


void Result::saveAs(const std::string &filename)
{
  if (!canvas) {
    log_error("No canvas. Result cannot be saved.");
    return;
  }

  std::string fileextension = "";
  size_t dot = filename.rfind('.');
  if (dot != string::npos && dot > 0) {
    fileextension = filename.substr(dot+1);
  }

  if (fileextension == "root" || fileextension == "C") {
    canvas->SaveAs(filename.c_str());
    return;
  } else if (fileextension == "png" ||
             fileextension == "eps" ||
             fileextension == "pdf" ||
             fileextension == "jpg")
  {
    canvas->Print(filename.c_str());
    return;
  } else if (fileextension != ".txt") {
    log_info("Unknown file extension (%s), saving as text file.",
             fileextension.c_str());
  }

  if (!flux || !fraction || !angle) {
    log_error("Missing graphs, cannot save text file.");
    return;
  }

  if (flux->GetN() != fraction->GetN() || flux->GetN() != angle->GetN()) {
    log_error("The number of points in the graphs don't match. "
              "Cannot save as text.");
    return;
  }

  std::ofstream outfile(filename.c_str());
  if (!outfile.is_open()) {
    log_error("Failed to open output file %s", filename.c_str());
    return;
  }

  outfile << "# <E low/keV> <E high/keV> <flux> <flux error> "
          << "<fraction> <sigma fraction low> <sigma fraction high> "
          << "<angle> <sigma angle>\n";

  int npoints = flux->GetN();
  for (int i = 0; i < npoints; ++i) {
    double elow = flux->GetX()[i] - flux->GetEXlow()[i];
    double ehigh = flux->GetX()[i] + flux->GetEXhigh()[i];
    
    // note: for flux and angle EYlow and EYhigh are the same
    outfile << elow << " " << ehigh << " "
            << flux->GetY()[i] << " " << flux->GetEYlow()[i] << " "
            << fraction->GetY()[i] << " " << fraction->GetEYlow()[i] << " "
            << fraction->GetEYhigh()[i] << " "
            << angle->GetY()[i] << " " << angle->GetEYlow()[i] << '\n';
  }

  outfile << std::flush;
}


void gpst(const char *name, bool plot_clean, bool large, bool debug)
{
  // store in global variable
  do_debug = debug;

  std::cout << "\x1b[1m\n==========\n";
  gpst_version();
  std::cout << "==========\n\x1b[0m" << std::endl;

  /* Sensitivity Information */

  int NumBins = 340;
  double eSensitivityMin = 2.;
  double deltaESensitivity = 0.2;
  string sensitivityFile = "";
  double sensitivityScale = 1.;
  string backgroundFile = "";
  double modulation_factor = 0.52;
  bool have_modulationfactor_keyword = false;
  string modulationFile = "";

#define MaxB 100
  double bEnergy[MaxB];
  double bPhase[MaxB];
    
#define MaxC  20
#define MaxP  500
  ComponentInfo component[MaxC];
  int    nP[MaxC];
  double dE[MaxC][MaxP];
  double dI[MaxC][MaxP];
  double dPi[MaxC][MaxP];
  double dChi[MaxC][MaxP];

  double dQ[MaxC][MaxP];
  double dU[MaxC][MaxP];

  bool have_phase = false;
  int nPPhi = 0;
  double dPhi[MaxP];
  double dPhiI[MaxP];
  double dPhiPi[MaxP];
  double dPhiChi[MaxP];

  log_info("Opening input file: '%s'.", name);

  stack<ifstream*> filestack;
  stack<string> filenamestack;
  stack<int> linecountstack;
  ifstream *fp_in = new ifstream(name);
  string lastfilename = make_filename("", name);
    
  if (!fp_in->good()) {
    delete fp_in;
    log_fatal("Failed to open input file '%s'.", name);
  }

  // flux scaling defaults to the RXTE ASM 2-12keV fluxes
  double scaleEmin = 2.;
  double scaleEmax = 12.;
  double aFlux = 0.;
  bool have_reference_range = false;
  double time = -1;
  int nBins = 0;
  int nPhaseBins = 0;
  int nComp = 0;
  double Hymin = 0;
  double Hymax = 0;
  double Hpmin = 0;
  double Hpmax = 0;
  double Hchimin = 0;
  double Hchimax = 0;
  double deltachi = 0;
  double userXmin = NAN;
  double userXmax = NAN;
  int fluxNdiv = 510;
  int piNdiv = 510;
  int chiNdiv = 510;
  int xNdiv = 510;

  vector<Panel> panels;
  panels.push_back(Flux);
  panels.push_back(Fraction);
  panels.push_back(Angle);

  vector<LegendConfig> legends;
  int model_legend = -1;
  int data_legend = -1;
  int mdp_legend = -1;

  vector<ExtraModel> extraModels;
  ExtraModel currentExtraModel;

  bool showmdp = false;
  bool customizeMdp = false;
  Color_t mdpColor = kBlack;
  Style_t mdpStyle = kSolid;
  Width_t mdpWidth = 2;
  string mdpLabel = "MDP";

  string userFluxLabel = "";
  string fluxLabel = "E f_{E} [erg cm^{-2} s^{-1}]";
  bool fluxMoreLabels = false;
  string userFracLabel = "";
  string polFracLabel = "Pol. Fraction";
  string userDirLabel = "";
  string polDirLabel = "Pol. Direction";
  string userXLabel = "";

  FluxUnit fluxunit = CGS;
  FractionUnit fractionunit = Frac;

  enum ModelOrSimulation {
    CustomizeModel,
    CustomizeSimulation
  };

  ModelOrSimulation modelOrSimulation = CustomizeModel;
  
  string datalabel = "Measurement";
  Color_t datacolor = kBlack;
  Width_t datawidth = 1;
  
  string simulation_name;
  Color_t linecolor = kBlack;
  Style_t linestyle = kSolid;
  Width_t linewidth = 3;
  
  ParserState pState = ExpectVersionInfo;
  int pFlags = NoParserFlags;
  int phaseFlags = 0;
  int linecount = 0;
  int read_version = -1;
  string array_key = "";
  vector<string> array;
  int currentComponentLine = 0;
  string currentComponentName = "";
  int currentComponentColor = -1;
  int currentComponentStyle = -1;
  Width_t currentComponentWidth = 2;
  int currentComponentLegend = -1;
  double ontime = 1;

  while (fp_in) {
    while (!fp_in->eof()) {
      string line;
      getline(*fp_in, line);
      ++linecount;
      
      strings::strip(line);
      if (line.empty()) continue;
      
      // comments start with a ;
      if (line[0] == ';') continue;
      
      string key = "", value = "";
      bool is_array = false;
      bool assignment_complete = false;
      
      if (pState == ExpectVersionInfo) {
        if (!read_version_info(line, read_version, linecount)) abort();
        pState = ExpectAssignment;
        if (read_version < 12) {
          log_warn("The way the 'legend' keyword is handled changed in Version "
                   "12.0 of XST. Results may not look as expected.");
          log_warn("The simplest fix is to move the 'legend' keyword to the top "
                   "of your input file (if you have a legend).");
        }
      } else if (pState == ExpectAssignment) {
        if (line[0] == '@') {
          static const string INCLUDE = "@include";
          if (line.substr(0, INCLUDE.length()) == INCLUDE) {
            line = line.substr(INCLUDE.length());
            line = line.substr(0, line.find(';'));
            strings::strip(line);

            string fn = make_filename(lastfilename, line);
            
            filestack.push(fp_in);
            filenamestack.push(lastfilename);
            lastfilename = fn;
            log_info("Opening input file: '%s'.", fn.c_str());
            fp_in = new ifstream(fn.c_str());
            if (!fp_in->good()) {
              delete fp_in;
              while (!filestack.empty()) {
                delete filestack.top();
                filestack.pop();
              }
              log_fatal("On line %d of input file: Failed to open file '%s'.",
                        linecount, fn.c_str());
            }
            linecountstack.push(linecount);
            linecount = 0;
          } else {
            log_fatal("On line %d of input file: Found invalid @-statement: %s",
                      linecount, line.c_str());
          }
        } else {
          if (!strings::split_assignment(line, key, value)) {
            log_fatal("On line %d of input file: Expecting assignment.",
                      linecount);
          }
          strings::strip(key);
          strings::strip(value);
          
          if (key.empty() || value.empty()) {
            log_fatal("On line %d of input file: Expecting assignment.",
                      linecount);
          }
          
          if (value[0] == '[') {
            size_t array_end = value.find(']');
            array.clear();
            strings::split(value.substr(1, array_end == string::npos ?
                                        string::npos : array_end-1),
                           ' ', array);
            array_key = key;
            is_array = true;
            if (array_end != string::npos) {
              string after_array = value.substr(array_end + 1);
              strings::strip(after_array);
              if (!after_array.empty() && after_array[0] != ';')
                log_warn("On line %d of input file: Invalid text after array.",
                         linecount);
              assignment_complete = true;
            } else {
              pState = ExpectArray;
            }
          } else {
            assignment_complete = true;
          }
        }
      } else if (pState == ExpectArray) {
        size_t array_end = line.find(']');
        if (array_end != 0) {
          strings::split(line.substr(0, array_end == string::npos ?
                                     string::npos : array_end-1), ' ', array);
        }
        if (array_end != string::npos) {
          string after_array = line.substr(array_end + 1);
          strings::strip(after_array);
          if (!after_array.empty() && after_array[0] != ';')
            log_warn("On line %d of input file: Invalid text after array.",
                     linecount);
          assignment_complete = true;
          pState = ExpectAssignment;
        }
        key = array_key;
        is_array = true;
      }
        
      if (assignment_complete) {
        double double_value = 0.;
        bool have_double = false;
          
        if (!is_array) {
          if (read_version >= 9 && value[0] == '"') {
            string parsed_value = "";
            bool have_mask = false;
            size_t value_iter;
            for (value_iter = 1; value_iter < value.size(); ++value_iter) {
              if (value[value_iter] == '"') {
                if (have_mask) {
                  parsed_value += value[value_iter];
                  have_mask = false;
                } else {
                  ++value_iter;
                  break;
                }
              } else if (value[value_iter] == '\\') {
                if (have_mask) {
                  parsed_value += value[value_iter];
                  have_mask = false;
                } else {
                  have_mask = true;
                  continue;
                }
              } else {
                if (have_mask) {
                  log_warn("On line %d of input file: found \\ followed by %c. "
                           "Please use \\\\ if you want to have a single back "
                           "slash.", linecount, int(value[value_iter]));
                  parsed_value += '\\';
                  have_mask = false;
                }
                parsed_value += value[value_iter];
              }
            }
            value = value.substr(value_iter);
            strings::strip(value);
            if (!value.empty() && value[0] != ';') {
              log_warn("On line %d of input file: Invalid text after end of "
                       "string parameter.", linecount);
            }
            value = parsed_value;
          } else {
            try {
              double_value = strings::to<double>(value);
              have_double = true;
            } catch (...) {
              if (read_version >= 9 && !is_special_value(value))
                log_warn("On line %d of input file: Found unquoted string (%s). "
                         "Since version 9 of GPST string parameters must be "
                         "enclosed by double quotes(\"...\").",
                         linecount, value.c_str());
            }
          }
          log_debug("Key '%s' ==> Value '%s'", key.c_str(), value.c_str());
        } else if (is_array && do_debug) {
          stringstream str;
          str << "[ ";
          for (vector<string>::const_iterator it = array.begin();
               it != array.end(); ++it)
          {
            str << *it << " ";
          }
          str << "]";
          log_debug("Key '%s' ==> Value %s", key.c_str(), str.str().c_str());
        }
          
#define GPST_ASSERT_NOT_ARRAY(variable)                                  \
        do {                                                            \
          if (is_array) {                                               \
            log_fatal("On line %d of input file: %s cannot be an array.", \
                      linecount, #variable);                            \
          }                                                             \
        } while(false)
          
#define GPST_ASSERT_ARRAY(variable)                                      \
        do {                                                            \
          if (!is_array) {                                              \
            log_fatal("On line %d of input file: %s must be an array.", \
                      linecount, #variable);                            \
          }                                                             \
        } while(false)
          
#define GPST_ASSERT_PFLAG(flag, key)                                     \
        do {                                                            \
          if (!(pFlags & flag)) {                                       \
            log_fatal("On line %d of input file: Unexpected key '%s'.", \
                      linecount, #key);                                 \
          }                                                             \
        } while(false)
          
#define GPST_DISABLE_PFLAG(flag) pFlags &= ~flag
          
#define GPST_USE_DOUBLE(key, to_variable)                                \
        do {                                                            \
          if (have_double) {                                            \
            to_variable = double_value;                                 \
          } else {                                                      \
            log_fatal("On line %d of input file: '%s' is supposed to be " \
                      "of type double but conversion failed.",          \
                      linecount, #key);                                 \
          }                                                             \
        } while (false)
          
#define GPST_CONVERT(key, to_type, to_variable, from)                    \
        do {                                                            \
          try {                                                         \
            to_variable = strings::to<to_type>(from);                   \
          } catch (...) {                                               \
            log_fatal("On line %d of input file: '%s' is supposed to be " \
                      "of type %s but conversion failed.",              \
                      linecount, #key, #to_type);                       \
          }                                                             \
        } while (false)
        
#define GPST_CHECK_LENGTH(key, length)                                   \
        do {                                                            \
          if (pFlags & HaveComponentPoints) {                           \
            if (array.size() != (size_t)length) {                       \
              log_fatal("On line %d of input file: Number of components " \
                        "in key '%s' does not match (expected %d, got " \
                        "%zu).", linecount, #key, length, array.size()); \
            }                                                           \
          } else {                                                      \
            length = array.size();                                      \
            pFlags |= HaveComponentPoints;                              \
          }                                                             \
        } while(false)
          
#define GPST_CONVERT_ARRAY(key, to_type, to_variable, from)              \
        do {                                                            \
          for (size_t i = 0; i < from.size(); ++i) {                    \
            try {                                                       \
              to_variable[i] = strings::to<to_type>(from[i]);           \
            } catch (...) {                                             \
              log_fatal("On line %d of input file: %s is supposed to be " \
                        "an array of type %s but conversion failed.",   \
                        linecount, #key, #to_type);                     \
            }                                                           \
          }                                                             \
        } while (false)
          
#define GPST_REQUIRE_MIN_VERSION(keyword, vmin)                          \
        do {                                                            \
          if (read_version < vmin) {                                    \
            log_fatal("On line %d of input file: Keyword '%s' only "    \
                      "supported in file format version %d or later.",  \
                      linecount, #keyword, vmin);                       \
          }                                                             \
        } while(false)
          
        if (key == "panels") {
          GPST_REQUIRE_MIN_VERSION(panels, 13);
          if (is_array) {
            panels.resize(array.size());
            for (size_t i = 0; i < array.size(); ++i) {
              try {
                panels[i] = interpretPanel(array[i]);
              } catch (...) {
                log_fatal("On line %d of input file: Expected panel name "
                          "('flux', 'frac', 'angle') but got '%s'",
                          linecount, array[i].c_str());
              }
              for (size_t j = 0; j < i; ++j) {
                if (panels[i] == panels[j]) {
                  log_fatal("On line %d of input file: Panel %s occured more "
                            "than once", linecount, array[i].c_str());
                }
              }
            }
          } else {
            panels.resize(1);
            try {
              panels[0] = interpretPanel(value);
            } catch (...) {
              log_fatal("On line %d of input file: Expected panel name ('flux', "
                        "'frac', 'angle') but got '%s'",
                        linecount, value.c_str());
            }
          }
        } else if (key == "name") {
          GPST_ASSERT_NOT_ARRAY(name);
          simulation_name = value;
          model_legend = legends.size()-1;
          modelOrSimulation = CustomizeModel;
        } else if (key == "renormalize") {
          if (is_array) {
            if (read_version < 14) {
              log_fatal("On line %d of input file: Keyword 'renormalize' as "
                        "array requires version 14 or greater of GPST.",
                        linecount);
            }
            if (array.size() != 3) {
              log_fatal("On line %d of input file: Number of components in key "
                        "'renormalize' has to be 3.", linecount);
            }
            double renormalize_data[3];
            GPST_CONVERT_ARRAY(renormalize, double, renormalize_data, array);
            aFlux = renormalize_data[0];
            scaleEmin = renormalize_data[1];
            scaleEmax = renormalize_data[2];
            have_reference_range = true;
          } else {
            GPST_USE_DOUBLE(renormalize, aFlux);
          }
        } else if (key == "referencerange") {
          GPST_REQUIRE_MIN_VERSION(referencerange, 14);
          if (have_reference_range) {
            log_error("Keyword 'referencerange' already specified, or energy "
                      "range given in 'renormalize' keyword. Ignoring second "
                      "occurrence.");
          } else {
            GPST_ASSERT_ARRAY(referencerange);
            if (array.size() != 2) {
              log_fatal("On line %d of input file: Number of components in key "
                        "'referencerange' has to be 2.", linecount);
            }
            double referencerange_data[4];
            GPST_CONVERT_ARRAY(referencerange, double, referencerange_data, array);
            scaleEmin = referencerange_data[0];
            scaleEmax = referencerange_data[1];
          }
        } else if (key == "fluxunit") {
          GPST_REQUIRE_MIN_VERSION(fluxunit, 11);
          GPST_ASSERT_NOT_ARRAY(fluxunit);
          try {
            fluxunit = interpretFluxUnit(value);
          } catch (std::runtime_error &e) {
            log_fatal("On line %d of input file: %s", linecount, e.what());
          }
        } else if (key == "fracunit") {
          GPST_REQUIRE_MIN_VERSION(fractionunit, 11);
          GPST_ASSERT_NOT_ARRAY(fractionunit);
          try {
            fractionunit = interpretFractionUnit(value);
          } catch (std::runtime_error &e) {
            log_fatal("On line %d of input file: %s", linecount, e.what());
          }
        } else if (key == "time") {
          GPST_ASSERT_NOT_ARRAY(time);
          GPST_USE_DOUBLE(time, time);
          time*=86400.;
        } else if (key == "ontime") {
          GPST_REQUIRE_MIN_VERSION(ontime, 14);
          GPST_ASSERT_NOT_ARRAY(ontime);
          GPST_USE_DOUBLE(ontime, ontime);
        } else if (key == "fluxmin") {
          GPST_ASSERT_NOT_ARRAY(fluxmin);
          GPST_USE_DOUBLE(fluxmin, Hymin);
        } else if (key == "fluxmax") {
          GPST_ASSERT_NOT_ARRAY(fluxmax);
          GPST_USE_DOUBLE(fluxmax, Hymax);
        } else if (key == "fracmin") {
          GPST_ASSERT_NOT_ARRAY(fracmin);
          GPST_USE_DOUBLE(fracmin, Hpmin);
        } else if (key == "fracmax") {
          GPST_ASSERT_NOT_ARRAY(fracmax);
          GPST_USE_DOUBLE(fracmax, Hpmax);
        } else if (key == "chimin") {
          GPST_ASSERT_NOT_ARRAY(chimin);
          GPST_USE_DOUBLE(chimin, Hchimin);
        } else if (key == "chimax") {
          GPST_ASSERT_NOT_ARRAY(chimax);
          GPST_USE_DOUBLE(chimax, Hchimax);
        } else if (key == "xmin") {
          GPST_REQUIRE_MIN_VERSION(xmin, 11);
          GPST_ASSERT_NOT_ARRAY(xmin);
          GPST_USE_DOUBLE(xmin, userXmin);
        } else if (key == "xmax") {
          GPST_REQUIRE_MIN_VERSION(xmax, 11);
          GPST_ASSERT_NOT_ARRAY(xmax);
          GPST_USE_DOUBLE(xmax, userXmax);
        } else if (key == "fluxndiv") {
          GPST_REQUIRE_MIN_VERSION(fluxndiv, 10);
          GPST_ASSERT_NOT_ARRAY(fluxndiv);
          GPST_CONVERT(fluxndiv, int, fluxNdiv, value);
        } else if (key == "fracndiv") {
          GPST_REQUIRE_MIN_VERSION(fracndiv, 10);
          GPST_ASSERT_NOT_ARRAY(fracndiv);
          GPST_CONVERT(fracndiv, int, piNdiv, value);
        } else if (key == "chindiv") {
          GPST_REQUIRE_MIN_VERSION(chindiv, 10);
          GPST_ASSERT_NOT_ARRAY(chindiv);
          GPST_CONVERT(chindiv, int, chiNdiv, value);
        } else if (key == "xndiv") {
          GPST_REQUIRE_MIN_VERSION(xndiv, 11);
          GPST_ASSERT_NOT_ARRAY(xndiv);
          GPST_CONVERT(xndiv, int, xNdiv, value);
        } else if (key == "fluxlabel") {
          GPST_REQUIRE_MIN_VERSION(fluxlabel, 11);
          GPST_ASSERT_NOT_ARRAY(fluxlabel);
          userFluxLabel = value;
        } else if (key == "fraclabel") {
          GPST_REQUIRE_MIN_VERSION(fraclabel, 11);
          GPST_ASSERT_NOT_ARRAY(fraclabel);
          userFracLabel = value;
        } else if (key == "chilabel") {
          GPST_REQUIRE_MIN_VERSION(chilabel, 11);
          GPST_ASSERT_NOT_ARRAY(chilabel);
          userDirLabel = value;
        } else if (key == "xlabel") {
          GPST_REQUIRE_MIN_VERSION(xlabel, 12);
          GPST_ASSERT_NOT_ARRAY(xlabel);
          userXLabel = value;
        } else if (key == "fluxmorelabels") {
          GPST_REQUIRE_MIN_VERSION(fluxmorelabels, 12);
          GPST_ASSERT_NOT_ARRAY(fluxmorelabels);
          try {
            fluxMoreLabels = interpretBoolean(value);
          } catch (...) {
            log_fatal("On line %d of input file: Expected boolean value but got "
                      "'%s'", linecount, value.c_str());
          }
        } else if (key == "legend") {
          GPST_REQUIRE_MIN_VERSION(legend, 3);
          GPST_ASSERT_ARRAY(legend);
          if (array.size() != 4) {
            log_fatal("On line %d of input file: Number of components in key "
                      "'legend' has to be 4.", linecount);
          }
          double legendxy[4];
          GPST_CONVERT_ARRAY(legend, double, legendxy, array);
          LegendConfig legendconfig = {
            legendxy[0],
            legendxy[1],
            legendxy[2],
            legendxy[3],
            Flux, /* default is top panel */
            1,
            0., /* use a panel-dependent default */
            1, /* use one column by default */
            NULL
          };
          legends.push_back(legendconfig);
        } else if (key == "legendsize") {
          GPST_REQUIRE_MIN_VERSION(legendsize, 10);
          if (legends.empty()) {
            log_error("Found keyword 'legendsize' without 'legend'. Ignoring.");
          } else {
            GPST_ASSERT_NOT_ARRAY(legendsize);
            GPST_USE_DOUBLE(legendsize, legends.back().legendfont);
          }
        } else if (key == "legendpanel") {
          GPST_REQUIRE_MIN_VERSION(legendpanel, 12);
          if (legends.empty()) {
            log_error("Found keyword 'legendpanel' without 'legend'. Ignoring.");
          } else {
            GPST_ASSERT_NOT_ARRAY(legendpanel);
            try {
              legends.back().panel = interpretPanel(value);
            } catch (...) {
              log_fatal("On line %d of input file: Expected panel name ('flux', "
                        "'frac', 'angle') value but got '%s'",
                        linecount, value.c_str());
            }
          }
        } else if (key == "legendcolumns") {
          GPST_REQUIRE_MIN_VERSION(legendcolumns, 13);
          if (legends.empty()) {
            log_error("Found keyword 'legendcolumns' without 'legend'. "
                      "Ignoring.");
          } else {
            GPST_ASSERT_NOT_ARRAY(legendcolumns);
            GPST_CONVERT(legendcolumns, int, legends.back().columns, value);
          }
        } else if (key == "deltachi") {
          GPST_ASSERT_NOT_ARRAY(deltachi);
          GPST_USE_DOUBLE(deltachi, deltachi);
        } else if (key == "energybins") {
          GPST_ASSERT_ARRAY(energybins);
          GPST_CONVERT_ARRAY(energybins, double, bEnergy, array);
          nBins = array.size();
        } else if (key == "phasebins") {
          GPST_ASSERT_ARRAY(phasebins);
          GPST_CONVERT_ARRAY(phasebins, double, bPhase, array);
          nPhaseBins = array.size();
        } else if (key == "component") {
          if (currentComponentLine > 0) {
            int cf = pFlags & (ExpectComponentFlux | ExpectComponentEnergies |
                               ExpectComponentPolarization |
                               ExpectComponentAngle);
            ComponentInfo ci;
            ci.line = currentComponentLine;
            ci.name = currentComponentName;
            ci.flags = cf;
            ci.color = currentComponentColor;
            ci.style = currentComponentStyle;
            ci.width = currentComponentWidth;
            ci.legend = currentComponentLegend;
            component[nComp-1] = ci;
          }
          if (!currentExtraModel.name.empty()) {
            extraModels.push_back(currentExtraModel);
            currentExtraModel.clear();
          }
          customizeMdp = false;
          GPST_ASSERT_NOT_ARRAY(component);
          ++nComp;
          currentComponentLine = linecount;
          currentComponentName = value;
          currentComponentColor = -1;
          currentComponentStyle = -1;
          currentComponentWidth = 2;
          currentComponentLegend = legends.size()-1;
          pFlags |= (ExpectComponentFlux | ExpectComponentEnergies |
                     ExpectComponentPolarization | ExpectComponentAngle);
          GPST_DISABLE_PFLAG(HaveComponentPoints);
          GPST_DISABLE_PFLAG(PreparePhase);
          GPST_DISABLE_PFLAG(PrepareExtraModel);
        } else if (key == "energy") {
          GPST_ASSERT_ARRAY(energy);
          if (pFlags & PrepareExtraModel) {
            GPST_ASSERT_PFLAG(ExpectModelEnergiesOrPhases, energy);
            currentExtraModel.energy.resize(array.size());
            GPST_CONVERT_ARRAY(energy, double, currentExtraModel.energy, array);
          } else {
            GPST_ASSERT_PFLAG(ExpectComponentEnergies, energy);
            GPST_CHECK_LENGTH(energy, nP[nComp-1]);
            GPST_CONVERT_ARRAY(energy, double, dE[nComp-1], array);
          }
          GPST_DISABLE_PFLAG(ExpectComponentEnergies);
        } else if (key == "flux") {
          GPST_ASSERT_PFLAG(ExpectComponentFlux, flux);
          GPST_ASSERT_ARRAY(flux);
          if (pFlags & PreparePhase) {
            GPST_CHECK_LENGTH(flux, nPPhi);
            GPST_CONVERT_ARRAY(flux, double, dPhiI, array);
            phaseFlags &= ~ExpectComponentFlux;
          } else if (pFlags & PrepareExtraModel) {
            currentExtraModel.flux.resize(array.size());
            GPST_CONVERT_ARRAY(flux, double, currentExtraModel.flux, array);
          } else {
            GPST_CHECK_LENGTH(flux, nP[nComp-1]);
            GPST_CONVERT_ARRAY(flux, double, dI[nComp-1], array);
          }
          GPST_DISABLE_PFLAG(ExpectComponentFlux);
        } else if (key == "fraction") {
          GPST_ASSERT_PFLAG(ExpectComponentPolarization, fraction);
          GPST_ASSERT_ARRAY(fraction);
          if (pFlags & PreparePhase) {
            GPST_CHECK_LENGTH(fraction, nPPhi);
            GPST_CONVERT_ARRAY(fraction, double, dPhiPi, array);
            phaseFlags &= ~ExpectComponentPolarization;
          } else if (pFlags & PrepareExtraModel) {
            currentExtraModel.fraction.resize(array.size());
            GPST_CONVERT_ARRAY(fraction, double, currentExtraModel.fraction,
                               array);
          } else {
            GPST_CHECK_LENGTH(fraction, nP[nComp-1]);
            GPST_CONVERT_ARRAY(fraction, double, dPi[nComp-1], array);
          }
          GPST_DISABLE_PFLAG(ExpectComponentPolarization);
        } else if (key == "angle") {
          GPST_ASSERT_PFLAG(ExpectComponentAngle, angle);
          GPST_ASSERT_ARRAY(angle);
          if (pFlags & PreparePhase) {
            GPST_CHECK_LENGTH(angle, nPPhi);
            GPST_CONVERT_ARRAY(fraction, double, dPhiChi, array);
            phaseFlags &= ~ExpectComponentAngle;
          } else if (pFlags & PrepareExtraModel) {
            currentExtraModel.angle.resize(array.size());
            GPST_CONVERT_ARRAY(angle, double, currentExtraModel.angle, array);
          } else {
            GPST_CHECK_LENGTH(angle, nP[nComp-1]);
            GPST_CONVERT_ARRAY(angle, double, dChi[nComp-1], array);
          }
          GPST_DISABLE_PFLAG(ExpectComponentAngle);
        } else if (key == "phase") {
          if (currentComponentLine > 0) {
            int cf = pFlags & (ExpectComponentFlux | ExpectComponentEnergies |
                               ExpectComponentPolarization |
                               ExpectComponentAngle);
            ComponentInfo ci;
            ci.line = currentComponentLine;
            ci.name = currentComponentName;
            ci.flags = cf;
            ci.color = currentComponentColor;
            ci.style = currentComponentStyle;
            ci.width = currentComponentWidth;
            ci.legend = currentComponentLegend;
            component[nComp-1] = ci;
            currentComponentLine = 0;
            currentComponentName = "";
            currentComponentColor = -1;
            currentComponentStyle = -1;
            currentComponentWidth = 2;
          }
          if (!currentExtraModel.name.empty()) {
            if (currentExtraModel.phase.size() > 0 ||
                currentExtraModel.energy.size() > 0)
            {
              extraModels.push_back(currentExtraModel);
              currentExtraModel.clear();
            } else {
              currentExtraModel.phase.resize(array.size());
              GPST_CONVERT_ARRAY(phase, double, currentExtraModel.phase, array);
              GPST_DISABLE_PFLAG(ExpectModelEnergiesOrPhases);
              continue;
            }
          }
          GPST_ASSERT_ARRAY(phase);
          GPST_CONVERT_ARRAY(phase, double, dPhi, array);
          phaseFlags = (PreparePhase | ExpectComponentFlux |
                        ExpectComponentPolarization | ExpectComponentAngle);
          pFlags |= phaseFlags;
          nPPhi = array.size();
          have_phase = true;
          GPST_DISABLE_PFLAG(PrepareExtraModel);
          GPST_DISABLE_PFLAG(ExpectComponentEnergies);
        } else if (key == "model") {
          GPST_REQUIRE_MIN_VERSION(model, 4);
          if (currentComponentLine > 0) {
            int cf = pFlags & (ExpectComponentFlux | ExpectComponentEnergies |
                               ExpectComponentPolarization |
                               ExpectComponentAngle);
            ComponentInfo ci;
            ci.line = currentComponentLine;
            ci.name = currentComponentName;
            ci.flags = cf;
            ci.color = currentComponentColor;
            ci.style = currentComponentStyle;
            ci.width = currentComponentWidth;
            ci.legend = currentComponentLegend;
            component[nComp-1] = ci;
            currentComponentLine = 0;
            currentComponentColor = -1;
            currentComponentStyle = -1;
            currentComponentWidth = 2;
          }
          customizeMdp = false;
          if (!currentExtraModel.name.empty()) {
            extraModels.push_back(currentExtraModel);
          }
          GPST_ASSERT_NOT_ARRAY(model);
          currentExtraModel.clear();
          currentExtraModel.line = linecount;
          currentExtraModel.name = value;
          currentExtraModel.legend = legends.size()-1;
          GPST_DISABLE_PFLAG(PreparePhase);
          GPST_DISABLE_PFLAG(ExpectComponentEnergies);
          pFlags |= (PrepareExtraModel | ExpectModelEnergiesOrPhases |
                     ExpectComponentFlux | ExpectComponentPolarization |
                     ExpectComponentAngle);
        } else if (key == "seed") {
          GPST_REQUIRE_MIN_VERSION(seed, 5);
          GPST_ASSERT_NOT_ARRAY(seed);
          GPST_CONVERT(seed, int, random_seed, value);
        } else if (key == "color") {
          GPST_REQUIRE_MIN_VERSION(color, 6);
          GPST_ASSERT_NOT_ARRAY(color);
          Color_t color;
          try {
            color = interpretColor(value);
          } catch (...) {
            log_fatal("On line %d of input file: Invalid color '%s'",
                      linecount, value.c_str());
          }
          if (currentComponentLine > 0) {
            currentComponentColor = color;
          } else if (!currentExtraModel.name.empty()) {
            currentExtraModel.color = color;
          } else if (customizeMdp) {
            mdpColor = color;
          } else if (modelOrSimulation == CustomizeSimulation) {
            datacolor = color;
          } else {
            linecolor = color;
          }
        } else if (key == "style") {
          GPST_REQUIRE_MIN_VERSION(style, 6);
          GPST_ASSERT_NOT_ARRAY(style);
          Style_t style;
          try {
            style = interpretStyle(value);
          } catch (...) {
            log_fatal("On line %d of input file: Invalid line style '%s'",
                      linecount, value.c_str());
          }
          if (currentComponentLine > 0) {
            currentComponentStyle = style;
          } else if (!currentExtraModel.name.empty()) {
            currentExtraModel.style = style;
          } else if (customizeMdp) {
            mdpStyle = style;
          } else {
            linestyle = style;
          }
        } else if (key == "width") {
          GPST_REQUIRE_MIN_VERSION(width, 6);
          GPST_ASSERT_NOT_ARRAY(width);
          Width_t width;
          GPST_CONVERT(width, Width_t, width, value);
          if (currentComponentLine > 0) {
            currentComponentWidth = width;
          } else if (!currentExtraModel.name.empty()) {
            currentExtraModel.width = width;
          } else if (customizeMdp) {
            mdpWidth = width;
          } else if (modelOrSimulation == CustomizeSimulation) {
            datawidth = width;
          }  else {
            linewidth = width;
          }
        } else if (key == "datalabel") {
          GPST_REQUIRE_MIN_VERSION(datalabel, 7);
          GPST_ASSERT_NOT_ARRAY(datalabel);
          datalabel = value;
          data_legend = legends.size()-1;
          modelOrSimulation = CustomizeSimulation;
        } else if (key == "mirror") {
          GPST_REQUIRE_MIN_VERSION(mirror, 8);
          GPST_ASSERT_NOT_ARRAY(mirror);
          log_fatal("On line %d of input file: Keyword 'mirror' is no longer "
                    "supported by GPST version 14 and later", linecount);
        } else if (key == "rings") {
          GPST_REQUIRE_MIN_VERSION(rings, 8);
          GPST_ASSERT_NOT_ARRAY(rings);
          log_fatal("On line %d of input file: Keyword 'rings' is no longer "
                    "supported by GPST version 14 and later", linecount);
        } else if (key == "showmdp") {
          GPST_REQUIRE_MIN_VERSION(showmdp, 10);
          GPST_ASSERT_NOT_ARRAY(showmdp);
          try {
            showmdp = interpretBoolean(value);
          } catch (...) {
            log_fatal("On line %d of input file: Expected boolean value but got "
                      "'%s'", linecount, value.c_str());
          }
          mdp_legend = legends.size()-1;
          customizeMdp = true;
          if (currentComponentLine > 0) {
            int cf = pFlags & (ExpectComponentFlux | ExpectComponentEnergies |
                               ExpectComponentPolarization | ExpectComponentAngle);
            ComponentInfo ci;
            ci.line = currentComponentLine;
            ci.name = currentComponentName;
            ci.flags = cf;
            ci.color = currentComponentColor;
            ci.style = currentComponentStyle;
            ci.width = currentComponentWidth;
            ci.legend = currentComponentLegend;
            component[nComp-1] = ci;
            currentComponentLine = 0;
          }
          if (!currentExtraModel.name.empty()) {
            extraModels.push_back(currentExtraModel);
            currentExtraModel.clear();
          }
        } else if (key == "label") {
          GPST_REQUIRE_MIN_VERSION(label, 9);
          GPST_ASSERT_NOT_ARRAY(label);
          if (!customizeMdp) {
            log_fatal("On line %d of input file: Keyword 'label' must follow "
                      "keyword 'showmdp'.", linecount);
          }
          mdpLabel = value;
        } else if (key == "large") {
          GPST_REQUIRE_MIN_VERSION(large, 10);
          GPST_ASSERT_NOT_ARRAY(large);
          try {
            large = interpretBoolean(value);
          } catch (...) {
            log_fatal("On line %d of input file: Expected boolean value but got "
                      "'%s'", linecount, value.c_str());
          }
        } else if (key == "clean") {
          GPST_REQUIRE_MIN_VERSION(clean, 10);
          GPST_ASSERT_NOT_ARRAY(clean);
          try {
            plot_clean = interpretBoolean(value);
          } catch (...) {
            log_fatal("On line %d of input file: Expected boolean value but got "
                      "'%s'", linecount, value.c_str());
          }
        } else if (key == "detectoremin") {
          GPST_REQUIRE_MIN_VERSION(detectoremin, 14);
          GPST_ASSERT_NOT_ARRAY(detectoremin);
          GPST_USE_DOUBLE(detectoremin, eSensitivityMin);
        } else if (key == "detectordeltae") {
          GPST_REQUIRE_MIN_VERSION(detectordeltae, 14);
          GPST_ASSERT_NOT_ARRAY(detectordeltae);
          GPST_USE_DOUBLE(detectordeltae, deltaESensitivity);
        } else if (key == "detectornbins") {
          GPST_REQUIRE_MIN_VERSION(detectornbins, 14);
          GPST_ASSERT_NOT_ARRAY(detectornbins);
          GPST_CONVERT(detectornbins, int, NumBins, value);
        } else if (key == "sensitivityfile") {
          GPST_REQUIRE_MIN_VERSION(sensitivityfile, 14);
          GPST_ASSERT_NOT_ARRAY(sensitivityfile);
          sensitivityFile = value;
        } else if (key == "sensitivityscale") {
          GPST_REQUIRE_MIN_VERSION(sensitivityscale, 14);
          GPST_ASSERT_NOT_ARRAY(sensitivityscale);
          GPST_USE_DOUBLE(sensitivityscale, sensitivityScale);
        } else if (key == "backgroundfile") {
          GPST_REQUIRE_MIN_VERSION(backgroundfile, 14);
          GPST_ASSERT_NOT_ARRAY(backgroundfile);
          backgroundFile = value;
        } else if (key == "modulationfactor") {
          GPST_REQUIRE_MIN_VERSION(modulationfactor, 14);
          GPST_ASSERT_NOT_ARRAY(modulationfactor);
          GPST_USE_DOUBLE(modulationfactor, modulation_factor);
          have_modulationfactor_keyword = true;
        } else if (key == "modulationfile") {
          GPST_REQUIRE_MIN_VERSION(modulationfile, 14);
          GPST_ASSERT_NOT_ARRAY(modulationfile);
          modulationFile = value;
        } else {
          log_fatal("On line %d of input file: Unknown key '%s'.", linecount,
                    key.c_str());
        }

#undef GPST_ASSERT_ARRAY
#undef GPST_ASSERT_NOT_ARRAY
#undef GPST_ASSERT_PFLAG
#undef GPST_DISABLE_PFLAG
#undef GPST_CONVERT
#undef GPST_CHECK_LENGTH
#undef GPST_CONVERT_ARRAY
    }
  }

    log_info("Closing input file '%s'", lastfilename.c_str());
    delete fp_in;
    if (!filestack.empty()) {
      fp_in = filestack.top();
      filestack.pop();
      lastfilename = filenamestack.top();
      filenamestack.pop();
      linecount = linecountstack.top();
      linecountstack.pop();
    } else {
      fp_in = NULL;
    }
  }

  const double mdpFactor = (fractionunit == Frac ? 4.29 : 429);

  double *nS = new double[NumBins];
  double *nBG = new double[NumBins];
  double *mu = new double[NumBins];

  loadSensitivity(nS, nBG, mu, NumBins, modulation_factor, eSensitivityMin,
                  deltaESensitivity, sensitivityFile, backgroundFile,
                  modulationFile, sensitivityScale);
  log_info("Loaded the sensitivity information");

  if (currentComponentLine > 0) {
    int cf = pFlags & (ExpectComponentFlux | ExpectComponentEnergies |
                       ExpectComponentPolarization | ExpectComponentAngle);
    ComponentInfo ci;
    ci.line = currentComponentLine;
    ci.name = currentComponentName;
    ci.flags = cf;
    ci.color = currentComponentColor;
    ci.style = currentComponentStyle;
    ci.width = currentComponentWidth;
    ci.legend = currentComponentLegend;
    component[nComp-1] = ci;
  }
  if (!currentExtraModel.name.empty()) {
    extraModels.push_back(currentExtraModel);
  }

  // apply on time fraction
  time *= ontime;

  // file is read, now it's time for sanity checks and summary output
  if (time < 0) {
    log_fatal("On input: time parameter not specified");
  }

  // check that energy and phase bins are consistent
  if (nPhaseBins > 0) {
    if (nBins != 2) {
      log_fatal("Phase bins have been specified. In this case exactly 2 energy "
                "bins must be given, but I got %d bins.", nBins);
    }
    if (!have_phase) {
      log_fatal("Phase bins have been specified but no phase-dependent "
                "polarization model.");
    }
  }

  bool fatal_error = false;
  for (int i = 0; i < nComp; ++i) {
    if (component[i].flags & ExpectComponentFlux) {
      log_error("Component '%s' defined on line %d of input file: no flux "
                "defined.", component[i].name.c_str(), component[i].line);
      fatal_error = true;
    }
    if (!have_phase) {
      if (component[i].flags & ExpectComponentPolarization) {
        log_error("Component '%s' defined on line %d of input file: no "
                  "polarization fraction defined and this is no phase-dependent "
                  "simulation.", component[i].name.c_str(), component[i].line);
        fatal_error = true;
      }
      if (component[i].flags & ExpectComponentAngle) {
        log_error("Component '%s' defined on line %d of input file: no "
                  "polarization angle defined and this is no phase-dependent "
                  "simulation.", component[i].name.c_str(), component[i].line);
        fatal_error = true;
      }
    }
  }

  if (have_phase) {
    if (phaseFlags & ExpectComponentFlux) {
      log_error("No phase-dependent flux defined.");
      fatal_error = true;
    }
    if (phaseFlags & ExpectComponentPolarization) {
      log_error("No phase-dependent polarization fraction defined.");
      fatal_error = true;
    }
    if (phaseFlags & ExpectComponentAngle) {
      log_error("No phase-dependent polarization angle defined.");
      fatal_error = true;
    }
  }

  if (fatal_error) log_fatal("Critical error(s) in input file.");

  if (have_modulationfactor_keyword && !modulationFile.empty()) {
    log_fatal("Input file contains both 'modulationfactor' and 'modulationfile'"
              " keywords. Please use only one of them.");
  }

  log_info("Observation time %.1lf day(s) (on-time)", time/86400);
  log_info("Number of spectral components: %d", nComp);
  log_info("Number of energy bins: %d", nBins-1);

  for (int i = 0; i < nBins; i++) {
    log_debug("Energy %d = %lgkeV", i+1, bEnergy[i]);
  }

  double modelEmin = numeric_limits<double>::max();
  double modelEmax = 0;

  for (int i = 0; i < nComp; i++) {
    log_info("Component #%d has %d data points.", i+1, nP[i]);
        
    for (int j=0; j < nP[i]; j++)
    {
      log_debug("Energy  #%d = %lgkeV", j+1, dE[i][j]);
      if (!have_phase) {
        modelEmin = std::min(modelEmin, dE[i][j]);
        modelEmax = std::max(modelEmax, dE[i][j]);
      }
      log_debug("Flux    #%d = %lgcm-2 s-1 keV-1", j+1, dI[i][j]);
      log_debug("Pol-Deg #%d = %lg%%", j+1, dPi[i][j]*100);
      log_debug("Pol-Dir #%d = %lgdegree", j+1, dChi[i][j]);
      dChi[i][j] += deltachi;
      if (dChi[i][j] > 90.) dChi[i][j] -= 180.;
      else if (dChi[i][j] < -90.) dChi[i][j] += 180.;
      dChi[i][j] *= M_PI/180.;
      dQ[i][j] = dI[i][j]*dPi[i][j]*cos(2.*dChi[i][j]);
      dU[i][j] = dI[i][j]*dPi[i][j]*sin(2.*dChi[i][j]);
    }
  }

  double modelXmin = have_phase ? dPhi[0] : modelEmin;
  double modelXmax = have_phase ? dPhi[nPPhi-1] : modelEmax;

    
  for (vector<ExtraModel>::iterator it = extraModels.begin();
       it != extraModels.end(); ++it)
  {
    if (it->energy.empty() && it->phase.empty()) {
      log_fatal("Model '%s' defined on line %d of input file: neither energy "
                "nor phase points defined.", it->name.c_str(), it->line);
    }

    if (it->energy.size() > 0 && it->phase.size() > 0) {
      log_fatal("Model '%s' defined on line %d of input file: both energy "
                "and phase points defined.", it->name.c_str(), it->line);
    }

    if (it->energy.size() > 0 && have_phase) {
      log_fatal("Model '%s' defined on line %d of input file: energy points "
                "given but this is a phase-dependent simulation.",
                it->name.c_str(), it->line);
    } else if (it->phase.size() && !have_phase) {
      log_fatal("Model '%s' defined on line %d of input file: phase points "
                "given but this is a energy-dependent simulation.",
                it->name.c_str(), it->line);
    }

    size_t npoints = max(it->energy.size(), it->phase.size());

    if (it->flux.empty() && it->fraction.empty() && it->angle.empty()) {
      log_fatal("Model '%s' defined on line %d of input file: empty model. "
                "Neither flux, polarization fraction, nor polarization angle "
                "given.", it->name.c_str(), it->line);
    }

    if (it->flux.size() > 0 && it->flux.size() != npoints) {
      log_fatal("Model '%s' defined on line %d of input file: number of flux "
                "points does not match.", it->name.c_str(), it->line);
    }

    if (it->fraction.size() > 0 && it->fraction.size() != npoints) {
      log_fatal("Model '%s' defined on line %d of input file: number of "
                "polarization fraction points does not match.",
                it->name.c_str(), it->line);
    }

    if (it->angle.size() > 0 && it->angle.size() != npoints) {
      log_fatal("Model '%s' defined on line %d of input file: number of "
                "polarization angle points does not match (expected %zu, got "
                "%zu).", it->name.c_str(), it->line, npoints, it->angle.size());
    }

    for (size_t i = 0; i < it->angle.size(); ++i) {
      it->angle[i] += deltachi;
      if (it->angle[i] > 90.) it->angle[i] -= 180.;
      else if (it->angle[i] < -90.) it->angle[i] += 180.;
    }

    if (fractionunit == Percent) {
      for (size_t i = 0; i < it->fraction.size(); ++i)
        it->fraction[i] *= 100;
    }
  }

  double xmin = have_phase ? bPhase[0] : bEnergy[0];
  double xmax = have_phase ? bPhase[nPhaseBins-1] : bEnergy[nBins-1];
    
  style();
  Dimensions d;
  if (large) {
    d.frameHeight = 240;
    d.frameWidth = 564;
    d.leftMargin = 125;
    d.rightMargin = 40;
    d.bottomMargin = 85;
    d.topMargin = 22;
    d.textHeight = 32;
    d.yoffset = 22;
    gStyle->SetFrameLineWidth(2);
  } else {
    d.frameHeight = 240;
    d.frameWidth = 564;
    d.leftMargin = 88;
    d.rightMargin = 20;
    d.bottomMargin = 58;
    d.topMargin = 22;
    d.textHeight = 22;
    d.yoffset = 19;
  }
  d.width = d.frameWidth + d.leftMargin + d.rightMargin;
  d.height = panels.size() * d.frameHeight + d.bottomMargin + d.topMargin;

  vector<PanelConfig> panelConfig;
  {
    double top = 1;
    for (size_t i = 1; i <= panels.size(); ++i) {
      panelConfig.push_back(makePanelConfig(i, panels.size(), top, d, large));
    }
  }
    
  TCanvas *c00[4];
  char fname[50];
  for (int i=0;i<1;i++)
  {
    sprintf(fname,"c00%1d",i);
    c00[i] = new TCanvas(fname, fname, d.width, d.height);
    c00[i]->SetWindowSize(d.width + (d.width - c00[i]->GetWw()),
                          d.height + (d.height - c00[i]->GetWh()));
    c00[i]->SetMargin(0, 0, 0, 0);
    c00[i]->SetFixedAspectRatio();
  }
  
  TGraph *sumFlux;
  TGraph *grFlux[MaxC];
  TGraph *grPi[MaxC];
  TGraph *grChi[MaxC];

  TMultiGraph *mgrFlux = new TMultiGraph();
  TMultiGraph *mgrPi = new TMultiGraph();
  TMultiGraph *mgrChi = new TMultiGraph();

  TGraph *tPi;
  TGraph *tChi;

  double ymin=99.;
  double ymax=0.;
  double Pimax=0.;
  double chimin=+90.;
  double chimax=-90.;

    
  /* go over all bins and sum up information over 0.2 keV bins */
  int points=0;
    
  vector<double> x(MaxBins), x2(MaxBins), y2(MaxBins);
  vector<double> yI(MaxBins), yQ(MaxBins), yU(MaxBins), yPi(MaxBins), yChi(MaxBins), yPi2(MaxBins), yChi2(MaxBins);
  

  for (int i=0; i<nBins-1; i++) {
    for (double energy = bEnergy[i]; energy <= bEnergy[i+1]; energy += deltaESensitivity) {
      double eC=energy + 0.5 * deltaESensitivity;

      double isum=0;
      double qsum=0;
      double usum=0;

      for (int j = 0; j < nComp; j++) {
        double I, Q, U;
        interpolate_stokes(dE[j], dI[j], dQ[j], dU[j], nP[j], eC, I, Q, U);
        isum+=I;
        qsum+=Q;
        usum+=U;
      }
            
      x[points]=eC;
      yI[points]=isum;
      yQ[points]=qsum;
      yU[points]=usum;

//      log_debug("E = %g ==> I = %g, Q = %g, U = %g", eC, isum, qsum, usum);

      if (yI[points]>0.) {
        yPi[points]=sqrt(yQ[points]*yQ[points]+yU[points]*yU[points])/yI[points];
        yChi[points]=getChi(qsum,usum);
        if (yPi[points]>Pimax) Pimax=yPi[points];
        if (yChi[points]>chimax) chimax=yChi[points];
        if (yChi[points]<chimin) chimin=yChi[points];                
      }
            
      points++;
    }
    log_debug("Energy bin [%g -- %g] ==> points %d",
              bEnergy[i], bEnergy[i+1], points);
  }

    
  if (aFlux > 0.) {
    double average1=0.;
    double average2=0.;
    for (int i = 0; i < NumBins; ++i) {
      double eC =
        eSensitivityMin + 0.5*deltaESensitivity + i * deltaESensitivity;
      if (eC < scaleEmin) continue;
      if (eC > scaleEmax) break;
      
      double isum = 0;
      for (int j = 0; j < nComp; ++j) {
        isum += get_table(dE[j], dI[j], nP[j], eC);
      }
      
      if (isum > 0.) {
        /* we divide by eC here to get the flux - rather than the nu-Fnu flux */
        average1 += isum/eC;
        average2 += Crab(eC)/eC; /* do the same for the Crab */
      }
    }
    
    log_info("The fluxes you entered correspond to: %.3lg mCrab", average1/average2*1000.);
    log_info("Renormalizing flux to %.3lg mCrab", aFlux);

    double scale = aFlux/(average1/average2)/1000.;
    log_info("Scaling factor: %.3f", scale);
 
    for (int i = 0; i < points; i++)
    {
      yI[i] *= scale;
      yQ[i] *= scale;
      yU[i] *= scale;
    }
        
    for (int i=0; i<nComp; i++) 
      for (int j=0;j<nP[i];j++)
        dI[i][j] *= scale;
  }

  if (have_phase) {
    for (int i = 0; i < nPPhi; ++i) {
      dPhiChi[i] += deltachi;
      if (dPhiChi[i] > 90) dPhiChi[i] -= 180;
      else if (dPhiChi[i] < -90) dPhiChi[i] += 180;
    }

    // integrate total flux from lower to upper energy bin
    // remember: a phase-dependent model requires exactly 1 energy bin
    double fluxsum = 0;
    for (int i = 0; i < points; ++i) {
      double eC2 = math::sqr(x[i]);
      fluxsum += deltaESensitivity*yI[i]/eC2;
    }
    fluxsum *= math::sqr(exp(0.5*(log(bEnergy[0])+log(bEnergy[1]))))/(bEnergy[1]-bEnergy[0]);
    double displayflux[MaxP];
    for (int i = 0; i < nPPhi; ++i) displayflux[i] = dPhiI[i] * fluxsum;


    for (vector<ExtraModel>::iterator it = extraModels.begin();
         it != extraModels.end(); ++it)
    {
      for (size_t i = 0; i < it->flux.size(); ++i) {
        it->flux[i] *= fluxsum;
      }
    }

    sumFlux = new TGraph(nPPhi, dPhi, displayflux);
    tPi = new TGraph(nPPhi, dPhi, dPhiPi);
    tChi = new TGraph(nPPhi, dPhi, dPhiChi);
    Pimax = *max_element(dPhiPi, dPhiPi+nPPhi);

    for (int i = 0; i < nPPhi; ++i) {
      dPhiChi[i] *= M_PI/180;
    }
  } else {
    log_debug("Number of points for graph: %d", points);
    sumFlux = new TGraph(points, &x[0], &yI[0]);
    tPi  = new TGraph(points, &x[0], &yPi[0]);
    tChi = new TGraph(points, &x[0], &yChi[0]);
  }

  for (size_t i = 0; i < legends.size(); ++i) {
    for (legends[i].panel_index = 0;
         legends[i].panel_index < panels.size();
         ++legends[i].panel_index)
    {
      if (panels[legends[i].panel_index] == legends[i].panel) break;
    }
    ++legends[i].panel_index;  // because in root sub-panels are counted from 1
    if (legends[i].panel_index > panels.size()) {
      log_error("%zuth legend should be drawn on a panel that is not shown", i+1);
      legends[i].panel_index = 0;
    }
    legends[i].legend = new TLegend(legends[i].x1, legends[i].y1,
                                    legends[i].x2, legends[i].y2);
    legends[i].legend->SetBorderSize(0);
    legends[i].legend->SetFillStyle(0);
    legends[i].legend->SetNColumns(legends[i].columns);
    log_debug("Legend %zu size %f", i, legends[i].legendfont);
    if (legends[i].legendfont > 0) {
      legends[i].legend->SetTextSize(legends[i].legendfont);
    } else {
      if (legends[i].panel_index > 0) {
        legends[i].legend->SetTextSize(panelConfig[legends[i].panel_index-1].fontsize);
      }  // else we don't care because the legend won't be shown anyways.
    }
    legends[i].legend->SetTextFont(42);
    legends[i].legend->SetName(TString::Format("legend%zu", i));
  }

  while (sumFlux->GetX()[0] < modelXmin)
    sumFlux->RemovePoint(0);
  while (sumFlux->GetX()[sumFlux->GetN()-1] > modelXmax)
    sumFlux->RemovePoint(sumFlux->GetN()-1);

  if (fluxunit == CGSe9) {
    for (int i = 0; i < sumFlux->GetN(); ++i) {
      double px, py;
      sumFlux->GetPoint(i, px, py);
      sumFlux->SetPoint(i, px, 1e9 * py);
    }
  }
  sumFlux->SetLineWidth(linewidth);            
  sumFlux->SetLineColor(linecolor);
  sumFlux->SetMarkerColor(linecolor);
  sumFlux->SetLineStyle(linestyle);
  mgrFlux->Add(sumFlux);

  if (!simulation_name.empty() && model_legend >= 0) {
    legends[model_legend].legend->AddEntry(sumFlux, simulation_name.c_str(), "L");
  }
    
  while (tPi->GetX()[0] < modelXmin) tPi->RemovePoint(0);
  while (tPi->GetX()[tPi->GetN()-1] > modelXmax) tPi->RemovePoint(tPi->GetN()-1);
  if (fractionunit == Percent) {
    for (int i = 0; i < tPi->GetN(); ++i) {
      double px, py;
      tPi->GetPoint(i, px, py);
      tPi->SetPoint(i, px, 100*py);
    }
  }
  tPi->SetLineWidth(linewidth);            
  tPi->SetLineColor(linecolor);
  tPi->SetMarkerColor(linecolor);
  tPi->SetLineStyle(linestyle);
  mgrPi->Add(tPi);
    
  while (tChi->GetX()[0] < modelXmin) tChi->RemovePoint(0);
  while (tChi->GetX()[tChi->GetN()-1] > modelXmax) tChi->RemovePoint(tChi->GetN()-1);
  tChi->SetLineWidth(linewidth);            
  tChi->SetLineColor(linecolor);
  tChi->SetMarkerColor(linecolor);
  tChi->SetLineStyle(linestyle);
  mgrChi->Add(tChi);
    
  for (int i=0;i<points;i++) 
  {
    if ((yI[i]<ymin)&&(yI[i]>0.)) ymin=yI[i];
    if (yI[i]>ymax) ymax=yI[i];
  }
    
  ymin/=3.4;
  ymax*=1.2;
  if (!have_phase) {
    xmin/=1.1;
    xmax*=1.1;
  }

  if (!have_phase && extraModels.empty()) {
    int style_count = 0;
    for (int i=0; i<nComp; i++) {
      for (int j=0;j<nP[i];j++) {
        x2[j]=dE[i][j];
        y2[j]=dI[i][j];
        yPi2[j]=dPi[i][j];
        yChi2[j]=dChi[i][j]*180./PI;     
        
        if (yChi2[j]>chimax) chimax=yChi2[j];
        if (yChi2[j]<chimin) chimin=yChi2[j];
      }

      if (component[i].style == -1 || component[i].color == -1)
        ++style_count;

      grFlux[i] = new TGraph(nP[i], &x2[0], &y2[0]);
      grFlux[i]->SetLineWidth(component[i].width);    
      if (fluxunit == CGSe9) {
        for (int j = 0; j < grFlux[i]->GetN(); ++j) {
          double px, py;
          grFlux[i]->GetPoint(j, px, py);
          grFlux[i]->SetPoint(j, px, 1e9 * py);
        }
      }

      grPi[i] = new TGraph(nP[i], &x2[0], &yPi2[0]);
      grPi[i]->SetLineWidth(component[i].width);
      if (fractionunit == Percent) {
        for (int j = 0; j < grPi[i]->GetN(); ++j) {
          double px, py;
          grPi[i]->GetPoint(j, px, py);
          grPi[i]->SetPoint(j, px, 100*py);
        }
      }

      grChi[i] = new TGraph(nP[i], &x2[0], &yChi2[0]);
      grChi[i]->SetLineWidth(component[i].width);

      if (component[i].style >= 0) {
        grFlux[i]->SetLineStyle(component[i].style);
        grPi[i]->SetLineStyle(component[i].style);
        grChi[i]->SetLineStyle(component[i].style);
      } else {
        grFlux[i]->SetLineStyle(style_count);
        grPi[i]->SetLineStyle(style_count);
        grChi[i]->SetLineStyle(style_count);
      }
      if (component[i].color >= 0) {
        grFlux[i]->SetLineColor(component[i].color);
        grFlux[i]->SetMarkerColor(component[i].color);
        grPi[i]->SetLineColor(component[i].color);
        grPi[i]->SetMarkerColor(component[i].color);
        grChi[i]->SetLineColor(component[i].color);
        grChi[i]->SetMarkerColor(component[i].color);
      } else {
        grFlux[i]->SetLineColor(style_count+1);
        grFlux[i]->SetMarkerColor(style_count+1);
        grPi[i]->SetLineColor(style_count+1);
        grPi[i]->SetMarkerColor(style_count+1);
        grChi[i]->SetLineColor(style_count+1);
        grChi[i]->SetMarkerColor(style_count+1);
      }

      if (nComp>1) {
        mgrFlux->Add(grFlux[i]);
        mgrPi->Add(grPi[i]);
        mgrChi->Add(grChi[i]);
        if (component[i].legend >= 0) {
          legends[component[i].legend].legend
            ->AddEntry(grFlux[i], component[i].name.c_str(), "L");
        }
      }
    }
  } else if (!extraModels.empty()) {
    int i = 0;
    int style_count = 0;
    for (vector<ExtraModel>::const_iterator it = extraModels.begin();
         it != extraModels.end(); ++it, ++i)
    {
      const std::vector<double> &xPoints = have_phase ? it->phase : it->energy;

      if (it->color == -1 || it->style == -1) ++style_count;

      Width_t width = it->width;
      Color_t color = it->color == -1 ? style_count+1 : it->color;
      Style_t style = it->style == -1 ? style_count : it->style;

      TGraph *gRef = NULL;
      if (!it->flux.empty()) {
        grFlux[i] = new TGraph(xPoints.size(), &xPoints[0], &it->flux[0]);
        grFlux[i]->SetLineWidth(width);
        grFlux[i]->SetLineStyle(style);
        grFlux[i]->SetLineColor(color);   
        grFlux[i]->SetMarkerColor(color);   
        mgrFlux->Add(grFlux[i]);
        gRef = grFlux[i];
      }

      if (!it->fraction.empty()) {
        grPi[i] = new TGraph(xPoints.size(), &xPoints[0], &it->fraction[0]);
        grPi[i]->SetLineWidth(width);    
        grPi[i]->SetLineStyle(style);
        grPi[i]->SetLineColor(color);   
        grPi[i]->SetMarkerColor(color);   
        mgrPi->Add(grPi[i]);  
        gRef = grPi[i];
      }

      if (!it->angle.empty()) {
        grChi[i] = new TGraph(xPoints.size(), &xPoints[0], &it->angle[0]);
        grChi[i]->SetLineWidth(width);    
        grChi[i]->SetLineStyle(style);
        grChi[i]->SetLineColor(color);   
        grChi[i]->SetMarkerColor(color);   
        mgrChi->Add(grChi[i]);  
        gRef = grChi[i];
      }

      if ((it->legend >= 0) && gRef)
        legends[it->legend].legend->AddEntry(gRef, it->name.c_str(), "L");
    }
  }


  double xc[MaxBins], xem[MaxBins], xep[MaxBins], ycFlux[MaxBins],
    yeFlux[MaxBins], ycPi[MaxBins], yePiLow[MaxBins], yePiHigh[MaxBins],
    ycChi[MaxBins], yeChi[MaxBins], ycMdp[MaxBins];
    

  int kd=0;
  if (have_phase) {
    // first, we compute the statistics in the energy range, without scaling
    double nSignal=0;
    double nBackground=0.;
    double fluxInt = 0;
    int intP = 0;
    double debug_crab_rate = 0;
    double debug_signal_rate = 0;
    double modulationFactor = 0;
    for (double energy=bEnergy[0]; energy<=bEnergy[1];
         energy += deltaESensitivity, ++intP)
    {
      double eC = energy + 0.5*deltaESensitivity;
      double crab = Crab(eC);
      double weight = time*get_table(&x[0], &yI[0], points, eC)/crab;
      
      // calculate bin
      int bin = (energy-eSensitivityMin)/deltaESensitivity;
      if (bin<0) { log_fatal("bin < 0: %lg", energy); }
      if (bin>=NumBins) { log_fatal("bin >= NumBins: %lg", energy); }
        
      nSignal += weight * nS[bin];
      nBackground += time*nBG[bin];
      double this_flux =
        deltaESensitivity * get_table(&x[0], &yI[0], points, eC)/(eC*eC);
      fluxInt += this_flux;
      modulationFactor += this_flux * mu[bin];

      // if (eC > 2 && eC <= 12) {
        debug_crab_rate += nS[bin];
        debug_signal_rate += weight*nS[bin]/time;
        log_debug("Energy %.1fkeV: Weight = %g, nS = %g", eC, weight/time, nS[bin]);
      // }
    }
    modulationFactor /= fluxInt;

    log_info("Integral source rate: %.3fHz", nSignal/time);
    log_info("Integral crab rate: %.3fHz", debug_crab_rate);
    log_info("Integral background rate: %.3fHz", nBackground/time);
    log_info("Modulation factor: %.3f", modulationFactor);

    double eC = exp(0.5*(log(bEnergy[0]) + log(bEnergy[1])));
    fluxInt *= eC*eC/(deltaESensitivity*intP);

    for (int i = 0; i < nPhaseBins-1; ++i) {
      xc[kd] = 0.5*(bPhase[i] + bPhase[i+1]);
      double binwidth = (bPhase[i+1] - bPhase[i]);
      xep[kd] = xem[kd] = 0.5*binwidth;

      double binfraction = binwidth/(xmax-xmin);

      double qsum, usum, isum;
      linearIntegralQUI(bPhase[i], bPhase[i+1], nPPhi, dPhi, dPhiI, dPhiPi,
                        dPhiChi, qsum, usum, isum);
      double pi = sqrt(qsum*qsum + usum*usum)/isum;
      double chi = getChi(qsum, usum) * M_PI/180;
      double flux = isum/(bPhase[i+1]-bPhase[i]);

      double thisSignal = nSignal*flux*binfraction;
      double thisBackground = nBackground*binfraction;

      log_debug("Bin %d, Phase %.2lf: Nsrc = %.1lf, Nbg = %.1lf",
                i+1, xc[kd], thisSignal, thisBackground);

      double FluxM, PiM, ChiM, FluxE, PiELow, PiEHigh, ChiE;
      getValuesAndErrors(thisSignal, thisBackground, modulationFactor, pi, chi,
                         FluxM, PiM, ChiM, FluxE, PiELow, PiEHigh, ChiE);
      FluxM /= binfraction;
      FluxE /= binfraction;
      double mdp = mdpFactor * sqrt(thisSignal + thisBackground) /
                     (modulationFactor * thisSignal);

      log_debug("Bin %d, Phase %.2lf: aT = %.2lf, aM = %.2lf, mu = %.2lf, "
                "chi = %.2lf, ChiM = %.2lf, ChiE = %.2lf", i+1, xc[kd], pi, PiM,
                modulationFactor, chi*180/M_PI, ChiM*180/M_PI, ChiE*180/M_PI);

      double eFlux = sqrt(thisSignal + 2*thisBackground)/thisSignal;
      log_debug("sigmaFlux = %.2lf%%", 100*eFlux);
      ycFlux[kd] = getRNG().Gaus(fluxInt*flux, fluxInt*flux * eFlux);
      yeFlux[kd] = ycFlux[kd] * eFlux;
      ycPi[kd]=PiM;
      ycChi[kd]=ChiM*180./PI;
      ycMdp[kd]=mdp;
      yePiLow[kd]=PiELow;
      yePiHigh[kd]=PiEHigh;
      yeChi[kd]=ChiE*180./PI;
      if (PiM > 2*PiELow) kd++;
    }
  } else {
    double nSignalTotal = 0;
    log_debug("Energy dependent model with %d bins", nBins);
    for (int i=0;i<nBins-1;i++)
    {
      double eC=exp((log(bEnergy[i])+log(bEnergy[i+1]))/2.);
      xc[kd]=eC;
      xep[kd]=bEnergy[i+1]-eC;
      xem[kd]=eC-bEnergy[i];
      
      // Now, let's get the theoretical pol degree and pol. direction at the
      // centers of all bins
      double tPiV  = get_table(&x[0], &yPi[0], points, eC);
      double tChiV = get_table(&x[0], &yChi[0], points, eC);
      
      // and now we compute the statistics in that bin
      double nSignal=0;
      double nBackground=0.;
      double fluxInt=0;
      double modulationFactor = 0;
      int intP = 0;
      for (double energy=bEnergy[i];
           energy<=bEnergy[i+1];
           energy += deltaESensitivity, ++intP)
      {
        double eC2=energy + 0.5 * deltaESensitivity;
        double crab=Crab(eC2);
        double weight = time*get_table(&x[0], &yI[0], points, eC2)/crab;
        //log_debug("eC2 = %g, w = %g", eC2, get_table(&x[0], &yI[0], points, eC2));
        
        // calculate bin
        int bin = int(floor((eC2-eSensitivityMin)/deltaESensitivity));
        if (bin<0) { log_fatal("bin < 0: %lg", energy); }
        if (bin>=NumBins) { log_fatal("bin >= NumBins: %lg", energy); }

        nSignal += weight * nS[bin];
        nBackground += time*nBG[bin];
        double this_flux =
          deltaESensitivity * get_table(&x[0], &yI[0], points, eC2)/(eC2*eC2);
        fluxInt += this_flux;
        modulationFactor += this_flux * mu[bin];

        //log_debug("Energy %.1fkeV: Weight = %g, Crab = %g, nS = %g", eC2, weight/time, crab, nS[bin]);
      }
      
      modulationFactor /= fluxInt;
      nSignalTotal += nSignal;
      log_debug("eC: %.2f -- nSignal: %.2f", eC, nSignal);
      
      fluxInt *= eC*eC/(deltaESensitivity*intP);

      ycPi[kd]=tPiV;
      ycChi[kd]=tChiV;
      
      // let's get "measured values" and error bars on these quantities
      double FluxM,PiM,ChiM,FluxE,PiELow,PiEHigh,ChiE;
      
      tChiV *= PI/180.;
      getValuesAndErrors(nSignal, nBackground, modulationFactor, tPiV, tChiV,
                         FluxM, PiM, ChiM, FluxE, PiELow, PiEHigh, ChiE);
      double mdp = mdpFactor * sqrt(nSignal + nBackground) / (modulationFactor * nSignal);

      double eFlux = sqrt(nSignal + 2*nBackground)/nSignal;
      ycFlux[kd]=getRNG().Gaus(fluxInt, fluxInt * eFlux);
      yeFlux[kd]=ycFlux[kd] * eFlux;
      ycPi[kd]=PiM;
      ycChi[kd]=ChiM*180./PI;
      ycMdp[kd] = mdp;
      yePiLow[kd]=PiELow;
      yePiHigh[kd]=PiEHigh;
      yeChi[kd]=ChiE*180./PI;
      if (PiM > 2*PiELow) kd++;
    }
    log_info("Integral source rate: %.3fHz", nSignalTotal/time);
  }    

  TGraphAsymmErrors *bFlux;
  TGraphAsymmErrors *bPi;
  TGraphAsymmErrors *bChi;
  TGraphAsymmErrors *bMdp = NULL;
  
  if (fluxunit == CGSe9) {
    for (int i = 0; i < kd; ++i) {
      ycFlux[i] *= 1e9;
      yeFlux[i] *= 1e9;
    }
  }
  bFlux = new TGraphAsymmErrors(kd,xc,ycFlux,xem,xep,yeFlux,yeFlux);
  bFlux->SetLineColor(datacolor);
  bFlux->SetMarkerColor(datacolor);
  if (large)
    bFlux->SetLineWidth(2);
  else
    bFlux->SetLineWidth(datawidth);
  bFlux->SetMarkerStyle(24);
  if (data_legend >= 0)
    legends[data_legend].legend->AddEntry(bFlux, datalabel.c_str(), "LP");

  if (fractionunit == Percent) {
    for (int i = 0; i < kd; ++i) {
      ycPi[i] *= 100;
      yePiLow[i] *= 100;
      yePiHigh[i] *= 100;
    }
  }
  bPi = new TGraphAsymmErrors(kd,xc,ycPi,xem,xep,yePiLow,yePiHigh);
  bPi->SetLineColor(datacolor);
  bPi->SetMarkerColor(datacolor);
  if (large)
    bPi->SetLineWidth(2);
  else
    bPi->SetLineWidth(datawidth);
  bPi->SetMarkerStyle(24);    
    
  bChi = new TGraphAsymmErrors(kd,xc,ycChi,xem,xep,yeChi,yeChi);
  bChi->SetLineColor(datacolor);
  bChi->SetMarkerColor(datacolor);
  if (large)
    bChi->SetLineWidth(2);
  else
    bChi->SetLineWidth(datawidth);
  bChi->SetMarkerStyle(24);

  if (showmdp) {
    bMdp = new TGraphAsymmErrors(kd, xc, ycMdp, xem, xep);
    bMdp->SetLineWidth(mdpWidth);
    bMdp->SetLineStyle(mdpStyle);
    bMdp->SetLineColor(mdpColor);
    bMdp->SetMarkerStyle(kDot);
    bMdp->SetMarkerColor(mdpColor);
    if ((mdp_legend >= 0) && !mdpLabel.empty())
      legends[mdp_legend].legend->AddEntry(bMdp, mdpLabel.c_str(), "L");

    const double mdpPrintScale = (fractionunit == Frac ? 100 : 1);

    int coutPrecision = cout.precision();
    ios::fmtflags coutFlags = cout.flags();
    cout.setf(ios::fixed, ios::floatfield);
    std::cout << "\n";
    log_info("MDP data:");
    if (have_phase) {
      std::cout << "      Phase       | MDP [%] \n"
                << "----------------------------\n";
      for (int i = 0; i < kd; ++i) {
        std::cout << std::setprecision(3) << std::setw(6) << xc[i]-xem[i]
                  << " ... " << std::setw(6) << xc[i]+xep[i] << " | "
                  << std::setprecision(2) << std::setw(5)
                  << mdpPrintScale * ycMdp[i] << "\n";
      }
    } else {
      std::cout << "  Energy  [MeV]   | MDP [%] \n"
                << "----------------------------\n";
      for (int i = 0; i < kd; ++i) {
        std::cout << std::setprecision(3) << std::setw(6) << (xc[i]-xem[i])*1e-3
                  << " ... " << std::setw(6) << (xc[i]+xep[i])*1e-3 << " | "
                  << std::setprecision(2) << std::setw(5)
                  << mdpPrintScale * ycMdp[i] << "\n";
      }
    }
    std::cout << std::endl;
    cout.precision(coutPrecision);
    cout.flags(coutFlags);
  }
  
  if (fluxunit == CGSe9) {
    ymin *= 1e9;
    ymax *= 1e9;
  }

  if (Hymin>0.) ymin=Hymin;
  if (Hymax>0.) ymax=Hymax;

  if (!math::isnan(userXmin)) xmin = userXmin;
  if (!math::isnan(userXmax)) xmax = userXmax;

  c00[0]->Divide(1, panels.size());

  TH1F *hFrame = NULL;
  vector<Panel>::const_iterator panel_iter =
    std::find(panels.begin(), panels.end(), Flux);

  if (panel_iter != panels.end()) {
    log_info("Now drawing intensities");
    int panel = panel_iter - panels.begin() + 1;

    c00[0]->cd(panel);
    gPad->SetPad(0, panelConfig[panel-1].bottom, 1, panelConfig[panel-1].top);
    c00[0]->Update();
    c00[0]->cd(panel);
    double leftMargin = double(d.leftMargin)/d.width;
    double rightMargin = double(d.rightMargin)/d.width;
    gPad->SetMargin(leftMargin, rightMargin,
                    panelConfig[panel-1].bottomMargin,
                    panelConfig[panel-1].topMargin);
    gPad->SetBorderMode(0);
    if (!have_phase) {
      gPad->SetLogx();
      gPad->SetLogy();
    }
    gPad->SetTickx(1);

    if (fluxunit == CGSe9) fluxLabel = "E f_{E} [10^{-9} erg cm^{-2} s^{-1}]";
    if (!userFluxLabel.empty()) fluxLabel = userFluxLabel;
    hFrame = new TH1F("hFrame1", (";Energy [keV];" + fluxLabel).c_str(), 100, xmin, xmax);
    if (!userXLabel.empty()) {
      hFrame->SetXTitle(userXLabel.c_str());
    } else if (have_phase) {
      hFrame->SetXTitle("Phase");
    }
    hFrame->SetMinimum(ymin);
    hFrame->SetMaximum(ymax);
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetXaxis()->SetMoreLogLabels();
    hFrame->GetXaxis()->SetNoExponent();
    hFrame->SetStats(false);
    hFrame->GetXaxis()->SetTickLength(0.05);
    hFrame->GetXaxis()->SetNdivisions(xNdiv);
    hFrame->GetXaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetXaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetTickLength(0.02);
    hFrame->GetYaxis()->SetNdivisions(fluxNdiv);
    hFrame->GetYaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    if (have_phase) {
      hFrame->GetYaxis()->SetLabelOffset(0.011);
    } else {
      hFrame->GetYaxis()->SetLabelOffset(0.);
    }
    if (large) {
      hFrame->GetXaxis()->SetTitleOffset(1.0);
    } else {
      hFrame->GetXaxis()->SetTitleOffset(1.1);
    }
    if (fluxMoreLabels) {
      hFrame->GetYaxis()->SetMoreLogLabels();
      hFrame->GetYaxis()->SetNoExponent();
    }

    double titleoffset = calculateTitleOffset(d.yoffset, xmin, xmax,
                                              panelConfig[panel-1].fontsize,
                                              !have_phase, d);
    hFrame->GetYaxis()->SetTitleOffset(titleoffset);

    hFrame->DrawCopy("AXIS");
    
    mgrFlux->Draw("L");
    bFlux->Draw("P");
    
    c00[0]->cd();
    c00[0]->Update();        
  }

  panel_iter = std::find(panels.begin(), panels.end(), Fraction);
  if (panel_iter != panels.end()) {
    log_info("Now drawing polarization fraction");
    int panel = panel_iter - panels.begin() + 1;

    c00[0]->cd(panel);
    gPad->SetPad(0, panelConfig[panel-1].bottom, 1, panelConfig[panel-1].top);
    c00[0]->Update();
    c00[0]->cd(panel);
    gPad->SetMargin(double(d.leftMargin)/d.width, double(d.rightMargin)/d.width,
                    panelConfig[panel-1].bottomMargin,
                    panelConfig[panel-1].topMargin);
    gPad->SetBorderMode(0);
    if (!have_phase) gPad->SetLogx();
    gPad->SetLogy(false);
    gPad->SetTickx(1);
    
    if (hFrame) hFrame->Delete();
    
    double pmin = 0.001;
    Pimax *= 1.4;
    
    if (fractionunit == Percent) {
      pmin *= 100;
      Pimax *= 100;
    }    
    
    if (Hpmin>0.) pmin=Hpmin;
    if (Hpmax>0.) Pimax=Hpmax;
    
    if (fractionunit == Percent) polFracLabel = "Pol. Fraction [%]";
    if (!userFracLabel.empty()) polFracLabel = userFracLabel;
    hFrame = new TH1F("hFrame2", (";Energy [keV];" + polFracLabel).c_str(), 100, xmin, xmax);
    if (!userXLabel.empty()) {
      hFrame->SetXTitle(userXLabel.c_str());
    } else if (have_phase) {
      hFrame->SetXTitle("Phase");
    }
    hFrame->SetMinimum(pmin);
    hFrame->SetMaximum(Pimax);
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetXaxis()->SetMoreLogLabels();
    hFrame->GetXaxis()->SetNoExponent();
    hFrame->SetStats(false);
    hFrame->GetXaxis()->SetTickLength(0.05);
    hFrame->GetXaxis()->SetNdivisions(xNdiv);
    hFrame->GetXaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetXaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetTickLength(0.02);
    hFrame->GetYaxis()->SetNdivisions(piNdiv);
    hFrame->GetYaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    if (large) {
      hFrame->GetXaxis()->SetTitleOffset(1.0);
    } else {
      hFrame->GetXaxis()->SetTitleOffset(1.1);
    }

    double titleoffset = calculateTitleOffset(d.yoffset, xmin, xmax,
                                              panelConfig[panel-1].fontsize,
                                              !have_phase, d);
    hFrame->GetYaxis()->SetTitleOffset(titleoffset);
    
    hFrame->DrawCopy("AXIS");

    mgrPi->Draw("LP");
    bPi->Draw("P");
    if (bMdp) bMdp->Draw("PZ");
    
    double avgPolError = 0;
    for (int i = 0; i < bPi->GetN(); ++i) {
      avgPolError += bPi->GetEYlow()[i] + bPi->GetEYhigh()[i];
    }
    avgPolError /= 2*bPi->GetN();
    log_info("Average error on the polarization fraction: %.2f%%",
             (fractionunit == Frac ? 100 : 1)*avgPolError);

    c00[0]->cd();
    c00[0]->Update();
  }

  panel_iter = std::find(panels.begin(), panels.end(), Angle);
  if (panel_iter != panels.end()) {
    log_info("Now drawing polarization angle");
    int panel = panel_iter - panels.begin() + 1;

    c00[0]->cd(panel);
    gPad->SetPad(0, panelConfig[panel-1].bottom, 1, panelConfig[panel-1].top);
    c00[0]->Update();
    c00[0]->cd(panel);
    gPad->SetMargin(double(d.leftMargin)/d.width, double(d.rightMargin)/d.width,
                    panelConfig[panel-1].bottomMargin,
                    panelConfig[panel-1].topMargin);
    gPad->SetBorderMode(0);
    if (!have_phase) gPad->SetLogx();
    gPad->SetLogy(false);
    gPad->SetTickx(1);

    if (hFrame) hFrame->Delete();
    
    chimin-=6.;
    chimax+=6.;
    if (chimin<-90) chimin=-90.;
    if (chimax> 90) chimax= 90.;
    
    if (Hchimin != 0.) chimin=Hchimin;
    if (Hchimax  > 0.) chimax=Hchimax;

    if (!userDirLabel.empty()) polDirLabel = userDirLabel;
    hFrame = new TH1F("hFrame3", (";Energy [keV];" + polDirLabel).c_str(), 100, xmin, xmax);
    if (!userXLabel.empty()) {
      hFrame->SetXTitle(userXLabel.c_str());
    } else if (have_phase) {
      hFrame->SetXTitle("Phase");
    }
    hFrame->SetMinimum(chimin);
    hFrame->SetMaximum(chimax);
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetXaxis()->SetMoreLogLabels();
    hFrame->GetXaxis()->SetNoExponent();
    hFrame->SetStats(false);
    hFrame->GetXaxis()->SetTickLength(0.05);
    hFrame->GetXaxis()->SetNdivisions(xNdiv);
    hFrame->GetXaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetXaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetTickLength(0.02);
    hFrame->GetYaxis()->SetNdivisions(chiNdiv);
    hFrame->GetYaxis()->SetTitleSize(panelConfig[panel-1].fontsize);
    hFrame->GetYaxis()->SetLabelSize(panelConfig[panel-1].fontsize);
    if (large) {
      hFrame->GetXaxis()->SetTitleOffset(1.0);
    } else {
      hFrame->GetXaxis()->SetTitleOffset(1.1);
    }

    double titleoffset = calculateTitleOffset(d.yoffset, xmin, xmax,
                                              panelConfig[panel-1].fontsize,
                                              !have_phase, d);
    hFrame->GetYaxis()->SetTitleOffset(titleoffset);
    
    hFrame->DrawCopy("AXIS");

    mgrChi->Draw("LP");
    bChi->Draw("P");

    if (!plot_clean) {
      time_t t = ::time(NULL);
      char date[100];
      size_t len = strftime(date, 100, "%m/%d/%Y", localtime(&t));
      if (len == 0) {
        log_error("Cannot print date because it is too long ?!?!?!?");
      } else {
        TPaveText *label = new TPaveText(0.7, 0, 0.9995, 0.04, "BR NDC");
        label->SetFillStyle(0);
        label->SetBorderSize(0);
        label->SetTextFont(42);
        label->SetTextSize(0.04);
        label->SetTextAlign(31);
        
        string str = string("Created ") + date + " with GPST " GPST_VERSION_STRING;
        label->AddText(1, 0, str.c_str());
        label->Draw();
      }
    }

    c00[0]->cd();
    c00[0]->Update();
  }

  if (!legends.empty()) {
    log_info("Now drawing legends");
    for (vector<LegendConfig>::iterator it = legends.begin();
         it != legends.end(); ++it)
    {
      if (it->panel_index > 0) {
        c00[0]->cd(it->panel_index);
        it->legend->Draw();
      }
    }
  }

  c00[0]->cd();
  c00[0]->Update();

  result = Result(name, c00[0], bFlux, bPi, bChi);
  log_info("\x1b[1m=== Done ===\x1b[0m");
  log_info("You can now access and save the result through the variable 'result'.");
  result.info();

  delete[] nS;
  delete[] nBG;
  delete[] mu;
}


double Crab(double energy)
{
  return 10.17*pow(energy,-2.15)*energy*energy*1.6E-9; // ergs cm-2 s-1
}

double getChi(double q, double u)
{
  double res = atan(-u/q)/2.;
  if(q < 0.) res += PI/2.;
  else if (u > 0.) res += PI;
  res*=-1;
  res *= 180./PI;
  if (res<(-90.)) res+=180.;
  return res;
}

void getValuesAndErrors(double Nsrc, double Nbg, double mu, double aT,
                        double phiT, double &fluxM, double &aM, double &phiM,
                        double &fluxE, double &aELow, double &aEHigh,
                        double &phiE)
{
  TRandom3 &rng = getRNG();
  fluxM=rng.Gaus(Nsrc,sqrt(Nsrc));

  //	double ratio = N n*n / 4. / (sigma*sigma);
  double ratio = (fluxM*fluxM)/(fluxM+Nbg)/4.; 
	
  // get sigma for at deltaPhi=0.
  int steps=25000;
    
  double maxProb = 0.;
  double aMax = 0.;
  double sum1=0.;
  double delta1 = 1./(double)steps;

  for (double aLoop=delta1/2.; aLoop<10.; aLoop+=delta1) {
    double p = prob(ratio, aT, mu, aLoop, 0.);
    sum1 += p;
    if (p > maxProb) {
      maxProb = p;
      aMax = aLoop;
    }
  }
  log_debug("aMax = %g, maxProb = %g", aMax, maxProb);
    
  double sum1b=0;
  double a1=aMax;
  double a2=a1;

  while (sum1b<0.67*sum1)
  {
    if ((prob(ratio,aT,mu,a1-delta1,0.) > prob(ratio,aT,mu,a2+delta1,0.)) && (a1>delta1))
    {
      a1-=delta1;
      sum1b += prob(ratio,aT,mu,a1,0.);
    }
    else
    {
      a2+=delta1;
      sum1b += prob(ratio,aT,mu,a2,0.);
    }
  }
  double sigma1 = (a2-a1)/2.;

  aM=rng.Gaus(aMax, sigma1);
  aELow = aMax-a1;
  aEHigh = a2-aMax;
    
  double sum2=0;
  double delta2 = PI/(double)steps;

  for (double deltaPhi=-M_PI_2+delta2/2.; deltaPhi<M_PI_2; deltaPhi+=delta2)
    sum2 += prob(ratio,aT,mu,aMax,deltaPhi);

  double sum2b=0;
  double p1=phiT;
  double p2=phiT;
  
  while (sum2b<0.67*sum2)
  {
    if ((prob(ratio,aM,mu,aM,p1-delta2-phiT) > prob(ratio,aM,mu,aM,p2+delta2-phiT)) && (p1>-M_PI_2))
    {
      p1-=delta2;
      sum2b += prob(ratio,aM,mu,aM,p1-phiT);
    }
    else
    {
      p2+=delta2;
      sum2b += prob(ratio,aM,mu,aM,p2-phiT);
    }
  };
  double sigma2 = (p2-p1)/2.;

  phiM=rng.Gaus(phiT,sigma2);
  fluxE=sqrt(fluxM);
  phiE=sigma2;
}

double prob(double ratio,double aT,double mu,double aM,double deltaPhi)
{
  return ratio*aM*mu*mu / PI * exp(-ratio*(aT*aT*mu*mu+aM*aM*mu*mu-2.*aT*mu*aM*mu*cos(2*deltaPhi)));
}

void loadSensitivity(double *ns, double *nbg, double *mu, int numbins,
                     double modulation_factor, double start,
                     double deltaE, std::string sensitivityFile,
                     std::string backgroundFile, std::string modulationFile,
                     double scale)
{
  FILE *fp1,*fp2,*fp3=NULL;
  double n1,n2;
  double n3 = modulation_factor;
  double end=start;

  string prefix = "";
  const char *directory = getenv("GPST_DATA_DIR");
  if (directory && strlen(directory) > 0) {
    prefix = string(directory) + "/";
  }

  if (prefix.length() > 0) {
    log_info("Loading sensitivity files from directory %s", prefix.c_str());
  } else {
    log_info("Loading sensitivity files from current directory");
  }

  if (sensitivityFile.empty()) {
    log_fatal("No sensitivity file specified (Keyword 'sensitivityfile' "
              "missing).");
  }
  log_info("Using the sensitivity file '%s'", sensitivityFile.c_str());
  sensitivityFile = prefix + sensitivityFile;

  if (backgroundFile.empty()) {
    log_fatal("No background file specified (Keyword 'backgroundfile' "
              "missing).");
  }
  log_info("Using the background rate file '%s'", backgroundFile.c_str());
  backgroundFile = prefix + backgroundFile;

  if (!(fp1 = fopen(sensitivityFile.c_str(), "r"))) {
    log_fatal("Could not open the file %s.", sensitivityFile.c_str());
  }

  if (!(fp2 = fopen(backgroundFile.c_str(), "r"))) {
    log_fatal("Could not open the file %s.", backgroundFile.c_str());
  }

  if (!modulationFile.empty()) {
    log_info("Using modulation factors from file '%s'.", modulationFile.c_str());
    modulationFile = prefix + modulationFile;
    if (!(fp3 = fopen(modulationFile.c_str(), "r"))) {
      log_fatal("Could not open the file %s.", modulationFile.c_str());
    }
  } else {
    log_info("Using modulation factor %.3f", modulation_factor);
  }

  if (scale != 1.) {
    log_info("Scaling signal rate by %.2f", scale);
  }
  
  int count=0;
    
  while(count<numbins)// (flag)
  {
    if (fscanf(fp1,"%lf",&n1) != 1) {
      log_fatal("Failed to read from file %s.", sensitivityFile.c_str());
    }
    if (fscanf(fp2,"%lf",&n2) != 1) {
      log_fatal("Failed to read from file %s.", backgroundFile.c_str());
    }
    if (fp3) {
      if (fscanf(fp3,"%lf",&n3) != 1) {
        log_fatal("Failed to read from file %s.", modulationFile.c_str());
      }
    }
    // printf("%d %f %f\n",count,n1,n2);
    ns[count]=n1*scale;
    nbg[count]=n2;
    mu[count]=n3;
    ++count;
    end+=deltaE;
  }
  fclose(fp1);
  fclose(fp2);
  if (fp3) fclose(fp3);
  log_info("Read %d data records for energies from %lg to %lgkeV",
           count, start, end);
}

inline double linear_interpolate(double x1, double y1, double x2, double y2,
                                 double x)
{
  return y1 + (y2-y1)/(x2-x1) * (x-x1);
}

void interpolate_stokes(double *xs, double *Is, double *Qs, double *Us,
                        int numbins, double x, double &I, double &Q, double &U)
{
  if ((x<xs[0])||(x>xs[numbins-1]))
  {
    I = Q = U = 0;
    return;
  }

  int i;
  for (i = 0; i < numbins-1; ++i)
    if (x < xs[i]) break;

  double pi_left = sqrt(math::sqr(Qs[i-1]) + math::sqr(Us[i-1]))/Is[i-1];
  double chi_left = 0.5*atan2(Us[i-1], Qs[i-1]);
  double pi_right = sqrt(math::sqr(Qs[i]) + math::sqr(Us[i]))/Is[i];
  double chi_right = 0.5*atan2(Us[i], Qs[i]);

  double pi = linear_interpolate(xs[i-1], pi_left, xs[i], pi_right, x);
  double chi = linear_interpolate(xs[i-1], chi_left, xs[i], chi_right, x);

  I = get_table(xs, Is, numbins, x);
  Q = I * pi * cos(2*chi);
  U = I * pi * sin(2*chi);
}

double get_table(double *t_x, double *t_v, int numbins, double x)
{
  int i;
  double rest,res;
  double delta;
	
  if ((x<t_x[0])||(x>t_x[numbins-1]))
  {
    //log_debug("get_table: x = %g, min = %g, max = %g, points = %d", x, t_x[0], t_x[numbins-1], numbins);
    return 0.;
  }
	
  for (i=0; i<numbins-1; i++) {
    if (x == t_x[i]) return t_v[i];
    if (x<t_x[i]) break;
  }
	
  delta = log(t_x[i])-log(t_x[i-1]);
  rest  = (log(x)-log(t_x[i-1]))/delta;
  if ((t_v[i-1]>0)&&(t_v[i]>0))
    res   = exp((1.-rest)*log(t_v[i-1]) + rest*log(t_v[i]));
  else if ((t_v[i-1]<0)&&(t_v[i]<0))
    res   = -exp((1.-rest)*log(-t_v[i-1]) + rest*log(-t_v[i]));
  else
    res   = (1.-rest)*t_v[i-1] + rest*t_v[i];

  return res;
}

void style()
{
  TStyle *mystyle = new TStyle("mystyle","my own style");
    
  mystyle->SetOptDate(0);
  mystyle->SetCanvasColor(0);
  mystyle->SetCanvasBorderMode(0);
  mystyle->SetCanvasBorderSize(1);
  mystyle->SetCanvasDefH(500);
  mystyle->SetCanvasDefW(500);
  mystyle->SetFrameBorderMode(0);
  mystyle->SetFrameBorderSize(1);
    
  mystyle->SetPadColor(0);
  mystyle->SetPadBorderMode(1);
  mystyle->SetPadBorderSize(1);
	
  mystyle->SetHistLineWidth(2);
  mystyle->SetHistLineColor(1);
  mystyle->SetLabelSize(.08, "x");
  mystyle->SetLabelOffset(0.01, "x");
  mystyle->SetLabelSize(.08, "y");
  mystyle->SetLabelOffset(0.01, "y");
    
  mystyle->SetTitleColor(1);
  mystyle->SetStatColor(0);

  mystyle->SetPadTopMargin(0.04);
  mystyle->SetPadBottomMargin(0.2);
  mystyle->SetPadLeftMargin(0.1);
    
  mystyle->SetTitleYOffset(0.80);
  mystyle->SetTitleXOffset(1.1);
  mystyle->SetTitleSize(0.08, "X");
  mystyle->SetTitleSize(0.08, "Y");
  mystyle->cd();

  gROOT->ForceStyle();
}


bool read_version_info(const std::string &line, int &version, int linecount)
{
  if (line.substr(0, 4) != "GPST") {
    log_error("On line %d of input file: First line must be version info line "
              "(GPST %d).", linecount, GPST_VERSION);
    return false;
  }

  size_t strpos = 4;
  for ( ; strpos != line.length() && isspace(line[strpos]); ++strpos) ;

  if (strpos == line.length()) {
    log_error("On line %d of input file: First line must be version info line "
              "(GPST %d).", linecount, GPST_VERSION);
    return false;
  }

  string version_string;
  for ( ; strpos != line.length() && isdigit(line[strpos]); ++strpos)
    version_string.push_back(line[strpos]);

  if (version_string.empty()) {
    log_error("On line %d of input file: First line must be version info line "
              "(GPST %d).", linecount, GPST_VERSION);
    return false;
  }

  for ( ; strpos != line.length(); ++strpos) {
    if (!isspace(line[strpos])) {
      log_error("On line %d of input file: First line must be version info line "
                "(GPST %d).", linecount, GPST_VERSION);
      return false;
    }
  }

  // by construction, version_string is an integer
  version = atoi(version_string.c_str());
  if (version > GPST_VERSION) {
    log_error("File version is %d, but this version of the Gamma-ray "
              "Polarimetry Simulation Toolkit can only read version %d and "
              "older.", version, GPST_VERSION);
    return false;
  }

  return true;
}


void linearIntegralQUI(double xMin, double xMax, int nPoints, double *x,
                       double *flux, double *pi, double *chi,
                       double &qsum, double &usum, double &isum)
{
  qsum = usum = isum = 0;

  int p = 0;
  while (p < nPoints && x[p] <= xMin) ++p;
  if (p >= nPoints) return;
  if (p > 0) --p;

  while (x[p] < xMax) {
    if (p >= nPoints-1) break;
    double x1 = x[p];
    double flux1 = flux[p];
    double q1 = flux[p] * pi[p] * cos(2 * chi[p]);
    double u1 = flux[p] * pi[p] * sin(2 * chi[p]);

    double x2 = x[p+1];
    double flux2 = flux[p+1];
    double q2 = flux[p+1] * pi[p+1] * cos(2 * chi[p+1]);
    double u2 = flux[p+1] * pi[p+1] * sin(2 * chi[p+1]);

    if (x1 < xMin) {
      flux1 += (xMin - x1) * (flux2-flux1)/(x2-x1);
      q1 += (xMin - x1) * (q2-q1)/(x2-x1);
      u1 += (xMin - x1) * (u2-u1)/(x2-x1);
      x1 = xMin;
    }

    if (x2 > xMax) {
      flux2 += (xMax - x2) * (flux2-flux1)/(x2-x1);
      q2 += (xMax - x1) * (q2-q1)/(x2-x1);
      u2 += (xMax - x1) * (u2-u1)/(x2-x1);
      x2 = xMax;
    }

    isum += 0.5 * (flux1 + flux2) * (x2 - x1);
    qsum += 0.5 * (q1 + q2) * (x2 - x1);
    usum += 0.5 * (u1 + u2) * (x2 - x1);

     ++p;
  }
}


namespace strings {

  void strip(string &s)
  {
    // do nothing on empty strings
    if (s.empty()) return;

    string::iterator p;

    // remove from beginning --> find the first character that is not whitespace
    for (p = s.begin(); (p != s.end()) && isspace(*p); ++p);

    s.erase(s.begin(), p);

    // remove from end --> find the first of the trailing c characters
    for (p = s.end(); (p != s.begin()) && isspace(*(--p)););

    if (!isspace(*p)) ++p;

    s.erase(p, s.end());
  }

  void strip(string &s, char c)
  {
    // do nothing on empty strings
    if (s.empty()) return;

    string::iterator p;

    // remove from beginning --> find the first character that is not whitespace
    for (p = s.begin(); (p != s.end()) && *p == c; ++p);

    s.erase(s.begin(), p);

    // remove from end --> find the first of the trailing c characters
    for (p = s.end(); (p != s.begin()) && *(--p) == c;);

    if (*p != c) ++p;

    s.erase(p, s.end());
  }

  void split(const string &str, char delim, vector<string> &v)
  {
    if (str.empty()) return;

    string::size_type i = 0;
    string::size_type j = str.find(delim);

    if (j == string::npos) {
      v.push_back(str);
      return;
    }
    
    while (j != string::npos) {
      if (j != i)
        v.push_back(str.substr(i, j-i));
      i = ++j;
      j = str.find(delim, j);
      
      if (j == string::npos && i < str.length()) {
	v.push_back(str.substr(i));
      }
    }
  }

  bool split_assignment(const string &str, string &key, string &value)
  {
    size_t equal = str.find('=');
    if (equal == string::npos) return false;
    key = str.substr(0, equal);
    value = str.substr(equal+1);
    return true;
  }

}


std::string prepare_special_value(const std::string &value)
{
  size_t semicolon = value.find(';');
  std::string res = value.substr(0, semicolon);
  strings::strip(res);
  return res;
}


typedef std::map<std::string, Style_t> StyleMap_t;

StyleMap_t make_style_map()
{
  StyleMap_t style_map;
  GPST_ADD_SPECIAL_VALUE(style_map, "solid", 1);
  GPST_ADD_SPECIAL_VALUE(style_map, "short_dashed", 2);
  GPST_ADD_SPECIAL_VALUE(style_map, "dotted", 3);
  GPST_ADD_SPECIAL_VALUE(style_map, "short_dash_dotted", 4);
  GPST_ADD_SPECIAL_VALUE(style_map, "dash_dotted", 5);
  GPST_ADD_SPECIAL_VALUE(style_map, "dash_3dot", 6);
  GPST_ADD_SPECIAL_VALUE(style_map, "dashed", 7);
  GPST_ADD_SPECIAL_VALUE(style_map, "dash_dot_dot", 8);
  GPST_ADD_SPECIAL_VALUE(style_map, "long_dashed", 9);
  GPST_ADD_SPECIAL_VALUE(style_map, "long_dash_dotted", 10);
  return style_map;
}

Style_t interpretStyle(const std::string &value)
{
  try {
    return strings::to<Style_t>(value);
  } catch (...) {
    static const StyleMap_t style_map = make_style_map();
    StyleMap_t::const_iterator iter = style_map.find(prepare_special_value(value));
    if (iter == style_map.end()) throw;
    return iter->second;
  }
}


typedef std::map<std::string, Color_t> ColorMap_t;

ColorMap_t make_color_map()
{
  ColorMap_t color_map;
  GPST_ADD_SPECIAL_VALUE(color_map, "black", kBlack);
  GPST_ADD_SPECIAL_VALUE(color_map, "kBlack", kBlack);
  GPST_ADD_SPECIAL_VALUE(color_map, "gray", kGray);
  GPST_ADD_SPECIAL_VALUE(color_map, "kGray", kGray);
  GPST_ADD_SPECIAL_VALUE(color_map, "white", kWhite);
  GPST_ADD_SPECIAL_VALUE(color_map, "kWhite", kWhite);
  GPST_ADD_SPECIAL_VALUE(color_map, "blue", kBlue);
  GPST_ADD_SPECIAL_VALUE(color_map, "kBlue", kBlue);
  GPST_ADD_SPECIAL_VALUE(color_map, "violet", kViolet);
  GPST_ADD_SPECIAL_VALUE(color_map, "kViolet", kViolet);
  GPST_ADD_SPECIAL_VALUE(color_map, "magenta", kMagenta);
  GPST_ADD_SPECIAL_VALUE(color_map, "kMagenta", kMagenta);
  GPST_ADD_SPECIAL_VALUE(color_map, "pink", kPink);
  GPST_ADD_SPECIAL_VALUE(color_map, "kPink", kPink);
  GPST_ADD_SPECIAL_VALUE(color_map, "red", kRed);
  GPST_ADD_SPECIAL_VALUE(color_map, "kRed", kRed);
  GPST_ADD_SPECIAL_VALUE(color_map, "orange", kOrange);
  GPST_ADD_SPECIAL_VALUE(color_map, "kOrange", kOrange);
  GPST_ADD_SPECIAL_VALUE(color_map, "yellow", kYellow);
  GPST_ADD_SPECIAL_VALUE(color_map, "kYellow", kYellow);
  GPST_ADD_SPECIAL_VALUE(color_map, "spring", kSpring);
  GPST_ADD_SPECIAL_VALUE(color_map, "kSpring", kSpring);
  GPST_ADD_SPECIAL_VALUE(color_map, "green", kGreen);
  GPST_ADD_SPECIAL_VALUE(color_map, "kGreen", kGreen);
  GPST_ADD_SPECIAL_VALUE(color_map, "teal", kTeal);
  GPST_ADD_SPECIAL_VALUE(color_map, "kTeal", kTeal);
  GPST_ADD_SPECIAL_VALUE(color_map, "cyan", kCyan);
  GPST_ADD_SPECIAL_VALUE(color_map, "kCyan", kCyan);
  GPST_ADD_SPECIAL_VALUE(color_map, "azure", kAzure);
  GPST_ADD_SPECIAL_VALUE(color_map, "kAzure", kAzure);
  return color_map;
}

Color_t interpretColor(const std::string &value)
{
  try {
    return strings::to<Color_t>(value);
  } catch (...) {
    static const ColorMap_t color_map = make_color_map();
    size_t signpos = value.find_first_of("+-");
    string color_string = value.substr(0, signpos);
    strings::strip(color_string);
    ColorMap_t::const_iterator iter = color_map.find(prepare_special_value(color_string));
    if (iter == color_map.end()) throw;
    Color_t color = iter->second;
    if (signpos != string::npos) {
      int sign = value[signpos] == '+' ? 1 : -1;
      string shift_string = value.substr(signpos+1);
      int shift = strings::to<int>(shift_string);
      color += sign*shift;
    }
    return color;
  }
}

bool interpretBoolean(const std::string &value)
{
  static bool boolSetup = false;
  if (!boolSetup) {
    SPECIAL_VALUES.insert("true");
    SPECIAL_VALUES.insert("false");
  }

  std::string sv = prepare_special_value(value);
  if (sv == "true") return true;
  if (sv == "false") return false;

  return strings::to<int>(value);
}


typedef std::map<std::string, FluxUnit> FluxUnitMap;

FluxUnitMap makeFluxUnitMap()
{
  FluxUnitMap fluxunit_map;
  GPST_ADD_SPECIAL_VALUE(fluxunit_map, "cgs", CGS);
  GPST_ADD_SPECIAL_VALUE(fluxunit_map, "cgse9", CGSe9);
  return fluxunit_map;
}

FluxUnit interpretFluxUnit(const std::string &value)
{
  static const FluxUnitMap fluxunit_map = makeFluxUnitMap();
  FluxUnitMap::const_iterator iter = fluxunit_map.find(prepare_special_value(value));
  if (iter == fluxunit_map.end())
    throw std::runtime_error("Invalid flux unit: '" + value + "'");
  return iter->second;
}


typedef std::map<std::string, FractionUnit> FractionUnitMap;

FractionUnitMap makeFractionUnitMap()
{
  FractionUnitMap fractionunit_map;
  GPST_ADD_SPECIAL_VALUE(fractionunit_map, "fraction", Frac);
  GPST_ADD_SPECIAL_VALUE(fractionunit_map, "percent", Percent);
  return fractionunit_map;
}

FractionUnit interpretFractionUnit(const std::string &value)
{
  static const FractionUnitMap fractionunit_map = makeFractionUnitMap();
  FractionUnitMap::const_iterator iter = fractionunit_map.find(prepare_special_value(value));
  if (iter == fractionunit_map.end())
    throw std::runtime_error("Invalid unit for polarization fraction: '" + value + "'");
  return iter->second;
}


typedef std::map<std::string, Panel> PanelNameMap;

PanelNameMap makePanelNameMap()
{
  PanelNameMap panelname_map;
  GPST_ADD_SPECIAL_VALUE(panelname_map, "flux", Flux);
  GPST_ADD_SPECIAL_VALUE(panelname_map, "frac", Fraction);
  GPST_ADD_SPECIAL_VALUE(panelname_map, "angle", Angle);
  return panelname_map;
}

Panel interpretPanel(const std::string &value)
{
  static const PanelNameMap panelname_map = makePanelNameMap();
  PanelNameMap::const_iterator iter = panelname_map.find(prepare_special_value(value));
  if (iter == panelname_map.end()) {
    throw std::runtime_error("Invalid panel specification: '" + value + "'");
  }
  return iter->second;
}


string make_filename(const std::string &base, std::string input)
{
  strings::strip(input, '"');
  
  log_debug("make_filename input: %s", input.c_str());

  if (input.empty()) return input;
  if (input[0] == '/') return input;

  if (input[0] == '$') {
    size_t end_of_var =
      input.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                              "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "1234567890_", 1);
    string varname = input.substr(1, end_of_var-1);

    log_debug("make_filename varname: %s", varname.c_str());

    const char *envvar = getenv(varname.c_str());

    if (!envvar) return input;

    log_debug("make_filename envvar: %s", envvar);

    return envvar + input.substr(end_of_var);
  }

  size_t slash = base.rfind('/');

  if (slash == 0 || slash == string::npos) return input;
  
  string base_dir = base.substr(0, slash);

  return base_dir + "/" + input;
}


PanelConfig makePanelConfig(unsigned index, unsigned nPads, double &top,
                            const Dimensions &d, bool large)
{
  PanelConfig pc;

  pc.height = d.frameHeight;
  if (index == nPads) pc.height += d.bottomMargin;
  if (index == 1) pc.height += d.topMargin;

  pc.fontsize = double(d.textHeight)/pc.height;

  pc.top = top;

  if (index < nPads)
    pc.bottom = top - double(pc.height)/d.height;
  else
    pc.bottom = 0;

  top = pc.bottom;

  if (index == nPads)
    pc.bottomMargin = double(d.bottomMargin)/pc.height;
  else
    pc.bottomMargin = 0;

  if (index == 1)
    pc.topMargin = double(d.topMargin)/pc.height;
  else
    pc.topMargin = 0;

// the way root draws the axis title makes it almost impossible to convert
// the result of the title offset to pixels, so I'll use empirical values
  if (large) {
    if (index == nPads)
      pc.yoffset = 0.72;
    else if (index == 1)
      pc.yoffset = 0.61;
    else
      pc.yoffset = 0.51;
  } else {
    if (index == nPads)
      pc.yoffset = 0.83;
    else if (index == 1)
      pc.yoffset = 0.80;
    else
      pc.yoffset = 0.72;
  }

  return pc;
}


double calculateTitleOffset(int pixels, double xmin, double xmax,
                            double fontsize, bool logx, const Dimensions &d)
{
  if (logx) {
    xmin = log10(xmin);
    xmax = log10(xmax);
  }

  double leftMargin = double(d.leftMargin)/d.width;
  double rightMargin = double(d.rightMargin)/d.width;
  double framewidth = 1 - (leftMargin + rightMargin);

  double dx = xmax - xmin;
  double scale = dx/framewidth;
  double padmin = xmin - leftMargin * scale;
  double padmax = xmax + rightMargin * scale;
  double padwidth = padmax - padmin;

  double titleoffset =
    ((xmin - padmin)/padwidth - double(pixels)/d.width)/(1.6 * fontsize);

  log_debug("Title offset = %g", titleoffset);

  return titleoffset;
}


void save()
{
  gInterpreter->ProcessLine("result.save();");
}


void saveAs(const std::string &filename)
{
  gInterpreter->ProcessLine(("result.saveAs(" + filename + ");").c_str());
}


void gpst_version()
{
  cout << "This is the Gamma-ray Polarimetry Simulation Toolkit version "
       << GPST_VERSION_STRING << " (" << GPST_RELEASE_DATE << ')' << endl;
}

#endif // !__CINT__
