// ROOT styles (gROOT->SetStyle) and experiment labels (TLatex), shared by
// graphfit.cpp (included inside its anonymous namespace, after axisstyle.hpp)
// and render.cpp. ROOT and standard headers belong in the translation units
// (TCanvas, TLatex, TList, TPad, TROOT, TStyle, <algorithm>, <cctype>, <string>,
// <vector>). The including file provides check() and integer().
//
// Nothing here runs unless a request selects a style or a label: without them
// gStyle is never replaced and no extra object is drawn.

// Name of the selected TStyle; empty = ROOT's default gStyle, untouched.
std::string &selected_style() { static std::string name; return name; }
bool style_selected() { return !selected_style().empty(); }

std::string lowercase_ascii(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){ return std::tolower(c); });
  return s;
}

// "CMS" TStyle, built from ROOT TStyle setters with values that follow the
// public CMS publication style conventions (white canvas/pad/frame without
// borders, no statistics or title box, Helvetica (font 42) text, large axis
// titles and labels, ticks on all four frame sides). ROOT has no built-in CMS
// style; this one is created once per process and registered in gROOT's style
// list, so it is selected like ROOT's own styles.
TStyle *cms_style() {
  if(TStyle *existing=gROOT->GetStyle("CMS")) return existing;
  auto *s=new TStyle("CMS","CMS publication style conventions");
  // Canvas, pads and frame: white, no borders, no grid.
  s->SetCanvasBorderMode(0); s->SetCanvasColor(kWhite);
  s->SetPadBorderMode(0); s->SetPadColor(kWhite); s->SetPadGridX(false); s->SetPadGridY(false);
  s->SetGridColor(0); s->SetGridStyle(3); s->SetGridWidth(1);
  s->SetFrameBorderMode(0); s->SetFrameBorderSize(1); s->SetFrameFillColor(0); s->SetFrameFillStyle(0);
  s->SetFrameLineColor(1); s->SetFrameLineStyle(1); s->SetFrameLineWidth(1);
  // Histograms, markers, error bars and functions.
  s->SetHistLineColor(1); s->SetHistLineStyle(0); s->SetHistLineWidth(1);
  s->SetEndErrorSize(2); s->SetMarkerStyle(20);
  s->SetFuncColor(2); s->SetFuncStyle(1); s->SetFuncWidth(2);
  // No date/file/statistics/title boxes; box formats when a box is requested.
  s->SetOptDate(0); s->SetOptFile(0); s->SetOptStat(0); s->SetOptFit(0); s->SetOptTitle(0);
  s->SetFitFormat("5.4g"); s->SetStatFormat("6.4g");
  s->SetStatColor(kWhite); s->SetStatFont(42); s->SetStatFontSize(0.025); s->SetStatTextColor(1);
  s->SetStatBorderSize(1); s->SetStatH(0.1); s->SetStatW(0.15);
  // Pad margins (fractions of the pad): room for the large axis titles and the
  // label/luminosity line above the frame.
  s->SetPadTopMargin(0.08); s->SetPadBottomMargin(0.13);
  s->SetPadLeftMargin(0.16); s->SetPadRightMargin(0.04);
  // Plot title box (only drawn with CanvasStyle(opt_title=true)).
  s->SetTitleFont(42); s->SetTitleColor(1); s->SetTitleTextColor(1); s->SetTitleFillColor(10);
  s->SetTitleFontSize(0.05); s->SetTitleBorderSize(0);
  // Axis titles and labels.
  s->SetTitleColor(1,"XYZ"); s->SetTitleFont(42,"XYZ"); s->SetTitleSize(0.06,"XYZ");
  s->SetTitleXOffset(0.9); s->SetTitleYOffset(1.25);
  s->SetLabelColor(1,"XYZ"); s->SetLabelFont(42,"XYZ"); s->SetLabelOffset(0.007,"XYZ");
  s->SetLabelSize(0.05,"XYZ");
  // Axes: ticks on both sides of the frame.
  s->SetAxisColor(1,"XYZ"); s->SetStripDecimals(kTRUE); s->SetTickLength(0.03,"XYZ");
  s->SetNdivisions(510,"XYZ"); s->SetPadTickX(1); s->SetPadTickY(1);
  s->SetOptLogx(0); s->SetOptLogy(0); s->SetOptLogz(0);
  // Text and legends.
  s->SetTextFont(42); s->SetTextSize(0.05);
  s->SetLegendBorderSize(0); s->SetLegendFillColor(0); s->SetLegendFont(42);
  // PostScript/PDF output.
  s->SetPaperSize(20.,20.);
  return s;
}

// Names accepted by select_root_style: gROOT's style list plus "CMS".
std::vector<std::string> available_styles() {
  cms_style();
  std::vector<std::string> names;
  for(TObject *o:*gROOT->GetListOfStyles()) names.push_back(o->GetName());
  std::sort(names.begin(),names.end());
  names.erase(std::unique(names.begin(),names.end()),names.end());
  return names;
}

// gROOT->SetStyle(name) before any canvas or histogram exists, so every object
// takes its attributes from the style (no ForceStyle: nothing is read from a
// file, and explicit setters applied later must keep their values). Names are
// matched case-insensitively against gROOT's style list. The package keeps its
// opt-in statistics/fit boxes and title box off until a request asks for them.
void select_root_style(const std::string &name) {
  check(!name.empty() && name.size()<=64,"invalid ROOT style name");
  for(unsigned char c:name) check(std::isalnum(c) || c=='_' || c=='-',"invalid ROOT style name");
  if(lowercase_ascii(name)=="cms") cms_style();
  TStyle *found=nullptr;
  for(TObject *o:*gROOT->GetListOfStyles())
    if(lowercase_ascii(o->GetName())==lowercase_ascii(name)) { found=static_cast<TStyle*>(o); break; }
  if(!found) {
    std::string list;
    for(const auto &n:available_styles()) list+=(list.empty() ? "" : ", ")+n;
    throw std::runtime_error("unknown ROOT style \""+name+"\" (available: "+list+")");
  }
  gROOT->SetStyle(found->GetName());
  check(gStyle==found,"ROOT did not select style "+name);
  selected_style()=found->GetName();
  gStyle->SetOptStat(0); gStyle->SetOptFit(0); gStyle->SetOptTitle(0);
}

// A request may name its style in several rows (one per figure panel); all must agree.
template<class Rows> void select_requested_style(const Rows &rows,const std::string &tag) {
  std::string name;
  for(const auto &row:rows) {
    if(row.empty() || row[0]!=tag) continue;
    check(row.size()==2 && !row[1].empty(),"invalid "+tag+" row");
    check(name.empty() || lowercase_ascii(name)==lowercase_ascii(row[1]),
          "pads request different ROOT styles (the style is process-global)");
    name=row[1];
  }
  if(!name.empty()) select_root_style(name);
}

// Experiment label of one pad (EXPLABEL row): kind ("CMS" or "ATLAS"),
// position, extra text after/below the experiment name, text above the frame
// on the right (CMS luminosity/energy), and further lines below the label.
struct ExperimentLabelSpec {
  bool active=false;
  std::string kind,position,extra,right;
  std::vector<std::string> lines;
};
void label_text(const std::string &s) {
  check(s.find('\0')==std::string::npos && s.find('\r')==std::string::npos && s.size()<=512,"invalid label text");
}
// EXPLABEL kind position extra right nlines line...  ("-" = empty text)
ExperimentLabelSpec read_experiment_label(const std::vector<std::string> &f) {
  check(f.size()>=6 && f[0]=="EXPLABEL","invalid EXPLABEL");
  ExperimentLabelSpec l;
  l.active=true; l.kind=f[1]; l.position=f[2];
  l.extra=f[3]=="-" ? "" : f[3]; l.right=f[4]=="-" ? "" : f[4];
  check(l.kind=="CMS" || l.kind=="ATLAS","invalid experiment label kind");
  check(l.position=="inside_left" || l.position=="inside_right" || l.position=="outside","invalid label position");
  const int n=integer(f[5]);
  check(n>=0 && n<=8 && f.size()==size_t(6+n),"invalid EXPLABEL line count");
  for(int i=0;i<n;++i) l.lines.push_back(f[6+i]);
  label_text(l.extra); label_text(l.right);
  for(const auto &s:l.lines) label_text(s);
  return l;
}

// Text geometry of a pad that shares the canvas height with another pad
// (fit/pull, comparison). ROOT text sizes are fractions of the smaller pad
// dimension; k converts a one-pad size (a pad of the same width and the full
// canvas height) into this pad's size, hr is the pad's share of the canvas height.
struct PadText { double k=1,hr=1; };
PadText two_pad_text(TVirtualPad *pad) {
  const double w=pad->GetWw()*pad->GetAbsWNDC(),h=pad->GetWh()*pad->GetAbsHNDC(),H=pad->GetWh();
  if(w<=0 || h<=0 || H<=0) return {};
  return {std::min(w,H)/std::min(w,h),h/H};
}

// Pads that receive a label. owner identifies the request scene, so a scene
// applied to several pads (fit/pull) labels only its first (upper) pad.
// two_pad marks the upper pad of a two-pad layout (sizes as in a one-pad plot);
// fixed_top marks a top margin set by the request (CanvasStyle), never grown.
struct LabelSlot { TVirtualPad *pad; ExperimentLabelSpec label; int owner; bool two_pad,fixed_top; };
std::vector<LabelSlot> &label_slots() { static std::vector<LabelSlot> slots; return slots; }
std::vector<TVirtualPad*> &two_pad_uppers() { static std::vector<TVirtualPad*> pads; return pads; }
int next_label_owner() { static int owner=0; return ++owner; }
bool reserve_label_room(TVirtualPad *pad);
void register_label(TVirtualPad *pad,const ExperimentLabelSpec &l,int owner,bool fixed_top=false) {
  if(!l.active || !pad) return;
  for(const auto &s:label_slots()) if(s.owner==owner) return;
  const auto &uppers=two_pad_uppers();
  const bool two_pad=std::find(uppers.begin(),uppers.end(),pad)!=uppers.end();
  label_slots().push_back({pad,l,owner,two_pad,fixed_top});
  reserve_label_room(pad);
}
// May be called before or after the pad's label is registered.
void mark_two_pad(TVirtualPad *pad) {
  two_pad_uppers().push_back(pad);
  for(auto &s:label_slots()) if(s.pad==pad) s.two_pad=true;
  reserve_label_room(pad);
}

// NDC width of a text in pad, measured by ROOT.
double text_width_ndc(TVirtualPad *pad,const std::string &s,int font,double size) {
  if(s.empty()) return 0;
  TLatex probe; probe.SetTextFont(font); probe.SetTextSize(size); probe.SetTitle(s.c_str());
  UInt_t w=0,h=0; probe.GetBoundingBox(w,h);
  const double pixels=pad->GetWw()*pad->GetAbsWNDC();
  return pixels>0 ? w/pixels : 0;
}
// NDC width of a word space at text size size (TLatex sizes are fractions of
// the smaller pad dimension).
double text_gap_ndc(TVirtualPad *pad,double size) {
  const double w=pad->GetWw()*pad->GetAbsWNDC(),h=pad->GetWh()*pad->GetAbsHNDC();
  return w>0 ? 0.28*size*std::min(w,h)/w : 0.01;
}

// NDC height of a text's TLatex box (ascent plus descent) in pad, measured by
// ROOT; the pad must be gPad.
double text_height_ndc(TVirtualPad *pad,const std::string &s,int font,double size) {
  if(s.empty()) return 0;
  TLatex probe; probe.SetTextFont(font); probe.SetTextSize(size); probe.SetTitle(s.c_str());
  UInt_t w=0,h=0; probe.GetBoundingBox(w,h);
  const double pixels=pad->GetWh()*pad->GetAbsHNDC();
  return pixels>0 ? h/pixels : 0;
}

// One text of a label in pad NDC. above: drawn above the frame, aligned on its
// baseline (vertical alignment 1), so the text spans y to y + its box height.
struct LabelText { double x,y; std::string s; int font; double size; int align; bool above; };

// Texts of one label, positioned from the pad's current margins (gPad = pad).
// CMS: bold "CMS" (font 61) with the extra text in italics (font 52), sizes
// relative to the top margin as in CMS publications; the right text
// (luminosity and energy, font 42) sits above the frame on the right.
// ATLAS: bold-italic "ATLAS" (font 72) followed by regular text (font 42),
// further lines (energy/luminosity) below, at the style's text size. Outside
// the frame, ATLAS puts "ATLAS <text>" on the left and the further lines,
// joined, on the right of one row above the frame; when the row is too narrow
// for both, the lines are stacked under "ATLAS <text>", all left-aligned.
std::vector<LabelText> label_layout(const LabelSlot &slot) {
  TVirtualPad *pad=slot.pad;
  const auto &L=slot.label;
  const double l=pad->GetLeftMargin(),r=pad->GetRightMargin(),t=pad->GetTopMargin(),b=pad->GetBottomMargin();
  const double fw=1-l-r,fh=1-t-b;
  const PadText g=slot.two_pad ? two_pad_text(pad) : PadText();
  // Text sizes are fractions of the smaller pad dimension; to_y converts one into a height fraction.
  const double pw=pad->GetWw()*pad->GetAbsWNDC(),ph=pad->GetWh()*pad->GetAbsHNDC();
  const double to_y=pw>0 && ph>0 ? std::min(pw,ph)/ph : 1;
  std::vector<LabelText> out;
  auto put=[&](double x,double y,const std::string &s,int font,double size,int align,bool above) {
    if(!s.empty()) out.push_back({x,y,s,font,size,align,above});
  };
  const bool left=L.position!="inside_right";
  if(L.kind=="CMS") {
    // Scaled with the (one-pad) top margin, capped for unusually tall margins.
    const double tm=std::min(t*g.hr,0.1),cms=0.75*tm*g.k,extra=0.76*cms,small=0.6*tm*g.k;
    const double above=1-t+0.2*tm/g.hr;
    put(1-r,above,L.right,42,small,31,true);
    double x=left ? l+0.045*fw : 1-r-0.045*fw, y=1-t-0.035*fh;
    if(L.position=="outside") {
      put(l,above,"CMS",61,cms,11,true);
      put(l+text_width_ndc(pad,"CMS",61,cms)+text_gap_ndc(pad,cms),above,L.extra,52,extra,11,true);
      for(const auto &line:L.lines) { put(x,y,line,42,small,13,false); y-=1.2*small*to_y; }
      return out;
    }
    const int align=left ? 13 : 33;
    put(x,y,"CMS",61,cms,align,false); y-=1.2*cms*to_y;
    if(!L.extra.empty()) { put(x,y,L.extra,52,extra,align,false); y-=1.2*extra*to_y; }
    for(const auto &line:L.lines) { put(x,y,line,42,small,align,false); y-=1.2*small*to_y; }
    return out;
  }
  const double size=gStyle->GetTextSize()*g.k;
  const double gap=text_gap_ndc(pad,size);
  if(L.position=="outside") {
    // Baseline a quarter text size above the frame.
    const double base=1-t+0.25*size*to_y;
    const double head=text_width_ndc(pad,"ATLAS",72,size)+(L.extra.empty() ? 0 : gap+text_width_ndc(pad,L.extra,42,size));
    std::string joined;
    for(const auto &line:L.lines) joined+=(joined.empty() ? "" : ", ")+line;
    const bool one_row=joined.empty() || l+head+3*gap<=1-r-text_width_ndc(pad,joined,42,size);
    const size_t rows=one_row ? 0 : L.lines.size();
    const double y=base+rows*1.2*size*to_y;
    put(l,y,"ATLAS",72,size,11,true);
    put(l+text_width_ndc(pad,"ATLAS",72,size)+gap,y,L.extra,42,size,11,true);
    if(one_row) put(1-r,y,joined,42,size,31,true);
    else for(size_t i=0;i<rows;++i) put(l,base+(rows-1-i)*1.2*size*to_y,L.lines[i],42,size,11,true);
    return out;
  }
  double y=1-t-0.045*fh;
  if(left) {
    const double x=l+0.05*fw;
    put(x,y,"ATLAS",72,size,13,false);
    put(x+text_width_ndc(pad,"ATLAS",72,size)+gap,y,L.extra,42,size,13,false);
    for(const auto &line:L.lines) { y-=1.2*size*to_y; put(x,y,line,42,size,13,false); }
  } else {
    const double x=1-r-0.05*fw;
    const double w=text_width_ndc(pad,L.extra,42,size);
    put(x,y,L.extra,42,size,33,false);
    put(L.extra.empty() ? x : x-w-gap,y,"ATLAS",72,size,33,false);
    for(const auto &line:L.lines) { y-=1.2*size*to_y; put(x,y,line,42,size,33,false); }
  }
  return out;
}

// How far the label's texts above the frame reach past the pad top, keeping a
// white edge of a tenth of their size (at least 3 pixels); <= 0 when they fit.
double label_overflow(const LabelSlot &slot,const std::vector<LabelText> &texts) {
  TVirtualPad *pad=slot.pad;
  const double pw=pad->GetWw()*pad->GetAbsWNDC(),ph=pad->GetWh()*pad->GetAbsHNDC();
  const double to_y=pw>0 && ph>0 ? std::min(pw,ph)/ph : 1;
  double over=-1;
  for(const auto &x:texts) if(x.above) {
    const double edge=std::max(0.1*x.size*to_y,ph>0 ? 3/ph : 0.0);
    over=std::max(over,x.y+text_height_ndc(pad,x.s,x.font,x.size)+edge-1);
  }
  return over;
}

// Grow the top margins of pad's labels (not fixed by the request) until the
// texts above the frame fit inside the pad, e.g. the outside ATLAS label under
// the ATLAS style's 0.05 top margin. Runs while the plot is being set up (the
// margins are final before ROOT paints frames, palettes, legends and boxes);
// layouts that widen a margin later run it again. Returns whether a margin changed.
bool reserve_label_room(TVirtualPad *pad) {
  bool changed=false;
  TVirtualPad *previous=gPad;
  for(const auto &slot:label_slots()) {
    if(slot.pad!=pad || slot.fixed_top) continue;
    pad->cd();
    // CMS sizes follow the top margin, so a few rounds may be needed.
    for(int round=0;round<12;++round) {
      const double over=label_overflow(slot,label_layout(slot));
      const double t=pad->GetTopMargin(),limit=std::min(0.35,1-pad->GetBottomMargin()-0.3);
      const double next=std::min(limit,t+over);
      if(over<=1e-6 || next<=t+1e-6) break;
      pad->SetTopMargin(next); changed=true;
    }
  }
  if(previous) previous->cd();
  return changed;
}

// Draw one label with ROOT TLatex in its pad, positioned from the pad's final
// margins (label_layout). Texts above the frame that would still end past the
// pad top (a fixed top margin too small for them) move down together.
void draw_experiment_label(const LabelSlot &slot) {
  TVirtualPad *pad=slot.pad;
  pad->cd();
  const auto texts=label_layout(slot);
  const double shift=std::max(0.0,label_overflow(slot,texts));
  TLatex text; text.SetNDC(); text.SetTextColor(kBlack); text.SetTextAngle(0);
  for(const auto &x:texts) {
    text.SetTextFont(x.font); text.SetTextSize(x.size); text.SetTextAlign(x.align);
    text.DrawLatexNDC(x.x,x.above ? x.y-shift : x.y,x.s.c_str());
  }
}

// Two-pad layouts (fit/pull, comparison) under a selected style: the style's
// one-pad text sizes, offsets and x tick lengths are converted so that text
// and distances keep their size in a one-pad plot of the same width and the
// full canvas height. TGaxis places a y title 1.6*size*offset pad widths from
// the axis and an x title 1.3*size*offset pad heights below it; x label
// offsets and x tick lengths are pad-height fractions.
void style_scaled_axis(TAxis *a,const char *which,TVirtualPad *pad,bool horizontal) {
  if(!a || !pad) return;
  const PadText g=two_pad_text(pad);
  a->SetTitleSize(gStyle->GetTitleSize(which)*g.k);
  a->SetLabelSize(gStyle->GetLabelSize(which)*g.k);
  if(horizontal) {
    a->SetTitleOffset(gStyle->GetTitleOffset(which)/(g.k*g.hr));
    a->SetLabelOffset(gStyle->GetLabelOffset(which)/g.hr);
    a->SetTickLength(gStyle->GetTickLength(which)/g.hr);
  } else {
    a->SetTitleOffset(gStyle->GetTitleOffset(which)/g.k);
  }
}

// Last step before saving: draw the registered labels, then run the
// ROOTPLOTS_AXIS_DUMP hook. Without labels the canvas is not touched.
void finish_canvas(TVirtualPad &canvas) {
  if(!label_slots().empty()) {
    for(const auto &slot:label_slots()) { draw_experiment_label(slot); slot.pad->Modified(); }
    canvas.cd(); canvas.Modified(); canvas.Update();
  }
  axis_dump(canvas);
}
