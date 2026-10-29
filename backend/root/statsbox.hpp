// Included inside graphfit.cpp's anonymous namespace, after draw_scene and before
// the histogram/graph painters. ROOT headers belong in the translation unit.
//
// Statistics and fit boxes are ROOT TPaveStats objects created by ROOT's own
// painters (THistPainter/TGraphPainter) from gStyle OptStat/OptFit. Each box is
// created by one pad update and then frozen: detached from its owner and drawn
// as an ordinary pad primitive. This keeps per-pad options independent of the
// process-global gStyle and allows ROOT-like placement next to legends.
std::string formula(const std::string &model);
int nparameters(const std::string &model);

struct StatsBox {
  bool active=false,headroom=false;
  int optstat=0,optfit=0;
  std::string corner="auto",name="Histogram";
  double x1=0,y1=0,x2=0,y2=0;
};

bool root_digits(int value,int ndigits) {
  if(value<0) return false;
  for(int i=0;i<ndigits;++i) { if(value%10>2) return false; value/=10; }
  return value==0;
}

StatsBox read_stats_box(const std::vector<std::vector<std::string>> &rows,size_t &at) {
  StatsBox b;
  if(at>=rows.size() || rows[at].empty() || rows[at][0]!="STATSBOX") return b;
  const auto &f=rows[at++];
  check(f.size()==10,"invalid STATSBOX");
  b.optstat=integer(f[1]); b.optfit=integer(f[2]); b.corner=f[3];
  check(root_digits(b.optstat,9) && root_digits(b.optfit,4),"invalid statistics box options");
  check(b.corner=="auto" || b.corner=="topright" || b.corner=="topleft" || b.corner=="bottomright" ||
        b.corner=="bottomleft" || b.corner=="manual","invalid statistics box position");
  b.x1=real(f[4]); b.y1=real(f[5]); b.x2=real(f[6]); b.y2=real(f[7]);
  if(b.corner=="manual")
    check(0<=b.x1 && b.x1<b.x2 && b.x2<=1 && 0<=b.y1 && b.y1<b.y2 && b.y2<=1,"invalid statistics box NDC");
  text_field(f[8]);
  if(!f[8].empty()) b.name=f[8];
  b.headroom=bit(f[9]);
  b.active=b.optstat>0 || b.optfit>0;
  return b;
}

// Optional histogram bookkeeping after the visible bins: entries, flow bins
// (ROOT indices), and unbinned moments in TH1/TH2/TH3::GetStats order. Without
// moments ROOT recomputes them from bin centres; without entries ROOT's
// ResetStats sets the effective entries.
void read_hist_stats(const std::vector<std::vector<std::string>> &rows,size_t &at,TH1 &hist) {
  if(at>=rows.size() || rows[at].empty() || rows[at][0]!="HSTATS") { hist.ResetStats(); return; }
  const auto &f=rows[at++];
  const int ndim=hist.GetDimension();
  check(f.size()>=4,"invalid HSTATS");
  const bool known=f[1]!="-";
  const double entries=known ? real(f[1]) : 0.0;
  check(entries>=0,"invalid histogram entries");
  const int nflow=integer(f[2]),nm=integer(f[3]);
  const int expected=ndim==1 ? 4 : ndim==2 ? 7 : 11;
  check((nm==0 || nm==expected) && f.size()==static_cast<size_t>(4+nm),"invalid HSTATS moments");
  check(nflow>=0 && nflow<=hist.GetNcells(),"invalid flow-bin count");
  const int nx=hist.GetNbinsX(),ny=hist.GetNbinsY(),nz=hist.GetNbinsZ();
  for(int k=0;k<nflow;++k) {
    check(at<rows.size() && rows[at].size()==static_cast<size_t>(ndim+3) && rows[at][0]=="FLOW","invalid FLOW");
    const auto &r=rows[at++];
    const int ix=integer(r[1]),iy=ndim>1 ? integer(r[2]) : 0,iz=ndim>2 ? integer(r[3]) : 0;
    check(ix>=0 && ix<=nx+1 && iy>=0 && (ndim<2 || iy<=ny+1) && iz>=0 && (ndim<3 || iz<=nz+1),
          "FLOW index out of range");
    const bool flow=ix==0 || ix==nx+1 || (ndim>1 && (iy==0 || iy==ny+1)) || (ndim>2 && (iz==0 || iz==nz+1));
    check(flow,"FLOW must address an under/overflow bin");
    const double w=real(r[ndim+1]),w2=real(r[ndim+2]);
    check(w2>=0,"negative flow-bin variance");
    const int bin=hist.GetBin(ix,iy,iz);
    hist.SetBinContent(bin,w);
    if(hist.GetSumw2N()) hist.SetBinError(bin,std::sqrt(w2));
  }
  if(nm) {
    double s[TH1::kNstat]={};
    for(int i=0;i<nm;++i) s[i]=real(f[4+i]);
    check(s[1]>=0,"invalid histogram moments");
    hist.PutStats(s);
  } else hist.ResetStats();
  if(known) hist.SetEntries(entries);
}

// FITFUNC (model, range, parameters) and optional FITERR (chi2, ndf, errors,
// fixed flags) rows describe a fitted TF1 whose values ROOT prints in a fit box.
std::unique_ptr<TF1> read_fit_function(const std::vector<std::vector<std::string>> &rows,size_t &at) {
  static int counter=0;
  if(at>=rows.size() || rows[at].empty() || rows[at][0]!="FITFUNC") return nullptr;
  const auto &function_row=rows[at++];
  check(function_row.size()>=6,"invalid FITFUNC");
  const std::string model=function_row[1];
  const double fit_min=real(function_row[2]),fit_max=real(function_row[3]);
  const int npar=integer(function_row[4]);
  check(fit_min<fit_max && npar==nparameters(model) &&
        function_row.size()==static_cast<size_t>(5+npar),"invalid FITFUNC model/range");
  auto attached=std::make_unique<TF1>(("rootplots_attached_fit_"+std::to_string(++counter)).c_str(),
                                      formula(model).c_str(),fit_min,fit_max);
  check(attached->IsValid() && attached->GetNpar()==npar,"ROOT rejected attached function");
  for(int i=0;i<npar;++i) attached->SetParameter(i,real(function_row[5+i]));
  if(model=="erfc_edge") check(attached->GetParameter(3)>0,"edge width must be positive");
  if(model=="double_erf") check(attached->GetParameter(4)>0 && attached->GetParameter(5)>0,"edge widths must be positive");
  if(at<rows.size() && !rows[at].empty() && rows[at][0]=="FITERR") {
    const auto &e=rows[at++];
    check(e.size()==static_cast<size_t>(3+2*npar),"invalid FITERR");
    const double chi2=real(e[1]); const int ndf=integer(e[2]);
    check(chi2>=0 && ndf>=0,"invalid fit chi2/ndf");
    for(int i=0;i<npar;++i) {
      const double error=real(e[3+i]);
      check(error>=0,"invalid fit parameter error");
      if(bit(e[3+npar+i])) attached->FixParameter(i,attached->GetParameter(i));
      attached->SetParError(i,error);
    }
    attached->SetChisquare(chi2); attached->SetNDF(ndf);
  }
  return attached;
}

struct FitTarget { std::unique_ptr<TF1> function; int target=0; std::string label; };
FitTarget read_fit_target(const std::vector<std::vector<std::string>> &rows,size_t &at) {
  FitTarget out;
  if(at>=rows.size() || rows[at].empty() || rows[at][0]!="FITTARGET") return out;
  const auto &f=rows[at++];
  check(f.size()==3,"invalid FITTARGET");
  out.target=integer(f[1]); text_field(f[2]); out.label=f[2];
  check(out.target>=1,"invalid fit target");
  out.function=read_fit_function(rows,at);
  check(out.function!=nullptr,"FITTARGET needs FITFUNC");
  return out;
}

bool option_has_palette(const std::string &option,int dimension) {
  if(dimension<2) return false;
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  return upper.find('Z')!=std::string::npos;
}

// Flat (2-D view) palette: the 3-D painters place their palette themselves.
bool flat_palette(const std::string &option,int dimension) {
  if(!option_has_palette(option,dimension)) return false;
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  for(const char *view:{"SURF","LEGO","BOX","ISO","TRI","POL","CYL","SPH","PSR","GL"})
    if(upper.find(view)!=std::string::npos) return false;
  return true;
}
// Palette labels as TGaxis paints them (linear: THLimitsFinder divisions with
// the bin-width decimals, a "#times10^{n}" exponent above the axis beyond the
// axis max digits; logarithmic: decade labels), measured with ROOT's TLatex.
// Widths are fractions of the pad width. An estimate: the margin only needs to
// be large enough, so the label set errs on the wide side.
struct PaletteText { double labels=0,exponent=0,exponent_height=0; };
std::string palette_number(double v,int decimals) {
  char buffer[64];
  std::snprintf(buffer,sizeof buffer,"%.*f",decimals,v);
  std::string s=buffer;
  if(s.find_first_not_of("-0.")==std::string::npos) s=s.substr(s[0]=='-' ? 1 : 0);
  return s[0]=='-' ? "#minus"+s.substr(1) : s;
}
PaletteText palette_text(TVirtualPad &pad,TH1 *h) {
  TAxis *z=h->GetZaxis();
  double lo=h->GetMinimum(),hi=h->GetMaximum();
  if(!std::isfinite(lo) || !std::isfinite(hi)) return {};
  std::vector<std::string> labels;
  std::string exponent;
  if(pad.GetLogz()) {
    if(hi<=0) hi=1;
    if(lo<=0 || lo>=hi) lo=std::min(1.0,1e-3*hi);
    const int top=int(std::ceil(std::log10(hi))),bottom=std::max(int(std::floor(std::log10(lo))),top-60);
    for(int e=bottom;e<=top;++e)
      labels.push_back(e==0 ? "1" : e==1 ? "10" : z->GetNoExponent() ? palette_number(std::pow(10.0,e),std::max(0,-e))
                                                                      : "10^{"+std::to_string(e)+"}");
    if(top-bottom<=2) labels.push_back("5#times10^{"+std::to_string(top-1)+"}");
  } else {
    if(!(hi>lo)) hi=lo+1;
    int ndiv=std::abs(z->GetNdivisions())%100;
    if(ndiv<1) ndiv=10;
    double low=lo,high=hi,width=hi-lo; int bins=1;
    THLimitsFinder::Optimize(lo,hi,ndiv,low,high,bins,width,"");
    if(!(width>0) || bins<1 || bins>1000) { low=lo; high=hi; bins=1; width=hi-lo; }
    const int digits=z->GetMaxDigits()>0 ? z->GetMaxDigits() : gStyle->GetAxisMaxDigits();
    double magnitude=std::max(std::abs(low),std::abs(high));
    int power=0;
    if(!z->GetNoExponent() && magnitude>0) {
      if(width<std::pow(10.0,-digits) && magnitude<1) {
        power=-3*int(std::ceil(-std::log10(magnitude)/3-1e-9));
      } else if(magnitude>=1 && int(std::log10(magnitude)+1e-9)+1>digits) {
        while(true) { ++power; magnitude/=10; if(power%3==0 && magnitude<=std::pow(10.0,digits-1)) break; }
      }
    }
    const double scale=std::pow(10.0,-power),step=width*scale;
    int decimals=0;
    while(decimals<digits+3) {
      const double shifted=step*std::pow(10.0,decimals);
      if(std::abs(shifted-std::round(shifted))<=1e-6*std::max(1.0,shifted)) break;
      ++decimals;
    }
    for(int i=0;i<=bins;++i) labels.push_back(palette_number((low+i*width)*scale,decimals));
    if(power) exponent="#times10^{"+std::to_string(power)+"}";
  }
  PaletteText out;
  const int font=z->GetLabelFont(); const double size=z->GetLabelSize();
  for(const auto &label:labels) out.labels=std::max(out.labels,text_width_ndc(&pad,label,font,size));
  if(!exponent.empty()) {
    TLatex probe; probe.SetTextFont(font); probe.SetTextSize(size); probe.SetTitle(exponent.c_str());
    UInt_t w=0,ht=0; probe.GetBoundingBox(w,ht);
    const double pw=pad.GetWw()*pad.GetAbsWNDC(),ph=pad.GetWh()*pad.GetAbsHNDC();
    if(pw>0 && ph>0) { out.exponent=w/pw; out.exponent_height=ht/ph; }
  }
  return out;
}

// THistPainter::PaintPalette puts the palette (TPaletteAxis) between 0.005 and
// 0.05 pad widths right of the frame; TGaxis paints its labels at the z label
// offset right of it, a "#times10^{n}" exponent above it, and the z title
// 1.6*size*offset pad widths from it. The room this needs grows with the z
// label and title sizes: a selected style such as CMS (title 0.06, labels 0.05)
// needs a wider right margin than ROOT's default (0.035).
// Callers apply the explicit z-axis and canvas styles first. Unless the request
// fixed them (CanvasStyle right/top margin in canvas), the right margin is
// widened until labels, exponent and title end inside the pad, and the top
// margin until the exponent does (small top margins such as ATLAS's 0.05);
// unless the z title offset was fixed (offset=false), it is chosen clear of
// the labels. A titled palette keeps at least the earlier 0.2 margin and 1.3
// offset, so ROOT's default style keeps its look.
void palette_title_room(TVirtualPad &pad,TH1 *h,bool palette,const std::vector<std::string> &canvas,bool offset) {
  // The layout may have widened the right margin since the label was registered.
  if(!palette || !h) { reserve_label_room(&pad); return; }
  TVirtualPad *previous=gPad;
  pad.cd();
  TAxis *z=h->GetZaxis();
  const PaletteText text=palette_text(pad,h);
  double extent=std::max(z->GetLabelOffset()+text.labels,text.exponent),least=0;
  if(z->GetTitle()[0]) {
    const PaletteTitle t=palette_title_geometry(&pad,z);
    if(offset && t.unit>0) {
      // Title centre: past the labels, a gap of 0.3 title sizes, and half the title height.
      const double clear=z->GetLabelOffset()+text.labels+0.3*t.size+0.5*t.thickness;
      z->SetTitleOffset(std::max(1.3,clear/t.unit));
    }
    extent=std::max(extent,t.unit*z->GetTitleOffset()+0.5*t.thickness);
    least=0.2;
  }
  bool moved=false;
  const double right=std::min(0.5,std::max(least,0.05+extent+0.025));
  if(canvas[CV_RIGHT]=="-" && right>pad.GetRightMargin()) { pad.SetRightMargin(right); moved=true; }
  // TGaxis puts the exponent 0.1 label sizes above the palette top (precision-3 fonts: pixels).
  const double ph=pad.GetWh()*pad.GetAbsHNDC(),lift=0.1*z->GetLabelSize()/(z->GetLabelFont()%10==3 && ph>0 ? ph : 1);
  const double top=std::min(0.3,lift+text.exponent_height+0.015);
  if(text.exponent_height>0 && canvas[CV_TOP]=="-" && top>pad.GetTopMargin()) { pad.SetTopMargin(top); moved=true; }
  // An experiment label above the frame needs its room with the final side margins.
  if(reserve_label_room(&pad)) moved=true;
  // A palette painted before (a pad update) keeps its old place: let ROOT redo it.
  if(moved) if(TObject *old=h->GetListOfFunctions()->FindObject("palette")) {
    h->GetListOfFunctions()->Remove(old); delete old;
  }
  pad.Modified();
  if(previous) previous->cd();
}

void stats_begin(const StatsBox &b) {
  gStyle->SetOptStat(b.optstat); gStyle->SetOptFit(b.optfit);
}

TList *stats_owner_functions(TObject *owner) {
  if(auto *h=dynamic_cast<TH1*>(owner)) return h->GetListOfFunctions();
  if(auto *g=dynamic_cast<TGraph*>(owner)) return g->GetListOfFunctions();
  return nullptr;
}

// Test/debug hook: ROOTPLOTS_STATS_DUMP=<file> appends every frozen box with
// its printed lines and ROOT's full-precision values.
void stats_dump(TPaveStats &st,TObject *owner) {
  const char *path=std::getenv("ROOTPLOTS_STATS_DUMP");
  if(!path || !*path) return;
  std::ofstream out(path,std::ios::app);
  check(bool(out),"cannot open ROOTPLOTS_STATS_DUMP");
  out<<std::setprecision(17)<<"BOX\t"<<st.GetX1NDC()<<'\t'<<st.GetY1NDC()<<'\t'
     <<st.GetX2NDC()<<'\t'<<st.GetY2NDC()<<'\t'<<st.GetOptStat()<<'\t'<<st.GetOptFit()<<'\n';
  for(TObject *line:*st.GetListOfLines()) out<<"LINE\t"<<line->GetTitle()<<'\n';
  if(auto *h=dynamic_cast<TH1*>(owner)) {
    out<<"VALUE\tentries\t"<<h->GetEntries()<<'\n';
    for(int axis=1;axis<=h->GetDimension() || (axis==2 && h->InheritsFrom("TProfile"));++axis)
      out<<"VALUE\tmean"<<axis<<'\t'<<h->GetMean(axis)<<"\nVALUE\tstd"<<axis<<'\t'<<h->GetStdDev(axis)<<'\n';
  }
  TF1 *fit=nullptr;
  if(TList *list=stats_owner_functions(owner))
    for(TObject *o:*list) if(o->InheritsFrom(TF1::Class())) { fit=static_cast<TF1*>(o); break; }
  if(fit) {
    out<<"VALUE\tchi2\t"<<fit->GetChisquare()<<"\nVALUE\tndf\t"<<fit->GetNDF()
       <<"\nVALUE\tprob\t"<<fit->GetProb()<<'\n';
    for(int i=0;i<fit->GetNpar();++i)
      out<<"PARAM\t"<<fit->GetParName(i)<<'\t'<<fit->GetParameter(i)<<'\t'<<fit->GetParError(i)<<'\n';
  }
  out<<"END\n";
  check(bool(out),"failed writing ROOTPLOTS_STATS_DUMP");
}

// NDC width at which every "label = value" line fits at ROOT's height-limited
// text size (0.85 of the line pitch). ROOT shrinks the text of narrow boxes
// instead; widening keeps many-line or "±" boxes readable.
double stats_required_width(TPaveStats &st,TPad &pad,double height) {
  TList *lines=st.GetListOfLines();
  const int nlines=lines ? lines->GetSize() : 0;
  const double pixels=pad.GetWw()*pad.GetAbsWNDC();
  if(nlines==0 || pixels<=0) return 0;
  pad.cd();
  TLatex probe; probe.SetNDC(); probe.SetTextFont(st.GetTextFont());
  probe.SetTextSize(0.85*height/nlines);
  auto measure=[&](std::string text) {
    while(!text.empty() && text.back()==' ') text.pop_back();
    while(!text.empty() && text.front()==' ') text.erase(text.begin());
    if(text.empty()) return 0.0;
    probe.SetTitle(text.c_str());
    UInt_t w=0,h=0; probe.GetBoundingBox(w,h);
    return static_cast<double>(w);
  };
  double widest=0;
  for(TObject *line:*lines) {
    const std::string text=line->GetTitle();
    if(text.find('|')!=std::string::npos) continue; // 2-D flow table, drawn in columns
    const auto split=text.find('=');
    const double need=split==std::string::npos ? measure(text) :
      measure(text.substr(0,split))+measure(text.substr(split+1))+0.6*probe.GetTextSize()*pad.GetWh()*pad.GetAbsHNDC();
    widest=std::max(widest,need);
  }
  return widest/pixels/0.9;
}

struct FrozenStats {
  std::vector<std::unique_ptr<TPaveStats>> boxes;
  StatsRect rect;
};

// Paint the pad once so ROOT creates the boxes of the given owners, then detach,
// place, and redraw them as pad primitives. Later boxes (overlays drawn with
// SAMES, or a graph fit next to a histogram) are stacked below the first one.
FrozenStats freeze_stats(TPad &pad,const StatsBox &b,const std::vector<TObject*> &owners,bool palette) {
  FrozenStats out;
  if(!b.active) return out;
  pad.cd(); pad.Modified(); pad.Update();
  std::vector<std::pair<TPaveStats*,TObject*>> found;
  for(TObject *owner:owners) {
    TList *list=owner ? stats_owner_functions(owner) : nullptr;
    if(!list) continue;
    TPaveStats *st=nullptr;
    for(TObject *o:*list) if(o->InheritsFrom(TPaveStats::Class())) { st=static_cast<TPaveStats*>(o); break; }
    if(!st) continue;
    list->Remove(st);
    if(auto *h=dynamic_cast<TH1*>(owner)) h->SetStats(false);
    else owner->SetBit(TGraph::kNoStats);
    st->SetParent(nullptr); st->ResetBit(TObject::kCanDelete);
    out.boxes.emplace_back(st);
    found.emplace_back(st,owner);
  }
  const double lm=pad.GetLeftMargin(),rm=pad.GetRightMargin(),bm=pad.GetBottomMargin(),tm=pad.GetTopMargin();
  const bool bottom=b.corner=="bottomright" || b.corner=="bottomleft";
  double right=0,top=0,low=0;
  for(size_t k=0;k<found.size();++k) {
    TPaveStats *st=found[k].first;
    double w=st->GetX2NDC()-st->GetX1NDC(),h=st->GetY2NDC()-st->GetY1NDC(),x1=st->GetX1NDC(),y1=st->GetY1NDC();
    if(b.corner!="manual") {
      const double need=std::min(0.6,stats_required_width(*st,pad,h));
      if(need>w) { x1-=need-w; w=need; } // keep ROOT's right edge
    }
    if(k==0) {
      if(b.corner=="manual") { x1=b.x1; y1=b.y1; w=b.x2-b.x1; h=b.y2-b.y1; }
      else if(b.corner=="auto") { if(palette) x1=1-rm-0.01-w; }
      else {
        x1=b.corner.find("left")!=std::string::npos ? lm+0.02 : 1-rm-0.02-w;
        y1=bottom ? bm+0.02 : 1-tm-0.02-h;
      }
      right=x1+w; top=y1+h; low=y1;
    } else {
      x1=right-w;
      if(bottom) { y1=top; top+=h; } else { y1=low-h; low=y1; }
    }
    st->SetX1NDC(x1); st->SetY1NDC(y1); st->SetX2NDC(x1+w); st->SetY2NDC(y1+h);
    if(found.size()>1) {
      if(auto *hist=dynamic_cast<TH1*>(found[k].second)) st->SetTextColor(hist->GetLineColor());
    }
    st->Draw();
    stats_dump(*st,found[k].second);
    if(!out.rect.on) out.rect={true,x1,y1,x1+w,y1+h};
    else out.rect={true,std::min(out.rect.x1,x1),std::min(out.rect.y1,y1),
                   std::max(out.rect.x2,x1+w),std::max(out.rect.y2,y1+h)};
  }
  gStyle->SetOptStat(0); gStyle->SetOptFit(0);
  pad.Modified();
  return out;
}

// Highest data point under the boxes as NDC (pad user coordinates are log10 on
// log axes); -inf when no data lie in the boxes' horizontal span.
double stats_data_top(TPad &pad,const StatsRect &rect,const std::vector<TObject*> &data);
// Test/debug hook (ROOTPLOTS_STATS_DUMP): final 1-D frame y range in pad
// coordinates, highest data point under the boxes and lowest box edge (NDC).
void stats_dump_frame(TPad &pad,const StatsBox &b,const StatsRect &rect,const std::vector<TObject*> &data) {
  const char *path=std::getenv("ROOTPLOTS_STATS_DUMP");
  if(!path || !*path) return;
  std::ofstream out(path,std::ios::app);
  check(bool(out),"cannot open ROOTPLOTS_STATS_DUMP");
  out<<std::setprecision(17)<<"FRAME\t"<<int(b.headroom)<<'\t'<<pad.GetUymin()<<'\t'<<pad.GetUymax()<<'\t'
     <<stats_data_top(pad,rect,data)<<'\t'<<rect.y1<<'\n';
  check(bool(out),"failed writing ROOTPLOTS_STATS_DUMP");
}

// Opt-in (stats_headroom) remedy ROOT users apply when a box covers data:
// raise the y maximum of the 1-D frame until every drawn bin, point, error bar
// and function value in the horizontal span of the boxes lies below them. The
// box positions and the y minimum are unchanged; nothing changes when the data
// already clear the boxes or the boxes sit in the lower half of the frame.
void stats_headroom(TPad &pad,const StatsBox &b,const StatsRect &rect,TH1 *frame,
                    const std::vector<TObject*> &data) {
  if(!rect.on || !frame || frame->GetDimension()!=1) return;
  pad.cd(); pad.Modified(); pad.Update();
  const double fb=pad.GetBottomMargin(),ft=1-pad.GetTopMargin();
  // Only boxes in the upper part of the frame can be cleared by a higher maximum.
  const bool upper=b.corner=="auto" || b.corner=="topright" || b.corner=="topleft" ||
                   (b.corner=="manual" && rect.y1+rect.y2>fb+ft);
  if(b.headroom && upper) {
    const double uy1=pad.GetUymin(),uy2=pad.GetUymax();
    const double high=stats_data_top(pad,rect,data);
    const double target=rect.y1-0.02; // NDC gap below the lowest box
    if(std::isfinite(high) && high>target && target>fb+0.1*(ft-fb) && uy2>uy1) {
      // Pad coordinate of the data top, then the maximum that maps it to target.
      const double y=uy1+(high-fb)/(ft-fb)*(uy2-uy1);
      const double top=uy1+(y-uy1)*(ft-fb)/(target-fb);
      const bool logy=pad.GetLogy();
      frame->SetMinimum(logy ? std::pow(10.0,uy1) : uy1);
      frame->SetMaximum(logy ? std::pow(10.0,top) : top);
      pad.Modified(); pad.Update();
    }
  }
  stats_dump_frame(pad,b,rect,data);
}

double stats_data_top(TPad &pad,const StatsRect &rect,const std::vector<TObject*> &data) {
  const double fl=pad.GetLeftMargin(),fr=1-pad.GetRightMargin(),fb=pad.GetBottomMargin(),ft=1-pad.GetTopMargin();
  const double ux1=pad.GetUxmin(),ux2=pad.GetUxmax(),uy1=pad.GetUymin(),uy2=pad.GetUymax();
  const bool logx=pad.GetLogx(),logy=pad.GetLogy();
  if(!rect.on || !(fr>fl && ft>fb && ux2>ux1 && uy2>uy1)) return -INFINITY;
  // Box span in pad coordinates (log10 on log axes).
  const double bx1=ux1+(rect.x1-fl)/(fr-fl)*(ux2-ux1),bx2=ux1+(rect.x2-fl)/(fr-fl)*(ux2-ux1);
  auto px=[&](double x){ return logx ? (x>0 ? std::log10(x) : -INFINITY) : x; };
  double high=-INFINITY;
  auto consider=[&](double xlo,double xhi,double y) {
    if(!std::isfinite(y) || (logy && y<=0)) return;
    if(px(xhi)<bx1 || px(xlo)>bx2 || px(xhi)<ux1 || px(xlo)>ux2) return;
    high=std::max(high,logy ? std::log10(y) : y);
  };
  auto function=[&](TF1 *f) {
    double lo=0,hi=0; f->GetRange(lo,hi);
    for(int i=0;i<=200;++i) { const double x=lo+(hi-lo)*i/200.0; consider(x,x,f->Eval(x)); }
  };
  auto functions=[&](TList *list) {
    if(list) for(TObject *o:*list) if(o->InheritsFrom(TF1::Class())) function(static_cast<TF1*>(o));
  };
  for(TObject *o:data) {
    if(!o) continue;
    if(auto *h=dynamic_cast<TH1*>(o)) {
      if(h->GetDimension()!=1) continue;
      // Error bars count only when drawn: an E option, or ROOT's default for Sumw2 histograms.
      std::string option=h->GetDrawOption();
      std::transform(option.begin(),option.end(),option.begin(),[](unsigned char c){return std::toupper(c);});
      for(const std::string same:{"SAMES","SAME"})
        for(size_t k;(k=option.find(same))!=std::string::npos;) option.erase(k,same.size());
      option.erase(std::remove(option.begin(),option.end(),' '),option.end());
      const bool errors=option.find('E')!=std::string::npos || (option.empty() && h->GetSumw2N()>0);
      for(int i=1;i<=h->GetNbinsX();++i)
        consider(h->GetXaxis()->GetBinLowEdge(i),h->GetXaxis()->GetBinUpEdge(i),
                 h->GetBinContent(i)+(errors ? h->GetBinError(i) : 0.0));
      functions(h->GetListOfFunctions());
    } else if(auto *g=dynamic_cast<TGraph*>(o)) {
      for(int i=0;i<g->GetN();++i) {
        const double ex1=std::max(0.0,g->GetErrorXlow(i)),ex2=std::max(0.0,g->GetErrorXhigh(i));
        consider(g->GetX()[i]-ex1,g->GetX()[i]+ex2,g->GetY()[i]+std::max(0.0,g->GetErrorYhigh(i)));
      }
      functions(g->GetListOfFunctions());
    } else if(auto *f=dynamic_cast<TF1*>(o)) function(f);
  }
  return std::isfinite(high) ? fb+(high-uy1)/(uy2-uy1)*(ft-fb) : -INFINITY;
}
