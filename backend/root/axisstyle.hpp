// ROOT axis, pad-title and axis-dump helpers shared by graphfit.cpp (included
// inside its anonymous namespace) and render.cpp. ROOT and standard headers
// belong in the translation units (TAxis, TAxisModLab, TColor, TF1, TGaxis,
// TGraph2D, TH1, THStack, THashList, TLatex, TList, TPad, TPaveText, TROOT, TStyle,
// TText, <fstream>, <iomanip>, <string>, <vector>). The including file provides
// check(), integer(), real(), bit() and color().
//
// Every value maps to one ROOT setter; "-" means "not set": ROOT's (or the
// painter's documented) default stays in place.

// AXIS row scalar fields, in protocol order (see src/scene.jl axis_scalars).
enum AxisFieldIndex {
  AX_TITLE_FONT, AX_TITLE_SIZE, AX_TITLE_OFFSET, AX_LABEL_FONT, AX_LABEL_SIZE,
  AX_LABEL_OFFSET, AX_TICK_LENGTH, AX_DIVISIONS, AX_MAX_DIGITS, AX_CENTER_TITLE,
  AX_OPTIMIZE, AX_NO_EXPONENT, AX_MORE_LOG_LABELS, AX_DECIMALS, AX_CENTER_LABELS,
  AX_LABEL_COLOR, AX_TITLE_COLOR, AX_AXIS_COLOR, AX_TICKS, AX_LABELS_OPTION,
  AX_TIME_DISPLAY, AX_TIME_OFFSET, AX_EXPONENT_X, AX_EXPONENT_Y, AX_FIELDS
};

struct LabelChangeSpec {
  int index=0,align=-1,font=-1; double angle=-1,size=-1;
  std::string hex="-"; bool has_text=false; std::string text;
};
struct AxisSpec {
  std::vector<std::string> v=std::vector<std::string>(AX_FIELDS,"-");
  bool has_time_format=false;
  std::string time_format;
  std::vector<std::pair<int,std::string>> bins;
  std::vector<LabelChangeSpec> changes;
};

bool axis_font_code(int f) { return f>=10 && f<=159 && f%10<=3; }
void axis_text(const std::string &s) {
  check(s.find('\0')==std::string::npos && s.find('\r')==std::string::npos,"invalid axis text");
}

// Scalars of one AXIS/ZSTYLE/PAXIS row, starting at f[first].
void axis_scalars(AxisSpec &a,const std::vector<std::string> &f,size_t first) {
  check(f.size()==first+AX_FIELDS,"invalid axis field count");
  a.v.assign(f.begin()+first,f.end());
  const auto &v=a.v;
  auto set=[&](int i){ return v[i]!="-"; };
  for(int i:{AX_TITLE_FONT,AX_LABEL_FONT}) if(set(i)) check(axis_font_code(integer(v[i])),"invalid axis font");
  for(int i:{AX_TITLE_SIZE,AX_TITLE_OFFSET,AX_LABEL_SIZE}) if(set(i)) check(real(v[i])>0,"invalid axis size/offset");
  if(set(AX_LABEL_OFFSET)) real(v[AX_LABEL_OFFSET]);
  if(set(AX_TICK_LENGTH)) check(real(v[AX_TICK_LENGTH])>=0 && real(v[AX_TICK_LENGTH])<=1,"invalid tick length");
  if(set(AX_DIVISIONS)) check(integer(v[AX_DIVISIONS])>=1 && integer(v[AX_DIVISIONS])<=999999,"invalid divisions");
  if(set(AX_MAX_DIGITS)) check(integer(v[AX_MAX_DIGITS])>=1 && integer(v[AX_MAX_DIGITS])<=10,"invalid max digits");
  for(int i:{AX_CENTER_TITLE,AX_OPTIMIZE,AX_NO_EXPONENT,AX_MORE_LOG_LABELS,AX_DECIMALS,AX_CENTER_LABELS,AX_TIME_DISPLAY})
    if(set(i)) bit(v[i]);
  for(int i:{AX_LABEL_COLOR,AX_TITLE_COLOR,AX_AXIS_COLOR}) if(set(i)) color(v[i]);
  if(set(AX_TICKS)) check(v[AX_TICKS]=="plus" || v[AX_TICKS]=="minus" || v[AX_TICKS]=="both","invalid axis ticks");
  if(set(AX_LABELS_OPTION)) {
    const auto &o=v[AX_LABELS_OPTION];
    check(o=="h" || o=="v" || o=="u" || o=="d","invalid labels option");
  }
  if(set(AX_TIME_OFFSET)) real(v[AX_TIME_OFFSET]);
  check(set(AX_EXPONENT_X)==set(AX_EXPONENT_Y),"exponent offset needs two values");
  if(set(AX_EXPONENT_X)) for(int i:{AX_EXPONENT_X,AX_EXPONENT_Y}) check(std::abs(real(v[i]))<=1,"invalid exponent offset");
}

// AXTIME / AXBIN / AXCHANGE rows following an axis row; false when f is not one.
bool axis_extra(AxisSpec &a,const std::vector<std::string> &f,const std::string &name) {
  if(f.size()<2 || f[1]!=name) return false;
  if(f[0]=="AXTIME") {
    check(f.size()==3 && !a.has_time_format && !f[2].empty() && f[2].size()<=256,"invalid AXTIME");
    for(unsigned char c:f[2]) check(c>=0x20 && c<=0x7e,"invalid time format character");
    a.time_format=f[2]; a.has_time_format=true; return true;
  }
  if(f[0]=="AXBIN") {
    check(f.size()==4 && a.bins.size()<10000,"invalid AXBIN");
    int bin=integer(f[2]); check(bin>=1 && bin<=1000000 && !f[3].empty(),"invalid bin label");
    for(const auto &b:a.bins) check(b.first!=bin,"duplicate bin label");
    axis_text(f[3]); a.bins.emplace_back(bin,f[3]); return true;
  }
  if(f[0]=="AXCHANGE") {
    check(f.size()==10 && a.changes.size()<1000,"invalid AXCHANGE");
    LabelChangeSpec c;
    c.index=integer(f[2]); c.angle=real(f[3]); c.size=real(f[4]); c.align=integer(f[5]);
    c.hex=f[6]; c.font=integer(f[7]); c.has_text=bit(f[8]); c.text=f[9];
    check(c.index!=0 && std::abs(c.index)<=1000,"invalid label number");
    check(c.size==-1 || c.size>=0,"invalid label size");
    check(c.align==-1 || (c.align/10>=1 && c.align/10<=3 && c.align%10>=1 && c.align%10<=3),"invalid label alignment");
    check(c.font==-1 || axis_font_code(c.font),"invalid label font");
    if(c.hex!="-") color(c.hex);
    if(c.has_text) { check(!c.text.empty(),"empty label text"); axis_text(c.text); }
    a.changes.push_back(c); return true;
  }
  return false;
}

// Apply the requested attributes with the ROOT setter of each field.
// direction ("x"/"y"/"z") selects TGaxis::SetExponentOffset's axis.
void apply_axis_spec(TAxis *axis,const AxisSpec &a,const char *direction) {
  if(!axis) return;
  const auto &v=a.v;
  auto set=[&](int i){ return v[i]!="-"; };
  if(set(AX_TITLE_FONT)) axis->SetTitleFont(integer(v[AX_TITLE_FONT]));
  if(set(AX_TITLE_SIZE)) axis->SetTitleSize(real(v[AX_TITLE_SIZE]));
  if(set(AX_TITLE_OFFSET)) axis->SetTitleOffset(real(v[AX_TITLE_OFFSET]));
  if(set(AX_LABEL_FONT)) axis->SetLabelFont(integer(v[AX_LABEL_FONT]));
  if(set(AX_LABEL_SIZE)) axis->SetLabelSize(real(v[AX_LABEL_SIZE]));
  if(set(AX_LABEL_OFFSET)) axis->SetLabelOffset(real(v[AX_LABEL_OFFSET]));
  if(set(AX_TICK_LENGTH)) axis->SetTickLength(real(v[AX_TICK_LENGTH]));
  if(set(AX_DIVISIONS) || set(AX_OPTIMIZE)) {
    // TAttAxis::SetNdivisions(n, optim); without explicit divisions keep the current n.
    int n=set(AX_DIVISIONS) ? integer(v[AX_DIVISIONS]) : std::abs(axis->GetNdivisions())%1000000;
    axis->SetNdivisions(n,set(AX_OPTIMIZE) ? bit(v[AX_OPTIMIZE]) : kTRUE);
  }
  if(set(AX_MAX_DIGITS)) axis->SetMaxDigits(integer(v[AX_MAX_DIGITS])); // Must follow SetNdivisions.
  if(set(AX_CENTER_TITLE)) axis->CenterTitle(bit(v[AX_CENTER_TITLE]));
  if(set(AX_NO_EXPONENT)) axis->SetNoExponent(bit(v[AX_NO_EXPONENT]));
  if(set(AX_MORE_LOG_LABELS)) axis->SetMoreLogLabels(bit(v[AX_MORE_LOG_LABELS]));
  if(set(AX_DECIMALS)) axis->SetDecimals(bit(v[AX_DECIMALS]));
  if(set(AX_CENTER_LABELS)) axis->CenterLabels(bit(v[AX_CENTER_LABELS]));
  if(set(AX_LABEL_COLOR)) axis->SetLabelColor(TColor::GetColor(v[AX_LABEL_COLOR].c_str()));
  if(set(AX_TITLE_COLOR)) axis->SetTitleColor(TColor::GetColor(v[AX_TITLE_COLOR].c_str()));
  if(set(AX_AXIS_COLOR)) axis->SetAxisColor(TColor::GetColor(v[AX_AXIS_COLOR].c_str()));
  if(set(AX_TICKS)) axis->SetTicks(v[AX_TICKS]=="plus" ? "+" : v[AX_TICKS]=="minus" ? "-" : "+-");
  if(set(AX_TIME_DISPLAY)) axis->SetTimeDisplay(bit(v[AX_TIME_DISPLAY]));
  if(a.has_time_format) axis->SetTimeFormat(a.time_format.c_str());
  if(set(AX_TIME_OFFSET)) axis->SetTimeOffset(real(v[AX_TIME_OFFSET]),"gmt"); // appends %F<offset> to the format
  if(set(AX_EXPONENT_X) && (std::string(direction)=="x" || std::string(direction)=="y"))
    TGaxis::SetExponentOffset(real(v[AX_EXPONENT_X]),real(v[AX_EXPONENT_Y]),direction);
  for(const auto &b:a.bins) {
    check(b.first<=axis->GetNbins(),"bin label index exceeds the axis bins");
    axis->SetBinLabel(b.first,b.second.c_str());
  }
  if(set(AX_LABELS_OPTION)) {
    check(!a.bins.empty(),"labels option needs bin labels");
    axis->LabelsOption(v[AX_LABELS_OPTION].c_str());
  }
  for(const auto &c:a.changes)
    axis->ChangeLabel(c.index,c.angle,c.size,c.align,c.hex=="-" ? -1 : TColor::GetColor(c.hex.c_str()),
                      c.font,c.has_text ? TString(c.text.c_str()) : TString(""));
}

// A DrawFrame frame has 1000 bins, so a frame whose axis carries bin labels
// (TAxis::SetBinLabel) takes the binning of the histogram it frames; the
// requested visible range is kept.
void frame_bins(TH1 *frame,const TH1 *source,double xmin,double xmax) {
  if(!frame || !source || frame==source) return;
  const TAxis *a=source->GetXaxis();
  if(a->GetXbins()->GetSize()>0) frame->SetBins(a->GetNbins(),a->GetXbins()->GetArray());
  else frame->SetBins(a->GetNbins(),a->GetXmin(),a->GetXmax());
  frame->GetXaxis()->SetRangeUser(xmin,xmax);
}

// CANVAS row fields (see src/scene.jl CanvasStyle): margins, ticks, global error/PS
// settings, grids, then the ROOT title box (gStyle title attributes).
enum CanvasFieldIndex {
  CV_LEFT, CV_RIGHT, CV_BOTTOM, CV_TOP, CV_TICKX, CV_TICKY, CV_FRAME_WIDTH, CV_END_ERROR,
  CV_LINE_SCALE, CV_GRIDX, CV_GRIDY, CV_OPT_TITLE, CV_TITLE_X, CV_TITLE_Y, CV_TITLE_W,
  CV_TITLE_H, CV_TITLE_ALIGN, CV_TITLE_FONT, CV_TITLE_FONT_SIZE, CV_TITLE_BORDER,
  CV_TITLE_FILL, CV_TITLE_STYLE, CV_TITLE_TEXT_COLOR, CV_FIELDS
};
void canvas_title_fields(const std::vector<std::string> &c) {
  check(c.size()==CV_FIELDS,"invalid canvas field count");
  auto set=[&](int i){ return c[i]!="-"; };
  for(int i:{CV_GRIDX,CV_GRIDY,CV_OPT_TITLE}) if(set(i)) bit(c[i]);
  const bool title=set(CV_OPT_TITLE) && bit(c[CV_OPT_TITLE]);
  for(int i=CV_TITLE_X;i<CV_FIELDS;++i) check(!set(i) || title,"title attributes need opt_title");
  for(int i:{CV_TITLE_X,CV_TITLE_Y,CV_TITLE_W,CV_TITLE_H}) if(set(i)) check(real(c[i])>=0 && real(c[i])<=1,"invalid title position/size");
  if(set(CV_TITLE_ALIGN)) { int a=integer(c[CV_TITLE_ALIGN]); check(a/10>=1 && a/10<=3 && a%10>=1 && a%10<=3,"invalid title alignment"); }
  if(set(CV_TITLE_FONT)) check(axis_font_code(integer(c[CV_TITLE_FONT])),"invalid title font");
  if(set(CV_TITLE_FONT_SIZE)) check(real(c[CV_TITLE_FONT_SIZE])>0,"invalid title font size");
  if(set(CV_TITLE_BORDER)) check(integer(c[CV_TITLE_BORDER])>=0 && integer(c[CV_TITLE_BORDER])<=100,"invalid title border");
  if(set(CV_TITLE_STYLE)) { int s=integer(c[CV_TITLE_STYLE]); check(s==0 || s==1001 || (s>=3000 && s<=4100),"invalid title fill style"); }
  for(int i:{CV_TITLE_FILL,CV_TITLE_TEXT_COLOR}) if(set(i)) color(c[i]);
}
void apply_pad_grid(TVirtualPad &pad,const std::vector<std::string> &c) {
  if(c.size()!=CV_FIELDS) return;
  if(c[CV_GRIDX]!="-") pad.SetGridx(bit(c[CV_GRIDX]));
  if(c[CV_GRIDY]!="-") pad.SetGridy(bit(c[CV_GRIDY]));
}

// ROOT's own plot title: THistPainter creates the TPaveText "title" from gStyle
// (OptTitle, TitleX/Y/W/H, TitleAlign, TitleFont, TitleFontSize, TitleBorderSize,
// TitleFillColor, TitleStyle, TitleTextColor) when the pad is painted. gStyle is
// process-global, so this pad's values are set only while ROOT creates the box;
// the frame histogram then gets TH1::kNoTitle so later repaints (with OptTitle 0
// again) keep the box unchanged. The title goes through TH1::SetTitle, so ROOT's
// "title;x title;y title" form also sets the axis titles. Returns false when
// the ROOT title is not requested (the caller then draws its own heading).
bool root_title_on(const std::vector<std::string> &c) { return c.size()==CV_FIELDS && c[CV_OPT_TITLE]=="1"; }
bool root_title(TVirtualPad &pad,TH1 *frame,const std::string &title,const std::vector<std::string> &c) {
  if(!root_title_on(c)) return false;
  if(!frame || title.empty()) return true;
  frame->SetTitle(title.c_str());
  const Int_t opt=gStyle->GetOptTitle(); const Float_t x=gStyle->GetTitleX(),y=gStyle->GetTitleY(),
    w=gStyle->GetTitleW(),h=gStyle->GetTitleH(),font_size=gStyle->GetTitleFontSize();
  const Int_t align=gStyle->GetTitleAlign(); const Style_t font=gStyle->GetTitleFont(""),style=gStyle->GetTitleStyle();
  const Width_t border=gStyle->GetTitleBorderSize(); const Color_t fill=gStyle->GetTitleFillColor(),text=gStyle->GetTitleTextColor();
  auto set=[&](int i){ return c[i]!="-"; };
  if(set(CV_TITLE_X)) gStyle->SetTitleX(real(c[CV_TITLE_X]));
  if(set(CV_TITLE_Y)) gStyle->SetTitleY(real(c[CV_TITLE_Y]));
  if(set(CV_TITLE_W)) gStyle->SetTitleW(real(c[CV_TITLE_W]));
  if(set(CV_TITLE_H)) gStyle->SetTitleH(real(c[CV_TITLE_H]));
  if(set(CV_TITLE_ALIGN)) gStyle->SetTitleAlign(integer(c[CV_TITLE_ALIGN]));
  if(set(CV_TITLE_FONT)) gStyle->SetTitleFont(integer(c[CV_TITLE_FONT]),"t");
  if(set(CV_TITLE_FONT_SIZE)) gStyle->SetTitleFontSize(real(c[CV_TITLE_FONT_SIZE]));
  if(set(CV_TITLE_BORDER)) gStyle->SetTitleBorderSize(integer(c[CV_TITLE_BORDER]));
  if(set(CV_TITLE_FILL)) gStyle->SetTitleFillColor(TColor::GetColor(c[CV_TITLE_FILL].c_str()));
  if(set(CV_TITLE_STYLE)) gStyle->SetTitleStyle(integer(c[CV_TITLE_STYLE]));
  if(set(CV_TITLE_TEXT_COLOR)) gStyle->SetTitleTextColor(TColor::GetColor(c[CV_TITLE_TEXT_COLOR].c_str()));
  gStyle->SetOptTitle(1);
  pad.cd(); pad.Modified(); pad.Update();
  frame->SetBit(TH1::kNoTitle);
  gStyle->SetOptTitle(opt); gStyle->SetTitleX(x); gStyle->SetTitleY(y); gStyle->SetTitleW(w); gStyle->SetTitleH(h);
  gStyle->SetTitleAlign(align); gStyle->SetTitleFont(font,"t"); gStyle->SetTitleFontSize(font_size);
  gStyle->SetTitleBorderSize(border); gStyle->SetTitleFillColor(fill); gStyle->SetTitleStyle(style);
  gStyle->SetTitleTextColor(text);
  pad.Modified();
  return true;
}

// Geometry of a palette's z title as TGaxis paints it, in fractions of the pad
// width: the rotated title is centred unit*offset right of the palette axis and
// is thickness wide; size is the title text size in the same unit. ROOT text
// sizes are fractions of the smaller pad dimension (precision-3 fonts: pixels).
// TLatex measures the title in the current pad (gPad must be pad).
struct PaletteTitle { double unit=0,thickness=0,size=0; };
PaletteTitle palette_title_geometry(TVirtualPad *pad,TAxis *z) {
  PaletteTitle t;
  const double w=pad->GetWw()*pad->GetAbsWNDC(),h=pad->GetWh()*pad->GetAbsHNDC();
  if(w<=0 || h<=0) return t;
  const int font=z->GetTitleFont(); const double size=z->GetTitleSize();
  const bool pixels=font%10==3;
  t.unit=1.6*(pixels ? size/w : size);
  t.size=pixels ? size/w : size*std::min(w,h)/w;
  TLatex probe; probe.SetTextFont(font); probe.SetTextSize(size); probe.SetTitle(z->GetTitle());
  UInt_t bw=0,bh=0; probe.GetBoundingBox(bw,bh);
  t.thickness=bh/w;
  return t;
}

// Test/debug hook: ROOTPLOTS_AXIS_DUMP=<file> appends, just before the image is
// saved, every pad's ROOT state (log/grid/tick flags, margins, place in the
// canvas, title box), the attributes of the pad's frame axes as ROOT reports
// them, and the place of a palette and its z title.
std::string axis_dump_color(int index) {
  TColor *c=gROOT->GetColor(index);
  return c ? std::string(c->AsHexString()) : std::to_string(index);
}
TH1 *axis_dump_frame(TVirtualPad *pad) {
  for(TObject *o:*pad->GetListOfPrimitives()) {
    if(auto *h=dynamic_cast<TH1*>(o)) return h;
    if(auto *s=dynamic_cast<THStack*>(o)) return s->GetHistogram();
    if(auto *f=dynamic_cast<TF1*>(o)) return f->GetHistogram();
    if(auto *g=dynamic_cast<TGraph2D*>(o)) return g->GetHistogram("empty");
  }
  return nullptr;
}
void axis_dump_axis(std::ofstream &out,const char *name,TAxis *a) {
  std::string labels;
  if(a->TestBit(TAxis::kLabelsHori)) labels+="h";
  if(a->TestBit(TAxis::kLabelsVert)) labels+="v";
  if(a->TestBit(TAxis::kLabelsDown)) labels+="d";
  if(a->TestBit(TAxis::kLabelsUp)) labels+="u";
  THashList *bins=a->GetLabels();
  out<<"AXIS\t"<<name<<'\t'<<a->GetTitle()<<'\t'<<a->GetNdivisions()<<'\t'<<a->GetMaxDigits()<<'\t'
     <<a->GetNoExponent()<<'\t'<<a->GetMoreLogLabels()<<'\t'<<a->GetDecimals()<<'\t'<<a->GetCenterTitle()<<'\t'
     <<a->GetCenterLabels()<<'\t'<<(a->TestBit(TAxis::kTickPlus) ? "+" : "")<<(a->TestBit(TAxis::kTickMinus) ? "-" : "")<<'\t'
     <<a->GetTimeDisplay()<<'\t'<<a->GetTimeFormat()<<'\t'<<axis_dump_color(a->GetLabelColor())<<'\t'
     <<axis_dump_color(a->GetTitleColor())<<'\t'<<axis_dump_color(a->GetAxisColor())<<'\t'<<a->GetLabelSize()<<'\t'
     <<a->GetTitleSize()<<'\t'<<a->GetTitleOffset()<<'\t'<<a->GetLabelOffset()<<'\t'<<a->GetTickLength()<<'\t'
     <<a->GetLabelFont()<<'\t'<<a->GetTitleFont()<<'\t'<<(labels.empty() ? "-" : labels)<<'\t'
     <<(bins ? bins->GetSize() : 0)<<'\t'<<a->GetNbins()<<'\t'<<a->GetFirst()<<'\t'<<a->GetLast()<<'\n';
  if(bins) for(TObject *o:*bins) out<<"BINLABEL\t"<<name<<'\t'<<o->GetUniqueID()<<'\t'<<o->GetName()<<'\n';
  if(TList *mods=a->GetModifiedLabels()) for(TObject *o:*mods) {
    auto *m=static_cast<TAxisModLab*>(o);
    out<<"MODLABEL\t"<<name<<'\t'<<m->GetLabNum()<<'\t'<<m->GetAngle()<<'\t'<<m->GetSize()<<'\t'<<m->GetAlign()<<'\t'
       <<(m->GetColor()<0 ? std::string("-") : axis_dump_color(m->GetColor()))<<'\t'<<m->GetFont()<<'\t'<<m->GetText().Data()<<'\n';
  }
}
void axis_dump_pad(std::ofstream &out,TVirtualPad *pad) {
  out<<"PAD\t"<<pad->GetName()<<'\t'<<pad->GetLogx()<<'\t'<<pad->GetLogy()<<'\t'<<pad->GetLogz()<<'\t'
     <<pad->GetGridx()<<'\t'<<pad->GetGridy()<<'\t'<<pad->GetTickx()<<'\t'<<pad->GetTicky()<<'\t'
     <<pad->GetLeftMargin()<<'\t'<<pad->GetRightMargin()<<'\t'<<pad->GetBottomMargin()<<'\t'<<pad->GetTopMargin()<<'\t'
     <<pad->GetAbsXlowNDC()<<'\t'<<pad->GetAbsYlowNDC()<<'\t'<<pad->GetAbsWNDC()<<'\t'<<pad->GetAbsHNDC()<<'\n';
  for(TObject *o:*pad->GetListOfPrimitives()) {
    auto *t=dynamic_cast<TPaveText*>(o);
    if(!t || std::string(t->GetName())!="title") continue;
    TText *line=t->GetLine(0);
    out<<"TITLE\t"<<(line ? line->GetTitle() : "")<<'\t'<<t->GetX1NDC()<<'\t'<<t->GetY1NDC()<<'\t'<<t->GetX2NDC()<<'\t'
       <<t->GetY2NDC()<<'\t'<<axis_dump_color(t->GetFillColor())<<'\t'<<t->GetFillStyle()<<'\t'<<t->GetBorderSize()<<'\t'
       <<axis_dump_color(t->GetTextColor())<<'\t'<<t->GetTextFont()<<'\n';
  }
  // TLatex: text, font, size, position, alignment, NDC bit, and the box width
  // and height (ascent plus descent) as pad NDC fractions, measured by ROOT.
  for(TObject *o:*pad->GetListOfPrimitives()) if(auto *t=dynamic_cast<TLatex*>(o)) {
    TVirtualPad *previous=gPad; pad->cd();
    UInt_t w=0,h=0; t->GetBoundingBox(w,h);
    if(previous) previous->cd();
    const double pw=pad->GetWw()*pad->GetAbsWNDC(),ph=pad->GetWh()*pad->GetAbsHNDC();
    out<<"LATEX\t"<<t->GetTitle()<<'\t'<<t->GetTextFont()<<'\t'<<t->GetTextSize()<<'\t'<<t->GetX()<<'\t'
       <<t->GetY()<<'\t'<<t->GetTextAlign()<<'\t'<<t->TestBit(TLatex::kTextNDC)<<'\t'
       <<(pw>0 ? w/pw : 0)<<'\t'<<(ph>0 ? h/ph : 0)<<'\n';
  }
  if(TH1 *frame=axis_dump_frame(pad)) {
    axis_dump_axis(out,"x",frame->GetXaxis()); axis_dump_axis(out,"y",frame->GetYaxis());
    if(frame->GetDimension()>1 || frame->GetZaxis()->GetTitle()[0]) axis_dump_axis(out,"z",frame->GetZaxis());
    // Palette box (NDC) and the horizontal extent of its z title (equal to the box's right edge without a title).
    TObject *palette=frame->GetListOfFunctions()->FindObject("palette");
    if(palette && palette->InheritsFrom("TPave")) {
      auto *box=static_cast<TPave*>(palette);
      TAxis *z=frame->GetZaxis();
      double left=box->GetX2NDC(),right=left;
      if(z->GetTitle()[0]) {
        TVirtualPad *previous=gPad; pad->cd();
        const PaletteTitle t=palette_title_geometry(pad,z);
        if(previous) previous->cd();
        const double centre=box->GetX2NDC()+t.unit*z->GetTitleOffset();
        left=centre-0.5*t.thickness; right=centre+0.5*t.thickness;
      }
      out<<"PALETTE\t"<<box->GetX1NDC()<<'\t'<<box->GetY1NDC()<<'\t'<<box->GetX2NDC()<<'\t'<<box->GetY2NDC()<<'\t'
         <<left<<'\t'<<right<<'\n';
    }
  }
  for(TObject *o:*pad->GetListOfPrimitives()) if(auto *sub=dynamic_cast<TPad*>(o)) axis_dump_pad(out,sub);
  out<<"ENDPAD\t"<<pad->GetName()<<'\n';
}
void axis_dump(TVirtualPad &canvas) {
  const char *path=std::getenv("ROOTPLOTS_AXIS_DUMP");
  if(!path || !*path) return;
  std::ofstream out(path,std::ios::app);
  check(bool(out),"cannot open ROOTPLOTS_AXIS_DUMP");
  out<<std::setprecision(17)<<"CANVAS\t"<<canvas.GetName()<<'\t'<<gStyle->GetOptTitle()<<'\n';
  out<<"GSTYLE\t"<<gStyle->GetName()<<'\t'<<gStyle->GetPadLeftMargin()<<'\t'<<gStyle->GetPadRightMargin()<<'\t'
     <<gStyle->GetPadBottomMargin()<<'\t'<<gStyle->GetPadTopMargin()<<'\t'<<gStyle->GetLabelFont("X")<<'\t'
     <<gStyle->GetLabelSize("X")<<'\t'<<gStyle->GetTitleFont("X")<<'\t'<<gStyle->GetTitleSize("X")<<'\t'
     <<gStyle->GetPadTickX()<<'\t'<<gStyle->GetPadTickY()<<'\t'<<gStyle->GetOptStat()<<'\n';
  axis_dump_pad(out,&canvas);
  out<<"END\n";
  check(bool(out),"failed writing ROOTPLOTS_AXIS_DUMP");
}
