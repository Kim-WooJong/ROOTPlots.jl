#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TList.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <THStack.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TH2Poly.h>
#include <TArrow.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TGaxis.h>
#include <TGraphPolar.h>
#include <TGraphQQ.h>
#include <TSpline.h>
#include <TPie.h>
#include <TEfficiency.h>
#include <TImage.h>


#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool ok, const std::string &why) {
  if (!ok) throw std::runtime_error(why);
}
std::vector<std::string> fields(const std::string &line) {
  std::vector<std::string> out;
  size_t start=0;
  while (true) {
    auto pos=line.find('\t',start);
    out.push_back(line.substr(start,pos==std::string::npos ? pos : pos-start));
    if (pos==std::string::npos) break;
    start=pos+1;
  }
  return out;
}
double real(const std::string &s, bool finite=true) {
  check(!s.empty(),"empty number");
  char *end=nullptr;
  double v=std::strtod(s.c_str(),&end);
  check(end==s.c_str()+s.size() && (!finite || std::isfinite(v)),"invalid number");
  return v;
}
int integer(const std::string &s) {
  check(!s.empty(),"empty integer");
  size_t end=0;
  int v=std::stoi(s,&end);
  check(end==s.size(),"invalid integer");
  return v;
}
bool bit(const std::string &s) {
  check(s=="0" || s=="1","invalid boolean");
  return s=="1";
}
void color(const std::string &s) {
  check(s.size()==7 && s[0]=='#',"invalid color");
  for(size_t i=1;i<7;++i) check(std::isxdigit(static_cast<unsigned char>(s[i]))!=0,"invalid color");
}
struct Point { double x,y,exl,exh,eyl,eyh; };
struct Series {
  std::string kind,draw,errors,label,hex;
  int marker,line_width,line_style;
  double marker_size;
  std::vector<Point> points;
  std::string fill_hex="-";
  int fill_style=1001;
  double fill_alpha=0.35;
};
bool graph_raw(const Series &s) { return s.draw.rfind("ROOT:",0)==0; }
void check_graph_option(const std::string &option) {
  check(option.size()<=127,"ROOT graph draw option exceeds 127 bytes");
  int depth=0;
  for(unsigned char c:option) {
    check(c>=0x20 && c<=0x7e && c!='\"' && c!='\\' && c!=';',"invalid ROOT graph option");
    check(std::toupper(c)!='A',"graph frame owns axes; A/SAME is not allowed");
    if(c=='(') ++depth;
    if(c==')') --depth;
    check(depth>=0,"unbalanced ROOT graph option");
  }
  check(depth==0,"unbalanced ROOT graph option");
}
bool visible_errors(const Series &s) {
  if(s.kind=="plain") return false;
  if(!graph_raw(s)) return s.errors!="none";
  auto option=s.draw.substr(5);
  std::transform(option.begin(),option.end(),option.begin(),[](unsigned char c){return std::toupper(c);});
  return option.find('X')==std::string::npos;
}
void validate_series_log(const Series &s,bool logx,bool logy) {
  for(const auto &p:s.points) {
    check(!logx || (p.x>0 && (!visible_errors(s) || p.x-p.exl>0)),"log x needs positive visible graph extents");
    check(!logy || (p.y>0 && (!visible_errors(s) || p.y-p.eyl>0)),"log y needs positive visible graph extents");
  }
}
bool filled_graph(const Series &s) {
  if(s.errors=="band" || s.errors=="smoothband" || s.errors=="boxes" ||
     s.draw=="fill" || s.draw=="bar" || s.draw=="bar_offset") return true;
  if(!graph_raw(s)) return false;
  auto option=s.draw.substr(5);
  std::transform(option.begin(),option.end(),option.begin(),[](unsigned char c){return std::toupper(c);});
  return option.find_first_of("FB2345")!=std::string::npos;
}
Series read_series(const std::vector<std::vector<std::string>> &lines, size_t &at, int id) {
  check(at<lines.size(),"missing SERIES");
  const auto &f=lines[at++];
  check((f.size()==12 || f.size()==15) && f[0]=="SERIES" && integer(f[1])==id,"invalid SERIES");
  int n=integer(f[2]);
  check(n>0 && n<=1000000,"invalid point count");
  Series s{f[3],f[4],f[5],f[6],f[7],integer(f[8]),integer(f[10]),integer(f[11]),real(f[9]),{}};
  check(s.kind=="plain" || s.kind=="symmetric" || s.kind=="asymmetric","invalid graph kind");
  if(graph_raw(s)) check_graph_option(s.draw.substr(5));
  else check(s.draw=="points" || s.draw=="line" || s.draw=="linepoints" ||
             s.draw=="curve" || s.draw=="curvepoints" || s.draw=="fill" ||
             s.draw=="bar" || s.draw=="bar_offset","invalid draw mode");
  check(s.errors=="bars" || s.errors=="none" || s.errors=="band" ||
        s.errors=="smoothband" || s.errors=="boxes" || s.errors=="caps" ||
        s.errors=="nocaps" || s.errors=="arrows" || s.errors=="filledarrows","invalid error mode");
  check(s.kind!="plain" || s.errors=="bars" || s.errors=="none","error display needs graph errors");
  if(f.size()==15) {
    s.fill_hex=f[12]; s.fill_style=integer(f[13]); s.fill_alpha=real(f[14]);
    if(s.fill_hex!="-") color(s.fill_hex);
    check(s.fill_style>=0 && s.fill_style<=9999 && s.fill_alpha>=0 && s.fill_alpha<=1,"invalid fill style");
  }
  color(s.hex);
  check(s.marker>=1 && s.marker<=50 && s.marker_size>0 && s.line_width>0 &&
        s.line_style>=1 && s.line_style<=10,"invalid graph style");
  for(int i=1;i<=n;++i) {
    check(at<lines.size(),"missing POINT");
    const auto &p=lines[at++];
    check(p.size()==9 && p[0]=="POINT" && integer(p[1])==id && integer(p[2])==i,
          "invalid POINT id/index");
    Point q{real(p[3]),real(p[4]),real(p[5]),real(p[6]),real(p[7]),real(p[8])};
    check(q.exl>=0 && q.exh>=0 && q.eyl>=0 && q.eyh>=0 &&
          std::isfinite(q.x-q.exl) && std::isfinite(q.x+q.exh) &&
          std::isfinite(q.y-q.eyl) && std::isfinite(q.y+q.eyh),"invalid error extent");
    if(s.kind=="plain")
      check(q.exl==0 && q.exh==0 && q.eyl==0 && q.eyh==0,"plain graph has errors");
    if(s.kind=="symmetric")
      check(q.exl==q.exh && q.eyl==q.eyh,"symmetric graph has asymmetric errors");
    s.points.push_back(q);
  }
  return s;
}
std::unique_ptr<TGraph> make_graph(const Series &s) {
  std::unique_ptr<TGraph> g;
  int n=static_cast<int>(s.points.size());
  if(s.kind=="asymmetric") g=std::make_unique<TGraphAsymmErrors>(n);
  else if(s.kind=="symmetric") g=std::make_unique<TGraphErrors>(n);
  else g=std::make_unique<TGraph>(n);
  for(int i=0;i<n;++i) {
    const auto &p=s.points[i];
    g->SetPoint(i,p.x,p.y);
    if(auto *a=dynamic_cast<TGraphAsymmErrors*>(g.get()))
      a->SetPointError(i,p.exl,p.exh,p.eyl,p.eyh);
    else if(auto *e=dynamic_cast<TGraphErrors*>(g.get()))
      e->SetPointError(i,p.exl,p.eyl);
  }
  g->SetMarkerColor(TColor::GetColor(s.hex.c_str()));
  g->SetLineColor(TColor::GetColor(s.hex.c_str()));
  g->SetMarkerStyle(s.marker); g->SetMarkerSize(s.marker_size);
  g->SetLineWidth(s.line_width); g->SetLineStyle(s.line_style);
  if(filled_graph(s)) {
    g->SetFillColorAlpha(TColor::GetColor((s.fill_hex=="-" ? s.hex : s.fill_hex).c_str()),s.fill_alpha);
    g->SetFillStyle(s.fill_style);
  }
  return g;
}
std::string draw_option(const Series &s) {
  if(graph_raw(s)) return s.draw.substr(5)+" SAME";
  std::string opt=s.draw=="line" ? "L" : s.draw=="linepoints" ? "LP" :
    s.draw=="curve" ? "C" : s.draw=="curvepoints" ? "CP" :
    s.draw=="fill" ? "F" : s.draw=="bar" ? "B" : s.draw=="bar_offset" ? "B1" : "P";
  if(s.kind!="plain") {
    if(s.errors=="bars") opt+="E";
    else if(s.errors=="none") opt+="X";
    else if(s.errors=="band") opt+="3";
    else if(s.errors=="smoothband") opt+="4";
    else if(s.errors=="boxes") opt+="2";
    else if(s.errors=="caps") opt+="[]";
    else if(s.errors=="nocaps") opt+="Z";
    else if(s.errors=="arrows") opt+=">";
    else if(s.errors=="filledarrows") opt+="|>";
  }
  return opt+" SAME";
}
std::string legend_option(const Series &s) {
  if(filled_graph(s)) return "f";
  if(s.draw=="line" || s.draw=="curve") return "l";
  if(s.draw=="linepoints" || s.draw=="curvepoints") return visible_errors(s) ? "lep" : "lp";
  return visible_errors(s) ? "lep" : "p";
}
std::vector<std::vector<std::string>> read_input(const std::string &path) {
  std::ifstream in(path);
  check(in.good(),"cannot open input");
  std::string line;
  check(static_cast<bool>(std::getline(in,line)) && line=="ROOTPlotsGraph\t1","unsupported graph format");
  std::vector<std::vector<std::string>> out;
  while(std::getline(in,line)) {
    check(!line.empty(),"empty record");
    out.push_back(fields(line));
  }
  check(!out.empty(),"empty graph request");
  return out;
}
struct LineSpec {
  bool ndc; double x1,y1,x2,y2; std::string hex,label; int width,style;
};
struct TextSpec {
  bool ndc; double x,y,size,angle; int font,align; std::string hex,text;
};
struct LegendItem { std::string target,label,option; int index; };
struct Scene {
  bool active=false,legend=false,explicit_entries=false;
  std::vector<std::string> canvas=std::vector<std::string>(9,"-");
  std::vector<std::string> axes[2]={std::vector<std::string>(10,"-"),std::vector<std::string>(10,"-")};
  std::vector<LineSpec> lines;
  std::vector<TextSpec> texts;
  std::vector<LegendItem> entries;
  double x1=0,y1=0,x2=0,y2=0,text_size=0.04;
  int font=42,border=0,columns=1;
  std::string background="-";
};
bool font_code(int f) { return f>=10 && f<=159 && f%10<=3; }
void text_field(const std::string &s) {
  check(s.find('\0')==std::string::npos && s.find('\r')==std::string::npos,"invalid text field");
}
bool legend_option(const std::string &s) {
  return s=="l" || s=="p" || s=="pe" || s=="ep" || s=="lp" || s=="lep" || s=="f";
}
Scene read_scene(const std::vector<std::vector<std::string>> &rows,size_t &at) {
  Scene s;
  if(at>=rows.size() || rows[at][0]!="SCENE") return s;
  check(rows[at].size()==2 && rows[at][1]=="1","invalid SCENE version"); ++at;
  s.active=true;
  check(at<rows.size() && rows[at].size()==10 && rows[at][0]=="CANVAS","invalid CANVAS");
  s.canvas.assign(rows[at].begin()+1,rows[at].end()); ++at;
  for(int i=0;i<4;++i) if(s.canvas[i]!="-") {
    double v=real(s.canvas[i]); check(v>=0 && v<1,"invalid canvas margin");
  }
  auto margin=[&](int i){return s.canvas[i]=="-" ? 0.1 : real(s.canvas[i]);};
  check(margin(0)+margin(1)<1 && margin(2)+margin(3)<1,"overlapping canvas margins");
  for(int i=4;i<6;++i) if(s.canvas[i]!="-") bit(s.canvas[i]);
  if(s.canvas[6]!="-") check(integer(s.canvas[6])>=1 && integer(s.canvas[6])<=100,"invalid frame width");
  if(s.canvas[7]!="-") check(real(s.canvas[7])>=0,"invalid error cap size");
  if(s.canvas[8]!="-") check(real(s.canvas[8])>0,"invalid PS line scale");
  for(int a=0;a<2;++a) {
    check(at<rows.size() && rows[at].size()==12 && rows[at][0]=="AXIS" &&
          rows[at][1]==(a==0 ? "x" : "y"),"invalid AXIS");
    auto &v=s.axes[a]; v.assign(rows[at].begin()+2,rows[at].end()); ++at;
    for(int i:{0,3}) if(v[i]!="-") check(font_code(integer(v[i])),"invalid axis font");
    for(int i:{1,2,4}) if(v[i]!="-") check(real(v[i])>0,"invalid axis size/offset");
    if(v[5]!="-") real(v[5]);
    if(v[6]!="-") check(real(v[6])>=0 && real(v[6])<=1,"invalid tick length");
    if(v[7]!="-") check(integer(v[7])>=1 && integer(v[7])<=999999,"invalid divisions");
    if(v[8]!="-") check(integer(v[8])>=1 && integer(v[8])<=10,"invalid max digits");
    if(v[9]!="-") bit(v[9]);
  }
  while(at<rows.size() && rows[at][0]=="LINE") {
    const auto &v=rows[at++];
    check(v.size()==11 && integer(v[1])==static_cast<int>(s.lines.size())+1 && s.lines.size()<1000,"invalid LINE");
    LineSpec l{bit(v[2]),real(v[3]),real(v[4]),real(v[5]),real(v[6]),v[7],v[10],integer(v[8]),integer(v[9])};
    color(l.hex); text_field(l.label);
    check(l.width>=1 && l.width<=100 && l.style>=1 && l.style<=10,"invalid line style");
    s.lines.push_back(l);
  }
  while(at<rows.size() && rows[at][0]=="TEXT") {
    const auto &v=rows[at++];
    check(v.size()==11 && integer(v[1])==static_cast<int>(s.texts.size())+1 && s.texts.size()<1000,"invalid TEXT");
    TextSpec t{bit(v[2]),real(v[3]),real(v[4]),real(v[6]),real(v[9]),integer(v[5]),integer(v[8]),v[7],v[10]};
    color(t.hex); text_field(t.text);
    check(font_code(t.font) && t.size>0 && t.align/10>=1 && t.align/10<=3 &&
          t.align%10>=1 && t.align%10<=3,"invalid text style");
    s.texts.push_back(t);
  }
  if(at<rows.size() && rows[at][0]=="LEGEND") {
    const auto &v=rows[at++];
    check(v.size()==11,"invalid LEGEND");
    s.legend=true; s.x1=real(v[1]); s.y1=real(v[2]); s.x2=real(v[3]); s.y2=real(v[4]);
    s.font=integer(v[5]); s.text_size=real(v[6]); s.border=integer(v[7]);
    s.background=v[8]; s.columns=integer(v[9]); s.explicit_entries=bit(v[10]);
    check(0<=s.x1 && s.x1<s.x2 && s.x2<=1 && 0<=s.y1 && s.y1<s.y2 && s.y2<=1,"invalid legend bounds");
    check(font_code(s.font) && s.text_size>0 && s.border>=0 && s.border<=100 &&
          s.columns>=1 && s.columns<=16,"invalid legend style");
    if(s.background!="-") color(s.background);
    while(at<rows.size() && rows[at][0]=="LENTRY") {
      const auto &e=rows[at++];
      check(e.size()==5 && s.explicit_entries && s.entries.size()<1000,"invalid LENTRY");
      check((e[1]=="series" || e[1]=="line") && legend_option(e[3]),"invalid legend target/option");
      text_field(e[4]); s.entries.push_back({e[1],e[4],e[3],integer(e[2])});
    }
  }
  check(at<rows.size() && rows[at].size()==1 && rows[at][0]=="ENDSCENE","missing ENDSCENE"); ++at;
  return s;
}
void validate_scene(const Scene &s,bool logx,bool logy,size_t nseries) {
  auto coordinate=[&](double x,double y,bool ndc) {
    check(ndc ? x>=0 && x<=1 && y>=0 && y<=1 : (!logx || x>0) && (!logy || y>0),"invalid annotation coordinates");
  };
  for(const auto &l:s.lines) { coordinate(l.x1,l.y1,l.ndc); coordinate(l.x2,l.y2,l.ndc); }
  for(const auto &t:s.texts) coordinate(t.x,t.y,t.ndc);
  for(const auto &e:s.entries)
    check(e.index>=1 && static_cast<size_t>(e.index)<=(e.target=="line" ? s.lines.size() : nseries),"invalid legend reference");
}
void apply_canvas(TPad &canvas,const Scene &s) {
  const auto &v=s.canvas;
  if(v[0]!="-") canvas.SetLeftMargin(real(v[0]));
  if(v[1]!="-") canvas.SetRightMargin(real(v[1]));
  if(v[2]!="-") canvas.SetBottomMargin(real(v[2]));
  if(v[3]!="-") canvas.SetTopMargin(real(v[3]));
  if(v[4]!="-") canvas.SetTickx(bit(v[4]));
  if(v[5]!="-") canvas.SetTicky(bit(v[5]));
  if(v[6]!="-") { canvas.SetFrameLineWidth(integer(v[6])); gStyle->SetFrameLineWidth(integer(v[6])); }
  if(v[7]!="-") gStyle->SetEndErrorSize(real(v[7]));
  if(v[8]!="-") gStyle->SetLineScalePS(real(v[8]));
}
void apply_axis(TAxis *axis,const std::vector<std::string> &v) {
  if(v[0]!="-") axis->SetTitleFont(integer(v[0]));
  if(v[1]!="-") axis->SetTitleSize(real(v[1]));
  if(v[2]!="-") axis->SetTitleOffset(real(v[2]));
  if(v[3]!="-") axis->SetLabelFont(integer(v[3]));
  if(v[4]!="-") axis->SetLabelSize(real(v[4]));
  if(v[5]!="-") axis->SetLabelOffset(real(v[5]));
  if(v[6]!="-") axis->SetTickLength(real(v[6]));
  if(v[7]!="-") axis->SetNdivisions(integer(v[7]));
  if(v[8]!="-") axis->SetMaxDigits(integer(v[8])); // Must follow SetNdivisions.
  if(v[9]!="-") axis->CenterTitle(bit(v[9]));
}
void apply_axes(TH1 *frame,const Scene &s) {
  apply_axis(frame->GetXaxis(),s.axes[0]); apply_axis(frame->GetYaxis(),s.axes[1]);
}
struct SceneObjects {
  std::vector<std::unique_ptr<TLine>> lines;
  std::vector<std::unique_ptr<TLatex>> texts;
  std::unique_ptr<TLegend> legend;
};
SceneObjects draw_scene(const Scene &s,const std::string &corner,
                       const std::vector<TObject*> &objects,
                       const std::vector<std::string> &labels,
                       const std::vector<std::string> &options) {
  SceneObjects out;
  for(const auto &l:s.lines) {
    auto line=std::make_unique<TLine>(l.x1,l.y1,l.x2,l.y2);
    line->SetNDC(l.ndc); line->SetLineColor(TColor::GetColor(l.hex.c_str()));
    line->SetLineWidth(l.width); line->SetLineStyle(l.style); line->Draw();
    out.lines.push_back(std::move(line));
  }
  for(const auto &t:s.texts) {
    auto text=std::make_unique<TLatex>(t.x,t.y,t.text.c_str());
    text->SetNDC(t.ndc); text->SetTextFont(t.font); text->SetTextSize(t.size);
    text->SetTextColor(TColor::GetColor(t.hex.c_str()));
    text->SetTextAlign(t.align); text->SetTextAngle(t.angle); text->Draw();
    out.texts.push_back(std::move(text));
  }
  if(s.legend || corner!="none") {
    double x1=corner.find("left")!=std::string::npos ? 0.13 : 0.68;
    double y1=corner.find("bottom")!=std::string::npos ? 0.16 : 0.73;
    out.legend=s.legend ? std::make_unique<TLegend>(s.x1,s.y1,s.x2,s.y2) :
                         std::make_unique<TLegend>(x1,y1,x1+0.20,y1+0.16);
    out.legend->SetBorderSize(s.legend ? s.border : 0); out.legend->SetFillStyle(0);
    if(s.legend) {
      out.legend->SetTextFont(s.font); out.legend->SetTextSize(s.text_size);
      out.legend->SetNColumns(s.columns);
      if(s.background!="-") { out.legend->SetFillColor(TColor::GetColor(s.background.c_str())); out.legend->SetFillStyle(1001); }
    }
    if(s.explicit_entries) {
      for(const auto &e:s.entries) {
        TObject *object=e.target=="line" ? static_cast<TObject*>(out.lines.at(e.index-1).get()) : objects.at(e.index-1);
        out.legend->AddEntry(object,e.label.c_str(),e.option.c_str());
      }
    } else {
      for(size_t i=0;i<objects.size();++i) if(!labels[i].empty())
        out.legend->AddEntry(objects[i],labels[i].c_str(),options[i].c_str());
      for(size_t i=0;i<s.lines.size();++i) if(!s.lines[i].label.empty())
        out.legend->AddEntry(out.lines[i].get(),s.lines[i].label.c_str(),"l");
    }
    out.legend->Draw();
  }
  return out;
}
void draw(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==10 && f[0]=="DRAW","invalid DRAW");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  bool logx=bit(f[3]),logy=bit(f[4]);
  double xmin=real(f[5]),xmax=real(f[6]),ymin=real(f[7]),ymax=real(f[8]);
  check(xmin<xmax && ymin<ymax && (!logx || xmin>0) && (!logy || ymin>0),"invalid axes");
  std::string corner=f[9];
  check(corner=="topright" || corner=="topleft" || corner=="bottomright" ||
        corner=="bottomleft" || corner=="none","invalid legend");
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="LABELS","invalid LABELS");
  auto title=rows[at][1],xlabel=rows[at][2],ylabel=rows[at][3]; ++at;
  auto scene=read_scene(rows,at);
  std::vector<Series> series;
  while(at<rows.size() && rows[at][0]=="SERIES")
    series.push_back(read_series(rows,at,static_cast<int>(series.size())+1));
  check(!series.empty() && at+1==rows.size() && rows[at].size()==2 &&
        rows[at][0]=="END" && integer(rows[at][1])==static_cast<int>(series.size()),
        "invalid END/count");
  validate_scene(scene,logx,logy,series.size());
  for(const auto &s:series) validate_series_log(s,logx,logy);

  // Graphs, frame and legend remain alive until both SaveAs calls finish.
  TCanvas canvas("rootplots_graph","",width,height);
  apply_canvas(canvas,scene);
  canvas.SetLogx(logx); canvas.SetLogy(logy);
  auto *frame=canvas.DrawFrame(xmin,ymin,xmax,ymax);
  frame->GetXaxis()->SetTitle(xlabel.c_str());
  frame->GetYaxis()->SetTitle(ylabel.c_str());
  apply_axes(frame,scene);
  std::vector<std::unique_ptr<TGraph>> graphs;
  std::vector<TObject*> objects;
  std::vector<std::string> labels,options;
  for(const auto &s:series) {
    graphs.push_back(make_graph(s));
    graphs.back()->Draw(draw_option(s).c_str());
    objects.push_back(graphs.back().get()); labels.push_back(s.label);
    options.push_back(legend_option(s));
  }
  auto decorations=draw_scene(scene,corner,objects,labels,options);
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!title.empty()) heading.DrawLatex(0.13,0.94,title.c_str());
  canvas.Modified(); canvas.Update();
  canvas.SaveAs((stem+".png").c_str());
  canvas.SaveAs((stem+".pdf").c_str());
}
std::string hist_draw_option(const std::string &mode) {
  if(mode=="hist") return "HIST";
  if(mode=="e") return "E";
  if(mode=="e0") return "E0";
  if(mode=="e1") return "E1";
  if(mode=="e2") return "E2";
  if(mode=="e3") return "E3";
  if(mode=="e4") return "E4";
  if(mode=="e5") return "E5";
  if(mode=="e6") return "E6";
  if(mode=="hist_e1") return "HIST E1";
  if(mode=="bar") return "HIST BAR";
  if(mode=="points") return "HIST P";
  if(mode=="points_zero") return "HIST P0";
  if(mode=="line") return "HIST L";
  if(mode=="curve") return "HIST C";
  if(mode=="text") return "HIST TEXT";
  throw std::runtime_error("unsupported native histogram draw option");
}
void check_root_option(const std::string &option, bool allow_same=false) {
  check(option.size()<=127,"ROOT draw option exceeds 127 bytes");
  int depth=0;
  for(unsigned char c:option) {
    check(c>=0x20 && c<=0x7e && c!='"' && c!='\\' && c!=';',
          "invalid ROOT draw option character");
    if(c=='(') ++depth;
    if(c==')') --depth;
    check(depth>=0,"unbalanced ROOT draw option parentheses");
  }
  check(depth==0,"unbalanced ROOT draw option parentheses");
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  for(const auto &part:{"GLLEGO","GLSURF","GLCOL","GLBOX","GLISO"})
    check(upper.find(part)==std::string::npos,
          "OpenGL histogram options cannot be exported by the batch renderer");
  check(allow_same || upper.find("SAME")==std::string::npos,
        "SAME needs an overlay canvas");
}
bool root_option_uses_fill(const std::string &option) {
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  if(upper.find("BAR")!=std::string::npos || upper=="B" || upper=="F" ||
     upper.find("LF")!=std::string::npos) return true;
  for(const auto &part:{"E2","E3","E4","E5","E6"})
    if(upper.find(part)!=std::string::npos) return true;
  return false;
}
std::string formula(const std::string &model);
int nparameters(const std::string &model);
struct FunctionDetach {
  TH1 &hist;
  std::unique_ptr<TF1> &function;
  ~FunctionDetach() {
    if(function) hist.GetListOfFunctions()->Remove(function.get());
  }
};
void draw_hist(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  bool raw=f[0]=="DRAWHISTRAW";
  check((raw && f.size()==14) || (!raw && f.size()==11 && f[0]=="DRAWHIST"),
        "invalid DRAWHIST");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  bool logx=bit(f[3]),logy=bit(f[4]);
  double xmin=real(f[5]),xmax=real(f[6]),ymin=real(f[7]),ymax=real(f[8]);
  bool explicit_x=raw && bit(f[11]),explicit_y=raw && bit(f[12]);
  check(xmin<xmax && ymin<ymax && (!logx || xmin>0) &&
        (!logy || !raw || !explicit_y || ymin>0),"invalid axes");
  std::string corner=f[9],mode=f[10];
  check(corner=="topright" || corner=="topleft" || corner=="bottomright" ||
        corner=="bottomleft" || corner=="none","invalid legend");
  std::string option=raw ? mode : hist_draw_option(mode);
  bool stats=raw && bit(f[13]);
  if(raw) check_root_option(option);
  check(at<rows.size() && rows[at].size()==5 && rows[at][0]=="LABELS","invalid LABELS");
  auto title=rows[at][1],xlabel=rows[at][2],ylabel=rows[at][3],label=rows[at][4]; ++at;
  auto scene=read_scene(rows,at);
  check(at<rows.size() && (rows[at].size()==6 || rows[at].size()==9) && rows[at][0]=="HSTYLE","invalid HSTYLE");
  const auto &s=rows[at++];
  color(s[1]);
  int marker=integer(s[2]),line_width=integer(s[4]),line_style=integer(s[5]);
  double marker_size=real(s[3]);
  check(marker>=1 && marker<=50 && marker_size>0 && line_width>0 &&
        line_style>=1 && line_style<=10,"invalid histogram style");
  check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="HIST","missing HIST");
  int n=integer(rows[at++][1]);
  check(n>=1 && n<=1000000,"invalid histogram bin count");
  std::vector<double> edges;
  edges.reserve(n+1);
  for(int i=1;i<=n+1;++i) {
    check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="EDGE" &&
          integer(rows[at][1])==i,"invalid EDGE");
    double value=real(rows[at++][2]);
    check(edges.empty() || value>edges.back(),"non-increasing EDGE");
    edges.push_back(value);
  }
  check(!logx || edges.front()>0,"log x needs positive bin edges");
  TH1D hist("rootplots_draw_hist","",n,edges.data());
  hist.SetDirectory(nullptr);
  hist.Sumw2();
  for(int i=1;i<=n;++i) {
    check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="BIN" &&
          integer(rows[at][1])==i,"invalid BIN");
    double value=real(rows[at][2]),variance=real(rows[at][3]); ++at;
    check(variance>=0,"negative histogram variance");
    hist.SetBinContent(i,value);
    hist.SetBinError(i,std::sqrt(variance));
  }
  std::unique_ptr<TF1> attached;
  FunctionDetach detach{hist,attached}; // Also detach during parser/scene exceptions.
  if(raw && at<rows.size() && rows[at][0]=="FITFUNC") {
    const auto &function_row=rows[at++];
    check(function_row.size()>=6,"invalid FITFUNC");
    std::string model=function_row[1];
    double fit_min=real(function_row[2]),fit_max=real(function_row[3]);
    int npar=integer(function_row[4]);
    check(fit_min<fit_max && npar==nparameters(model) &&
          function_row.size()==static_cast<size_t>(5+npar),"invalid FITFUNC model/range");
    attached=std::make_unique<TF1>("rootplots_attached_fit",formula(model).c_str(),
                                   fit_min,fit_max);
    check(attached->IsValid() && attached->GetNpar()==npar,"ROOT rejected attached function");
    for(int i=0;i<npar;++i) attached->SetParameter(i,real(function_row[5+i]));
    if(model=="erfc_edge") check(attached->GetParameter(3)>0,"edge width must be positive");
    if(model=="double_erf") check(attached->GetParameter(4)>0 && attached->GetParameter(5)>0,"edge widths must be positive");
    hist.GetListOfFunctions()->Add(attached.get());
  }
  if(raw) {
    std::string upper=option;
    std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
    check(upper.find("FUNC")==std::string::npos || attached,
          "FUNC needs an attached ROOT function");
  }
  std::vector<Series> overlays;
  while(at<rows.size() && rows[at][0]=="SERIES")
    overlays.push_back(read_series(rows,at,static_cast<int>(overlays.size())+1));
  check(overlays.size()<=64 && (raw || overlays.empty()),"invalid histogram overlays");
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        integer(rows[at][1])==static_cast<int>(std::max(size_t(1),overlays.size())),"invalid histogram END/count");
  validate_scene(scene,logx,logy,1+overlays.size());
  hist.SetStats(stats);
  int root_color=TColor::GetColor(s[1].c_str());
  hist.SetLineColor(root_color); hist.SetLineWidth(line_width); hist.SetLineStyle(line_style);
  hist.SetMarkerColor(root_color); hist.SetMarkerStyle(marker); hist.SetMarkerSize(marker_size);
  bool filled=(mode=="bar" || mode=="e2" || mode=="e3" || mode=="e4" ||
               mode=="e5" || mode=="e6");
  if(mode=="bar") {
    hist.SetFillColor(root_color);
    hist.SetBarWidth(0.8);
    hist.SetBarOffset(0.1);
  }
  else if(filled) hist.SetFillColorAlpha(root_color,0.35);
  if(filled) hist.SetFillStyle(1001);
  if(raw && root_option_uses_fill(option)) {
    hist.SetFillColor(root_color); hist.SetFillStyle(1001);
  }
  if(s.size()==9) {
    color(s[6]); int fill=integer(s[7]); double alpha=real(s[8]);
    check(fill>=0 && fill<=9999 && alpha>=0 && alpha<=1,"invalid histogram fill style");
    hist.SetFillColorAlpha(TColor::GetColor(s[6].c_str()),alpha); hist.SetFillStyle(fill);
  }

  TCanvas canvas("rootplots_hist_canvas","",width,height);
  apply_canvas(canvas,scene);
  // SetBinContent does not preserve the original event count; omit Entries.
  if(stats) gStyle->SetOptStat(1100);
  canvas.SetLogx(logx); canvas.SetLogy(logy);
  if(raw) {
    hist.GetXaxis()->SetTitle(xlabel.c_str());
    hist.GetYaxis()->SetTitle(ylabel.c_str());
    if(explicit_x) hist.GetXaxis()->SetRangeUser(xmin,xmax);
    if(!explicit_y && !overlays.empty()) {
      double lo=INFINITY,hi=-INFINITY;
      auto extent=[&](double y,double low,double high) {
        if(logy) {
          for(double v:{y,low,high}) if(v>0) lo=std::min(lo,v);
        } else lo=std::min(lo,low);
        hi=std::max(hi,high);
      };
      for(int i=1;i<=n;++i) extent(hist.GetBinContent(i),hist.GetBinContent(i)-hist.GetBinError(i),hist.GetBinContent(i)+hist.GetBinError(i));
      for(const auto &s:overlays) for(const auto &p:s.points)
        extent(p.y,p.y-(visible_errors(s) ? p.eyl : 0),p.y+(visible_errors(s) ? p.eyh : 0));
      check(std::isfinite(lo) && std::isfinite(hi) && (!logy || hi>0),"invalid overlay y range");
      if(logy) { ymin=lo/2; ymax=hi*1.5; }
      else { lo=std::min(0.0,lo); hi=std::max(0.0,hi); double pad=hi==lo ? 0.5 : 0.05*(hi-lo); ymin=lo-pad; ymax=hi+pad; }
      check(std::isfinite(ymin) && std::isfinite(ymax) && ymin<ymax,"unrepresentable overlay y range");
    }
    if(explicit_y || !overlays.empty()) { hist.SetMinimum(ymin); hist.SetMaximum(ymax); }
    apply_axes(&hist,scene);
    hist.Draw(option.c_str());
  } else {
    auto *frame=canvas.DrawFrame(xmin,ymin,xmax,ymax);
    frame->GetXaxis()->SetTitle(xlabel.c_str());
    frame->GetYaxis()->SetTitle(ylabel.c_str());
    apply_axes(frame,scene);
    hist.Draw((option+" SAME").c_str());
  }
  std::vector<std::unique_ptr<TGraph>> graphs;
  std::vector<TObject*> objects{&hist};
  std::vector<std::string> labels{label},options{filled ? "f" : mode=="hist" || mode=="line" || mode=="curve" ? "l" : "lep"};
  for(const auto &s:overlays) {
    graphs.push_back(make_graph(s)); graphs.back()->Draw(draw_option(s).c_str());
    objects.push_back(graphs.back().get()); labels.push_back(s.label);
    options.push_back(legend_option(s));
  }
  auto decorations=draw_scene(scene,corner,objects,labels,options);
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!title.empty()) heading.DrawLatex(0.13,0.94,title.c_str());
  canvas.Modified(); canvas.Update();
  canvas.SaveAs((stem+".png").c_str());
  canvas.SaveAs((stem+".pdf").c_str());
}
void draw_hist_overlay(const std::vector<std::vector<std::string>> &rows,
                       const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==14 && f[0]=="DRAWHISTOVERLAY","invalid DRAWHISTOVERLAY");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  bool logx=bit(f[3]),logy=bit(f[4]);
  double xmin=real(f[5]),xmax=real(f[6]),ymin=real(f[7]),ymax=real(f[8]);
  check(xmin<xmax && ymin<ymax && (!logx || xmin>0),"invalid overlay limits");
  std::string corner=f[9];
  check(corner=="topright" || corner=="topleft" || corner=="bottomright" ||
        corner=="bottomleft" || corner=="none","invalid overlay legend");
  bool explicit_x=bit(f[10]),explicit_y=bit(f[11]);
  check(!explicit_y || !logy || ymin>0,"invalid log y limit");
  bool stats=bit(f[12]);
  int nlayer=integer(f[13]);
  check(nlayer>=2 && nlayer<=64,"invalid overlay layer count");
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="LABELS","invalid LABELS");
  auto title=rows[at][1],xlabel=rows[at][2],ylabel=rows[at][3]; ++at;
  auto scene=read_scene(rows,at);
  std::vector<std::unique_ptr<TH1D>> hists;
  std::vector<std::string> options,labels;
  for(int k=1;k<=nlayer;++k) {
    check(at<rows.size() && (rows[at].size()==9 || rows[at].size()==12) && rows[at][0]=="HLAYER" &&
          integer(rows[at][1])==k,"invalid HLAYER");
    const auto &s=rows[at++];
    check_root_option(s[2],k>1);
    std::string upper=s[2];
    std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
    check(k==1 || upper.find("SAME")!=std::string::npos,"overlay layer needs SAME");
    color(s[4]);
    int marker=integer(s[5]),line_width=integer(s[7]),line_style=integer(s[8]);
    double marker_size=real(s[6]);
    check(marker>=1 && marker<=50 && marker_size>0 && line_width>0 &&
          line_style>=1 && line_style<=10,"invalid overlay style");
    check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="HIST","missing HIST");
    int n=integer(rows[at++][1]);
    check(n>=1 && n<=1000000,"invalid overlay bin count");
    std::vector<double> edges;
    edges.reserve(n+1);
    for(int i=1;i<=n+1;++i) {
      check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="EDGE" &&
            integer(rows[at][1])==i,"invalid overlay EDGE");
      double edge=real(rows[at++][2]);
      check(edges.empty() || edge>edges.back(),"non-increasing overlay EDGE");
      edges.push_back(edge);
    }
    check(!logx || edges.front()>0,"log x needs positive overlay edges");
    auto hist=std::make_unique<TH1D>(("rootplots_overlay_"+std::to_string(k)).c_str(),"",n,edges.data());
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    for(int i=1;i<=n;++i) {
      check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="BIN" &&
            integer(rows[at][1])==i,"invalid overlay BIN");
      double content=real(rows[at][2]),variance=real(rows[at][3]); ++at;
      check(variance>=0,"negative overlay variance");
      hist->SetBinContent(i,content);
      hist->SetBinError(i,std::sqrt(variance));
    }
    hist->SetStats(stats);
    int root_color=TColor::GetColor(s[4].c_str());
    hist->SetLineColor(root_color); hist->SetLineWidth(line_width); hist->SetLineStyle(line_style);
    hist->SetMarkerColor(root_color); hist->SetMarkerStyle(marker); hist->SetMarkerSize(marker_size);
    if(root_option_uses_fill(s[2])) {
      hist->SetFillColor(root_color); hist->SetFillStyle(1001);
    }
    if(s.size()==12) {
      color(s[9]); int fill=integer(s[10]); double alpha=real(s[11]);
      check(fill>=0 && fill<=9999 && alpha>=0 && alpha<=1,"invalid overlay fill style");
      hist->SetFillColorAlpha(TColor::GetColor(s[9].c_str()),alpha); hist->SetFillStyle(fill);
    }
    options.push_back(s[2]); labels.push_back(s[3]); hists.push_back(std::move(hist));
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        integer(rows[at][1])==nlayer,"invalid overlay END/count");
  validate_scene(scene,logx,logy,hists.size());
  TCanvas canvas("rootplots_overlay_canvas","",width,height);
  apply_canvas(canvas,scene);
  // SetBinContent does not preserve the original event count; omit Entries.
  if(stats) gStyle->SetOptStat(1100);
  canvas.SetLogx(logx); canvas.SetLogy(logy);
  hists.front()->GetXaxis()->SetTitle(xlabel.c_str());
  hists.front()->GetYaxis()->SetTitle(ylabel.c_str());
  apply_axes(hists.front().get(),scene);
  if(explicit_x) hists.front()->GetXaxis()->SetRangeUser(xmin,xmax);
  if(!explicit_y) {
    double lo=INFINITY,hi=-INFINITY;
    for(const auto &hist:hists) for(int i=1;i<=hist->GetNbinsX();++i) {
      double y=hist->GetBinContent(i),error=hist->GetBinError(i);
      if(logy) {
        if(y-error>0) lo=std::min(lo,y-error);
        if(y>0) lo=std::min(lo,y);
        if(y+error>0) lo=std::min(lo,y+error);
      } else lo=std::min(lo,y-error);
      hi=std::max(hi,y+error);
    }
    if(logy) {
      check(std::isfinite(lo) && hi>0,"automatic overlay log y needs positive data");
      ymin=lo/2; ymax=hi*1.5;
    } else {
      lo=std::min(lo,0.0); hi=std::max(hi,0.0);
      double pad=hi==lo ? 0.5 : 0.05*(hi-lo);
      ymin=lo-pad; ymax=hi+pad;
    }
    check(std::isfinite(ymin) && std::isfinite(ymax) && ymin<ymax,
          "invalid automatic overlay y range");
  }
  hists.front()->SetMinimum(ymin); hists.front()->SetMaximum(ymax);
  for(int k=0;k<nlayer;++k) hists[k]->Draw(options[k].c_str());
  std::vector<TObject*> objects;
  for(const auto &hist:hists) objects.push_back(hist.get());
  auto decorations=draw_scene(scene,corner,objects,labels,std::vector<std::string>(hists.size(),"lep"));
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!title.empty()) heading.DrawLatex(0.13,0.94,title.c_str());
  canvas.Modified(); canvas.Update();
  canvas.SaveAs((stem+".png").c_str());
  canvas.SaveAs((stem+".pdf").c_str());
}
void draw_fit(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==12 && f[0]=="DRAWFIT","invalid DRAWFIT");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  bool logx=bit(f[3]),logy=bit(f[4]);
  double xmin=real(f[5]),xmax=real(f[6]),ymin=real(f[7]),ymax=real(f[8]);
  double pullmin=real(f[9]),pullmax=real(f[10]);
  check(xmin<xmax && ymin<ymax && pullmin<0 && pullmax>0 &&
        (!logx || xmin>0) && (!logy || ymin>0),"invalid fit plot axes");
  std::string corner=f[11];
  check(corner=="topright" || corner=="topleft" || corner=="bottomright" ||
        corner=="bottomleft" || corner=="none","invalid legend");
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="LABELS","invalid LABELS");
  auto title=rows[at][1],xlabel=rows[at][2],ylabel=rows[at][3]; ++at;
  std::vector<Series> series;
  for(int id=1;id<=3;++id) series.push_back(read_series(rows,at,id));
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        integer(rows[at][1])==3,"invalid fit plot END/count");
  check(series[0].draw=="line" && series[1].draw=="points" &&
        series[2].draw=="points" && series[2].kind=="plain" &&
        series[2].errors=="none","invalid fit plot series modes");

  TCanvas canvas("rootplots_fit_plot","",width,height);
  auto *top=new TPad("fit_top","",0,0.30,1,1);
  auto *bottom=new TPad("fit_bottom","",0,0,1,0.30);
  top->Draw(); bottom->Draw();
  for(auto *pad:{top,bottom}) {
    pad->SetLeftMargin(0.14); pad->SetRightMargin(0.05);
    pad->SetTicks(1,1); pad->SetLogx(logx);
  }
  top->SetLogy(logy);
  top->SetTopMargin(0.09); top->SetBottomMargin(0.02);
  bottom->SetTopMargin(0.03); bottom->SetBottomMargin(0.34);
  std::vector<std::unique_ptr<TGraph>> graphs;
  for(const auto &s:series) graphs.push_back(make_graph(s));

  top->cd();
  auto *upper=top->DrawFrame(xmin,ymin,xmax,ymax);
  upper->GetXaxis()->SetLabelSize(0); upper->GetXaxis()->SetTitleSize(0);
  upper->GetYaxis()->SetTitle(ylabel.c_str());
  upper->GetYaxis()->SetTitleSize(0.052); upper->GetYaxis()->SetLabelSize(0.043);
  upper->GetYaxis()->SetTitleOffset(1.15);
  graphs[0]->Draw(draw_option(series[0]).c_str());
  graphs[1]->Draw(draw_option(series[1]).c_str());
  std::unique_ptr<TLegend> legend;
  if(corner!="none") {
    double x1=corner.find("left")!=std::string::npos ? 0.17 : 0.67;
    double y1=corner.find("bottom")!=std::string::npos ? 0.17 : 0.72;
    legend=std::make_unique<TLegend>(x1,y1,x1+0.25,y1+0.16);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    if(!series[0].label.empty()) legend->AddEntry(graphs[0].get(),series[0].label.c_str(),"l");
    if(!series[1].label.empty()) legend->AddEntry(graphs[1].get(),series[1].label.c_str(),"pe");
    legend->Draw();
  }
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!title.empty()) heading.DrawLatex(0.14,0.95,title.c_str());

  bottom->cd();
  auto *lower=bottom->DrawFrame(xmin,pullmin,xmax,pullmax);
  lower->GetXaxis()->SetTitle(xlabel.c_str());
  lower->GetXaxis()->SetTitleSize(0.12); lower->GetXaxis()->SetLabelSize(0.10);
  lower->GetXaxis()->SetTitleOffset(1.05);
  lower->GetYaxis()->SetTitle("Pull");
  lower->GetYaxis()->SetTitleSize(0.11); lower->GetYaxis()->SetLabelSize(0.09);
  lower->GetYaxis()->SetTitleOffset(0.51); lower->GetYaxis()->SetNdivisions(505);
  auto zero=std::make_unique<TLine>(xmin,0,xmax,0);
  zero->SetLineStyle(2); zero->SetLineColor(kGray+2); zero->Draw("SAME");
  graphs[2]->Draw(draw_option(series[2]).c_str());

  canvas.cd(); canvas.Modified(); canvas.Update();
  canvas.SaveAs((stem+".png").c_str());
  canvas.SaveAs((stem+".pdf").c_str());
}
std::string formula(const std::string &model) {
  if(model=="pol0" || model=="pol1" || model=="pol2" || model=="gaus") return model;
  if(model=="exp") return "[0]*exp([1]*x)";
  if(model=="erfc_edge") return "[0]+0.5*[1]*TMath::Erfc((x-[2])/(TMath::Sqrt(2.0)*[3]))";
  if(model=="double_erf") return "[0]+0.5*[1]*(TMath::Erf((x-[2])/(TMath::Sqrt(2.0)*[4]))-TMath::Erf((x-[3])/(TMath::Sqrt(2.0)*[5])))";
  throw std::runtime_error("invalid fit model");
}
int nparameters(const std::string &model) {
  if(model=="pol0") return 1;
  if(model=="pol1" || model=="exp") return 2;
  if(model=="pol2" || model=="gaus") return 3;
  if(model=="erfc_edge") return 4;
  if(model=="double_erf") return 6;
  throw std::runtime_error("invalid fit model");
}
void fit(const std::vector<std::vector<std::string>> &rows, const std::string &outpath) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==8 && f[0]=="FIT","invalid FIT");
  std::string type=f[1],model=f[2],objective=f[3],xerrors=f[4];
  check(type=="graph" || type=="hist","invalid fit data type");
  check(objective=="chi2" || objective=="poisson" || objective=="weighted_poisson",
        "invalid objective");
  check(type=="hist" || objective=="chi2","graph likelihood fit unsupported");
  check(xerrors=="effective" || xerrors=="ignore","invalid xerror mode");
  double xmin=real(f[5]),xmax=real(f[6]);
  int npar=integer(f[7]);
  check(xmin<xmax && npar==nparameters(model),"invalid fit range/model parameters");
  TF1 function("rootplots_fit",formula(model).c_str(),xmin,xmax);
  check(function.IsValid() && function.GetNpar()==npar,"ROOT rejected fit model");
  for(int i=1;i<=npar;++i) {
    check(at<rows.size(),"missing PARAM");
    const auto &p=rows[at++];
    check(p.size()==6 && p[0]=="PARAM" && integer(p[1])==i,"invalid PARAM");
    double value=real(p[2]),lo=real(p[3],false),hi=real(p[4],false);
    bool fixed=bit(p[5]);
    check(lo<hi && lo<=value && value<=hi,"invalid parameter bound");
    bool width=(model=="erfc_edge" && i==4) || (model=="double_erf" && (i==5 || i==6));
    check(!width || (value>0 && (fixed || (std::isfinite(lo) && lo>0 && std::isfinite(hi)))),
          "edge widths need positive initial values and finite positive bounds, or must be fixed");
    function.SetParameter(i-1,value);
    if(fixed) function.FixParameter(i-1,value);
    else if(std::isfinite(lo) && std::isfinite(hi)) function.SetParLimits(i-1,lo,hi);
    else check(std::isinf(lo) && lo<0 && std::isinf(hi) && hi>0,"partial bound unsupported");
  }
  std::unique_ptr<TGraph> graph;
  std::unique_ptr<TH1D> hist;
  if(type=="graph") {
    auto s=read_series(rows,at,1);
    if(s.kind!="plain")
      for(const auto &p:s.points) check(p.eyl>0 && p.eyh>0,"fit requires positive y errors");
    graph=make_graph(s);
  } else {
    check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="HIST","missing HIST");
    int n=integer(rows[at++][1]);
    check(n>0 && n<=1000000,"invalid bin count");
    std::vector<double> edges;
    std::vector<double> values,variances;
    for(int i=1;i<=n;++i) {
      check(at<rows.size(),"missing BIN");
      const auto &b=rows[at++];
      check(b.size()==6 && b[0]=="BIN" && integer(b[1])==i,"invalid BIN");
      double left=real(b[2]),right=real(b[3]),value=real(b[4]),variance=real(b[5]);
      check(left<right && (edges.empty() || left==edges.back()) && variance>=0,
            "invalid bin geometry/variance");
      if(objective!="chi2") check(value>=0,"likelihood fit forbids negative bins");
      if(objective=="weighted_poisson")
        check(value==0 || variance>0,"positive weighted bin requires positive variance");
      if(objective=="poisson")
        check(value==std::round(value) && std::abs(variance-value)<=1e-10*std::max(1.0,value),
              "Poisson fit requires unweighted counts");
      if(edges.empty()) edges.push_back(left);
      edges.push_back(right); values.push_back(value); variances.push_back(variance);
    }
    hist=std::make_unique<TH1D>("rootplots_hist","",n,edges.data());
    hist->SetDirectory(nullptr); hist->Sumw2(kTRUE);
    for(int i=1;i<=n;++i) {
      hist->SetBinContent(i,values[i-1]);
      hist->SetBinError(i,std::sqrt(variances[i-1]));
    }
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        rows[at][1]=="1","invalid END/trailing records");
  std::string options="SQN0RFB"; // F uses the general fitter; B preserves predefined gaus parameters/bounds.
  if(objective=="poisson") options+="L";
  if(objective=="weighted_poisson") options+="WL";
  if(type=="graph" && xerrors=="ignore") options+="EX0";
  TFitResultPtr result=type=="graph" ? graph->Fit(&function,options.c_str()) :
                                     hist->Fit(&function,options.c_str());
  check(result.Get()!=nullptr,"ROOT fit did not produce a result");
  std::ofstream out(outpath);
  check(out.good(),"cannot write fit result");
  out<<std::setprecision(17);
  out<<"ROOTPlotsFitResult\t1\n";
  out<<"FIT\t"<<result->Status()<<"\t"<<int(result->IsValid())<<"\t"
     <<result->CovMatrixStatus()<<"\t"<<result->Chi2()<<"\t"<<result->Ndf()
     <<"\t"<<function.GetNumberFitPoints()<<"\t"<<result->Edm()<<"\t"
     <<result->MinimizerType()<<"\t"<<int(result->NormalizedErrors())<<"\n";
  for(int i=0;i<npar;++i)
    out<<"PARAM\t"<<i+1<<"\t"<<result->Parameter(i)<<"\t"<<result->ParError(i)
       <<"\t"<<int(result->IsParameterFixed(i))<<"\n";
  for(int i=0;i<npar;++i) {
    out<<"COV\t"<<i+1;
    for(int j=0;j<npar;++j) out<<"\t"<<result->CovMatrix(i,j);
    out<<"\n";
  }
  out<<"END\n";
  out.flush();
  check(out.good(),"incomplete fit result write");
}
struct Labels { std::string title,x,y,z; };
Labels read_labels(const std::vector<std::vector<std::string>> &rows, size_t &at) {
  check(at<rows.size() && rows[at].size()==5 && rows[at][0]=="LABELS","invalid LABELS");
  const auto &f=rows[at++];
  for(size_t i=1;i<f.size();++i)
    check(f[i].find('\0')==std::string::npos && f[i].find('\r')==std::string::npos,
          "invalid label control character");
  return {f[1],f[2],f[3],f[4]};
}
void palette(const std::string &name) {
  check(name=="bird" || name=="viridis","invalid palette");
  gStyle->SetPalette(name=="bird" ? kBird : kViridis);
  gStyle->SetNumberContours(255);
}
void save(TCanvas &canvas, const std::string &stem) {
  canvas.Modified(); canvas.Update();
  canvas.SaveAs((stem+".png").c_str());
  canvas.SaveAs((stem+".pdf").c_str());
}
void draw_grid(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  bool raw=f[0]=="DRAWGRIDRAW";
  check((raw && f.size()==12) || (!raw && f.size()==9 && f[0]=="DRAWGRID"),
        "invalid DRAWGRID");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  std::string dim=f[3], mode=f[4];
  check(raw ? (dim=="2d" || dim=="3d") :
        ((dim=="2d" && (mode=="heatmap" || mode=="contour")) ||
         (dim=="3d" && mode=="surface")),"invalid grid mode");
  if(raw) check_root_option(mode);
  palette(f[5]);
  int contours=integer(f[6]);
  check(contours>=1 && contours<=255,"invalid contour count");
  double theta=real(f[7]),phi=real(f[8]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid view angles");
  bool logx=raw && bit(f[9]),logy=raw && bit(f[10]),logz=raw && bit(f[11]);
  auto labels=read_labels(rows,at);
  check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="GRID","missing GRID");
  int nx=integer(rows[at][1]),ny=integer(rows[at][2]); ++at;
  check(nx>=1 && ny>=1 && nx<=2048 && ny<=2048 &&
        static_cast<long long>(nx)*ny<=1000000,"invalid grid dimensions");
  if(dim=="3d" && !raw) check(nx>=2 && ny>=2,"3-D surface needs at least two cells on each axis");
  std::vector<double> xedges,yedges;
  for(int i=1;i<=nx+1;++i) {
    check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="XEDGE" &&
          integer(rows[at][1])==i,"invalid XEDGE");
    double value=real(rows[at++][2]);
    check(xedges.empty() || value>xedges.back(),"non-increasing XEDGE");
    xedges.push_back(value);
  }
  for(int i=1;i<=ny+1;++i) {
    check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="YEDGE" &&
          integer(rows[at][1])==i,"invalid YEDGE");
    double value=real(rows[at++][2]);
    check(yedges.empty() || value>yedges.back(),"non-increasing YEDGE");
    yedges.push_back(value);
  }
  TH2D hist("rootplots_grid","",nx,xedges.data(),ny,yedges.data());
  hist.SetDirectory(nullptr);
  double low=INFINITY,high=-INFINITY;
  for(int iy=1;iy<=ny;++iy) for(int ix=1;ix<=nx;++ix) {
    check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="CELL" &&
          integer(rows[at][1])==ix && integer(rows[at][2])==iy,"invalid CELL");
    double value=real(rows[at++][3]);
    hist.SetBinContent(ix,iy,value);
    low=std::min(low,value); high=std::max(high,value);
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        rows[at][1]=="1","invalid END/trailing records");
  hist.GetXaxis()->SetTitle(labels.x.c_str());
  hist.GetYaxis()->SetTitle(labels.y.c_str());
  hist.GetZaxis()->SetTitle(labels.z.c_str());
  // ROOT suppresses zero-valued cells/surfaces when the display minimum is
  // nonnegative. Keep the bin contents unchanged and widen only the palette.
  if(low==0 && high==0) {
    hist.SetMinimum(-1.0);
    hist.SetMaximum(1.0);
  } else if(dim=="2d" && mode=="heatmap" && low==0) {
    double epsilon=std::max(high*1e-12,1e-300);
    hist.SetMinimum(-epsilon);
  }
  TCanvas canvas("rootplots_grid_canvas","",width,height);
  canvas.SetLogx(logx); canvas.SetLogy(logy); canvas.SetLogz(logz);
  canvas.SetLeftMargin(0.14); canvas.SetBottomMargin(0.13);
  if(dim=="2d") canvas.SetRightMargin(raw || mode=="heatmap" ? 0.18 : 0.08);
  else { canvas.SetRightMargin(0.22); canvas.SetTheta(theta); canvas.SetPhi(phi); }
  std::string upper=mode;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  hist.SetContour(raw ? (upper.find("CONT")!=std::string::npos ? contours :
                         upper.find("SURF")!=std::string::npos ? 80 : 255) :
                  mode=="contour" ? contours : dim=="3d" ? 80 : 255);
  hist.Draw(raw ? mode.c_str() : mode=="heatmap" ? "COLZ" :
            mode=="contour" ? "CONT3" : "SURF2Z");
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!labels.title.empty()) heading.DrawLatex(0.13,0.94,labels.title.c_str());
  save(canvas,stem);
}
void draw_grid_overlay(const std::vector<std::vector<std::string>> &rows,
                       const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==9 && f[0]=="DRAWGRIDOVERLAY","invalid DRAWGRIDOVERLAY");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  palette(f[3]);
  int contours=integer(f[4]);
  check(contours>=1 && contours<=255,"invalid contour count");
  bool logx=bit(f[5]),logy=bit(f[6]),logz=bit(f[7]);
  int nlayer=integer(f[8]);
  check(nlayer>=2 && nlayer<=32,"invalid grid overlay layer count");
  auto labels=read_labels(rows,at);
  std::vector<std::unique_ptr<TH2D>> grids;
  std::vector<std::string> options;
  for(int k=1;k<=nlayer;++k) {
    check(at<rows.size() && rows[at].size()==8 && rows[at][0]=="GLAYER" &&
          integer(rows[at][1])==k,"invalid GLAYER");
    const auto &s=rows[at++];
    check_root_option(s[2],k>1);
    std::string upper=s[2];
    std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
    check(k==1 || upper.find("SAME")!=std::string::npos,"grid overlay layer needs SAME");
    color(s[3]);
    int marker=integer(s[4]),line_width=integer(s[6]),line_style=integer(s[7]);
    double marker_size=real(s[5]);
    check(marker>=1 && marker<=50 && marker_size>0 && line_width>0 &&
          line_style>=1 && line_style<=10,"invalid grid overlay style");
    check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="GRID","missing GRID");
    int nx=integer(rows[at][1]),ny=integer(rows[at][2]); ++at;
    check(nx>=1 && ny>=1 && nx<=2048 && ny<=2048 &&
          static_cast<long long>(nx)*ny<=1000000,"invalid grid dimensions");
    auto read_edges=[&](const std::string &tag,int n) {
      std::vector<double> edges;
      edges.reserve(n+1);
      for(int i=1;i<=n+1;++i) {
        check(at<rows.size() && rows[at].size()==3 && rows[at][0]==tag &&
              integer(rows[at][1])==i,"invalid grid overlay edge");
        double edge=real(rows[at++][2]);
        check(edges.empty() || edge>edges.back(),"non-increasing grid overlay edge");
        edges.push_back(edge);
      }
      return edges;
    };
    auto xedges=read_edges("XEDGE",nx);
    auto yedges=read_edges("YEDGE",ny);
    check((!logx || xedges.front()>0) && (!logy || yedges.front()>0),
          "log axes need positive overlay edges");
    auto hist=std::make_unique<TH2D>(("rootplots_grid_overlay_"+std::to_string(k)).c_str(),
                                     "",nx,xedges.data(),ny,yedges.data());
    hist->SetDirectory(nullptr);
    for(int iy=1;iy<=ny;++iy) for(int ix=1;ix<=nx;++ix) {
      check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="CELL" &&
            integer(rows[at][1])==ix && integer(rows[at][2])==iy,
            "invalid grid overlay cell");
      hist->SetBinContent(ix,iy,real(rows[at++][3]));
    }
    int root_color=TColor::GetColor(s[3].c_str());
    hist->SetLineColor(root_color); hist->SetLineWidth(line_width);
    hist->SetLineStyle(line_style); hist->SetMarkerColor(root_color);
    hist->SetMarkerStyle(marker); hist->SetMarkerSize(marker_size);
    hist->SetStats(false);
    hist->SetContour(upper.find("CONT")!=std::string::npos ? contours :
                     upper.find("SURF")!=std::string::npos ? 80 : 255);
    options.push_back(s[2]); grids.push_back(std::move(hist));
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        integer(rows[at][1])==nlayer,"invalid grid overlay END/count");
  grids.front()->GetXaxis()->SetTitle(labels.x.c_str());
  grids.front()->GetYaxis()->SetTitle(labels.y.c_str());
  grids.front()->GetZaxis()->SetTitle(labels.z.c_str());
  TCanvas canvas("rootplots_grid_overlay_canvas","",width,height);
  canvas.SetLeftMargin(0.14); canvas.SetRightMargin(0.18); canvas.SetBottomMargin(0.13);
  canvas.SetLogx(logx); canvas.SetLogy(logy); canvas.SetLogz(logz);
  for(int k=0;k<nlayer;++k) grids[k]->Draw(options[k].c_str());
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!labels.title.empty()) heading.DrawLatex(0.13,0.94,labels.title.c_str());
  save(canvas,stem);
}
void draw_volume(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==10 && f[0]=="DRAWVOLUME","invalid DRAWVOLUME");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  std::string option=f[3];
  check_root_option(option);
  palette(f[4]);
  double theta=real(f[5]),phi=real(f[6]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid view angles");
  bool logx=bit(f[7]),logy=bit(f[8]),logz=bit(f[9]);
  auto labels=read_labels(rows,at);
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="VOLUME","missing VOLUME");
  int nx=integer(rows[at][1]),ny=integer(rows[at][2]),nz=integer(rows[at][3]); ++at;
  check(nx>=1 && ny>=1 && nz>=1 && nx<=256 && ny<=256 && nz<=256 &&
        static_cast<long long>(nx)*ny*nz<=1000000,"invalid volume dimensions");
  auto read_edges=[&](const std::string &tag,int n) {
    std::vector<double> edges;
    edges.reserve(n+1);
    for(int i=1;i<=n+1;++i) {
      check(at<rows.size() && rows[at].size()==3 && rows[at][0]==tag &&
            integer(rows[at][1])==i,"invalid volume edge");
      double value=real(rows[at++][2]);
      check(edges.empty() || value>edges.back(),"non-increasing volume edge");
      edges.push_back(value);
    }
    return edges;
  };
  auto xedges=read_edges("XEDGE",nx);
  auto yedges=read_edges("YEDGE",ny);
  auto zedges=read_edges("ZEDGE",nz);
  check((!logx || xedges.front()>0) && (!logy || yedges.front()>0) &&
        (!logz || zedges.front()>0),"log axes need positive volume edges");
  TH3D hist("rootplots_volume","",nx,xedges.data(),ny,yedges.data(),nz,zedges.data());
  hist.SetDirectory(nullptr);
  for(int iz=1;iz<=nz;++iz) for(int iy=1;iy<=ny;++iy) for(int ix=1;ix<=nx;++ix) {
    check(at<rows.size() && rows[at].size()==5 && rows[at][0]=="CELL3" &&
          integer(rows[at][1])==ix && integer(rows[at][2])==iy &&
          integer(rows[at][3])==iz,"invalid volume cell");
    hist.SetBinContent(ix,iy,iz,real(rows[at++][4]));
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        rows[at][1]=="1","invalid volume END/trailing records");
  hist.SetStats(false);
  hist.GetXaxis()->SetTitle(labels.x.c_str());
  hist.GetYaxis()->SetTitle(labels.y.c_str());
  hist.GetZaxis()->SetTitle(labels.z.c_str());
  hist.SetContour(80);
  TCanvas canvas("rootplots_volume_canvas","",width,height);
  canvas.SetLeftMargin(0.13); canvas.SetRightMargin(0.18); canvas.SetBottomMargin(0.13);
  canvas.SetTheta(theta); canvas.SetPhi(phi);
  canvas.SetLogx(logx); canvas.SetLogy(logy); canvas.SetLogz(logz);
  hist.Draw(option.c_str());
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.035);
  if(!labels.title.empty()) heading.DrawLatex(0.13,0.94,labels.title.c_str());
  save(canvas,stem);
}
struct FunctionDrawing {
  std::unique_ptr<TF1> function;
  SceneObjects decorations;
  std::unique_ptr<TLatex> heading;
};
FunctionDrawing paint_function(const std::vector<std::vector<std::string>> &rows,TPad &pad) {
  size_t at=0;
  const auto &f=rows[at++];
  check(f.size()==10 && f[0]=="DRAWFUNCTION","invalid DRAWFUNCTION");
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  bool logx=bit(f[3]),logy=bit(f[4]);
  double xmin=real(f[5]),xmax=real(f[6]),ymin=real(f[7]),ymax=real(f[8]);
  check(xmin<xmax && ymin<ymax && (!logx || xmin>0) && (!logy || ymin>0),"invalid function axes");
  std::string corner=f[9];
  check(corner=="topright" || corner=="topleft" || corner=="bottomright" ||
        corner=="bottomleft" || corner=="none","invalid function legend");
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="LABELS","invalid function labels");
  auto title=rows[at][1],xlabel=rows[at][2],ylabel=rows[at][3]; ++at;
  auto scene=read_scene(rows,at);
  check(at<rows.size() && rows[at].size()>=11 && rows[at][0]=="FUNCTION","missing FUNCTION");
  const auto &v=rows[at++];
  double lo=real(v[2]),hi=real(v[3]); int npx=integer(v[4]),npar=integer(v[5]);
  check(lo<hi && npx>=2 && npx<=1000000 && npar==nparameters(v[1]) &&
        v.size()==static_cast<size_t>(10+npar),"invalid function model/range/resolution");
  color(v[7]); int line_width=integer(v[8]),line_style=integer(v[9]);
  check(line_width>0 && line_style>=1 && line_style<=10,"invalid function style");
  text_field(v[6]);
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && rows[at][1]=="1","invalid function END");
  validate_scene(scene,logx,logy,1);
  FunctionDrawing out;
  out.function=std::make_unique<TF1>("rootplots_model_curve",formula(v[1]).c_str(),lo,hi);
  check(out.function->IsValid() && out.function->GetNpar()==npar,"ROOT rejected model curve");
  for(int i=0;i<npar;++i) out.function->SetParameter(i,real(v[10+i]));
  if(v[1]=="gaus") check(out.function->GetParameter(2)!=0,"Gaussian width cannot be zero");
  if(v[1]=="erfc_edge") check(out.function->GetParameter(3)>0,"edge width must be positive");
  if(v[1]=="double_erf") check(out.function->GetParameter(4)>0 && out.function->GetParameter(5)>0,"edge widths must be positive");
  for(int i=0;i<npx;++i) {
    double x=lo+(hi-lo)*(static_cast<double>(i)/(npx-1));
    double y=out.function->Eval(x);
    check(std::isfinite(x) && std::isfinite(y) && (!logx || x>0) && (!logy || y>0),"invalid visible function/log values");
  }
  out.function->SetNpx(npx);
  out.function->SetLineColor(TColor::GetColor(v[7].c_str()));
  out.function->SetLineWidth(line_width); out.function->SetLineStyle(line_style);
  pad.cd(); apply_canvas(pad,scene); pad.SetLogx(logx); pad.SetLogy(logy);
  auto *frame=pad.DrawFrame(xmin,ymin,xmax,ymax);
  frame->GetXaxis()->SetTitle(xlabel.c_str()); frame->GetYaxis()->SetTitle(ylabel.c_str());
  apply_axes(frame,scene);
  out.function->Draw("SAME");
  out.decorations=draw_scene(scene,corner,{out.function.get()},{v[6]},{"l"});
  out.heading=std::make_unique<TLatex>();
  out.heading->SetNDC(); out.heading->SetTextSize(0.035);
  if(!title.empty()) out.heading->DrawLatex(0.13,0.94,title.c_str());
  return out;
}
void draw_function(const std::vector<std::vector<std::string>> &rows,const std::string &stem) {
  check(!rows.empty() && rows[0].size()==10,"invalid DRAWFUNCTION");
  TCanvas canvas("rootplots_function_canvas","",integer(rows[0][1]),integer(rows[0][2]));
  auto objects=paint_function(rows,canvas);
  save(canvas,stem);
}
struct XYZDrawing {
  std::unique_ptr<TGraph2D> graph;
  std::unique_ptr<TLatex> heading;
};
XYZDrawing paint_xyz(const std::vector<std::vector<std::string>> &rows,TPad &pad) {
  size_t at=0;
  const auto &f=rows[at++];
  check((f.size()==8 || f.size()==20) && f[0]=="DRAWXYZ","invalid DRAWXYZ");
  bool extended=f.size()==20;
  int width=integer(f[1]),height=integer(f[2]);
  check(width>0 && height>0 && width<=10000 && height<=10000,"invalid canvas size");
  color(f[3]);
  int marker=integer(f[4]); double marker_size=real(f[5]);
  check(marker>=1 && marker<=50 && marker_size>0,"invalid marker style");
  double theta=real(f[6]),phi=real(f[7]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid view angles");
  std::string option=extended ? f[8] : "P0",kind=extended ? f[16] : "plain";
  check(kind=="plain" || kind=="symmetric","invalid XYZ error kind");
  int npx=extended ? integer(f[10]) : 40,npy=extended ? integer(f[11]) : 40;
  int contours=extended ? integer(f[12]) : 10;
  check(npx>=4 && npx<=500 && npy>=4 && npy<=500 && contours>=1 && contours<=255,"invalid XYZ grid resolution");
  bool interpolate=extended && bit(f[17]);
  if(extended) {
    check(!bit(f[13]) && !bit(f[14]) && !bit(f[15]),"XYZ does not support logarithmic axes");
    palette(f[9]);
  }
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  check_root_option(option);
  bool triangles=false;
  std::istringstream tokens(upper); std::string token; int ntoken=0;
  while(tokens>>token) {
    ++ntoken;
    const std::vector<std::string> allowed={"P","P0","PCOL","LINE","ERR","TRI","TRIW","TRI1","TRI2",
      "SURF","SURFZ","SURF0","SURF0Z","SURF1","SURF1Z","SURF2","SURF2Z","SURF3","SURF3Z",
      "SURF4","SURF4Z","SURF5","SURF5Z","SURF6","SURF6Z","SURF7","SURF7Z",
      "CONT","CONTZ","CONT0","CONT0Z","CONT1","CONT1Z","CONT2","CONT2Z","CONT3","CONT3Z",
      "CONT4","CONT4Z","CONT5","CONT5Z","COL","COLZ"};
    check(std::find(allowed.begin(),allowed.end(),token)!=allowed.end(),"unsupported TGraph2D option");
    check(token!="ERR" || kind=="symmetric","ERR needs TGraph2DErrors");
    triangles=triangles || token.find("TRI")!=std::string::npos || token.find("SURF")!=std::string::npos ||
      token.find("CONT")!=std::string::npos || token.find("COL")!=std::string::npos;
  }
  check(ntoken>0 && (!triangles || interpolate),"XYZ interpolation needs explicit opt-in");
  auto labels=read_labels(rows,at);
  check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="XYZ","missing XYZ");
  int n=integer(rows[at++][1]);
  check(n>=1 && n<=1000000,"invalid XYZ point count");
  XYZDrawing out;
  if(kind=="symmetric") out.graph=std::make_unique<TGraph2DErrors>(n);
  else out.graph=std::make_unique<TGraph2D>(n);
  auto &graph=*out.graph;
  std::vector<std::pair<double,double>> xy;
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  for(int i=1;i<=n;++i) {
    check(at<rows.size() && rows[at].size()==(kind=="symmetric" ? 8 : 5) &&
          rows[at][0]==(kind=="symmetric" ? "POINT3E" : "POINT3") &&
          integer(rows[at][1])==i,"invalid POINT3");
    const auto &p=rows[at++];
    double x=real(p[2]),y=real(p[3]),z=real(p[4]);
    graph.SetPoint(i-1,x,y,z);
    if(kind=="symmetric") {
      double ex=real(p[5]),ey=real(p[6]),ez=real(p[7]);
      check(ex>=0 && ey>=0 && ez>=0 && std::isfinite(x-ex) && std::isfinite(x+ex) &&
            std::isfinite(y-ey) && std::isfinite(y+ey) && std::isfinite(z-ez) && std::isfinite(z+ez),"invalid XYZ error extent");
      static_cast<TGraph2DErrors&>(graph).SetPointError(i-1,ex,ey,ez);
    }
    xy.emplace_back(x,y);
    xmin=std::min(xmin,x); xmax=std::max(xmax,x);
    ymin=std::min(ymin,y); ymax=std::max(ymax,y);
  }
  check(xmin<xmax && ymin<ymax,"XYZ rendering needs varying x and y coordinates");
  if(triangles) {
    check(n>=3 && std::isfinite(xmax-xmin) && std::isfinite(ymax-ymin),"invalid interpolation coordinate spans");
    auto sorted=xy; std::sort(sorted.begin(),sorted.end());
    check(std::adjacent_find(sorted.begin(),sorted.end())==sorted.end(),"interpolation needs unique xy samples");
    double dx=xmax-xmin,dy=ymax-ymin;
    for(auto &p:xy) {p.first=(p.first-xmin)/dx; p.second=(p.second-ymin)/dy;}
    int far=0; double distance=-1,area=0;
    for(int i=0;i<n;++i) {
      double d=std::pow(xy[i].first-xy[0].first,2)+std::pow(xy[i].second-xy[0].second,2);
      if(d>distance) {far=i;distance=d;}
    }
    for(const auto &p:xy) area=std::max(area,std::abs((xy[far].first-xy[0].first)*(p.second-xy[0].second)-
      (xy[far].second-xy[0].second)*(p.first-xy[0].first)));
    check(area>64*std::numeric_limits<double>::epsilon(),"interpolation needs non-collinear xy samples");
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        rows[at][1]=="1","invalid END/trailing records");
  graph.SetMarkerColor(TColor::GetColor(f[3].c_str()));
  graph.SetMarkerStyle(marker); graph.SetMarkerSize(marker_size);
  graph.SetLineColor(TColor::GetColor(f[3].c_str()));
  if(extended) {
    int line_width=integer(f[18]),line_style=integer(f[19]);
    check(line_width>0 && line_style>=1 && line_style<=10,"invalid XYZ line style");
    graph.SetLineWidth(line_width); graph.SetLineStyle(line_style);
  }
  graph.SetNpx(npx); graph.SetNpy(npy);
  pad.cd();
  pad.SetLeftMargin(0.14); pad.SetBottomMargin(0.13);
  pad.SetRightMargin(0.20); pad.SetTheta(theta); pad.SetPhi(phi);
  if(triangles) graph.GetHistogram()->SetContour(contours);
  graph.Draw(option.c_str());
  graph.GetXaxis()->SetTitle(labels.x.c_str());
  graph.GetYaxis()->SetTitle(labels.y.c_str());
  graph.GetZaxis()->SetTitle(labels.z.c_str());
  out.heading=std::make_unique<TLatex>();
  out.heading->SetNDC(); out.heading->SetTextSize(0.035);
  if(!labels.title.empty()) out.heading->DrawLatex(0.13,0.94,labels.title.c_str());
  return out;
}
void draw_xyz(const std::vector<std::vector<std::string>> &rows, const std::string &stem) {
  check(!rows.empty() && (rows[0].size()==8 || rows[0].size()==20),"invalid DRAWXYZ");
  TCanvas canvas("rootplots_xyz_canvas","",integer(rows[0][1]),integer(rows[0][2]));
  auto objects=paint_xyz(rows,canvas);
  save(canvas,stem);
}
#include "advanced_histograms.hpp"
#include "primitives.hpp"
#include "specialgraphs.hpp"
#include "layouts.hpp"
}
// Compare two images with ROOT's own image reader (used by visual regression tests).
// Prints "SIZE w1 h1 w2 h2" for different dimensions, otherwise
// "DIFF width height max_channel_diff mean_channel_diff changed_pixels" where a
// pixel counts as changed when any RGB channel differs by more than 32/255.
int image_diff(const char *left,const char *right) {
  gROOT->SetBatch(kTRUE);
  std::unique_ptr<TImage> a(TImage::Open(left)), b(TImage::Open(right));
  check(a && b && a->IsValid() && b->IsValid(),"cannot read image");
  const unsigned w=a->GetWidth(), h=a->GetHeight();
  if(w!=b->GetWidth() || h!=b->GetHeight()) {
    std::cout<<"SIZE\t"<<w<<'\t'<<h<<'\t'<<b->GetWidth()<<'\t'<<b->GetHeight()<<'\n';
    return 0;
  }
  const UInt_t *pa=a->GetArgbArray(), *pb=b->GetArgbArray();
  check(pa && pb,"cannot decode image pixels");
  unsigned maxdiff=0; unsigned long long changed=0; double total=0;
  for(size_t i=0;i<size_t(w)*h;++i) {
    unsigned local=0;
    for(int shift: {16,8,0}) {
      const int ca=int((pa[i]>>shift)&0xff), cb=int((pb[i]>>shift)&0xff);
      const unsigned d=unsigned(std::abs(ca-cb));
      total+=d; local=std::max(local,d);
    }
    maxdiff=std::max(maxdiff,local);
    if(local>32) ++changed;
  }
  std::cout<<"DIFF\t"<<w<<'\t'<<h<<'\t'<<maxdiff<<'\t'<<total/(3.0*w*h)<<'\t'<<changed<<'\n';
  return 0;
}
int main(int argc,char **argv) {
  try {
    if(argc==4 && std::string(argv[1])=="--image-diff") return image_diff(argv[2],argv[3]);
    check(argc==3,"usage: rootplots-graphfit input.tsv output | --image-diff a.png b.png");
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0); gStyle->SetOptTitle(0);
    auto rows=read_input(argv[1]);
    if(rows[0][0]=="DRAW") draw(rows,argv[2]);
    else if(rows[0][0]=="DRAWHIST" || rows[0][0]=="DRAWHISTRAW") draw_hist(rows,argv[2]);
    else if(rows[0][0]=="DRAWHISTOVERLAY") draw_hist_overlay(rows,argv[2]);
    else if(rows[0][0]=="DRAWFIT") draw_fit(rows,argv[2]);
    else if(rows[0][0]=="FIT") fit(rows,argv[2]);
    else if(rows[0][0]=="DRAWGRID" || rows[0][0]=="DRAWGRIDRAW") draw_grid(rows,argv[2]);
    else if(rows[0][0]=="DRAWGRIDOVERLAY") draw_grid_overlay(rows,argv[2]);
    else if(rows[0][0]=="DRAWVOLUME") draw_volume(rows,argv[2]);
    else if(rows[0][0]=="DRAWXYZ") draw_xyz(rows,argv[2]);
    else if(rows[0][0]=="DRAWFUNCTION") draw_function(rows,argv[2]);
    else if(rows[0][0]=="DRAWLAYOUT") draw_layout(rows,argv[2]);
    else if(!advanced_dispatch(rows,argv[2]) && !special_dispatch(rows,argv[2])) throw std::runtime_error("unknown operation");
    return 0;
  } catch(const std::exception &e) {
    std::cerr<<"ROOTPlots graphfit: "<<e.what()<<'\n';
    return 2;
  }
}
