// Included inside graphfit.cpp's anonymous namespace, before layouts.hpp.
// Native TGraphBentErrors / TGraphMultiErrors layers for basic XY panels.
struct ErrorLayer {
  std::vector<std::unique_ptr<TObject>> owned;
  TObject *primary=nullptr;
  std::string label,option;
  struct Extra { TObject *object; std::string label,option; };
  std::vector<Extra> extra;
};
void style_error_graph(TGraph &g,const std::string &hex,int marker,double marker_size,int width,int line_style) {
  color(hex);
  check(marker>=1 && marker<=50 && marker_size>0 && width>0 && line_style>=1 && line_style<=10,"invalid error-graph style");
  int c=TColor::GetColor(hex.c_str());
  g.SetMarkerColor(c); g.SetLineColor(c); g.SetMarkerStyle(marker); g.SetMarkerSize(marker_size);
  g.SetLineWidth(width); g.SetLineStyle(line_style);
}
bool error_draw_mode(const std::string &d) {
  return d=="points" || d=="line" || d=="linepoints" || d=="curve" || d=="curvepoints";
}
std::string error_base_option(const std::string &d) {
  return d=="line" ? "L" : d=="linepoints" ? "LP" : d=="curve" ? "C" : d=="curvepoints" ? "CP" : "P";
}
ErrorLayer paint_bent_layer(const std::vector<std::vector<std::string>> &rows,size_t &at,bool logx,bool logy) {
  check(at<rows.size() && rows[at].size()==13 && rows[at][0]=="BENT","invalid BENT");
  const auto f=rows[at++];
  int n=integer(f[1]); check(n>0 && n<=1000000,"invalid bent point count");
  // Reuse the TGraphAsymmErrors option/legend rules: the bent painter accepts the same codes.
  Series s{"asymmetric",f[2],f[3],f[4],f[5],integer(f[6]),integer(f[8]),integer(f[9]),real(f[7]),{}};
  check(error_draw_mode(s.draw),"invalid bent draw mode");
  check(s.errors=="bars" || s.errors=="none" || s.errors=="band" || s.errors=="smoothband" ||
        s.errors=="boxes" || s.errors=="caps" || s.errors=="nocaps" || s.errors=="arrows" ||
        s.errors=="filledarrows","invalid bent error mode");
  text_field(s.label);
  s.fill_hex=f[10]; s.fill_style=integer(f[11]); s.fill_alpha=real(f[12]);
  if(s.fill_hex!="-") color(s.fill_hex);
  check(s.fill_style>=0 && s.fill_style<=9999 && s.fill_alpha>=0 && s.fill_alpha<=1,"invalid bent fill style");
  auto g=std::make_unique<TGraphBentErrors>(n);
  for(int i=1;i<=n;++i) {
    check(at<rows.size() && rows[at].size()==12 && rows[at][0]=="BPOINT" && integer(rows[at][1])==i,"invalid BPOINT");
    const auto &p=rows[at++];
    double v[10]; for(int k=0;k<10;++k) v[k]=real(p[k+2]);
    check(v[2]>=0 && v[3]>=0 && v[4]>=0 && v[5]>=0,"bent errors must be non-negative");
    double ends[8]={v[0]-v[2],v[0]+v[3],v[1]-v[4],v[1]+v[5],v[1]+v[6],v[1]+v[7],v[0]+v[8],v[0]+v[9]};
    for(double e:ends) check(std::isfinite(e),"bent error end points overflow");
    if(s.errors!="none") {
      check(!logx || (v[0]>0 && ends[0]>0 && ends[6]>0 && ends[7]>0),"log x needs positive bent error end points");
      check(!logy || (v[1]>0 && ends[2]>0 && ends[4]>0 && ends[5]>0),"log y needs positive bent error end points");
    } else check((!logx || v[0]>0) && (!logy || v[1]>0),"log axes need positive points");
    g->SetPoint(i-1,v[0],v[1]);
    g->SetPointError(i-1,v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]);
  }
  style_error_graph(*g,s.hex,s.marker,s.marker_size,s.line_width,s.line_style);
  if(filled_graph(s)) {
    g->SetFillColorAlpha(TColor::GetColor((s.fill_hex=="-" ? s.hex : s.fill_hex).c_str()),s.fill_alpha);
    g->SetFillStyle(s.fill_style);
  }
  g->Draw(draw_option(s).c_str());
  ErrorLayer out;
  out.primary=g.get(); out.label=s.label; out.option=legend_option(s);
  out.owned.push_back(std::move(g));
  return out;
}
ErrorLayer paint_multi_layer(const std::vector<std::vector<std::string>> &rows,size_t &at,bool logx,bool logy) {
  check(at<rows.size() && rows[at].size()==12 && rows[at][0]=="MULTI","invalid MULTI");
  const auto f=rows[at++];
  int n=integer(f[1]),nsets=integer(f[2]);
  check(n>0 && n<=1000000 && nsets>=1 && nsets<=16,"invalid multi-error dimensions");
  std::string draw=f[3],sum=f[4],label=f[5];
  check(error_draw_mode(draw),"invalid multi-error draw mode");
  check(sum=="first" || sum=="quadrature" || sum=="linear","invalid multi-error sum mode");
  text_field(label);
  bool xerrors=bit(f[11]);
  struct Set { std::string display,label,hex,fill_hex; int width,line_style,fill_style; double alpha,xscale; };
  std::vector<Set> sets;
  for(int k=1;k<=nsets;++k) {
    check(at<rows.size() && rows[at].size()==11 && rows[at][0]=="ESET" && integer(rows[at][1])==k,"invalid ESET");
    const auto &e=rows[at++];
    Set set{e[2],e[3],e[4],e[7],integer(e[5]),integer(e[6]),integer(e[8]),real(e[9]),real(e[10])};
    const auto &d=set.display;
    check(d=="bars" || d=="nocaps" || d=="caps" || d=="arrows" || d=="filledarrows" || d=="boxes" ||
          d=="boxes_outline" || d=="band" || d=="smoothband" || d=="none","invalid error-set display");
    text_field(set.label); color(set.hex); color(set.fill_hex);
    check(set.width>0 && set.line_style>=1 && set.line_style<=10 && set.fill_style>=0 && set.fill_style<=9999 &&
          set.alpha>=0 && set.alpha<=1 && set.xscale>0 && set.xscale<=100,"invalid error-set style");
    sets.push_back(set);
  }
  auto g=std::make_unique<TGraphMultiErrors>(n,nsets);
  for(int i=1;i<=n;++i) {
    check(at<rows.size() && rows[at].size()==static_cast<size_t>(6+2*nsets) && rows[at][0]=="MPOINT" &&
          integer(rows[at][1])==i,"invalid MPOINT");
    const auto &p=rows[at++];
    double x=real(p[2]),y=real(p[3]),exl=real(p[4]),exh=real(p[5]);
    check(exl>=0 && exh>=0 && std::isfinite(x-exl) && std::isfinite(x+exh),"invalid multi-error x errors");
    check(xerrors || (exl==0 && exh==0),"x errors present but not declared");
    check((!logx || (x>0 && x-exl>0)) && (!logy || y>0),"log axes need positive points and x errors");
    g->SetPoint(i-1,x,y); g->SetPointEX(i-1,exl,exh);
    for(int k=0;k<nsets;++k) {
      double lo=real(p[6+2*k]),hi=real(p[7+2*k]);
      check(lo>=0 && hi>=0 && std::isfinite(y-lo) && std::isfinite(y+hi),"invalid multi-error y errors");
      check(!logy || sets[k].display=="none" || y-lo>0,"log y needs positive visible error end points");
      g->SetPointEY(i-1,k,lo,hi);
    }
  }
  g->SetSumErrorsMode(sum=="first" ? TGraphMultiErrors::kOnlyFirst :
                      sum=="quadrature" ? TGraphMultiErrors::kSquareSum : TGraphMultiErrors::kAbsSum);
  style_error_graph(*g,f[6],integer(f[7]),real(f[8]),integer(f[9]),integer(f[10]));
  // Block 0: markers/line and x errors; "S" selects the per-set attributes below.
  std::string option=error_base_option(draw)+" S";
  if(!xerrors) option+=" X";
  bool painter_draws_markers=false;
  ErrorLayer out;
  for(int k=0;k<nsets;++k) {
    const auto &set=sets[k];
    const auto &d=set.display;
    std::string block=d=="nocaps" ? "Z" : d=="caps" ? "[]" : d=="arrows" ? ">" : d=="filledarrows" ? "|>" :
                      d=="boxes" ? "2" : d=="boxes_outline" ? "5" : d=="band" ? "3" : d=="smoothband" ? "4" :
                      d=="none" ? "X" : "";
    if(set.xscale!=1) { std::ostringstream scale; scale<<std::setprecision(17)<<" s="<<set.xscale; block+=scale.str(); }
    option+=" ; "+block;
    painter_draws_markers|=!(d=="caps" || d=="boxes" || d=="boxes_outline" || d=="band" || d=="smoothband");
    auto *line=g->GetAttLine(k); auto *fill=g->GetAttFill(k);
    check(line!=nullptr && fill!=nullptr,"missing error-set attributes");
    int c=TColor::GetColor(set.hex.c_str());
    line->SetLineColor(c); line->SetLineWidth(set.width); line->SetLineStyle(set.line_style);
    fill->SetFillColorAlpha(TColor::GetColor(set.fill_hex.c_str()),set.alpha); fill->SetFillStyle(set.fill_style);
    if(!set.label.empty()) {
      // A small undrawn proxy carries this set's attributes into the legend.
      auto proxy=std::make_unique<TGraph>(1);
      proxy->SetLineColor(c); proxy->SetLineWidth(set.width); proxy->SetLineStyle(set.line_style);
      proxy->SetFillColorAlpha(TColor::GetColor(set.fill_hex.c_str()),set.alpha); proxy->SetFillStyle(set.fill_style);
      bool filled=d=="boxes" || d=="boxes_outline" || d=="band" || d=="smoothband";
      out.extra.push_back({proxy.get(),set.label,filled ? "f" : "l"});
      out.owned.push_back(std::move(proxy));
    }
  }
  g->Draw(option.c_str());
  out.primary=g.get(); out.label=label;
  bool has_points=draw=="points" || draw=="linepoints" || draw=="curvepoints";
  bool has_line=draw!="points";
  out.option=std::string(has_line ? "l" : "")+(has_points ? "p" : "");
  out.owned.push_back(std::move(g));
  if(!painter_draws_markers) {
    // ROOT's multi-error painter omits markers/lines when every set is drawn as
    // boxes, bands or caps; draw the points on top with the same attributes.
    auto points=std::make_unique<TGraph>(n);
    auto *multi=static_cast<TGraphMultiErrors*>(out.owned.back().get());
    for(int i=0;i<n;++i) points->SetPoint(i,multi->GetX()[i],multi->GetY()[i]);
    style_error_graph(*points,f[6],integer(f[7]),real(f[8]),integer(f[9]),integer(f[10]));
    points->Draw((error_base_option(draw)+" SAME").c_str());
    out.owned.push_back(std::move(points));
  }
  return out;
}
