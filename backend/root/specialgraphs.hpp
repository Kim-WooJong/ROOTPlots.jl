struct SpecialData {
  std::string kind,corner,title,xlabel,ylabel;
  int n,order,npx,width,height,marker,line_width,line_style,fill_style;
  double xmin,xmax,ymin,ymax,confidence,marker_size,alpha;
  bool logx,logy;
  std::string color,fill_color;
  Scene scene;
  std::vector<double> a,b,c,d;
  std::vector<std::string> labels,colors;
};
bool special_command(const std::string &command) { return command=="DRAWSPECIAL"; }
SpecialData read_special(const std::vector<std::vector<std::string>>& rows) {
  check(!rows.empty() && rows[0].size()==15 && (special_command(rows[0][0]) || rows[0][0]=="SPECIALPOINTS"),"invalid special graph request");
  auto f=rows[0]; SpecialData s; size_t at=1;
  s.width=integer(f[1]); s.height=integer(f[2]); s.kind=f[3]; s.n=integer(f[4]); s.logx=bit(f[5]); s.logy=bit(f[6]);
  s.xmin=real(f[7]); s.xmax=real(f[8]); s.ymin=real(f[9]); s.ymax=real(f[10]); s.corner=f[11];
  s.confidence=real(f[12]); s.order=integer(f[13]); s.npx=integer(f[14]);
  check(s.width>0 && s.width<=10000 && s.height>0 && s.height<=10000 && s.n>0 && s.n<=1000000 &&
        s.xmin<s.xmax && s.ymin<s.ymax && (!s.logx || s.xmin>0) && (!s.logy || s.ymin>0),"invalid special graph dimensions/ranges");
  check(s.kind=="polar" || s.kind=="qq" || s.kind=="spline" || s.kind=="pie" || s.kind=="efficiency","invalid special graph kind");
  check(s.corner=="none" || s.corner=="topleft" || s.corner=="topright" || s.corner=="bottomleft" || s.corner=="bottomright","invalid special legend");
  check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="LABELS","invalid special labels");
  s.title=rows[at][1]; s.xlabel=rows[at][2]; s.ylabel=rows[at++][3];
  s.scene=read_scene(rows,at); validate_scene(s.scene,s.logx,s.logy,1);
  check(at<rows.size() && rows[at].size()==9 && rows[at][0]=="STYLE","invalid special style");
  f=rows[at++]; s.color=f[1]; s.marker=integer(f[2]); s.marker_size=real(f[3]); s.line_width=integer(f[4]); s.line_style=integer(f[5]);
  s.fill_color=f[6]; s.fill_style=integer(f[7]); s.alpha=real(f[8]); color(s.color); color(s.fill_color);
  check(s.marker>=1 && s.marker<=50 && s.marker_size>0 && s.line_width>0 && s.line_style>=1 && s.line_style<=10 && s.fill_style>=0 && s.fill_style<=9999 && s.alpha>=0 && s.alpha<=1,"invalid special style");
  for(int i=0;i<s.n;++i) {
    check(at<rows.size(),"missing special data"); f=rows[at++];
    if(s.kind=="pie") {
      check(s.n<=1000 && f.size()==5 && f[0]=="SLICE" && integer(f[1])==i+1,"invalid pie slice");
      double v=real(f[2]); check(v>=0,"negative pie slice"); color(f[4]); text_field(f[3]);
      s.a.push_back(v); s.labels.push_back(f[3]); s.colors.push_back(f[4]);
    } else if(s.kind=="efficiency") {
      check(f.size()==6 && f[0]=="BIN" && integer(f[1])==i+1,"invalid efficiency bin");
      double left=real(f[2]),right=real(f[3]),passed=real(f[4]),total=real(f[5]);
      check(left<right && (i==0 || left==s.b.back()) && passed>=0 && passed<=total && total<=9007199254740992. &&
            passed==std::round(passed) && total==std::round(total),"invalid efficiency counts/edges");
      s.a.push_back(left); s.b.push_back(right); s.c.push_back(passed); s.d.push_back(total);
    } else {
      check(f.size()==static_cast<size_t>(s.kind=="polar" ? 6 : 4) && f[0]=="POINT" && integer(f[1])==i+1,"invalid special point");
      double x=real(f[2]),y=real(f[3]);
      if(s.kind=="spline") check(s.a.empty() || x>s.a.back(),"spline x must increase");
      s.a.push_back(x); s.b.push_back(y);
      if(s.kind=="polar") { double ex=real(f[4]),ey=real(f[5]); check(y>=0 && ex>=0 && ey>=0,"invalid polar radius/errors"); s.c.push_back(ex); s.d.push_back(ey); }
    }
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && rows[at][1]=="1","invalid special END");
  if(s.kind=="spline") check((s.order==3 || s.order==5) && s.n>=s.order && s.npx>=2 && s.npx<=1000000,"invalid spline order/samples");
  if(s.kind=="efficiency") check(s.confidence>0 && s.confidence<1,"invalid efficiency confidence");
  if(s.kind=="pie") { double sum=0; for(double v:s.a) sum+=v; check(std::isfinite(sum) && sum>0,"invalid pie total"); }
  if(s.kind=="polar" || s.kind=="qq") check(s.n>=2,"not enough graph samples");
  return s;
}
std::unique_ptr<TSpline> special_spline(SpecialData &s) {
  std::unique_ptr<TSpline> spline;
  if(s.order==3) spline=std::make_unique<TSpline3>("rootplots_spline",s.a.data(),s.b.data(),s.n);
  else spline=std::make_unique<TSpline5>("rootplots_spline",s.a.data(),s.b.data(),s.n);
  spline->SetNpx(s.npx); spline->SetLineColor(TColor::GetColor(s.color.c_str()));
  spline->SetLineWidth(s.line_width); spline->SetLineStyle(s.line_style);
  for(int i=0;i<s.npx;++i) check(std::isfinite(spline->Eval(s.a.front()+(s.a.back()-s.a.front())*i/(s.npx-1))),"nonfinite ROOT spline");
  return spline;
}
std::unique_ptr<TEfficiency> special_efficiency(const SpecialData &s) {
  std::vector<double> edges=s.a; edges.push_back(s.b.back());
  TH1D passed("rootplots_passed","",s.n,edges.data()),total("rootplots_total","",s.n,edges.data());
  passed.SetDirectory(nullptr); total.SetDirectory(nullptr);
  for(int i=0;i<s.n;++i) { passed.SetBinContent(i+1,s.c[i]); total.SetBinContent(i+1,s.d[i]); }
  check(TEfficiency::CheckConsistency(passed,total),"ROOT rejected efficiency inputs");
  auto efficiency=std::make_unique<TEfficiency>(passed,total); efficiency->SetDirectory(nullptr);
  efficiency->SetConfidenceLevel(s.confidence); efficiency->SetStatisticOption(TEfficiency::kFCP);
  return efficiency;
}
struct SpecialDrawing { std::vector<std::unique_ptr<TObject>> objects; SceneObjects decorations; };
SpecialDrawing paint_special(const std::vector<std::vector<std::string>>& rows,TPad& pad) {
  auto s=read_special(rows); SpecialDrawing out;
  pad.cd(); apply_canvas(pad,s.scene); pad.SetLogx(s.logx); pad.SetLogy(s.logy);
  TH1 *frame=nullptr; TGraph *graph=nullptr;
  if(s.kind!="pie" && s.kind!="polar") {
    frame=pad.DrawFrame(s.xmin,s.ymin,s.xmax,s.ymax);
    frame->GetXaxis()->SetTitle(s.xlabel.c_str()); frame->GetYaxis()->SetTitle(s.ylabel.c_str()); apply_axes(frame,s.scene);
  }
  if(s.kind=="polar") {
    auto polar=std::make_unique<TGraphPolar>(s.n,s.a.data(),s.b.data(),s.c.data(),s.d.data());
    graph=polar.get(); polar->Draw("AP"); out.objects.push_back(std::move(polar));
  } else if(s.kind=="qq") {
    // ROOT's constructor puts its first sample on Y; the Julia API names
    // samples by displayed axis, so reverse arguments deliberately.
    auto qq=std::make_unique<TGraphQQ>(s.n,s.b.data(),s.n,s.a.data()); graph=qq.get(); qq->Draw("P SAME"); out.objects.push_back(std::move(qq));
  } else if(s.kind=="spline") {
    auto spline=special_spline(s); spline->Draw("SAME"); out.objects.push_back(std::move(spline));
  } else if(s.kind=="efficiency") {
    auto efficiency=special_efficiency(s);
    std::unique_ptr<TGraphAsymmErrors> points(efficiency->CreateGraph()); check(points!=nullptr,"ROOT failed to create efficiency graph");
    graph=points.get(); points->Draw("P SAME"); out.objects.push_back(std::move(points)); out.objects.push_back(std::move(efficiency));
  } else {
    auto pie=std::make_unique<TPie>("rootplots_pie",s.title.c_str(),s.n,s.a.data());
    for(int i=0;i<s.n;++i) { pie->SetEntryLabel(i,s.labels[i].c_str()); pie->SetEntryFillColor(i,TColor::GetColor(s.colors[i].c_str())); }
    pie->Draw("3d"); out.objects.push_back(std::move(pie));
  }
  if(graph) {
    graph->SetMarkerColor(TColor::GetColor(s.color.c_str())); graph->SetLineColor(TColor::GetColor(s.color.c_str()));
    graph->SetMarkerStyle(s.marker); graph->SetMarkerSize(s.marker_size); graph->SetLineWidth(s.line_width); graph->SetLineStyle(s.line_style);
  }
  std::vector<TObject*> refs; for(const auto &object:out.objects) refs.push_back(object.get());
  // Only the displayed object participates in the one-series scene contract.
  refs.resize(1);
  out.decorations=draw_scene(s.scene,s.corner,refs,{""},{"l"});
  if(!root_title(pad,frame,s.title,s.scene.canvas) && !s.title.empty()) { auto heading=std::make_unique<TLatex>(0.13,0.94,s.title.c_str()); heading->SetNDC(); heading->SetTextSize(.035); heading->Draw(); out.objects.push_back(std::move(heading)); }
  return out;
}
bool special_dispatch(const std::vector<std::vector<std::string>>& rows,const std::string &output) {
  if(rows[0][0]=="SPECIALPOINTS") {
    auto s=read_special(rows); check(s.kind=="spline" || s.kind=="efficiency" || s.kind=="qq","unsupported point summary");
    std::ofstream stream(output); check(stream.good(),"cannot write special points"); stream<<std::setprecision(17)<<"ROOTPlotsSpecialPoints\t1\n";
    if(s.kind=="spline") {
      auto spline=special_spline(s);
      for(int i=0;i<s.npx;++i) { double x=s.a.front()+(s.a.back()-s.a.front())*i/(s.npx-1); stream<<"POINT\t"<<x<<'\t'<<spline->Eval(x)<<"\t0\t0\t0\t0\n"; }
    } else if(s.kind=="qq") {
      TGraphQQ graph(s.n,s.b.data(),s.n,s.a.data());
      for(int i=0;i<graph.GetN();++i) stream<<"POINT\t"<<graph.GetPointX(i)<<'\t'<<graph.GetPointY(i)<<"\t0\t0\t0\t0\n";
    } else {
      auto efficiency=special_efficiency(s);
      for(int i=0;i<s.n;++i) if(s.d[i]>0) { double x=(s.a[i]+s.b[i])/2;
        stream<<"POINT\t"<<x<<'\t'<<efficiency->GetEfficiency(i+1)<<'\t'<<x-s.a[i]<<'\t'<<s.b[i]-x<<'\t'<<efficiency->GetEfficiencyErrorLow(i+1)<<'\t'<<efficiency->GetEfficiencyErrorUp(i+1)<<'\n'; }
    }
    stream<<"END\n"; stream.flush(); check(stream.good(),"failed writing special points"); return true;
  }
  if(!special_command(rows[0][0])) return false;
  TCanvas canvas("rootplots_special","",integer(rows[0][1]),integer(rows[0][2]));
  auto objects=paint_special(rows,canvas); save(canvas,output); return true;
}
