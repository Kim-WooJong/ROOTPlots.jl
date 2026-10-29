// Included inside graphfit.cpp's anonymous namespace. ROOT headers belong in
// the translation unit, never here. All painters retain their primitives until
// the caller has exported the containing canvas.
struct AdvancedHeader {
  std::string command,option,palette_name,error,corner;
  int width,height,dim,contours,count;
  bool logx,logy,logz,explicit_x,explicit_y,stats;
  double theta,phi,xmin,xmax,ymin,ymax;
  Labels labels;
  Scene scene;
};

bool advanced_command(const std::string &command) {
  return command=="DRAWSTACK" || command=="DRAWPROFILE" || command=="DRAWPOLY";
}

AdvancedHeader advanced_header(const std::vector<std::vector<std::string>> &rows,size_t &at) {
  check(!rows.empty(),"empty advanced request");
  const auto &f=rows[at++];
  check(f.size()==22 && (advanced_command(f[0]) || f[0]=="PROFILESTATS"),"invalid advanced header");
  AdvancedHeader h;
  h.command=f[0]; h.width=integer(f[1]); h.height=integer(f[2]); h.dim=integer(f[3]);
  h.option=f[4]; h.logx=bit(f[5]); h.logy=bit(f[6]); h.logz=bit(f[7]);
  h.palette_name=f[8]; h.contours=integer(f[9]); h.theta=real(f[10]); h.phi=real(f[11]);
  h.explicit_x=bit(f[12]); h.xmin=real(f[13]); h.xmax=real(f[14]);
  h.explicit_y=bit(f[15]); h.ymin=real(f[16]); h.ymax=real(f[17]);
  h.stats=bit(f[18]); h.count=integer(f[19]); h.error=f[20]; h.corner=f[21];
  check(h.width>=1 && h.width<=10000 && h.height>=1 && h.height<=10000,"invalid advanced canvas size");
  check(h.dim>=1 && h.dim<=3,"invalid advanced dimension");
  check_root_option(h.option);
  check(h.palette_name=="bird" || h.palette_name=="viridis","invalid advanced palette");
  check(h.contours>=1 && h.contours<=255 && h.theta>=0 && h.theta<=180 && h.phi>=-180 && h.phi<=180,
        "invalid advanced view");
  check(h.xmin<h.xmax && h.ymin<h.ymax && (!h.explicit_x || !h.logx || h.xmin>0) &&
        (!h.explicit_y || !h.logy || h.ymin>0),"invalid advanced limits");
  check(h.corner=="topright" || h.corner=="topleft" || h.corner=="bottomright" ||
        h.corner=="bottomleft" || h.corner=="none","invalid advanced legend");
  check(at<rows.size() && rows[at].size()==5 && rows[at][0]=="LABELS","invalid advanced LABELS");
  for(size_t i=1;i<5;++i) text_field(rows[at][i]);
  h.labels={rows[at][1],rows[at][2],rows[at][3],rows[at][4]}; ++at;
  h.scene=read_scene(rows,at);
  return h;
}

struct AdvancedDrawing {
  std::vector<std::unique_ptr<TH1>> histograms;
  std::unique_ptr<THStack> stack;
  SceneObjects decorations;
  TH1 *profile=nullptr;
  int dimension=0;
};

std::vector<double> advanced_edges(const std::vector<std::vector<std::string>> &rows,
                                  size_t &at,int axis,int dim) {
  check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="PEDGES" &&
        integer(rows[at][1])==axis,"missing profile edges");
  int n=integer(rows[at++][2]);
  check(n>=2 && n<=(dim==3 ? 257 : 2049),"invalid profile axis size");
  std::vector<double> edges;
  for(int i=1;i<=n;++i) {
    check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="PEDGE" &&
          integer(rows[at][1])==axis && integer(rows[at][2])==i,"invalid PEDGE");
    double value=real(rows[at++][3]);
    check(edges.empty() || value>edges.back(),"profile edges must increase");
    edges.push_back(value);
  }
  return edges;
}

std::string advanced_error_option(const std::string &error) {
  if(error=="mean") return "";
  if(error=="spread") return "s";
  if(error=="integer") return "i";
  if(error=="gaussian") return "g";
  throw std::runtime_error("invalid profile error mode");
}

template<class P> void advanced_check_profile(P &profile) {
  double statistics[TH1::kNstat]={};
  profile.GetStats(statistics);
  for(double value:statistics) check(std::isfinite(value),"profile aggregate statistics overflow");
  for(int bin=0;bin<profile.GetNcells();++bin) {
    check(std::isfinite(profile.GetBinContent(bin)) && std::isfinite(profile.GetBinError(bin)) &&
          std::isfinite(profile.GetBinEntries(bin)) && std::isfinite(profile.GetBinEffectiveEntries(bin)),
          "profile moments overflow");
    check(std::isfinite(profile.GetSumw2()->At(bin)) && std::isfinite(profile.GetBinSumw2()->At(bin)),
          "profile squared moments overflow");
  }
}

void advanced_profile(const std::vector<std::vector<std::string>> &rows,size_t &at,
                      const AdvancedHeader &h,AdvancedDrawing &drawing) {
  check(h.count>=0 && h.count<=1000000,"invalid profile sample count");
  std::vector<std::vector<double>> edges;
  long long cells=1;
  for(int axis=1;axis<=h.dim;++axis) {
    edges.push_back(advanced_edges(rows,at,axis,h.dim)); cells*=edges.back().size()-1;
    check(cells<=1000000,"too many profile cells");
    bool logarithmic=axis==1 ? h.logx : axis==2 ? h.logy : h.logz;
    check(!logarithmic || edges.back().front()>0,"log profile axis needs positive edges");
  }
  std::string error=advanced_error_option(h.error);
  std::unique_ptr<TH1> owner;
  TProfile *p1=nullptr; TProfile2D *p2=nullptr; TProfile3D *p3=nullptr;
  if(h.dim==1) {
    auto p=std::make_unique<TProfile>("rootplots_profile_1d","",edges[0].size()-1,edges[0].data());
    p1=p.get(); p1->SetErrorOption(error.c_str()); owner=std::move(p);
  } else if(h.dim==2) {
    auto p=std::make_unique<TProfile2D>("rootplots_profile_2d","",edges[0].size()-1,edges[0].data(),
                                      edges[1].size()-1,edges[1].data());
    p2=p.get(); p2->SetErrorOption(error.c_str()); owner=std::move(p);
  } else {
    auto p=std::make_unique<TProfile3D>("rootplots_profile_3d","",edges[0].size()-1,edges[0].data(),
        edges[1].size()-1,edges[1].data(),edges[2].size()-1,edges[2].data());
    p3=p.get(); p3->SetErrorOption(error.c_str()); owner=std::move(p);
  }
  owner->SetDirectory(nullptr); owner->Sumw2(); owner->SetStats(h.stats);
  for(int i=1;i<=h.count;++i) {
    check(at<rows.size() && rows[at].size()==static_cast<size_t>(h.dim+4) && rows[at][0]=="SAMPLE" &&
          integer(rows[at][1])==i,"invalid profile SAMPLE");
    const auto &r=rows[at++];
    std::vector<double> v;
    for(int j=0;j<h.dim+1;++j) v.push_back(real(r[j+2]));
    double w=real(r[h.dim+3]);
    check(w>=0 && (h.error!="gaussian" || w>0) && std::isfinite(w*w),"invalid profile weight");
    for(double value:v) check(std::isfinite(w*value) && std::isfinite(w*value*value),"profile sample moments overflow");
    if(p1) p1->Fill(v[0],v[1],w);
    else if(p2) p2->Fill(v[0],v[1],v[2],w);
    else p3->Fill(v[0],v[1],v[2],v[3],w);
  }
  if(p1) advanced_check_profile(*p1);
  if(p2) advanced_check_profile(*p2);
  if(p3) advanced_check_profile(*p3);
  drawing.profile=owner.get(); drawing.dimension=h.dim;
  drawing.histograms.push_back(std::move(owner));
}

void advanced_stack(const std::vector<std::vector<std::string>> &rows,size_t &at,
                    const AdvancedHeader &h,AdvancedDrawing &drawing,
                    std::vector<std::string> &labels) {
  check(h.dim==1 && h.count>=1 && h.count<=64 && h.error=="-","invalid stack header");
  std::string upper=h.option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  bool unstacked=upper.find("NOSTACK")!=std::string::npos;
  drawing.stack=std::make_unique<THStack>("rootplots_stack","");
  std::vector<double> first_edges;
  for(int k=1;k<=h.count;++k) {
    check(at<rows.size() && rows[at].size()==10 && rows[at][0]=="HLAYER" &&
          integer(rows[at][1])==k,"invalid stack HLAYER");
    const auto style=rows[at++];
    text_field(style[2]); color(style[3]); color(style[8]);
    int marker=integer(style[4]),width=integer(style[6]),line=integer(style[7]),fill=integer(style[9]);
    double size=real(style[5]);
    check(marker>=1 && marker<=50 && size>0 && width>0 && line>=1 && line<=10 &&
          (fill==0 || fill==1001 || (fill>=3000 && fill<=3999)),"invalid stack style");
    check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="HIST","missing stack HIST");
    int n=integer(rows[at++][1]); check(n>=1 && n<=1000000,"invalid stack bin count");
    std::vector<double> edges;
    for(int i=1;i<=n+1;++i) {
      check(at<rows.size() && rows[at].size()==3 && rows[at][0]=="EDGE" &&
            integer(rows[at][1])==i,"invalid stack EDGE");
      double edge=real(rows[at++][2]);
      check(edges.empty() || edge>edges.back(),"stack edges must increase"); edges.push_back(edge);
    }
    check(!h.logx || edges.front()>0,"log stack x needs positive edges");
    if(k==1) first_edges=edges;
    check(unstacked || edges==first_edges,"stacked histograms need identical edges");
    auto hist=std::make_unique<TH1D>(("rootplots_stack_"+std::to_string(k)).c_str(),"",n,edges.data());
    hist->SetDirectory(nullptr); hist->Sumw2(); hist->SetStats(h.stats);
    for(int i=1;i<=n;++i) {
      check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="BIN" && integer(rows[at][1])==i,"invalid stack BIN");
      double value=real(rows[at][2]),variance=real(rows[at][3]); ++at;
      check(variance>=0 && (unstacked || value>=0),"negative stack bin or variance");
      hist->SetBinContent(i,value); hist->SetBinError(i,std::sqrt(variance));
    }
    int c=TColor::GetColor(style[3].c_str());
    hist->SetLineColor(c); hist->SetMarkerColor(c); hist->SetLineWidth(width); hist->SetLineStyle(line);
    hist->SetMarkerStyle(marker); hist->SetMarkerSize(size);
    hist->SetFillColor(TColor::GetColor(style[8].c_str())); hist->SetFillStyle(fill);
    drawing.stack->Add(hist.get()); labels.push_back(style[2]); drawing.histograms.push_back(std::move(hist));
  }
  if(!unstacked) {
    for(int bin=1;bin<=drawing.histograms.front()->GetNbinsX();++bin) {
      double total=0,variance=0;
      for(const auto &hist:drawing.histograms) {
        total+=hist->GetBinContent(bin); double error=hist->GetBinError(bin); variance+=error*error;
      }
      check(std::isfinite(total) && std::isfinite(variance),"stack sums overflow");
    }
  }
}

void advanced_poly(const std::vector<std::vector<std::string>> &rows,size_t &at,
                   const AdvancedHeader &h,AdvancedDrawing &drawing) {
  check(h.dim==2 && h.count>=1 && h.count<=10000 && h.error=="-","invalid polygon header");
  struct Polygon { std::vector<double> x,y; double value; };
  std::vector<Polygon> polygons;
  int vertices=0;
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  for(int k=1;k<=h.count;++k) {
    check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="POLYGON" && integer(rows[at][1])==k,"invalid POLYGON");
    int n=integer(rows[at][2]); double value=real(rows[at][3]); ++at;
    vertices+=n; check(n>=3 && n<=256 && vertices<=100000,"invalid polygon vertex count");
    Polygon polygon{{},{},value};
    for(int i=1;i<=n;++i) {
      check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="VERTEX" && integer(rows[at][1])==i,"invalid VERTEX");
      double x=real(rows[at][2]),y=real(rows[at][3]); ++at;
      polygon.x.push_back(x); polygon.y.push_back(y);
      xmin=std::min(xmin,x); xmax=std::max(xmax,x); ymin=std::min(ymin,y); ymax=std::max(ymax,y);
    }
    // Simple polygons only: convex or non-convex, no holes or self-intersections.
    check_simple_polygon(polygon.x,polygon.y);
    polygons.push_back(std::move(polygon));
  }
  check(xmin<xmax && ymin<ymax && (!h.logx || xmin>0) && (!h.logy || ymin>0),"invalid polygon range");
  auto poly=std::make_unique<TH2Poly>("rootplots_polygon","",xmin,xmax,ymin,ymax);
  poly->SetDirectory(nullptr); poly->SetStats(h.stats);
  for(const auto &p:polygons) {
    int bin=poly->AddBin(p.x.size(),p.x.data(),p.y.data());
    check(bin>0,"ROOT rejected polygon"); poly->SetBinContent(bin,p.value);
  }
  drawing.histograms.push_back(std::move(poly));
}

AdvancedDrawing advanced_parse(const std::vector<std::vector<std::string>> &rows,
                               AdvancedHeader &h,std::vector<std::string> &labels) {
  size_t at=0; h=advanced_header(rows,at);
  AdvancedDrawing drawing;
  if(h.command=="DRAWSTACK") advanced_stack(rows,at,h,drawing,labels);
  else if(h.command=="DRAWPOLY") advanced_poly(rows,at,h,drawing);
  else advanced_profile(rows,at,h,drawing);
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" &&
        integer(rows[at][1])==h.count,"invalid advanced END/count or trailing rows");
  validate_scene(h.scene,h.logx,h.logy,drawing.histograms.size());
  return drawing;
}

AdvancedDrawing paint_advanced(const std::vector<std::vector<std::string>> &rows,TPad &pad) {
  check(!rows.empty() && advanced_command(rows[0][0]),"advanced request is not drawable");
  AdvancedHeader h; std::vector<std::string> labels;
  auto drawing=advanced_parse(rows,h,labels);
  pad.cd(); apply_canvas(pad,h.scene);
  gStyle->SetOptStat(h.stats ? (drawing.profile ? 1110 : 1100) : 0);
  pad.SetLogx(h.logx); pad.SetLogy(h.logy); pad.SetLogz(h.logz);
  pad.SetTheta(h.theta); pad.SetPhi(h.phi);
  if(h.dim>1) palette(h.palette_name);
  TH1 *frame=nullptr;
  if(drawing.stack) {
    if(h.explicit_y) { drawing.stack->SetMinimum(h.ymin); drawing.stack->SetMaximum(h.ymax); }
    drawing.stack->Draw(h.option.c_str()); pad.Modified(); pad.Update();
    frame=drawing.stack->GetHistogram(); check(frame!=nullptr,"ROOT did not construct stack axes");
  } else {
    frame=drawing.histograms.front().get();
    frame->SetContour(h.contours);
    if(h.explicit_y && h.dim==1) { frame->SetMinimum(h.ymin); frame->SetMaximum(h.ymax); }
    if(h.dim>=2) pad.SetRightMargin(0.18);
    frame->Draw(h.option.c_str());
  }
  frame->GetXaxis()->SetTitle(h.labels.x.c_str()); frame->GetYaxis()->SetTitle(h.labels.y.c_str());
  if(h.dim>=2) frame->GetZaxis()->SetTitle(h.labels.z.c_str());
  apply_axes(frame,h.scene);
  if(h.explicit_x) frame->GetXaxis()->SetRangeUser(h.xmin,h.xmax);
  std::vector<TObject*> objects;
  for(const auto &hist:drawing.histograms) objects.push_back(hist.get());
  if(labels.empty()) labels.resize(objects.size());
  drawing.decorations=draw_scene(h.scene,h.corner,objects,labels,
                                std::vector<std::string>(objects.size(),drawing.stack ? "f" : "lep"));
  if(!h.labels.title.empty()) {
    auto heading=std::make_unique<TLatex>(0.13,0.94,h.labels.title.c_str());
    heading->SetNDC(); heading->SetTextSize(0.035); heading->Draw();
    drawing.decorations.texts.push_back(std::move(heading));
  }
  pad.Modified(); pad.Update();
  return drawing;
}

template<class P> void advanced_profile_summary(P &p,const std::string &path) {
  std::ofstream out(path); check(bool(out),"cannot create profile summary");
  out<<"ROOTPlotsProfile\t1\n"<<std::setprecision(17);
  int index=0;
  int ny=p.GetDimension()>=2 ? p.GetNbinsY() : 1;
  int nz=p.GetDimension()>=3 ? p.GetNbinsZ() : 1;
  for(int iz=1;iz<=nz;++iz) for(int iy=1;iy<=ny;++iy) for(int ix=1;ix<=p.GetNbinsX();++ix) {
    int bin=p.GetBin(ix,iy,iz);
    out<<"BIN\t"<<++index<<'\t'<<p.GetBinContent(bin)<<'\t'<<p.GetBinError(bin)<<'\t'
       <<p.GetBinEntries(bin)<<'\t'<<p.GetBinEffectiveEntries(bin)<<'\n';
  }
  check(bool(out),"failed writing profile summary");
}

bool advanced_dispatch(const std::vector<std::vector<std::string>> &rows,const std::string &stem) {
  if(rows.empty() || rows[0].empty()) return false;
  if(rows[0][0]=="PROFILESTATS") {
    AdvancedHeader h; std::vector<std::string> labels;
    auto drawing=advanced_parse(rows,h,labels);
    if(auto *p=dynamic_cast<TProfile*>(drawing.profile)) advanced_profile_summary(*p,stem);
    else if(auto *p=dynamic_cast<TProfile2D*>(drawing.profile)) advanced_profile_summary(*p,stem);
    else if(auto *p=dynamic_cast<TProfile3D*>(drawing.profile)) advanced_profile_summary(*p,stem);
    else throw std::runtime_error("missing native profile");
    return true;
  }
  if(!advanced_command(rows[0][0])) return false;
  check(rows[0].size()>=3,"invalid advanced canvas header");
  int width=integer(rows[0][1]),height=integer(rows[0][2]);
  check(width>=1 && width<=10000 && height>=1 && height<=10000,"invalid advanced canvas size");
  TCanvas canvas("rootplots_advanced_canvas","",width,height);
  auto drawing=paint_advanced(rows,canvas);
  save(canvas,stem);
  return true;
}
