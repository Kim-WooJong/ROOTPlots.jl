// Included inside graphfit.cpp's anonymous namespace, after the reusable painters.
struct LayoutDrawing {
  std::vector<std::unique_ptr<TObject>> objects;
  SceneObjects decorations;
};
std::vector<double> layout_edges(const std::vector<std::vector<std::string>>& rows,
                                 size_t &at,const std::string &tag,int n) {
  check(n>=1 && n<=2048,"invalid panel bin count");
  std::vector<double> edges;
  for(int i=1;i<=n+1;++i) {
    check(at<rows.size() && rows[at].size()==3 && rows[at][0]==tag && integer(rows[at][1])==i,"invalid panel edge");
    double value=real(rows[at++][2]);
    check(edges.empty() || value>edges.back(),"panel edges must increase");
    edges.push_back(value);
  }
  return edges;
}
LayoutDrawing paint_layout_basic(const std::vector<std::vector<std::string>>& rows,TPad &pad) {
  size_t at=0;
  check(!rows.empty() && rows[0].size()==13 && rows[0][0]=="PLOT","invalid PLOT");
  const auto &header=rows[at++];
  int dimension=integer(header[3]),nlayers=integer(header[11]);
  check(dimension>=1 && dimension<=3 && nlayers>=1 && nlayers<=64 && (dimension==1 || nlayers==1),"invalid panel layers/dimension");
  bool logx=bit(header[4]),logy=bit(header[5]),logz=bit(header[6]);
  double theta=real(header[7]),phi=real(header[8]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid panel view");
  int contours=integer(header[10]); check(contours>=1 && contours<=255,"invalid panel contours");
  if(dimension>1) palette(header[9]);
  std::string corner=header[12];
  check(corner=="none" || corner=="topleft" || corner=="topright" || corner=="bottomleft" || corner=="bottomright","invalid panel legend");
  auto labels=read_labels(rows,at);
  check(at<rows.size() && rows[at].size()==10 && rows[at][0]=="LIMITS","missing panel limits");
  auto limits=rows[at++];
  bool explicit_z=bit(limits[7]);
  double xmin=real(limits[2]),xmax=real(limits[3]),ymin=real(limits[5]),ymax=real(limits[6]),
         zmin=real(limits[8]),zmax=real(limits[9]);
  check(bit(limits[1]) && bit(limits[4]) && xmin<xmax && ymin<ymax && zmin<zmax &&
        (!logx || xmin>0) && (!logy || ymin>0) && (!logz || !explicit_z || zmin>0),"invalid panel limits");
  auto scene=read_scene(rows,at);
  validate_scene(scene,logx,logy,nlayers);
  check(at<rows.size() && rows[at].size()==11 && rows[at][0]=="ZSTYLE","invalid panel z axis style");
  auto zstyle=rows[at++];
  std::vector<std::string> zaxis(zstyle.begin()+1,zstyle.end());
  check(at<rows.size() && rows[at].size()>=2 && rows[at][0]=="LEVELS","missing panel contour levels");
  auto levels_row=rows[at++]; int nlevels=integer(levels_row[1]);
  check(nlevels>=0 && nlevels<=255 && levels_row.size()==static_cast<size_t>(nlevels+2),"invalid contour levels");
  std::vector<double> levels;
  for(int i=0;i<nlevels;++i) { double v=real(levels_row[2+i]); check(levels.empty() || v>levels.back(),"unordered contour levels"); levels.push_back(v); }
  pad.cd(); apply_canvas(pad,scene);
  pad.SetLogx(logx); pad.SetLogy(logy); pad.SetLogz(logz); pad.SetTheta(theta); pad.SetPhi(phi);
  LayoutDrawing result;
  std::vector<TObject*> objects;
  std::vector<std::string> object_labels,legend_options;
  TH1 *frame=nullptr;
  if(dimension==1) frame=pad.DrawFrame(xmin,ymin,xmax,ymax);
  for(int layer=0;layer<nlayers;++layer) {
    check(at<rows.size() && rows[at].size()>=3 && rows[at][0]=="LAYER","missing panel layer");
    auto item=rows[at++];
    if(item[1]=="FUNCTION") {
      check(dimension==1 && item.size()==3 && at<rows.size(),"invalid function layer");
      const auto &f=rows[at++];
      check(f.size()>=10 && f[0]=="FUNCTION","missing layer function");
      int npar=integer(f[5]),npx=integer(f[4]); double lo=real(f[2]),hi=real(f[3]);
      check(npar==nparameters(f[1]) && f.size()==static_cast<size_t>(10+npar) && lo<hi && npx>=2 && npx<=1000000,"invalid layer function");
      auto curve=std::make_unique<TF1>(("rootplots_layer_function_"+std::to_string(layer)).c_str(),formula(f[1]).c_str(),lo,hi);
      check(curve->IsValid(),"ROOT rejected layer function");
      for(int i=0;i<npar;++i) curve->SetParameter(i,real(f[10+i]));
      if(f[1]=="gaus") check(curve->GetParameter(2)!=0,"invalid Gaussian width");
      if(f[1]=="erfc_edge") check(curve->GetParameter(3)>0,"invalid edge width");
      if(f[1]=="double_erf") check(curve->GetParameter(4)>0 && curve->GetParameter(5)>0,"invalid edge widths");
      color(f[7]); int width=integer(f[8]),line=integer(f[9]); check(width>0 && line>=1 && line<=10,"invalid function style");
      curve->SetNpx(npx); curve->SetLineColor(TColor::GetColor(f[7].c_str())); curve->SetLineWidth(width); curve->SetLineStyle(line);
      curve->Draw("SAME"); objects.push_back(curve.get()); object_labels.push_back(item[2]); legend_options.push_back("l");
      result.objects.push_back(std::move(curve)); continue;
    }
    if(item[1]=="GRAPH") {
      check(dimension==1 && item.size()==3,"invalid graph layer");
      auto series=read_series(rows,at,1);
      validate_series_log(series,logx,logy);
      auto graph=make_graph(series);
      graph->Draw(draw_option(series).c_str());
      objects.push_back(graph.get()); object_labels.push_back(series.label); legend_options.push_back(legend_option(series));
      result.objects.push_back(std::move(graph));
      continue;
    }
    check(item.size()==4 && (item[1]=="HIST1" || item[1]=="HIST2" || item[1]=="HIST3"),"invalid histogram layer");
    int ndim=item[1]=="HIST1" ? 1 : item[1]=="HIST2" ? 2 : 3;
    check(ndim==dimension,"inconsistent histogram dimension");
    check_root_option(item[3]);
    check(at<rows.size() && rows[at].size()==9 && rows[at][0]=="STYLE","invalid panel histogram style");
    auto style=rows[at++]; color(style[1]); color(style[6]);
    int marker=integer(style[2]),width=integer(style[4]),line_style=integer(style[5]),fill_style=integer(style[7]);
    double marker_size=real(style[3]),alpha=real(style[8]);
    check(marker>=1 && marker<=50 && marker_size>0 && width>0 && line_style>=1 && line_style<=10 &&
          fill_style>=0 && fill_style<=9999 && alpha>=0 && alpha<=1,"invalid panel histogram style");
    int nx=0,ny=1,nz=1; bool errors=true;
    if(ndim==1) {
      check(at<rows.size() && rows[at].size()==2 && rows[at][0]=="HIST","missing panel HIST");
      nx=integer(rows[at++][1]);
    } else {
      check(at<rows.size() && rows[at].size()==static_cast<size_t>(ndim+2) && rows[at][0]=="GRID","missing panel GRID");
      auto grid=rows[at++]; nx=integer(grid[1]); ny=integer(grid[2]);
      if(ndim==3) nz=integer(grid[3]);
      errors=bit(grid.back());
    }
    check(nx>=1 && ny>=1 && nz>=1 && static_cast<long long>(nx)*ny*nz<=1000000,"invalid panel histogram size");
    auto xe=layout_edges(rows,at,ndim==1 ? "EDGE" : "XEDGE",nx);
    std::vector<double> ye,ze;
    if(ndim>1) ye=layout_edges(rows,at,"YEDGE",ny);
    if(ndim>2) ze=layout_edges(rows,at,"ZEDGE",nz);
    check((!logx || xe.front()>0) && (ndim<2 || !logy || ye.front()>0) && (ndim<3 || !logz || ze.front()>0),"log axes need positive bin edges");
    std::string name="rootplots_panel_hist_"+std::to_string(layer);
    std::unique_ptr<TH1> hist;
    if(ndim==1) hist=std::make_unique<TH1D>(name.c_str(),"",nx,xe.data());
    else if(ndim==2) hist=std::make_unique<TH2D>(name.c_str(),"",nx,xe.data(),ny,ye.data());
    else hist=std::make_unique<TH3D>(name.c_str(),"",nx,xe.data(),ny,ye.data(),nz,ze.data());
    hist->SetDirectory(nullptr); hist->SetStats(false); if(errors) hist->Sumw2();
    for(int iz=1;iz<=nz;++iz) for(int iy=1;iy<=ny;++iy) for(int ix=1;ix<=nx;++ix) {
      check(at<rows.size(),"missing panel bin");
      auto bin=rows[at++];
      check(bin.size()==static_cast<size_t>(ndim+3) && bin[0]==(ndim==1 ? "BIN" : "CELL") && integer(bin[1])==ix &&
            (ndim<2 || integer(bin[2])==iy) && (ndim<3 || integer(bin[3])==iz),"invalid panel bin index");
      double value=real(bin[ndim+1]),variance=real(bin[ndim+2]); check(variance>=0,"negative panel variance");
      int global=hist->GetBin(ix,ndim>1 ? iy : 0,ndim>2 ? iz : 0);
      hist->SetBinContent(global,value); if(errors) hist->SetBinError(global,std::sqrt(variance));
    }
    int c=TColor::GetColor(style[1].c_str());
    hist->SetLineColor(c); hist->SetLineWidth(width); hist->SetLineStyle(line_style);
    hist->SetMarkerColor(c); hist->SetMarkerStyle(marker); hist->SetMarkerSize(marker_size);
    hist->SetFillColorAlpha(TColor::GetColor(style[6].c_str()),alpha); hist->SetFillStyle(fill_style);
    hist->SetContour(nlevels ? nlevels : contours,nlevels ? levels.data() : nullptr);
    std::string option=item[3];
    if(ndim==1) option+=" SAME";
    else {
      hist->GetXaxis()->SetRangeUser(xmin,xmax); hist->GetYaxis()->SetRangeUser(ymin,ymax);
      if(explicit_z) {
        if(ndim==2) { hist->SetMinimum(zmin); hist->SetMaximum(zmax); }
        else hist->GetZaxis()->SetRangeUser(zmin,zmax);
      }
      if(!explicit_z && hist->GetMinimum()==0 && hist->GetMaximum()==0) { hist->SetMinimum(-1); hist->SetMaximum(1); }
      frame=hist.get();
    }
    hist->Draw(option.c_str());
    objects.push_back(hist.get()); object_labels.push_back(item[2]);
    legend_options.push_back(root_option_uses_fill(item[3]) ? "f" : "lep");
    result.objects.push_back(std::move(hist));
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && integer(rows[at][1])==nlayers,"invalid panel END");
  check(frame!=nullptr,"missing panel frame");
  frame->GetXaxis()->SetTitle(labels.x.c_str()); frame->GetYaxis()->SetTitle(labels.y.c_str());
  frame->GetZaxis()->SetTitle(labels.z.c_str()); apply_axes(frame,scene); apply_axis(frame->GetZaxis(),zaxis);
  result.decorations=draw_scene(scene,corner,objects,object_labels,legend_options);
  if(!labels.title.empty()) {
    auto title=std::make_unique<TLatex>(0.13,0.94,labels.title.c_str());
    title->SetNDC(); title->SetTextSize(0.035); title->Draw(); result.objects.push_back(std::move(title));
  }
  return result;
}

void draw_layout(const std::vector<std::vector<std::string>>& rows,const std::string &stem) {
  check(!rows.empty() && rows[0].size()==6 && rows[0][0]=="DRAWLAYOUT","invalid DRAWLAYOUT");
  auto header=rows[0];
  int width=integer(header[1]),height=integer(header[2]),nr=integer(header[3]),nc=integer(header[4]),count=integer(header[5]);
  check(width>=100 && width<=10000 && height>=100 && height<=10000 && nr>=1 && nr<=16 && nc>=1 && nc<=16 &&
        count>=1 && count<=64 && nr*nc>=count && nr*nc<=64,"invalid layout dimensions");
  std::vector<std::vector<std::vector<std::string>>> panels;
  std::vector<std::vector<std::vector<std::string>>> primitive_rows;
  size_t at=1;
  for(int i=1;i<=count;++i) {
    check(at<rows.size() && rows[at].size()==4 && rows[at][0]=="PANEL" && integer(rows[at][1])==i,"invalid PANEL");
    int n=integer(rows[at][2]),nprimitives=integer(rows[at++][3]);
    check(n>0 && n<=5000000 && nprimitives>=0 && nprimitives<=1000 && static_cast<size_t>(n+nprimitives)<=rows.size()-at,"invalid panel record count");
    panels.emplace_back(rows.begin()+at,rows.begin()+at+n); at+=n;
    primitive_rows.emplace_back(rows.begin()+at,rows.begin()+at+nprimitives); at+=nprimitives;
  }
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && integer(rows[at][1])==count,"invalid layout END");
  // Native objects outlive all pads and final export; destroy canvas first.
  std::vector<LayoutDrawing> basics;
  std::vector<AdvancedDrawing> advanced;
  std::vector<FunctionDrawing> functions;
  std::vector<XYZDrawing> xyz;
  std::vector<SpecialDrawing> special;
  std::vector<std::vector<std::unique_ptr<TObject>>> primitives;
  TCanvas canvas("rootplots_layout","",width,height);
  canvas.Divide(nc,nr,0.001,0.001);
  const int base_frame_width=gStyle->GetFrameLineWidth();
  for(int i=0;i<count;++i) {
    gStyle->SetFrameLineWidth(base_frame_width);
    auto *pad=dynamic_cast<TPad*>(canvas.cd(i+1)); check(pad!=nullptr,"missing layout pad");
    const auto &payload=panels[i];
    check(!payload.empty() && !payload[0].empty(),"empty panel");
    auto command=payload[0][0];
    if(command=="PLOT") basics.push_back(paint_layout_basic(payload,*pad));
    else if(advanced_command(command)) advanced.push_back(paint_advanced(payload,*pad));
    else if(command=="DRAWFUNCTION") functions.push_back(paint_function(payload,*pad));
    else if(command=="DRAWXYZ") xyz.push_back(paint_xyz(payload,*pad));
    else if(special_command(command)) special.push_back(paint_special(payload,*pad));
    else throw std::runtime_error("unsupported layout panel command");
    primitives.push_back(paint_primitives(primitive_rows[i],*pad));
    pad->Modified(); pad->Update();
  }
  gStyle->SetFrameLineWidth(base_frame_width);
  canvas.cd(); save(canvas,stem);
}
