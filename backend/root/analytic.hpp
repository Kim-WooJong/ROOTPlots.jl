// Included inside graphfit.cpp's anonymous namespace, before layouts.hpp.
// Native TF2/TF3 from fixed formulas; the Julia side only selects a model name.
std::string analytic_formula(const std::string &model) {
  if(model=="gaus2d") return "[0]*exp(-0.5*((x-[1])/[2])^2-0.5*((y-[3])/[4])^2)";
  if(model=="bigaus") return "[0]*exp(-0.5/(1-[5]*[5])*(((x-[1])/[2])^2+((y-[3])/[4])^2-2*[5]*((x-[1])/[2])*((y-[3])/[4])))";
  if(model=="pol1_2d") return "[0]+[1]*x+[2]*y";
  if(model=="pol2_2d") return "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y";
  if(model=="gaus3d") return "[0]*exp(-0.5*((x-[1])/[2])^2-0.5*((y-[3])/[4])^2-0.5*((z-[5])/[6])^2)";
  if(model=="pol1_3d") return "[0]+[1]*x+[2]*y+[3]*z";
  if(model=="pol2_3d") return "[0]+[1]*x+[2]*y+[3]*z+[4]*x*x+[5]*y*y+[6]*z*z+[7]*x*y+[8]*x*z+[9]*y*z";
  throw std::runtime_error("unsupported analytic model");
}
int analytic_npar(const std::string &model) {
  if(model=="gaus2d") return 5;
  if(model=="bigaus" || model=="pol2_2d") return 6;
  if(model=="pol1_2d") return 3;
  if(model=="gaus3d") return 7;
  if(model=="pol1_3d") return 4;
  if(model=="pol2_3d") return 10;
  throw std::runtime_error("unsupported analytic model");
}
void check_analytic_parameters(const std::string &model,const std::vector<double> &p) {
  if(model=="gaus2d" || model=="bigaus") check(p[2]>0 && p[4]>0,"Gaussian widths must be positive");
  if(model=="bigaus") check(std::abs(p[5])<1,"bigaus correlation must satisfy |rho| < 1");
  if(model=="gaus3d") check(p[2]>0 && p[4]>0 && p[6]>0,"Gaussian widths must be positive");
}
bool tf2_option_token(std::string t) {
  std::transform(t.begin(),t.end(),t.begin(),[](unsigned char c){return std::toupper(c);});
  auto rest=[&](const std::string &prefix,const std::string &digits) {
    if(t.rfind(prefix,0)!=0) return false;
    std::string r=t.substr(prefix.size());
    if(!r.empty() && digits.find(r[0])!=std::string::npos) r=r.substr(1);
    return r.empty() || r=="Z";
  };
  return t=="TEXT" || t=="Z" || rest("COL","012") || rest("CONT","01234") || rest("SURF","1234567") || rest("LEGO","1234");
}
struct AnalyticDrawing {
  std::unique_ptr<TF1> function;
  std::unique_ptr<TH3D> sampled;
  std::unique_ptr<TLatex> heading;
};
AnalyticDrawing paint_tf2(const std::vector<std::vector<std::string>> &rows,TPad &pad) {
  size_t at=0;
  check(!rows.empty() && rows[0].size()==11 && rows[0][0]=="DRAWTF2","invalid DRAWTF2");
  const auto f=rows[at++];
  std::string option=f[3];
  check_root_option(option);
  {
    std::istringstream tokens(option); std::string token; bool any=false;
    while(tokens>>token) { check(tf2_option_token(token),"unsupported TF2 option"); any=true; }
    check(any,"empty TF2 option");
  }
  palette(f[4]);
  int contours=integer(f[5]); check(contours>=1 && contours<=255,"invalid TF2 contours");
  double theta=real(f[6]),phi=real(f[7]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid TF2 view");
  bool logx=bit(f[8]),logy=bit(f[9]),logz=bit(f[10]);
  auto labels=read_labels(rows,at);
  auto scene=read_scene(rows,at);
  check(at<rows.size() && rows[at].size()>=12 && rows[at][0]=="TF2","missing TF2");
  const auto v=rows[at++];
  std::string model=v[1];
  int npar=analytic_npar(model);
  check(model!="gaus3d" && model!="pol1_3d" && model!="pol2_3d","TF2 needs a two-dimensional model");
  check(integer(v[8])==npar && v.size()==static_cast<size_t>(12+npar),"invalid TF2 parameter count");
  double xlo=real(v[2]),xhi=real(v[3]),ylo=real(v[4]),yhi=real(v[5]);
  int npx=integer(v[6]),npy=integer(v[7]);
  check(xlo<xhi && ylo<yhi && npx>=4 && npx<=1000 && npy>=4 && npy<=1000,"invalid TF2 range/resolution");
  check((!logx || xlo>0) && (!logy || ylo>0),"log axes need positive TF2 ranges");
  color(v[9]); int width=integer(v[10]),line_style=integer(v[11]);
  check(width>0 && line_style>=1 && line_style<=10,"invalid TF2 line style");
  std::vector<double> p; for(int i=0;i<npar;++i) p.push_back(real(v[12+i]));
  check_analytic_parameters(model,p);
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && rows[at][1]=="1","invalid TF2 END");
  AnalyticDrawing out;
  auto tf2=std::make_unique<TF2>("rootplots_tf2",analytic_formula(model).c_str(),xlo,xhi,ylo,yhi);
  check(tf2->IsValid() && tf2->GetNpar()==npar,"ROOT rejected TF2 model");
  for(int i=0;i<npar;++i) tf2->SetParameter(i,p[i]);
  tf2->SetNpx(npx); tf2->SetNpy(npy);
  for(int i=0;i<npx;++i) for(int j=0;j<npy;++j) {
    double x=xlo+(xhi-xlo)*i/(npx-1.0),y=ylo+(yhi-ylo)*j/(npy-1.0),value=tf2->Eval(x,y);
    check(std::isfinite(value) && (!logz || value>0),"TF2 values must be finite (and positive on log z)");
  }
  tf2->SetLineColor(TColor::GetColor(v[9].c_str())); tf2->SetLineWidth(width); tf2->SetLineStyle(line_style);
  std::string upper=option;
  std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return std::toupper(c);});
  bool three_d=upper.find("SURF")!=std::string::npos || upper.find("LEGO")!=std::string::npos;
  // Same color-level policy as Grid2D: explicit count for contours, fine levels otherwise.
  tf2->SetContour(upper.find("CONT")!=std::string::npos ? contours : three_d ? 80 : 255);
  pad.cd(); pad.SetLogx(logx); pad.SetLogy(logy); pad.SetLogz(logz);
  // A selected ROOT style keeps its margins on flat views; a palette still gets room.
  const bool style_margins=style_selected() && !three_d;
  if(!style_margins) { pad.SetLeftMargin(0.14); pad.SetBottomMargin(0.13); }
  if(three_d) { pad.SetRightMargin(0.22); pad.SetTheta(theta); pad.SetPhi(phi); }
  else if(upper.find('Z')!=std::string::npos) pad.SetRightMargin(0.18);
  else if(!style_margins) pad.SetRightMargin(0.08);
  apply_canvas(pad,scene); // explicit CanvasStyle values win over the layout defaults
  tf2->Draw(option.c_str());
  pad.Update();
  TH1 *h=tf2->GetHistogram();
  if(h) {
    h->GetXaxis()->SetTitle(labels.x.c_str()); h->GetYaxis()->SetTitle(labels.y.c_str());
    h->GetZaxis()->SetTitle(labels.z.c_str());
    if(three_d) { h->GetXaxis()->SetTitleOffset(1.5); h->GetYaxis()->SetTitleOffset(1.5); h->GetZaxis()->SetTitleOffset(1.3); }
    apply_axes(h,scene);
    if(!three_d) palette_title_room(pad,h,flat_palette(option,2),scene.canvas,scene.axes[2].v[AX_TITLE_OFFSET]=="-");
  }
  out.heading=std::make_unique<TLatex>();
  out.heading->SetNDC(); out.heading->SetTextSize(0.035);
  if(root_title(pad,h,labels.title,scene.canvas)) tf2->SetTitle(labels.title.c_str());
  else if(!labels.title.empty()) out.heading->DrawLatex(0.13,0.94,labels.title.c_str());
  out.function=std::move(tf2);
  return out;
}
AnalyticDrawing paint_tf3(const std::vector<std::vector<std::string>> &rows,TPad &pad) {
  size_t at=0;
  check(!rows.empty() && rows[0].size()==7 && rows[0][0]=="DRAWTF3","invalid DRAWTF3");
  const auto f=rows[at++];
  std::string mode=f[3];
  check(mode=="iso" || mode=="box","TF3 draw must be iso or box");
  palette(f[4]);
  double theta=real(f[5]),phi=real(f[6]);
  check(theta>=0 && theta<=180 && phi>=-180 && phi<=180,"invalid TF3 view");
  auto labels=read_labels(rows,at);
  auto scene=read_scene(rows,at);
  check(at<rows.size() && rows[at].size()>=17 && rows[at][0]=="TF3","missing TF3");
  const auto v=rows[at++];
  std::string model=v[1];
  int npar=analytic_npar(model);
  check(model=="gaus3d" || model=="pol1_3d" || model=="pol2_3d","TF3 needs a three-dimensional model");
  check(integer(v[12])==npar && v.size()==static_cast<size_t>(16+npar),"invalid TF3 parameter count");
  double lo[3],hi[3]; int np[3];
  for(int k=0;k<3;++k) { lo[k]=real(v[2+2*k]); hi[k]=real(v[3+2*k]); np[k]=integer(v[8+k]);
    check(lo[k]<hi[k] && np[k]>=2 && np[k]<=100,"invalid TF3 range/sampling"); }
  double level=real(v[11]);
  color(v[13]); int width=integer(v[14]),line_style=integer(v[15]);
  check(width>0 && line_style>=1 && line_style<=10,"invalid TF3 line style");
  std::vector<double> p; for(int i=0;i<npar;++i) p.push_back(real(v[16+i]));
  check_analytic_parameters(model,p);
  check(at+1==rows.size() && rows[at].size()==2 && rows[at][0]=="END" && rows[at][1]=="1","invalid TF3 END");
  AnalyticDrawing out;
  // Parameter npar holds the iso level, so the implicit surface f - level = 0 is f = level.
  std::string formula="("+analytic_formula(model)+")-["+std::to_string(npar)+"]";
  auto tf3=std::make_unique<TF3>("rootplots_tf3",formula.c_str(),lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
  check(tf3->IsValid() && tf3->GetNpar()==npar+1,"ROOT rejected TF3 model");
  for(int i=0;i<npar;++i) tf3->SetParameter(i,p[i]);
  tf3->SetParameter(npar,level);
  int c=TColor::GetColor(v[13].c_str());
  tf3->SetLineColor(c); tf3->SetLineWidth(width); tf3->SetLineStyle(line_style); tf3->SetFillColor(c);
  pad.cd(); pad.SetLeftMargin(0.14); pad.SetBottomMargin(0.13); pad.SetTheta(theta); pad.SetPhi(phi);
  TH1 *frame=nullptr;
  if(mode=="iso") {
    pad.SetRightMargin(0.08);
    apply_canvas(pad,scene);
    tf3->Draw("");
    pad.Update();
    if(auto *h=tf3->GetHistogram()) {
      h->GetXaxis()->SetTitle(labels.x.c_str()); h->GetYaxis()->SetTitle(labels.y.c_str()); h->GetZaxis()->SetTitle(labels.z.c_str());
      h->GetXaxis()->SetTitleOffset(1.5); h->GetYaxis()->SetTitleOffset(1.5); h->GetZaxis()->SetTitleOffset(1.3);
      apply_axes(h,scene); frame=h;
    }
  } else {
    pad.SetRightMargin(0.22);
    apply_canvas(pad,scene);
    out.sampled=std::make_unique<TH3D>("rootplots_tf3_sampled","",np[0],lo[0],hi[0],np[1],lo[1],hi[1],np[2],lo[2],hi[2]);
    out.sampled->SetDirectory(nullptr); out.sampled->SetStats(false);
    for(int i=1;i<=np[0];++i) for(int j=1;j<=np[1];++j) for(int k=1;k<=np[2];++k) {
      double x=out.sampled->GetXaxis()->GetBinCenter(i),y=out.sampled->GetYaxis()->GetBinCenter(j),
             z=out.sampled->GetZaxis()->GetBinCenter(k);
      // Evaluate the model itself (level parameter temporarily zero).
      double value=tf3->Eval(x,y,z)+level;
      check(std::isfinite(value),"TF3 values must be finite");
      out.sampled->SetBinContent(i,j,k,value);
    }
    out.sampled->GetXaxis()->SetTitle(labels.x.c_str()); out.sampled->GetYaxis()->SetTitle(labels.y.c_str());
    out.sampled->GetZaxis()->SetTitle(labels.z.c_str());
    out.sampled->GetXaxis()->SetTitleOffset(1.5); out.sampled->GetYaxis()->SetTitleOffset(1.5);
    out.sampled->GetZaxis()->SetTitleOffset(1.3);
    // Bins at or below the iso level are not drawn, so both modes show the region f > level.
    out.sampled->SetMinimum(level);
    apply_axes(out.sampled.get(),scene); frame=out.sampled.get();
    out.sampled->Draw("BOX2Z");
  }
  out.heading=std::make_unique<TLatex>();
  out.heading->SetNDC(); out.heading->SetTextSize(0.035);
  if(root_title(pad,frame,labels.title,scene.canvas)) tf3->SetTitle(labels.title.c_str());
  else if(!labels.title.empty()) out.heading->DrawLatex(0.13,0.94,labels.title.c_str());
  out.function=std::move(tf3);
  return out;
}
