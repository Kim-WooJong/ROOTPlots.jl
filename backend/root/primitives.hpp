std::vector<std::unique_ptr<TObject>> paint_primitives(const std::vector<std::vector<std::string>>& rows,TPad& pad) {
  std::vector<std::unique_ptr<TObject>> objects;
  pad.cd();
  for(const auto &f:rows) {
    check(f.size()>=2 && f[0]=="PRIMITIVE","invalid primitive");
    std::unique_ptr<TObject> object;
    auto style_line=[&](TAttLine& line,const std::string &hex,int width) { color(hex); check(width>0,"invalid primitive width"); line.SetLineColor(TColor::GetColor(hex.c_str())); line.SetLineWidth(width); };
    auto style_fill=[&](TAttFill& fill,const std::string& hex,int style,double alpha) { color(hex); check(style>=0 && style<=9999 && alpha>=0 && alpha<=1,"invalid primitive fill"); fill.SetFillColorAlpha(TColor::GetColor(hex.c_str()),alpha); fill.SetFillStyle(style); };
    auto coordinate=[&](double x,double y) { check((!pad.GetLogx() || x>0) && (!pad.GetLogy() || y>0),"invalid log primitive coordinates"); };
    if(f[1]=="ARROW") {
      check(f.size()==10,"invalid ARROW");
      double x1=real(f[2]),y1=real(f[3]),x2=real(f[4]),y2=real(f[5]),size=real(f[8]);
      coordinate(x1,y1); coordinate(x2,y2); check(size>0 && size<=1,"invalid arrow size");
      check(f[9]==">" || f[9]=="<" || f[9]=="|>" || f[9]=="<|" || f[9]=="<>" || f[9]=="<|>" || f[9]=="-|","invalid arrow option");
      auto arrow=std::make_unique<TArrow>(x1,y1,x2,y2,size,f[9].c_str());
      style_line(*arrow,f[6],integer(f[7])); arrow->SetFillColor(TColor::GetColor(f[6].c_str())); arrow->Draw(); object=std::move(arrow);
    } else if(f[1]=="BOX") {
      check(f.size()==11,"invalid BOX");
      double x1=real(f[2]),y1=real(f[3]),x2=real(f[4]),y2=real(f[5]); coordinate(x1,y1); coordinate(x2,y2);
      check(x1<x2 && y1<y2,"invalid box bounds");
      auto box=std::make_unique<TBox>(x1,y1,x2,y2); style_line(*box,f[6],integer(f[7])); style_fill(*box,f[8],integer(f[9]),real(f[10])); box->Draw(); object=std::move(box);
    } else if(f[1]=="ELLIPSE") {
      check(f.size()==14,"invalid ELLIPSE");
      double x=real(f[2]),y=real(f[3]),rx=real(f[4]),ry=real(f[5]),a1=real(f[6]),a2=real(f[7]);
      check(rx>0 && ry>0 && a1<a2,"invalid ellipse bounds"); coordinate(x-rx,y-ry); coordinate(x+rx,y+ry);
      auto ellipse=std::make_unique<TEllipse>(x,y,rx,ry,a1,a2,real(f[8]));
      style_line(*ellipse,f[9],integer(f[10])); style_fill(*ellipse,f[11],integer(f[12]),real(f[13])); ellipse->Draw(); object=std::move(ellipse);
    } else if(f[1]=="AXIS") {
      check(f.size()==12,"invalid AXIS");
      double x1=real(f[2]),y1=real(f[3]),x2=real(f[4]),y2=real(f[5]),low=real(f[6]),high=real(f[7]);
      coordinate(x1,y1); coordinate(x2,y2); int divisions=integer(f[8]);
      check(low<high && divisions>=1 && divisions<=999999,"invalid secondary axis bounds");
      for(char c:f[11]) check(std::string("+-LURSMWN=><").find(c)!=std::string::npos,"unsupported secondary axis option");
      auto axis=std::make_unique<TGaxis>(x1,y1,x2,y2,low,high,divisions,f[11].c_str());
      style_line(*axis,f[10],1); axis->SetLabelColor(TColor::GetColor(f[10].c_str())); axis->SetTitle(f[9].c_str()); axis->Draw(); object=std::move(axis);
    } else throw std::runtime_error("unknown primitive kind");
    objects.push_back(std::move(object));
  }
  return objects;
}
