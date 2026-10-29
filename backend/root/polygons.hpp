// Included inside graphfit.cpp's anonymous namespace, before advanced_histograms.hpp.
// Simple-polygon validation shared by TH2Poly bins and polygon primitives.
// Mirrors the Julia checks: convex or non-convex polygons are accepted, holes
// and self-intersections are not.
double polygon_orient(double ax,double ay,double bx,double by,double cx,double cy) {
  return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
}
bool polygon_on_segment(double ax,double ay,double bx,double by,double px,double py) {
  return std::min(ax,bx)<=px && px<=std::max(ax,bx) && std::min(ay,by)<=py && py<=std::max(ay,by);
}
int polygon_sign(double v) { return (v>0)-(v<0); }
void check_simple_polygon(const std::vector<double> &x,const std::vector<double> &y) {
  const size_t n=x.size();
  check(n==y.size() && n>=3 && n<=256,"polygon needs 3:256 vertices");
  for(size_t i=0;i<n;++i) check(std::isfinite(x[i]) && std::isfinite(y[i]),"polygon vertices must be finite");
  double area=0;
  for(size_t i=0;i<n;++i) {
    size_t j=(i+1)%n;
    check(x[i]!=x[j] || y[i]!=y[j],"duplicate polygon neighbors");
    area+=x[i]*y[j]-x[j]*y[i];
  }
  check(std::isfinite(area),"polygon geometry overflows");
  check(area!=0,"polygon has zero area");
  auto intersects=[&](size_t a,size_t b,size_t c,size_t d) {
    double o1=polygon_orient(x[a],y[a],x[b],y[b],x[c],y[c]),o2=polygon_orient(x[a],y[a],x[b],y[b],x[d],y[d]);
    double o3=polygon_orient(x[c],y[c],x[d],y[d],x[a],y[a]),o4=polygon_orient(x[c],y[c],x[d],y[d],x[b],y[b]);
    check(std::isfinite(o1) && std::isfinite(o2) && std::isfinite(o3) && std::isfinite(o4),"polygon geometry overflows");
    if(polygon_sign(o1)*polygon_sign(o2)<0 && polygon_sign(o3)*polygon_sign(o4)<0) return true;
    if(o1==0 && polygon_on_segment(x[a],y[a],x[b],y[b],x[c],y[c])) return true;
    if(o2==0 && polygon_on_segment(x[a],y[a],x[b],y[b],x[d],y[d])) return true;
    if(o3==0 && polygon_on_segment(x[c],y[c],x[d],y[d],x[a],y[a])) return true;
    if(o4==0 && polygon_on_segment(x[c],y[c],x[d],y[d],x[b],y[b])) return true;
    return false;
  };
  for(size_t i=0;i<n;++i) {
    size_t a=i,b=(i+1)%n,c=(i+2)%n;
    bool folds=polygon_orient(x[a],y[a],x[b],y[b],x[c],y[c])==0 &&
               (x[b]-x[a])*(x[c]-x[b])+(y[b]-y[a])*(y[c]-y[b])<0;
    check(!folds,"polygon is self-intersecting");
    for(size_t j=i+2;j<n;++j) {
      if(i==0 && j==n-1) continue;
      check(!intersects(i,(i+1)%n,j,(j+1)%n),"polygon is self-intersecting");
    }
  }
}
