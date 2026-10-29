#include <TCanvas.h>
#include <TBox.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>


#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Rendering-only ROOT backend. Julia supplies all bin, error, and ratio numbers.
// Objects stay alive in Drawn until SaveAs finishes; canvases/pads are destroyed
// before Drawn. Histograms are detached from TDirectory. No .root file is written.
namespace {
struct Bin {
  int index;
  double left, right, data, data_lo, data_hi, pred, pred_lo, pred_hi;
  double ratio, ratio_lo, ratio_hi, band_lo, band_hi;
  bool valid;
};
struct Component {
  std::string label, color;
  std::vector<double> values;
};
struct Panel {
  int id, nbins;
  std::string xlabel, ylabel, annotation, lumi, uncertainty_label;
  std::vector<Bin> bins;
  std::vector<Component> components;
  bool logx = false, logy = false, ratio_logy = false, axis_seen = false;
  double ymin = 0, ymax = 0, ratio_min = 0.5, ratio_max = 1.5;
};
struct Figure {
  double ratio_min = 0.5, ratio_max = 1.5;
  std::string marker_color = "#000000", band_color = "#000000";
  int marker_style = 20, line_width = 1, stack_fill_style = 1001;
  int band_fill_style = 3345, legend_columns = 2;
  double marker_size = 0.8;
  int version = 1;
  bool settings_seen = false, style_seen = false;
  std::vector<Panel> panels;
};
struct Drawn {
  std::vector<std::unique_ptr<TH1D>> hists;
  std::unique_ptr<THStack> stack;
  std::unique_ptr<TGraphAsymmErrors> data, top_band, ratio_data, ratio_band;
  std::unique_ptr<TH1D> ratio_axis;
  std::unique_ptr<TLegend> legend;
  std::unique_ptr<TLine> unity;
  std::vector<std::unique_ptr<TLine>> clipped_errors;
  std::vector<std::unique_ptr<TBox>> clipped_bands;
  std::vector<std::unique_ptr<TBox>> clipped_ratio_bands;
};

std::vector<std::string> fields(const std::string &line) {
  std::vector<std::string> out;
  std::stringstream in(line);
  std::string field;
  while (std::getline(in, field, '\t')) out.push_back(field);
  return out;
}
void check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}
double number(const std::string &text) {
  char *end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  check(end != text.c_str() && end == text.c_str() + text.size() && std::isfinite(value),
        "invalid finite number: " + text);
  return value;
}
int integer(const std::string &text) {
  size_t read = 0;
  const int value = std::stoi(text, &read);
  check(read == text.size(), "invalid integer: " + text);
  return value;
}
bool boolean(const std::string &text) {
  check(text == "0" || text == "1", "invalid boolean: " + text);
  return text == "1";
}
bool hexcolor(const std::string &text) {
  if (text.size() != 7 || text[0] != '#') return false;
  return std::all_of(text.begin() + 1, text.end(),
                     [](unsigned char c) { return std::isxdigit(c) != 0; });
}
Figure parse(const std::string &path) {
  std::ifstream input(path);
  check(bool(input), "cannot open plot specification: " + path);
  std::string line;
  check(bool(std::getline(input, line)) &&
        (line == "ROOTPlots\t1" || line == "ROOTPlots\t2"), "unsupported plot schema");
  Figure figure;
  figure.version = line == "ROOTPlots\t2" ? 2 : 1;
  while (std::getline(input, line)) {
    const auto f = fields(line);
    check(!f.empty(), "blank specification row");
    if (f[0] == "SETTINGS") {
      check(f.size() == 3 && !figure.settings_seen && figure.panels.empty(),
            "malformed or duplicate SETTINGS row");
      figure.settings_seen = true;
      figure.ratio_min = number(f[1]);
      figure.ratio_max = number(f[2]);
      check(figure.ratio_min < figure.ratio_max, "invalid ratio range");
    } else if (f[0] == "STYLE") {
      check(figure.version == 2 && f.size() == 9 && !figure.style_seen && figure.panels.empty(),
            "malformed or duplicate STYLE row");
      figure.style_seen = true;
      figure.marker_style = integer(f[1]); figure.marker_size = number(f[2]);
      figure.marker_color = f[3]; figure.line_width = integer(f[4]);
      figure.stack_fill_style = integer(f[5]); figure.band_fill_style = integer(f[6]);
      figure.band_color = f[7]; figure.legend_columns = integer(f[8]);
      check(figure.marker_style >= 1 && figure.marker_style <= 50 && figure.marker_size > 0 &&
            figure.line_width > 0 && figure.legend_columns > 0 &&
            hexcolor(figure.marker_color) && hexcolor(figure.band_color),
            "invalid plot style");
      check((figure.stack_fill_style == 0 || figure.stack_fill_style == 1001 ||
             figure.stack_fill_style == 3345) &&
            (figure.band_fill_style == 0 || figure.band_fill_style == 1001 ||
             figure.band_fill_style == 3345), "unsupported fill style");
    } else if (f[0] == "PANEL") {
      check(f.size() == 8, "malformed PANEL row");
      Panel p{integer(f[1]), integer(f[7]), f[2], f[3], f[4], f[5], f[6], {}, {}};
      check(p.id == int(figure.panels.size()) + 1 && p.nbins > 0,
            "invalid panel id or bin count");
      p.ratio_min = figure.ratio_min; p.ratio_max = figure.ratio_max;
      figure.panels.push_back(std::move(p));
    } else if (f[0] == "AXIS") {
      check(figure.version == 2 && f.size() == 9 && !figure.panels.empty(),
            "malformed AXIS row");
      auto &p = figure.panels.back();
      check(integer(f[1]) == p.id && !p.axis_seen && p.bins.empty(),
            "invalid or duplicate AXIS row");
      p.axis_seen = true;
      p.logx = boolean(f[2]); p.logy = boolean(f[3]); p.ratio_logy = boolean(f[4]);
      p.ymin = number(f[5]); p.ymax = number(f[6]);
      p.ratio_min = number(f[7]); p.ratio_max = number(f[8]);
      check(p.ymin < p.ymax && p.ratio_min < p.ratio_max &&
            (!p.logy || p.ymin > 0) && (!p.ratio_logy || p.ratio_min > 0),
            "invalid axis limits");
    } else if (f[0] == "BIN") {
      check(f.size() == 17 && !figure.panels.empty(), "malformed BIN row");
      auto &p = figure.panels.back();
      check(integer(f[1]) == p.id, "BIN panel id disagrees");
      Bin b{integer(f[2]), number(f[3]), number(f[4]), number(f[5]),
            number(f[6]), number(f[7]), number(f[8]), number(f[9]),
            number(f[10]), number(f[11]), number(f[12]), number(f[13]),
            number(f[14]), number(f[15]), integer(f[16]) != 0};
      check(b.index == int(p.bins.size()) + 1 && b.left < b.right,
            "invalid bin index or edges");
      if (!p.bins.empty()) check(p.bins.back().right == b.left, "noncontiguous bins");
      p.bins.push_back(b);
    } else if (f[0] == "COMP") {
      check(!figure.panels.empty(), "COMP before PANEL");
      auto &p = figure.panels.back();
      check(f.size() == size_t(5 + p.nbins), "malformed COMP row");
      check(integer(f[1]) == p.id && integer(f[2]) == int(p.components.size()) + 1,
            "invalid component ordering");
      Component c{f[3], f[4], {}};
      for (int i = 0; i < p.nbins; ++i) {
        double value = number(f[5 + i]);
        check(value >= 0, "negative stack component");
        c.values.push_back(value);
      }
      p.components.push_back(std::move(c));
    } else {
      throw std::runtime_error("unknown specification row: " + f[0]);
    }
  }
  check(!figure.panels.empty(), "no panels");
  check(figure.settings_seen && (figure.version == 1 || figure.style_seen),
        "missing figure settings");
  for (const auto &p : figure.panels) {
    check(int(p.bins.size()) == p.nbins && !p.components.empty() &&
          (figure.version == 1 || p.axis_seen), "incomplete panel");
    if (p.logx) for (const auto &b : p.bins)
      check(b.left > 0 && b.right > 0, "logx needs positive bin edges");
  }
  return figure;
}

void draw_panel(const Panel &p, int col, int count, const Figure &figure,
                TCanvas &canvas, Drawn &d) {
  canvas.cd();
  const double x0 = double(col) / count, x1 = double(col + 1) / count;
  auto *top = new TPad(("top" + std::to_string(col)).c_str(), "", x0, 0.28, x1, 0.98);
  auto *bottom = new TPad(("bottom" + std::to_string(col)).c_str(), "", x0, 0.05, x1, 0.28);
  canvas.cd(); top->Draw(); bottom->Draw();
  top->SetTicks(1, 1); bottom->SetTicks(1, 1);
  top->SetLogx(p.logx); bottom->SetLogx(p.logx);
  top->SetLogy(p.logy); bottom->SetLogy(p.ratio_logy);
  top->SetLeftMargin(0.18); top->SetRightMargin(0.04);
  top->SetBottomMargin(0.012); top->SetTopMargin(0.20);
  bottom->SetLeftMargin(0.18); bottom->SetRightMargin(0.04);
  bottom->SetTopMargin(0.02); bottom->SetBottomMargin(0.35);

  std::vector<double> edges;
  edges.reserve(p.nbins + 1);
  for (const auto &b : p.bins) edges.push_back(b.left);
  edges.push_back(p.bins.back().right);
  d.stack = std::make_unique<THStack>(("stack" + std::to_string(col)).c_str(), "");
  for (size_t j = 0; j < p.components.size(); ++j) {
    const auto &c = p.components[j];
    auto h = std::make_unique<TH1D>(("h" + std::to_string(col) + "_" + std::to_string(j)).c_str(),
                                    "", p.nbins, edges.data());
    h->SetDirectory(nullptr);
    for (int i = 1; i <= p.nbins; ++i) h->SetBinContent(i, c.values[i-1]);
    h->SetFillColor(TColor::GetColor(c.color.c_str()));
    h->SetFillStyle(figure.stack_fill_style);
    h->SetLineColor(kBlack); h->SetLineWidth(figure.line_width);
    d.stack->Add(h.get(), "HIST");
    d.hists.push_back(std::move(h));
  }

  d.top_band = std::make_unique<TGraphAsymmErrors>();
  d.ratio_band = std::make_unique<TGraphAsymmErrors>();
  d.data = std::make_unique<TGraphAsymmErrors>();
  d.ratio_data = std::make_unique<TGraphAsymmErrors>();
  double ymax = 0;
  for (const auto &b : p.bins) {
    const double mid = (b.left + b.right) / 2, half = (b.right - b.left) / 2;
    ymax = std::max(ymax, std::max(b.pred + b.pred_hi, b.data + b.data_hi));
    if (!p.logy || b.pred > p.ymin) {
      const int t = d.top_band->GetN();
      d.top_band->SetPoint(t, mid, b.pred);
      d.top_band->SetPointError(t, half, half,
                                p.logy ? std::min(b.pred_lo, b.pred - p.ymin) : b.pred_lo,
                                b.pred_hi);
    } else if (b.pred + b.pred_hi > p.ymin) {
      auto box = std::make_unique<TBox>(b.left, p.ymin, b.right,
                                        std::min(p.ymax, b.pred + b.pred_hi));
      box->SetFillStyle(figure.band_fill_style);
      box->SetFillColor(TColor::GetColor(figure.band_color.c_str()));
      box->SetLineColor(TColor::GetColor(figure.band_color.c_str()));
      d.clipped_bands.push_back(std::move(box));
    }
    if (!p.logy || b.data > p.ymin) {
      const int q = d.data->GetN();
      d.data->SetPoint(q, mid, b.data);
      d.data->SetPointError(q, 0, 0,
                             p.logy ? std::min(b.data_lo, b.data - p.ymin) : b.data_lo,
                             b.data_hi);
    } else if (b.data + b.data_hi > p.ymin) {
      auto line = std::make_unique<TLine>(mid, p.ymin, mid,
                                          std::min(p.ymax, b.data + b.data_hi));
      line->SetLineWidth(figure.line_width);
      line->SetLineColor(TColor::GetColor(figure.marker_color.c_str()));
      d.clipped_errors.push_back(std::move(line));
    }
    if (b.valid) {
      if (!p.ratio_logy || 1 > p.ratio_min) {
        const int r = d.ratio_band->GetN();
        d.ratio_band->SetPoint(r, mid, 1);
        d.ratio_band->SetPointError(r, half, half,
                                    p.ratio_logy ? std::min(b.band_lo, 1 - p.ratio_min) : b.band_lo,
                                    b.band_hi);
      } else if (1 + b.band_hi > p.ratio_min) {
        auto box = std::make_unique<TBox>(b.left, p.ratio_min, b.right,
                                          std::min(p.ratio_max, 1 + b.band_hi));
        box->SetFillStyle(figure.band_fill_style);
        box->SetFillColor(TColor::GetColor(figure.band_color.c_str()));
        box->SetLineColor(TColor::GetColor(figure.band_color.c_str()));
        d.clipped_ratio_bands.push_back(std::move(box));
      }
      if (!p.ratio_logy || b.ratio > p.ratio_min) {
        const int s = d.ratio_data->GetN();
        d.ratio_data->SetPoint(s, mid, b.ratio);
        d.ratio_data->SetPointError(s, 0, 0,
                                    p.ratio_logy ? std::min(b.ratio_lo, b.ratio - p.ratio_min) : b.ratio_lo,
                                    b.ratio_hi);
      }
    }
  }
  for (auto *band : {d.top_band.get(), d.ratio_band.get()}) {
    band->SetFillStyle(figure.band_fill_style);
    band->SetFillColor(TColor::GetColor(figure.band_color.c_str()));
    band->SetLineColor(TColor::GetColor(figure.band_color.c_str()));
    band->SetLineWidth(figure.line_width);
  }
  for (auto *points : {d.data.get(), d.ratio_data.get()}) {
    points->SetMarkerStyle(figure.marker_style); points->SetMarkerSize(figure.marker_size);
    points->SetMarkerColor(TColor::GetColor(figure.marker_color.c_str()));
    points->SetLineColor(TColor::GetColor(figure.marker_color.c_str()));
    points->SetLineWidth(figure.line_width);
  }
  top->cd();
  d.stack->SetMinimum(p.axis_seen ? p.ymin : 0);
  d.stack->SetMaximum(p.axis_seen ? p.ymax : std::max(1.0, 1.46 * ymax));
  d.stack->Draw("HIST");
  d.stack->GetXaxis()->SetLabelSize(0);
  d.stack->GetXaxis()->SetTitleSize(0);
  d.stack->GetYaxis()->SetTitle(p.ylabel.c_str());
  d.stack->GetYaxis()->SetTitleOffset(1.55);
  d.stack->GetYaxis()->SetLabelSize(0.038);
  d.top_band->Draw("2 SAME");
  for (auto &box : d.clipped_bands) box->Draw("SAME");
  d.data->Draw("P E1 SAME");
  for (auto &line : d.clipped_errors) line->Draw("SAME");
  d.legend = std::make_unique<TLegend>(0.47, 0.80, 0.96, 0.94);
  d.legend->SetNColumns(figure.legend_columns); d.legend->SetBorderSize(0); d.legend->SetFillStyle(0);
  d.legend->SetTextSize(0.027);
  for (size_t j = 0; j < p.components.size(); ++j)
    d.legend->AddEntry(d.hists[j].get(), p.components[j].label.c_str(), "f");
  d.legend->AddEntry(d.top_band.get(), p.uncertainty_label.c_str(), "f");
  d.legend->AddEntry(d.data.get(), "Data", "pe");
  d.legend->Draw();
  TLatex text; text.SetNDC();
  // Caller-supplied free text only; no experiment label is added.
  text.SetTextFont(42); text.SetTextSize(0.034);
  if (!p.annotation.empty()) text.DrawLatex(0.18, 0.958, p.annotation.c_str());
  if (!p.lumi.empty()) {
    text.SetTextAlign(31); text.SetTextSize(0.029);
    text.DrawLatex(0.96, 0.966, p.lumi.c_str());
  }

  bottom->cd();
  d.ratio_axis = std::make_unique<TH1D>(("ratio" + std::to_string(col)).c_str(), "",
                                         p.nbins, edges.data());
  d.ratio_axis->SetDirectory(nullptr);
  d.ratio_axis->SetMinimum(p.ratio_min);
  d.ratio_axis->SetMaximum(p.ratio_max);
  d.ratio_axis->GetYaxis()->SetTitle("Data / pred.");
  d.ratio_axis->GetYaxis()->SetTitleSize(0.10);
  d.ratio_axis->GetYaxis()->SetTitleOffset(0.63);
  d.ratio_axis->GetYaxis()->SetLabelSize(0.09);
  d.ratio_axis->GetYaxis()->SetNdivisions(503);
  d.ratio_axis->GetXaxis()->SetTitle(p.xlabel.c_str());
  d.ratio_axis->GetXaxis()->SetTitleSize(0.12);
  d.ratio_axis->GetXaxis()->SetLabelSize(0.10);
  d.ratio_axis->GetXaxis()->SetTitleOffset(1.05);
  d.ratio_axis->Draw("AXIS");
  d.unity = std::make_unique<TLine>(edges.front(), 1, edges.back(), 1);
  d.unity->SetLineWidth(figure.line_width); d.unity->Draw("SAME");
  d.ratio_band->Draw("2 SAME");
  for (auto &box : d.clipped_ratio_bands) box->Draw("SAME");
  d.ratio_data->Draw("P E1 SAME");
  if (p.ratio_logy) for (const auto &b : p.bins) {
    if (!b.valid || b.ratio > p.ratio_min || b.ratio + b.ratio_hi <= p.ratio_min) continue;
    const double mid = (b.left + b.right) / 2;
    auto line = std::make_unique<TLine>(mid, p.ratio_min, mid,
                                        std::min(p.ratio_max, b.ratio + b.ratio_hi));
    line->SetLineWidth(figure.line_width);
    line->SetLineColor(TColor::GetColor(figure.marker_color.c_str()));
    line->Draw("SAME");
    d.clipped_errors.push_back(std::move(line));
  }
}
} // namespace

int main(int argc, char **argv) {
  try {
    check(argc == 3, "usage: rootplots-render spec.tsv output-stem");
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetEndErrorSize(4); gStyle->SetErrorX(0);
    auto figure = parse(argv[1]);
    // Declare drawings before canvas, so canvas/pads release references first.
    std::vector<std::unique_ptr<Drawn>> drawings;
    TCanvas canvas("rootplots", "ROOTPlots control distributions",
                   600 * int(figure.panels.size()), 720);
    for (size_t i = 0; i < figure.panels.size(); ++i) {
      auto drawing = std::make_unique<Drawn>();
      draw_panel(figure.panels[i], int(i), int(figure.panels.size()), figure, canvas, *drawing);
      drawings.push_back(std::move(drawing));
    }
    canvas.cd(); canvas.Modified(); canvas.Update();
    canvas.SaveAs((std::string(argv[2]) + ".png").c_str());
    canvas.SaveAs((std::string(argv[2]) + ".pdf").c_str());
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ROOTPlots ROOT backend: " << e.what() << '\n';
    return 1;
  }
}
