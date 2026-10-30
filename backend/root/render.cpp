#include <TAxis.h>
#include <TAxisModLab.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph2D.h>
#include <THashList.h>
#include <TPaveText.h>
#include <TText.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TError.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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
// Tab-separated fields; empty fields (including a trailing one) are kept.
std::vector<std::string> fields(const std::string &line) {
  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    const auto pos = line.find('\t', start);
    out.push_back(line.substr(start, pos == std::string::npos ? pos : pos - start));
    if (pos == std::string::npos) break;
    start = pos + 1;
  }
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
// Helper names used by the shared axis header.
double real(const std::string &text) { return number(text); }
bool bit(const std::string &text) { return boolean(text); }
void color(const std::string &text) { check(hexcolor(text), "invalid color: " + text); }
#include "axisstyle.hpp"
#include "experiment.hpp"
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
  // Optional PanelSpec axis styles (PCANVAS/PAXIS rows): shared x, upper y, ratio y.
  std::vector<std::string> canvas = std::vector<std::string>(CV_FIELDS, "-");
  AxisSpec xaxis, yaxis, ratio_axis;
  AxisSpec *current_axis = nullptr;
  std::string current_axis_name;
  // Optional CMS/ATLAS label of the upper pad (EXPLABEL row).
  ExperimentLabelSpec label;
};
struct Figure {
  double ratio_min = 0.5, ratio_max = 1.5;
  std::string marker_color = "#000000", band_color = "#000000";
  int marker_style = 20, line_width = 1, stack_fill_style = 1001;
  int band_fill_style = 3345, legend_columns = 2;
  double marker_size = 0.8;
  int version = 1;
  bool settings_seen = false, style_seen = false;
  std::string root_style; // ROOTSTYLE row: gROOT->SetStyle before the canvas exists
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
    } else if (f[0] == "ROOTSTYLE") {
      check(figure.version == 2 && f.size() == 2 && figure.root_style.empty() && figure.panels.empty(),
            "malformed or duplicate ROOTSTYLE row");
      figure.root_style = f[1];
    } else if (f[0] == "EXPLABEL") {
      check(!figure.panels.empty(), "EXPLABEL before PANEL");
      auto &p = figure.panels.back();
      check(!p.label.active && p.bins.empty(), "invalid or duplicate EXPLABEL row");
      p.label = read_experiment_label(f);
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
    } else if (f[0] == "PCANVAS") {
      check(f.size() == 2 + CV_FIELDS && !figure.panels.empty(), "malformed PCANVAS row");
      auto &p = figure.panels.back();
      check(integer(f[1]) == p.id && p.bins.empty(), "invalid PCANVAS row");
      p.canvas.assign(f.begin() + 2, f.end());
      canvas_title_fields(p.canvas);
      for (int i = 0; i < CV_FIELDS; ++i)
        check(p.canvas[i] == "-" || i == CV_TICKX || i == CV_TICKY || i == CV_GRIDX || i == CV_GRIDY,
              "comparison panels support grids and ticks only");
      for (int i : {CV_TICKX, CV_TICKY}) if (p.canvas[i] != "-") boolean(p.canvas[i]);
    } else if (f[0] == "PAXIS") {
      check(f.size() >= 3 && !figure.panels.empty(), "malformed PAXIS row");
      auto &p = figure.panels.back();
      check(integer(f[1]) == p.id && p.bins.empty(), "invalid PAXIS row");
      p.current_axis = f[2] == "x" ? &p.xaxis : f[2] == "y" ? &p.yaxis : f[2] == "ratio" ? &p.ratio_axis : nullptr;
      check(p.current_axis != nullptr, "invalid PAXIS axis");
      p.current_axis_name = f[2];
      axis_scalars(*p.current_axis, f, 3);
    } else if (f[0] == "AXTIME" || f[0] == "AXBIN" || f[0] == "AXCHANGE") {
      check(!figure.panels.empty() && figure.panels.back().current_axis, "axis detail row before PAXIS");
      auto &p = figure.panels.back();
      check(axis_extra(*p.current_axis, f, p.current_axis_name), "axis detail row for another axis");
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
  // With a selected ROOT style or an experiment label, the pads span the full
  // canvas height and take the style's margins, ticks and text sizes rescaled
  // to each pad's height (style_scaled_axis); the legend moves into the frame.
  const bool styled = style_selected() || p.label.active;
  const double split = styled ? 0.30 : 0.28, low = styled ? 0.0 : 0.05, high = styled ? 1.0 : 0.98;
  const double ft = high - split, fb = split - low;
  auto *top = new TPad(("top" + std::to_string(col)).c_str(), "", x0, split, x1, high);
  auto *bottom = new TPad(("bottom" + std::to_string(col)).c_str(), "", x0, low, x1, split);
  canvas.cd(); top->Draw(); bottom->Draw();
  if (!styled || !style_selected()) { top->SetTicks(1, 1); bottom->SetTicks(1, 1); }
  for (auto *pad : {top, bottom}) {
    if (p.canvas[CV_TICKX] != "-") pad->SetTickx(boolean(p.canvas[CV_TICKX]));
    if (p.canvas[CV_TICKY] != "-") pad->SetTicky(boolean(p.canvas[CV_TICKY]));
    apply_pad_grid(*pad, p.canvas);
  }
  top->SetLogx(p.logx); bottom->SetLogx(p.logx);
  top->SetLogy(p.logy); bottom->SetLogy(p.ratio_logy);
  if (styled) {
    if (!style_selected()) for (auto *pad : {top, bottom}) { pad->SetLeftMargin(0.18); pad->SetRightMargin(0.04); }
    top->SetTopMargin(gStyle->GetPadTopMargin() / ft); top->SetBottomMargin(0.02);
    bottom->SetTopMargin(0.03); bottom->SetBottomMargin(gStyle->GetPadBottomMargin() / fb);
  } else {
    top->SetLeftMargin(0.18); top->SetRightMargin(0.04);
    top->SetBottomMargin(0.012); top->SetTopMargin(0.20);
    bottom->SetLeftMargin(0.18); bottom->SetRightMargin(0.04);
    bottom->SetTopMargin(0.02); bottom->SetBottomMargin(0.35);
  }
  // The label (drawn last) may grow the top margin; register it before the
  // frame and legend are laid out from that margin.
  if (p.label.active) { mark_two_pad(top); register_label(top, p.label, p.id); }

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
  if (styled) {
    style_scaled_axis(d.stack->GetXaxis(), "X", top, true);
    style_scaled_axis(d.stack->GetYaxis(), "Y", top, false);
    // Event counts reach thousands; on the narrow comparison pad four-digit
    // labels would run into the y title, so ROOT shows an x10^n exponent instead
    // (yaxis=AxisStyle(max_digits=...) overrides).
    d.stack->GetYaxis()->SetMaxDigits(3);
  } else {
    d.stack->GetYaxis()->SetTitleOffset(1.55);
    d.stack->GetYaxis()->SetLabelSize(0.038);
  }
  // The shared x style also applies to the upper pad, whose labels stay hidden.
  apply_axis_spec(d.stack->GetXaxis(), p.xaxis, "x");
  d.stack->GetXaxis()->SetLabelSize(0); d.stack->GetXaxis()->SetTitleSize(0);
  apply_axis_spec(d.stack->GetYaxis(), p.yaxis, "y");
  d.top_band->Draw("2 SAME");
  for (auto &box : d.clipped_bands) box->Draw("SAME");
  d.data->Draw("P E1 SAME");
  for (auto &line : d.clipped_errors) line->Draw("SAME");
  int legend_columns = figure.legend_columns;
  if (styled) {
    // Inside the frame, top right, sized from the labels measured by ROOT; it
    // keeps clear of a label at the top left (fewer columns, then smaller text).
    const double l = top->GetLeftMargin(), r = top->GetRightMargin(), t = top->GetTopMargin();
    const double right = 1 - r - 0.03;
    const double left = p.label.active && p.label.position == "inside_left" ? l + 0.42 * (1 - l - r) : l + 0.03;
    double size = 0.032 * two_pad_text(top).k, widest = 0;
    std::vector<std::string> labels;
    for (const auto &c : p.components) labels.push_back(c.label);
    labels.push_back(p.uncertainty_label); labels.push_back("Data");
    for (const auto &label : labels) widest = std::max(widest, text_width_ndc(top, label, 42, size));
    const double column = (widest + 0.02) / 0.75; // TLegend symbol margin 0.25
    while (legend_columns > 1 && legend_columns * column > right - left) --legend_columns;
    const double width = std::min(legend_columns * column, right - left);
    if (legend_columns * column > width) size *= width / (legend_columns * column);
    const int rows = (int(labels.size()) + legend_columns - 1) / legend_columns;
    const double y2 = 1 - t - 0.03, y1 = y2 - rows * 1.35 * size;
    d.legend = std::make_unique<TLegend>(right - width, y1, right, y2);
    d.legend->SetTextSize(size); d.legend->SetTextFont(42); d.legend->SetMargin(0.25);
  } else {
    d.legend = std::make_unique<TLegend>(0.47, 0.80, 0.96, 0.94);
    d.legend->SetTextSize(0.027);
  }
  d.legend->SetNColumns(legend_columns); d.legend->SetBorderSize(0); d.legend->SetFillStyle(0);
  for (size_t j = 0; j < p.components.size(); ++j)
    d.legend->AddEntry(d.hists[j].get(), p.components[j].label.c_str(), "f");
  d.legend->AddEntry(d.top_band.get(), p.uncertainty_label.c_str(), "f");
  d.legend->AddEntry(d.data.get(), "Data", "pe");
  d.legend->Draw();
  TLatex text; text.SetNDC();
  if (p.label.active) {
    // The label (registered above) carries the annotation and the luminosity text.
  } else if (styled) {
    // Caller-supplied free text above the frame: annotation left, luminosity right.
    const double l = top->GetLeftMargin(), r = top->GetRightMargin(), t = top->GetTopMargin();
    const PadText g = two_pad_text(top);
    const double tm = std::min(t * g.hr, 0.1), y = 1 - t + 0.2 * tm / g.hr;
    text.SetTextFont(42); text.SetTextSize(0.6 * tm * g.k);
    if (!p.annotation.empty()) { text.SetTextAlign(11); text.DrawLatex(l, y, p.annotation.c_str()); }
    if (!p.lumi.empty()) { text.SetTextAlign(31); text.DrawLatex(1 - r, y, p.lumi.c_str()); }
  } else {
    // Caller-supplied free text only; no experiment label is added.
    text.SetTextFont(42); text.SetTextSize(0.034);
    if (!p.annotation.empty()) text.DrawLatex(0.18, 0.958, p.annotation.c_str());
    if (!p.lumi.empty()) {
      text.SetTextAlign(31); text.SetTextSize(0.029);
      text.DrawLatex(0.96, 0.966, p.lumi.c_str());
    }
  }

  bottom->cd();
  d.ratio_axis = std::make_unique<TH1D>(("ratio" + std::to_string(col)).c_str(), "",
                                         p.nbins, edges.data());
  d.ratio_axis->SetDirectory(nullptr);
  d.ratio_axis->SetMinimum(p.ratio_min);
  d.ratio_axis->SetMaximum(p.ratio_max);
  d.ratio_axis->GetYaxis()->SetTitle("Data / pred.");
  d.ratio_axis->GetXaxis()->SetTitle(p.xlabel.c_str());
  // The ratio pad's text takes the upper y axis's pixel sizes, and the two y
  // titles line up (finish_two_pad_axes).
  if (styled) style_scaled_axis(d.ratio_axis->GetXaxis(), "X", bottom, true);
  d.ratio_axis->GetYaxis()->SetNdivisions(503);
  apply_axis_spec(d.ratio_axis->GetXaxis(), p.xaxis, "x");
  apply_axis_spec(d.ratio_axis->GetYaxis(), p.ratio_axis, "y");
  finish_two_pad_axes(top, d.stack->GetYaxis(), bottom, d.ratio_axis->GetXaxis(), p.xaxis,
                      d.ratio_axis->GetYaxis(), p.ratio_axis, styled);
  bottom->cd();
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
// ROOT prints an Info message for every saved file. Show warnings and errors
// only, unless ROOTPLOTS_ROOT_VERBOSE=1 restores ROOT's default message level.
void quiet_root() {
  const char *verbose = std::getenv("ROOTPLOTS_ROOT_VERBOSE");
  if (!(verbose && std::string(verbose) == "1")) gErrorIgnoreLevel = kWarning;
}
} // namespace

int main(int argc, char **argv) {
  try {
    check(argc == 3, "usage: rootplots-render spec.tsv output-stem");
    gROOT->SetBatch(kTRUE);
    quiet_root();
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    auto figure = parse(argv[1]);
    // A selected style keeps its own error-bar settings.
    if (!figure.root_style.empty()) select_root_style(figure.root_style);
    else { gStyle->SetEndErrorSize(4); gStyle->SetErrorX(0); }
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
    finish_canvas(canvas);
    canvas.SaveAs((std::string(argv[2]) + ".png").c_str());
    canvas.SaveAs((std::string(argv[2]) + ".pdf").c_str());
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ROOTPlots ROOT backend: " << e.what() << '\n';
    return 1;
  }
}
