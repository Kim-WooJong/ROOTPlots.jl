// Output formats of one native rendering (shared by rootplots-graphfit and
// rootplots-render). The optional last command-line argument
// "--formats=png", "--formats=pdf" or "--formats=png,pdf" selects them; the
// legacy invocation without it writes both files. Every image route saves
// through save_outputs, so an unrequested format is never written.
#pragma once
#include <TCanvas.h>
#include <stdexcept>
#include <string>

struct OutputFormats { bool png=true, pdf=true; };

inline OutputFormats &output_formats() { static OutputFormats formats; return formats; }

// Parse "--formats=..."; anything else is a usage error (exit code 2 in main).
inline void set_output_formats(const std::string &argument) {
  const std::string prefix="--formats=";
  if(argument.compare(0,prefix.size(),prefix)!=0)
    throw std::runtime_error("unknown option (expected --formats=png, --formats=pdf or --formats=png,pdf)");
  const std::string value=argument.substr(prefix.size());
  OutputFormats formats{false,false};
  if(value=="png") formats.png=true;
  else if(value=="pdf") formats.pdf=true;
  else if(value=="png,pdf" || value=="pdf,png") formats.png=formats.pdf=true;
  else throw std::runtime_error("unsupported output formats (expected png, pdf or png,pdf)");
  output_formats()=formats;
}

// PNG is written before PDF, as in the legacy invocation.
inline void save_outputs(TCanvas &canvas,const std::string &stem) {
  if(output_formats().png) canvas.SaveAs((stem+".png").c_str());
  if(output_formats().pdf) canvas.SaveAs((stem+".pdf").c_str());
}
