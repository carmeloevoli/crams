#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "crams.h"

int main(int argc, char* argv[]) {
  bool quiet = false;
  std::string inifile;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "-q" || arg == "--quiet") {
      quiet = true;
    } else if (inifile.empty()) {
      inifile = arg;
    } else {
      throw std::runtime_error("unexpected argument: '" + arg + "'");
    }
  }

  log_startup_information(quiet);

  try {
    CRAMS::Input input;
    CRAMS::ParticleList injection;

    if (!inifile.empty()) {
      input.setSimname(inifile);
      input.readParamsFromFile(inifile);
      if (!quiet) input.print();
      injection.readParamsFromFile(inifile);
      if (!quiet) injection.print();
    } else {
      LOGI << "no input file provided, using default parameters";
    }

    CRAMS::Runner runner(input.inelasticModel(), input.fragmentationModel());
    runner.compute(injection, input, true, !quiet);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
