#include <iostream>
#include <vector>

#include "cgs.h"
#include "chi2.h"
#include "git_revision.h"
#include "input.h"
#include "logging.h"
#include "output.h"
#include "particle.h"
#include "utilities.h"

int main(int argc, char* argv[]) {
  log_startup_information();
  try {
    if (argc == 2) {
      CRAMS::Input input;
      input.setSimname(argv[1]);
      input.readParamsFromFile(argv[1]);
      input.print();

      CRAMS::ParticleList particleList;
      particleList.readParamsFromFile(argv[1]);
      particleList.print();

      CRAMS::Particles particles;
      auto list = particleList.getList();
      particles.reserve(list.size());
      for (auto it = list.rbegin(); it != list.rend(); ++it) {
        auto pid = it->first;
        auto nucleusParams = it->second;
        particles.emplace_back(pid, nucleusParams);
      }

      for (auto& particle : particles) {
        LOGI << "running : " << particle.getPid();
        particle.buildVectors(input);
        particle.buildGrammage(input);
        particle.buildLosses(input);
        particle.buildPrimarySource(input);
        particle.buildInelasticXsecs(input);
        particle.buildSecondarySource(input, particles);
        if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(particles);
        if (particle.getPid() == CRAMS::pbar) particle.buildAntiprotonSource(particles);
        // if (input.X_s > 0.) particle.buildGrammageAtSource(input, particles);
        particle.dump();
        particle.computeIntensity(input);
        particle.reset();
      }

      CRAMS::OutputManager outputManager(particles, input);
      outputManager.dumpSpectra();
      outputManager.dumpSpectraEkn();

      {
        std::ofstream fchi2("chi2_results.txt", std::ofstream::out);
        const std::pair<double, double> R_range = {10. * CRAMS::CGS::GeV, 1000. * CRAMS::CGS::GeV};
        std::vector<std::unique_ptr<CRAMS::Chi2>> chi2s;

        // ====================================================================
        // OPTION 1: Manual specification (flexible, for custom isotope selection)
        // ====================================================================
        // Example 1: Single element model for Hydrogen
        // chi2s.push_back(
        //    CRAMS::makeElementModel("H (manual)", particles, input.modulationPotential,
        //                   {CRAMS::H1, CRAMS::H2, CRAMS::H1_ter},
        //                   CRAMS::findDataFile("H"),
        //                   CRAMS::CGS::GeV,
        //                   1.0 / (CRAMS::CGS::GeV * CRAMS::CGS::m2 * CRAMS::CGS::sec * CRAMS::CGS::sr)));

        // ====================================================================
        // OPTION 2: Automatic isotope detection (simpler, for most use cases)
        // ====================================================================
        // Create element models by specifying only the atomic number (Z)
        // Automatically finds and includes all available isotopes
        
        // Single elements
        chi2s.push_back(CRAMS::get_Chi2("H", particles, input.modulationPotential, 1));
        chi2s.push_back(CRAMS::get_Chi2("He", particles, input.modulationPotential, 2));
        chi2s.push_back(CRAMS::get_Chi2("Be", particles, input.modulationPotential, 4));
        chi2s.push_back(CRAMS::get_Chi2("B", particles, input.modulationPotential, 5));
        chi2s.push_back(CRAMS::get_Chi2("C", particles, input.modulationPotential, 6));
        chi2s.push_back(CRAMS::get_Chi2("N", particles, input.modulationPotential, 7));
        chi2s.push_back(CRAMS::get_Chi2("O", particles, input.modulationPotential, 8));
        chi2s.push_back(CRAMS::get_Chi2("Ne", particles, input.modulationPotential, 10));
        chi2s.push_back(CRAMS::get_Chi2("Mg", particles, input.modulationPotential, 12));
        chi2s.push_back(CRAMS::get_Chi2("Si", particles, input.modulationPotential, 14));
        chi2s.push_back(CRAMS::get_Chi2("S", particles, input.modulationPotential, 16));
        chi2s.push_back(CRAMS::get_Chi2("Fe", particles, input.modulationPotential, 26));
        
        // Ratios
        chi2s.push_back(CRAMS::get_Chi2_ratio("HeO", particles, input.modulationPotential, 2, 8));
        chi2s.push_back(CRAMS::get_Chi2_ratio("BeB", particles, input.modulationPotential, 4, 5));
        chi2s.push_back(CRAMS::get_Chi2_ratio("BeC", particles, input.modulationPotential, 4, 6));
        chi2s.push_back(CRAMS::get_Chi2_ratio("BeO", particles, input.modulationPotential, 4, 8));
        chi2s.push_back(CRAMS::get_Chi2_ratio("BC", particles, input.modulationPotential, 5, 6));
        chi2s.push_back(CRAMS::get_Chi2_ratio("BO", particles, input.modulationPotential, 5, 8));
        chi2s.push_back(CRAMS::get_Chi2_ratio("CO", particles, input.modulationPotential, 6, 8));
        chi2s.push_back(CRAMS::get_Chi2_ratio("NeMg", particles, input.modulationPotential, 10, 12));
        chi2s.push_back(CRAMS::get_Chi2_ratio("SiMg", particles, input.modulationPotential, 14, 12));
        chi2s.push_back(CRAMS::get_Chi2_ratio("Be10/Be9", particles, input.modulationPotential, 
                                       {CRAMS::Be10}, 
                                       {CRAMS::Be9},
                                       CRAMS::findDataFile("Be10_Be9")));

        
        // DAMPE data ratios (kinetic energy per nucleon)
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_ratio("BC_DAMPE", particles, input.modulationPotential, 5, 6));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_ratio("BO_DAMPE", particles, input.modulationPotential, 5, 8));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_flux("H_DAMPE", particles, input.modulationPotential, 1));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_flux("He_DAMPE", particles, input.modulationPotential, 2));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_flux("C_DAMPE", particles, input.modulationPotential, 6));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_flux("O_DAMPE", particles, input.modulationPotential, 8));
        chi2s.push_back(CRAMS::get_Chi2_DAMPE_flux("Fe_DAMPE", particles, input.modulationPotential, 26));

        // Compute and print chi2 values for each model
        LOGI << "Computing chi2 values...";
        for (const auto& chi2model : chi2s) {
          try {
            double chi2_value = chi2model->computeChi2(R_range.first, R_range.second);
            LOGI << chi2model->getName() << ": " << chi2_value;
            fchi2 << chi2model->getName() << ", " << chi2_value << "\n";
          } catch (const std::exception& e) {
            LOGW << "Error computing chi2 for " << chi2model->getName() << ": " << e.what();
          }
        }
        fchi2.close();
      }
    } else {
      throw std::runtime_error("please provide an input file as './crams params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
