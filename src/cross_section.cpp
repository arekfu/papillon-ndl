/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est *
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <memory>

namespace pndl {

CrossSection::CrossSection(const ACE& ace, size_t i,
                           std::shared_ptr<EnergyGrid> E_grid, bool get_index)
    : energy_grid_(E_grid), values_(), index_(0) {
  uint32_t NE = ace.nxs(2);
  if (get_index) {
    index_ = ace.xss<uint32_t>(i) - 1;
    i++;
    NE = ace.xss<uint32_t>(i);
    i++;
  }

  values_ = ace.xss<float>(i, NE);

  if (energy_grid_->size() - index_ != values_.size()) {
    std::string mssg =
        "CrossSection::CrossSection: Different number of points in the energy "
        "grid and xs-values grid. Cross section begins at " +
        std::to_string(i) + " in XSS block.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  for (size_t l = 0; l < values_.size(); l++) {
    if (values_[l] < 0.) {
      std::string mssg =
          "CrossSection::CrossSection: Negative cross section found at "
          "element " +
          std::to_string(l) + " in cross section grid starting at " +
          std::to_string(i) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

CrossSection::CrossSection(const std::vector<double>& xs,
                           std::shared_ptr<EnergyGrid> E_grid, size_t index)
    : energy_grid_(E_grid),
      values_(xs.size(), 0.),
      index_(static_cast<uint32_t>(index)) {
  // Fill cross section values
  // DM: does std::copy work here?
  for (size_t e = 0; e < xs.size(); e++) {
    values_[e] = static_cast<float>(xs[e]);
  }

  if (index_ >= energy_grid_->size()) {
    std::string mssg =
        "CrossSection::CrossSection: Starting index is larger than size of the "
        "energy grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  for (size_t l = 0; l < values_.size(); l++) {
    if (values_[l] < 0.) {
      std::string mssg =
          "CrossSection::CrossSection: Negative cross section found at "
          "element " +
          std::to_string(l) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  if (energy_grid_->size() != values_.size()) {
    std::string mssg =
        "CrossSection::CrossSection: Different number of points in the energy "
        "grid and xs-values grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

CrossSection::CrossSection(const std::vector<float>& xs,
                           std::shared_ptr<EnergyGrid> E_grid, size_t index)
    : energy_grid_(E_grid),
      values_(xs.begin(), xs.end()),
      index_(static_cast<uint32_t>(index)) {
  if (index_ >= energy_grid_->size()) {
    std::string mssg =
        "CrossSection::CrossSection: Starting index is larger than size of the "
        "energy grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  for (size_t l = 0; l < values_.size(); l++) {
    if (values_[l] < 0.) {
      std::string mssg =
          "CrossSection::CrossSection: Nevative cross section found at "
          "element " +
          std::to_string(l) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  if (energy_grid_->size() + index_ != values_.size()) {
    std::string mssg =
        "CrossSection::CrossSection: Different number of points in the energy "
        "grid and xs-values grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

size_t CrossSection::size() const { return values_.size(); }

double CrossSection::xs(size_t i) const { return values_[i]; }

double CrossSection::energy(size_t i) const {
  return (*energy_grid_)[index_ + i];
}

uint32_t CrossSection::index() const { return index_; }

const std::vector<float>& CrossSection::xs() const { return values_; }

std::vector<float> CrossSection::energy() const {
  return {energy_grid_->grid().begin() + index_, energy_grid_->grid().end()};
}

}  // namespace pndl
