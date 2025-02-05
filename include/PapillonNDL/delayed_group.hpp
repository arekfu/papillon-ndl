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
 * de modification et de redistribution accordés par cette licence, il n'est
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
#ifndef PAPILLON_NDL_DELAYED_GROUP_H
#define PAPILLON_NDL_DELAYED_GROUP_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

// The delayed group numbers start at g = 1, and go up !
// g = 0 would be the prompt neutorns.

namespace pndl {

/**
 * @brief Contains data for a delayed neutron group.
 */
class DelayedGroup {
 public:
  /**
   * @param ace ACE file to take delayed neutron data from.
   * @param i Index to the beinning of the delayed group data
   *          in the XSS block.
   * @param g Delayed group index.
   */
  DelayedGroup(const ACE& ace, size_t i, size_t g);
  ~DelayedGroup() = default;

  /**
   * @brief Returns the decay constant for the group in inverse seconds.
   */
  double decay_constant() const { return decay_constant_; }

  /**
   * @brief Returns pointer to the Tabulated1D function for the probability
   *        of selecting the delayed group for a given energy.
   */
  std::shared_ptr<Tabulated1D> probability() const { return probability_; }

  /**
   * @brief Evaluates the probability of selecting the delayed group
   *        at incident energy E.
   * @param E Incident energy in MeV.
   */
  double probability(double E) const { return (*probability_)(E); }

  /**
   * @brief Samples and energy from the delayed group distribution.
   * @param E Incident energy.
   * @param rng Random number generation function.
   */
  double sample_energy(double E, std::function<double()> rng) const {
    return energy_->sample_energy(E, rng);
  }

  /**
   * @brief Returns a pointer to the EnergyLaw for the group.
   */
  std::shared_ptr<EnergyLaw> energy() const { return energy_; }

 private:
  double decay_constant_;
  std::shared_ptr<Tabulated1D> probability_;
  std::shared_ptr<EnergyLaw> energy_;
};

}  // namespace pndl

#endif
