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
#ifndef PAPILLON_NDL_EQUIPROBABLE_ANGLE_BINS_H
#define PAPILLON_NDL_EQUIPROBABLE_ANGLE_BINS_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_law.hpp>

namespace pndl {

/**
 * @brief Angular distribution represented as equiprobable cosine bins.
 */
class EquiprobableAngleBins : public AngleLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  EquiprobableAngleBins(const ACE& ace, size_t i);

  /**
   * @param bounds Vector of 33 bin bounds.
   */
  EquiprobableAngleBins(const std::vector<double>& bounds);
  ~EquiprobableAngleBins() = default;

  double sample_mu(double xi) const override final;

  double pdf(double mu) const override final;

  /**
   * @brief Returns the number of bin boundaries (number of bins + 1);
   */
  size_t size() const;

  /**
   * @brief Returns the vector with the bin boundaries.
   */
  const std::vector<double>& bin_bounds() const;

 private:
  std::vector<double> bounds_;

  static constexpr size_t NBOUNDS = 33;
  static constexpr size_t NBINS = 32;
  static constexpr double P_BIN = 1. / 32.;
};

}  // namespace pndl

#endif
