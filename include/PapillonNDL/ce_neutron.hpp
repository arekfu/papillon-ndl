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
#ifndef PAPILLON_NDL_CE_NEUTRON_H
#define PAPILLON_NDL_CE_NEUTRON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/delayed_group.hpp>
#include <PapillonNDL/reaction.hpp>
#include <array>

namespace pndl {

/**
 * @brief Holds all continuous energy neutron data for a single nuclide
 *        and at a single temperature.
 */
class CENeutron {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutron(const ACE& ace);

  /**
   * @param ace ACE file from which to take the new cross sections.
   * @param nuclide CENeutron containing another instance of the desired
   *                nuclide. Secondary distributions and fission data
   *                will be shared between the two data sets.
   */
  CENeutron(const ACE& ace, const CENeutron& nuclide);

  ~CENeutron() = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  uint32_t zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns true if the nuclide is fissile, and false otherwise.
   */
  bool fissile() const { return fissile_; }

  /**
   * @brief Returns the energy grid for the nuclide.
   */
  std::shared_ptr<EnergyGrid> energy_grid() const;
  // DM: unless you have specific reasons to return the pointer and expose this
  // internal implementation detail, I would probably return an EnergyGrid
  // const & here. It does not prevent you from sharing the internal pointers
  // between instances in copy constructors, but it leaves less room for
  // mischief. Note that the method is const, which I, as a user, read as a
  // promise that I cannot modify the object using this method. However, if I
  // do
  //   auto egrid = neutron.energy_grid();
  //   egrid.modify();  // totally existing method
  // then I have modified neutron, too. Returning an EnergyGrid const &
  // prevents users from doing that. This is what you do with DelayedGroups and
  // Reactions below, by the way.
  //
  // Also, I would perhaps put the definitions of this and the following
  // methods in the header file, to give the compiler better chances at
  // inlining them.

  /**
   * @brief Returns a pointer to the total CrossSection for the nuclide.
   */
  std::shared_ptr<CrossSection> total_cross_section() const;

  /**
   * @brief Returns a pointer to the elastic scattering CrossSection for the
   * nuclide.
   */
  std::shared_ptr<CrossSection> elastic_cross_section() const;

  /**
   * @brief Returns a pointer to the disappearance CrossSection for the nuclide.
   */
  std::shared_ptr<CrossSection> disappearance_cross_section() const;

  /**
   * @brief Returns a pointer to the photon production CrossSection for the
   * nuclide.
   */
  std::shared_ptr<CrossSection> photon_production_cross_section() const;

  /**
   * @brief Returns a pointer to the function for total nu.
   */
  std::shared_ptr<Function1D> nu_total() const { return nu_total_; }

  /**
   * @brief Returns a pointer to the function for prompt nu.
   */
  std::shared_ptr<Function1D> nu_prompt() const { return nu_prompt_; }

  /**
   * @brief Returns a pointer to the function for delayed nu.
   */
  std::shared_ptr<Function1D> nu_delayed() const { return nu_delayed_; }

  /**
   * @brief Returns a pointer to the AngleDistribution for elastic scattering.
   */
  std::shared_ptr<AngleDistribution> elastic_angle_distribution() const;

  /**
   * @brief Retrieves the index in the energy grid for an energy.
   * @param E Energy to find in the energy grid.
   */
  // DM: I would not clutter the API with methods such as this one and the
  // following ones. There is a general principle called Information Expert
  // that states that you should assign responsibilities to the class that has
  // the information required to fulfil it. This is also related to the
  // principle of Interface Segregation (the "I" in SOLID) that argues that you
  // need to keep interfaces as small and separated as possible.
  //
  // In this case, CENeutron already has an energy_grid() method. You don't
  // need the extra indirection and the compiler should inline and optimise
  // away the function calls anyway.
  size_t energy_grid_index(double E) const {
    return energy_grid_->get_lower_index(E);
  }

  /**
   * @brief Evaluates the total cross section at E using bisection search.
   * @param E Energy.
   */
  double total_xs(double E) const { return total_xs_->evaluate(E); }

  /**
   * @brief Evaluates the total cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double total_xs(double E, size_t i) const {
    return total_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the elastic scattering cross section at E using bisection
   * search.
   * @param E Energy.
   */
  double elastic_xs(double E) const { return elastic_xs_->evaluate(E); }

  /**
   * @brief Evaluates the elastic scattering cross section at energy E and index
   * i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double elastic_xs(double E, size_t i) const {
    return elastic_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the disappearance cross section at E using bisection
   * search.
   * @param E Energy.
   */
  double disappearance_xs(double E) const {
    return disappearance_xs_->evaluate(E);
  }

  /**
   * @brief Evaluates the disappearance cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double disappearance_xs(double E, size_t i) const {
    return disappearance_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the total fission cross section at E using bisection
   * search.
   * @param E Energy.
   */
  // DM: This method is somewhat of a special case. The methods above
  // (elastic_xs, etc.) may disappear, but this one contains some logic. I am
  // not sure what to do about this. Either live with the wrinkle, or perhaps
  // introduce a helper FissionCrossSection class with an evaluate() method
  // (the implementation of which would basically be the content of
  // fission_xs() below) and return that from here.
  double fission_xs(double E) const {
    if (!this->fissile_) return 0.;

    if (this->has_reaction(18)) return this->reaction_xs(18, E);

    double sigma_f = 0.;

    if (this->has_reaction(19)) sigma_f += this->reaction_xs(19, E);

    if (this->has_reaction(20)) sigma_f += this->reaction_xs(20, E);

    if (this->has_reaction(21)) sigma_f += this->reaction_xs(21, E);

    if (this->has_reaction(38)) sigma_f += this->reaction_xs(38, E);

    return sigma_f;
  }

  /**
   * @brief Evaluates the total fission cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double fission_xs(double E, size_t i) const {
    if (!this->fissile_) return 0.;

    if (this->has_reaction(18)) return this->reaction_xs(18, E, i);

    double sigma_f = 0.;

    if (this->has_reaction(19)) sigma_f += this->reaction_xs(19, E, i);

    if (this->has_reaction(20)) sigma_f += this->reaction_xs(20, E, i);

    if (this->has_reaction(21)) sigma_f += this->reaction_xs(21, E, i);

    if (this->has_reaction(38)) sigma_f += this->reaction_xs(38, E, i);

    return sigma_f;
  }

  /**
   * @brief Evaluates the photon production cross section at E using bisection
   * search.
   * @param E Energy.
   */
  double photon_production_xs(double E) const {
    // Need to check photon production XS exists, as this one is not
    // necessarily present.
    if (photon_production_xs_) return photon_production_xs_->evaluate(E);
    return 0.;
  }

  /**
   * @brief Evaluates the photon production cross section at energy E and index
   * i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double photon_production_xs(double E, size_t i) const {
    // Need to check photon production XS exists, as this one is not
    // necessarily present.
    if (photon_production_xs_) return photon_production_xs_->evaluate(E, i);
    return 0.;
  }

  /**
   * @brief Returns the total average number of fission neutron at energy E.
   * @param E Energy in MeV.
   */
  // DM: same remarks as in the case of fission_xs above
  double nu_total(double E) const {
    if (nu_total_)
      return (*nu_total_)(E);

    else if (nu_prompt_ && nu_delayed_)
      return (*nu_prompt_)(E) + (*nu_delayed_)(E);

    else if (nu_delayed_)
      return (*nu_delayed_)(E);

    return 0.;
  }

  /**
   * @brief Returns the average number of prompt fission neutron at energy E.
   * @param E Energy in MeV.
   */
  double nu_prompt(double E) const {
    if (nu_prompt_) return (*nu_prompt_)(E);

    // If no prompt, that means no delayed data either, so all
    // neutrons are treated as prompt.
    else if (nu_total_)
      return (*nu_total_)(E);

    return 0.;
  }

  /**
   * @brief Returns the average number of delayed fission neutron at energy E.
   * @param E Energy in MeV.
   */
  double nu_delayed(double E) const {
    if (nu_delayed_)
      return (*nu_delayed_)(E);

    else if (nu_total_ && nu_prompt_)
      return (*nu_total_)(E) - (*nu_prompt_)(E);

    return 0.;
  }

  /**
   * @brief Returns the number of delayed neutron groups.
   */
  size_t n_delayed_groups() const { return delayed_groups_.size(); }

  /**
   * @brief Returns the ith delayed group data.
   * @param i Index of the delayed group.
   */
  // DM: returning a DelayedGroup const & is fine. The integer index makes me a
  // bit uneasy. Either you expose the internal implementation detail that you
  // are using std::vector and return that, or maybe you can go the extra mile
  // and replace std::vector<DelayedGroup> with your class DelayedGroups, which
  // you can implement by either privately inheriting from std::vector or by
  // aggregating an std::vector. That hides the implementation detail at the
  // cost of some more boilerplate.
  const DelayedGroup& delayed_group(size_t i) const {
    return delayed_groups_[i];
  }

  /**
   * @brief Samples a scattering angle from the elastic scattering angular
   *        distribution.
   * @param E Incident energy.
   * @param rng Random number generation function.
   */
  double sample_elastic_angle(double E, std::function<double()> rng) const {
    return elastic_angle_->sample_angle(E, rng);
  }

  /**
   * @brief Returns a list of all MT reactions present for the nuclide.
   */
  // DM: same remark here about implementation details leaking out of the API,
  // but maybe that's fine.
  const std::vector<uint32_t>& mt_list() const { return mt_list_; }

  /**
   * @brief Checks to see if a nucldie has a given reaction.
   * @param mt MT reaction to search for.
   */
  bool has_reaction(uint32_t mt) const {
    return (mt > 891 || reaction_indices_[mt] < 0) ? false : true;
  }

  /**
   * @brief Retrieved a given MT reaction.
   * @param mt MT reaction to return.
   */
  const Reaction& reaction(uint32_t mt) const {
    if (mt > 891 || reaction_indices_[mt] < 0) {
      // DM: maybe replace the bounds check with a call to
      // reaction_indices_.size(), or by accessing it with at().
      std::string mssg = "CENeutron::reaction: MT = " + std::to_string(mt) +
                         " is not provided in ZAID = " + std::to_string(zaid_) +
                         ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    return reactions_[reaction_indices_[mt]];
  }

  /**
   * @brief Returns the cross section for a perscriped reaction at a
   *        provided energy. Uses bisection search.
   * @param mt MT value of the reaction.
   * @param E Energy to evaluate cross section at.
   */
  // DM: remove these two methods? see above
  double reaction_xs(uint32_t mt, double E) const {
    if (mt > 891) return 0.;

    auto indx = reaction_indices_[mt];

    return indx < 0 ? 0. : reactions_[indx].xs(E);
  }

  /**
   * @brief Returns the cross section for a perscriped reaction at a
   *        provided energy, and energy grid index.
   * @param mt MT value of the reaction.
   * @param E Energy to evaluate cross section at.
   * @param i Index to the energy grid for energy E.
   */
  double reaction_xs(uint32_t mt, double E, size_t i) const {
    if (mt > 891) return 0.;

    auto indx = reaction_indices_[mt];

    return indx < 0 ? 0. : reactions_[indx].xs(E, i);
  }

 private:
  uint32_t zaid_;
  double awr_;
  double temperature_;
  bool fissile_;

  std::shared_ptr<EnergyGrid> energy_grid_;

  std::shared_ptr<CrossSection> total_xs_;
  std::shared_ptr<CrossSection> disappearance_xs_;
  std::shared_ptr<CrossSection> elastic_xs_;
  std::shared_ptr<CrossSection> photon_production_xs_;
  // DM: this is a somewhat disruptive remark, and I know that we are not doing
  // that in PATMOS either, but I have always thought that it should be more
  // efficient to store cross sections as an array of structs, instead of a
  // struct of arrays like here. The reason is that when you have an AoS you
  // will load all the cross sections into the cache with just one load.
  // However, that will degrade the performance of cross section lookup while
  // iterating over energy (think Doppler broadening). So the best solution
  // overall is probably some weird AoSoA, but I'm too old for that.
  //
  // Anyway I think this comment goes to show that there may be a slightly
  // smaller class trying to get out here. You could bunch the energy grid and
  // the four cross sections into a smaller class and use that to play around
  // with different implementations of cross section lookups.

  std::shared_ptr<AngleDistribution> elastic_angle_;

  std::shared_ptr<Function1D> nu_total_;
  std::shared_ptr<Function1D> nu_prompt_;
  std::shared_ptr<Function1D> nu_delayed_;
  std::vector<DelayedGroup> delayed_groups_;

  std::vector<uint32_t> mt_list_;
  std::array<int32_t, 892> reaction_indices_;
  std::vector<Reaction> reactions_;
  // DM: have you tried replacing the double indexing with an
  // std::unordered_map<uint32_t,Reaction>?

  // Private helper methods
  void read_fission_data(const ACE& ace);
  std::shared_ptr<Function1D> read_nu(const ACE& ace, size_t i);
  std::shared_ptr<Function1D> read_polynomial_nu(const ACE& ace, size_t i);
  std::shared_ptr<Function1D> read_tabular_nu(const ACE& ace, size_t i);
};

}  // namespace pndl

#endif
