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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/discrete_photon.hpp>
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/watt.hpp>

namespace py = pybind11;

using namespace pndl;

// Trampoline class for abstract pndl::EnergyLaw
class PyEnergyLaw : public EnergyLaw {
 public:
  using EnergyLaw::EnergyLaw;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override {
    PYBIND11_OVERRIDE_PURE(double, EnergyLaw, sample_energy, E_in, rng);
  }

  double pdf(double E_in, double E_out) const override {
    PYBIND11_OVERRIDE_PURE(double, EnergyLaw, pdf, E_in, E_out);
  }
};

void init_EnergyLaw(py::module& m) {
  py::class_<EnergyLaw, PyEnergyLaw, std::shared_ptr<EnergyLaw>>(m, "EnergyLaw")
      .def(py::init<>())
      .def("sample_energy", &EnergyLaw::sample_energy)
      .def("pdf", &EnergyLaw::pdf);
}

void init_EquiprobableEnergyBins(py::module& m) {
  py::class_<EquiprobableEnergyBins, EnergyLaw,
             std::shared_ptr<EquiprobableEnergyBins>>(m,
                                                      "EquiprobableEnergyBins")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<const std::vector<double>&,
                    const std::vector<std::vector<double>>&>())
      .def("sample_energy", &EquiprobableEnergyBins::sample_energy)
      .def("pdf", &EquiprobableEnergyBins::pdf)
      .def("size", &EquiprobableEnergyBins::size)
      .def("incoming_energy", &EquiprobableEnergyBins::incoming_energy)
      .def("bin_bounds", &EquiprobableEnergyBins::bin_bounds);
}

void init_DiscretePhoton(py::module& m) {
  py::class_<DiscretePhoton, EnergyLaw, std::shared_ptr<DiscretePhoton>>(
      m, "DiscretePhoton")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<int, double, double>())
      .def("sample_energy", &DiscretePhoton::sample_energy)
      .def("pdf", &DiscretePhoton::pdf)
      .def("primary_indicator", &DiscretePhoton::primary_indicator)
      .def("photon_energy", &DiscretePhoton::photon_energy);
}

void init_LevelInelasticScatter(py::module& m) {
  py::class_<LevelInelasticScatter, EnergyLaw,
             std::shared_ptr<LevelInelasticScatter>>(m, "LevelInelasticScatter")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<double, double>())
      .def("sample_energy", &LevelInelasticScatter::sample_energy)
      .def("pdf", &LevelInelasticScatter::pdf)
      .def("C1", &LevelInelasticScatter::C1)
      .def("C2", &LevelInelasticScatter::C2);
}

void init_TabularEnergy(py::module& m) {
  py::class_<TabularEnergy, EnergyLaw, std::shared_ptr<TabularEnergy>>(
      m, "TabularEnergy")
      .def(py::init<const ACE&, size_t, size_t>())
      .def(py::init<const std::vector<double>&, const std::vector<PCTable>&>())
      .def("sample_energy", &TabularEnergy::sample_energy)
      .def("pdf", &TabularEnergy::pdf)
      .def("incoming_energy", &TabularEnergy::incoming_energy)
      .def("table", &TabularEnergy::table)
      .def("size", &TabularEnergy::size);
}

void init_GeneralEvaporation(py::module& m) {
  py::class_<GeneralEvaporation, EnergyLaw,
             std::shared_ptr<GeneralEvaporation>>(m, "GeneralEvaporation")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<std::shared_ptr<Tabulated1D>, const std::vector<double>&>())
      .def("sample_energy", &GeneralEvaporation::sample_energy)
      .def("pdf", &GeneralEvaporation::pdf)
      .def("temperature", &GeneralEvaporation::temperature)
      .def("bin_bounds", &GeneralEvaporation::bin_bounds);
}

void init_Evaporation(py::module& m) {
  py::class_<Evaporation, EnergyLaw, std::shared_ptr<Evaporation>>(
      m, "Evaporation")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<std::shared_ptr<Tabulated1D>, double>())
      .def("sample_energy", &Evaporation::sample_energy)
      .def("pdf", &Evaporation::pdf)
      .def("temperature", &Evaporation::temperature)
      .def("U", &Evaporation::U);
}

void init_Maxwellian(py::module& m) {
  py::class_<Maxwellian, EnergyLaw, std::shared_ptr<Maxwellian>>(m,
                                                                 "Maxwellian")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<std::shared_ptr<Tabulated1D>, double>())
      .def("sample_energy", &Maxwellian::sample_energy)
      .def("pdf", &Maxwellian::pdf)
      .def("temperature", &Maxwellian::temperature)
      .def("U", &Maxwellian::U);
}

void init_Watt(py::module& m) {
  py::class_<Watt, EnergyLaw, std::shared_ptr<Watt>>(m, "Watt")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<std::shared_ptr<Tabulated1D>, std::shared_ptr<Tabulated1D>,
                    double>())
      .def("sample_energy", &Watt::sample_energy)
      .def("pdf", &Watt::pdf)
      .def("a", &Watt::a)
      .def("b", &Watt::b)
      .def("U", &Watt::U);
}
