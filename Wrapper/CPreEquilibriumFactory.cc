/**
 * Copyright (c) 2026 Savva Savenkov, Artemii Novikov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "CPreEquilibriumFactory.hh"

#include "CoordinateMSTClustering.hh"

#include <memory>

using namespace cola;

std::unique_ptr<VFilter> CPreEquilibriumFactory::Create(const std::unordered_map<std::string, std::string>& paramMap) {
  if (paramMap.at("clustering_type") == "GMST") {
    if (auto it = paramMap.find("stat_exen_type"); it != paramMap.end()) {
      excitationEnergyType = std::stoi(it->second);
    }
    if (auto it = paramMap.find("consider_coulomb"); it != paramMap.end()) {
      repulsion = std::stoi(it->second);
    }
    if (auto it = paramMap.find("simulate_momentum"); it != paramMap.end()) {
      momentum = std::stoi(it->second);
    }
    return std::make_unique<CoordinateMSTClustering>(repulsion.value_or(false), momentum.value_or(true),
                                                     excitationEnergyType.value_or(7));
  } else {
    throw std::runtime_error("Clustering type is unrecognized");
  }
}