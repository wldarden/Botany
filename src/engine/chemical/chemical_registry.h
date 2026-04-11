// src/engine/chemical/chemical_registry.h
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/chemical/hormone/hormone.h"
#include "engine/genome.h"
#include <array>

namespace botany {

// All chemical IDs in the system. Add new chemicals here.
inline constexpr std::array<ChemicalID, 5> all_chemical_ids = {
    ChemicalID::Auxin,
    ChemicalID::Cytokinin,
    ChemicalID::Gibberellin,
    ChemicalID::Sugar,
    ChemicalID::Ethylene,
};

// Extract hormone transport params from a genome.
// Add new hormones here — one line per hormone.
inline std::array<HormoneTransportParams, 3> hormone_params(const Genome& g) {
    return {{
        {ChemicalID::Auxin, g.auxin_transport_rate, g.auxin_directional_bias, g.auxin_decay_rate},
        {ChemicalID::Cytokinin, g.cytokinin_transport_rate, g.cytokinin_directional_bias, g.cytokinin_decay_rate},
        {ChemicalID::Gibberellin, g.ga_transport_rate, g.ga_directional_bias, g.ga_decay_rate},
    }};
}

} // namespace botany
