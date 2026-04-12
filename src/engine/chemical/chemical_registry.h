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

// Chemicals that diffuse through the tree graph.
// One entry per chemical — rate and decay from the genome.
inline std::array<ChemicalDiffusionParams, 4> diffusion_params(const Genome& g) {
    return {{
        {ChemicalID::Auxin,       g.auxin_diffusion_rate,     g.auxin_decay_rate},
        {ChemicalID::Cytokinin,   g.cytokinin_diffusion_rate, g.cytokinin_decay_rate},
        {ChemicalID::Gibberellin, g.ga_diffusion_rate,        g.ga_decay_rate},
        {ChemicalID::Sugar,       g.sugar_diffusion_rate,     0.0f},  // sugar doesn't decay
    }};
}

} // namespace botany
