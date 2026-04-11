// src/engine/chemical/hormone/hormone.h
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

// Params extracted from a plant's genome for one hormone.
// Used by Node::transport_chemicals() to loop over all hormones generically.
struct HormoneTransportParams {
    ChemicalID id;
    float transport_rate;
    float directional_bias;
    float decay_rate;
};

} // namespace botany
