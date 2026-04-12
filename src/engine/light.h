#pragma once

#include <memory>
#include <vector>

namespace botany {

class Plant;
struct WorldParams;

// Compute light exposure for all leaves across all plants using a shadow map.
// Must be called before plant ticks so leaves have current light_exposure for photosynthesis.
void compute_light_exposure(const std::vector<std::unique_ptr<Plant>>& plants, const WorldParams& world);

// Single-plant convenience overload (wraps plant in a temporary vector).
void compute_light_exposure(Plant& plant, const WorldParams& world);

} // namespace botany
