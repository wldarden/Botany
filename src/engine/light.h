#pragma once

#include <memory>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

class Plant;
struct WorldParams;

struct ShadowMapCell {
    float x, z;         // world position of cell center (on ground plane)
    float coverage;     // total shadow coverage (0 = full sun, >=1 = fully shaded)
};

struct ShadowMapViz {
    float cell_size = 0.0f;
    std::vector<ShadowMapCell> cells;
};

// Compute light exposure for all leaves across all plants using a shadow map.
// Must be called before plant ticks so leaves have current light_exposure for photosynthesis.
// Optionally outputs visualization data.
void compute_light_exposure(const std::vector<std::unique_ptr<Plant>>& plants,
                            const WorldParams& world,
                            ShadowMapViz* viz_out = nullptr);

// Single-plant convenience overload.
void compute_light_exposure(Plant& plant, const WorldParams& world);

} // namespace botany
