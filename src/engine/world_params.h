// src/engine/world_params.h
#pragma once

namespace botany {

struct WorldParams {
    float light_level = 1.0f;
    int sugar_diffusion_iterations = 5;
};

inline WorldParams default_world_params() {
    return WorldParams{
        .light_level = 1.0f,
        .sugar_diffusion_iterations = 5,
    };
}

} // namespace botany
