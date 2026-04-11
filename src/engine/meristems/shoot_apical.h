// src/engine/meristems/shoot_apical.h
#pragma once

#include "engine/node.h"

namespace botany {

class ShootApicalMeristem : public Meristem {
public:
    ShootApicalMeristem() : Meristem(true) {}
    MeristemType type() const override { return MeristemType::APICAL; }
    bool is_tip() const override { return true; }
    void tick(Node& node, Plant& plant, const struct WorldParams& world) override;

    float target_internode_length = 0.0f;  // set per-internode; 0 = needs roll
    uint32_t phyllotaxis_index = 0;        // increments each internode; drives golden-angle spiral

private:
    // Extend the tip using sugar-scaled growth. Returns false if no sugar available.
    bool grow(Node& node, const struct Genome& g, const struct WorldParams& world, const glm::vec3& dir);
    // When internode is long enough, split into interior node + new tip.
    void split_internode(Node& node, Plant& plant, const struct Genome& g, const glm::vec3& dir);
};

} // namespace botany
