// src/engine/meristems/root_apical.h
#pragma once

#include "engine/node.h"

namespace botany {

class RootApicalMeristem : public Meristem {
public:
    RootApicalMeristem() : Meristem(true) {}
    MeristemType type() const override { return MeristemType::ROOT_APICAL; }
    bool is_tip() const override { return true; }
    void tick(Node& node, Plant& plant, const struct WorldParams& world) override;

    float target_internode_length = 0.0f;
};

} // namespace botany
