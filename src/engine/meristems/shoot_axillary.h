// src/engine/meristems/shoot_axillary.h
#pragma once

#include "engine/node.h"

namespace botany {

class ShootAxillaryMeristem : public Meristem {
public:
    ShootAxillaryMeristem() : Meristem(false) {}
    MeristemType type() const override { return MeristemType::AXILLARY; }
    bool is_tip() const override { return false; }
    void tick(Node& node, Plant& plant, const struct WorldParams& world) override;
};

} // namespace botany
