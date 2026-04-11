// src/engine/meristems/root_axillary.h
#pragma once

#include "engine/node.h"

namespace botany {

class RootAxillaryMeristem : public Meristem {
public:
    RootAxillaryMeristem() : Meristem(false) {}
    MeristemType type() const override { return MeristemType::ROOT_AXILLARY; }
    bool is_tip() const override { return false; }
    void tick(Node& node, Plant& plant, const struct WorldParams& world) override;
};

} // namespace botany
