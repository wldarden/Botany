#pragma once

#include "engine/node.h"

namespace botany {

class ShootApicalMeristem : public Meristem {
public:
    ShootApicalMeristem() : Meristem(true) {}
    MeristemType type() const override { return MeristemType::APICAL; }
    bool is_tip() const override { return true; }
    void tick(Node& node, Plant& plant) override;
};

class ShootAxillaryMeristem : public Meristem {
public:
    ShootAxillaryMeristem() : Meristem(false) {}
    MeristemType type() const override { return MeristemType::AXILLARY; }
    bool is_tip() const override { return false; }
    void tick(Node& node, Plant& plant) override;
};

class RootApicalMeristem : public Meristem {
public:
    RootApicalMeristem() : Meristem(true) {}
    MeristemType type() const override { return MeristemType::ROOT_APICAL; }
    bool is_tip() const override { return true; }
    void tick(Node& node, Plant& plant) override;
};

class RootAxillaryMeristem : public Meristem {
public:
    RootAxillaryMeristem() : Meristem(false) {}
    MeristemType type() const override { return MeristemType::ROOT_AXILLARY; }
    bool is_tip() const override { return false; }
    void tick(Node& node, Plant& plant) override;
};

} // namespace botany
