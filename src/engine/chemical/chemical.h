// src/engine/chemical/chemical.h
#pragma once

#include <cstdint>

namespace botany {

enum class ChemicalID : uint8_t {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
    Water,
};

enum class ChemicalCategory : uint8_t {
    Hormone,
    Resource,
    Volatile,
};

struct ChemicalDef {
    ChemicalID id;
    const char* name;
    ChemicalCategory category;
};

} // namespace botany
