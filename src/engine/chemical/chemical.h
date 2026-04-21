// src/engine/chemical/chemical.h
#pragma once

#include <cstdint>
#include <functional>

namespace botany {

enum class ChemicalID : uint8_t {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
    Water,
    Count,   // sentinel — keep last; used for fixed-size per-chemical arrays
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

// std::hash specialization — defined here (where ChemicalID lives) so any
// translation unit that includes chemical.h gets a consistent definition.
namespace std {
    template<>
    struct hash<botany::ChemicalID> {
        std::size_t operator()(botany::ChemicalID id) const noexcept {
            return static_cast<std::size_t>(id);
        }
    };
}
