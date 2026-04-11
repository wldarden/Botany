// src/engine/meristems/meristem.h
#pragma once
#include "engine/node/meristems/meristem_types.h"

namespace botany {

class Plant;
struct WorldParams;

void tick_meristems(Plant& plant, const WorldParams& world);

} // namespace botany
