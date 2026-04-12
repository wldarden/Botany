#pragma once

namespace botany {

class Plant;
struct WorldParams;

void compute_ethylene(Plant& plant, const WorldParams& world);

} // namespace botany
