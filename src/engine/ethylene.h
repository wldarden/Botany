#pragma once

namespace botany {

class Plant;
struct WorldParams;

void compute_ethylene(Plant& plant, const WorldParams& world);
void process_abscission(Plant& plant);

} // namespace botany
