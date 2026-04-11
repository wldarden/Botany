// src/engine/sugar.h
#pragma once

namespace botany {

class Plant;
struct Node;
struct Genome;
struct WorldParams;

float sugar_cap(const Node& node, const Genome& g);

void compute_light_exposure(Plant& plant, const WorldParams& world);
void prune_starved_nodes(Plant& plant, const WorldParams& world);
void transport_sugar(Plant& plant, const WorldParams& world);

} // namespace botany
