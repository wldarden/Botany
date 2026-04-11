// src/engine/sugar.h
#pragma once

namespace botany {

class Plant;
struct Node;
struct Genome;
struct WorldParams;

float sugar_cap(const Node& node, const Genome& g);

void produce_sugar(Plant& plant, const WorldParams& world);
void consume_sugar(Plant& plant);
void prune_starved_nodes(Plant& plant, const WorldParams& world);
void transport_sugar(Plant& plant, const WorldParams& world);

} // namespace botany
