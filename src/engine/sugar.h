// src/engine/sugar.h
#pragma once

namespace botany {

class Plant;
struct WorldParams;

void produce_sugar(Plant& plant, const WorldParams& world);
void consume_sugar(Plant& plant);
void diffuse_sugar(Plant& plant, const WorldParams& world);
void transport_sugar(Plant& plant, const WorldParams& world);

} // namespace botany
