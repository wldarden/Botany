#pragma once

namespace botany {
class Plant;
struct Genome;

void pin_transport(Plant& plant, const Genome& g);

// Explicit direct auxin transfer between every shoot-child and every
// root-child of the seed.  The seed is the formal junction between the two
// halves of the plant, but auxin flow THROUGH seed.local has historically
// been fragile — it depends on pin_transport phase ordering, per-tick
// timing, and transport_with_children flow signs — and the connection has
// broken several times during refactors.  This function adds a robust,
// explicit safety-net path that copies a diffusion-weighted fraction of
// auxin from shoot children directly into root children each tick,
// bypassing seed.local entirely.
//
// Called after pin_transport(), before tick_recursive().
void diffuse_auxin_across_seed_junction(Plant& plant, const Genome& g);

} // namespace botany
