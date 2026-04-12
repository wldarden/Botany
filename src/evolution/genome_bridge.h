#pragma once

#include "engine/genome.h"
#include <evolve/structured_genome.h>

namespace botany {

evolve::StructuredGenome build_genome_template(const Genome& g);
evolve::StructuredGenome to_structured(const Genome& g);
Genome from_structured(const evolve::StructuredGenome& sg);

} // namespace botany
