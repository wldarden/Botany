#pragma once

#include "engine/genome.h"
#include <evolve/structured_genome.h>

namespace botany {

// mutation_pct: mutation strength as fraction of each gene's valid range (e.g. 0.03 = 3%)
evolve::StructuredGenome build_genome_template(const Genome& g, float mutation_pct = 0.03f);
evolve::StructuredGenome to_structured(const Genome& g, float mutation_pct = 0.03f);
Genome from_structured(const evolve::StructuredGenome& sg);

} // namespace botany
