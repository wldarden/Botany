#pragma once

#include <string>
#include "engine/genome.h"

namespace botany {

// Parse a ".env"-style text genome file: one "key=value" per line.
// Unrecognized keys are silently ignored (so old files survive field additions).
// Missing file: returns default_genome().
Genome load_genome_file(const std::string& path);

// Write a genome as "key=value" lines.  Writes every field currently parsed by
// load_genome_file so round-trip is lossless for those fields.
// Returns true on success.
bool save_genome_file(const Genome& g, const std::string& path);

} // namespace botany
