#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/genome.h"

namespace botany {

bool can_merge(const Node& /*parent*/, const Node& /*child*/,
               const Genome& /*g*/, const CompressionParams& /*params*/) {
    return false; // Filled in by Task 2.
}

CompressionResult compress_plant(Plant& /*plant*/, const CompressionParams& /*params*/) {
    CompressionResult r;
    return r; // Filled in by Tasks 3-4.
}

} // namespace botany
