#include "serialization/plant_snapshot.h"
#include "engine/plant.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace botany {

// -- Binary helpers (file-local). Values are written host-endian; this format
// is not expected to be portable across differently-endian machines, same as
// the existing serializer.cpp convention.
namespace {
template<typename T>
void write_val(std::ostream& out, const T& v) {
    out.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template<typename T>
T read_val(std::istream& in) {
    T v;
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!in) throw std::runtime_error("plant_snapshot: unexpected EOF");
    return v;
}
} // namespace

void plant_snapshot_write_magic(std::ostream& out) {
    out.write(PLANT_SNAPSHOT_MAGIC, 4);
}

bool plant_snapshot_check_magic(std::istream& in) {
    char buf[4];
    in.read(buf, 4);
    if (!in) return false;
    return std::memcmp(buf, PLANT_SNAPSHOT_MAGIC, 4) == 0;
}

// Stubs — filled in by later tasks.
SaveResult save_plant_snapshot(const Plant&, uint64_t, const std::string&) {
    return SaveResult{false, "", "save_plant_snapshot not implemented yet"};
}

LoadedPlant load_plant_snapshot(const std::string&, const std::optional<Genome>&) {
    throw std::runtime_error("load_plant_snapshot not implemented yet");
}

} // namespace botany
