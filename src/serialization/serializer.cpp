#include "serialization/serializer.h"
#include "engine/engine.h"
#include "engine/node/tissues/leaf.h"

namespace botany {

template<typename T>
static void write_val(std::ostream& out, const T& val) {
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template<typename T>
static T read_val(std::istream& in) {
    T val;
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
    return val;
}

void save_recording_header(std::ostream& out, const Genome& genome, uint32_t num_ticks) {
    write_val(out, num_ticks);
    write_val(out, genome);
}

RecordingHeader load_recording_header(std::istream& in) {
    RecordingHeader header;
    header.num_ticks = read_val<uint32_t>(in);
    header.genome = read_val<Genome>(in);
    return header;
}

void save_tick(std::ostream& out, const Engine& engine, uint32_t plant_id) {
    const Plant& plant = engine.get_plant(plant_id);
    uint32_t tick = engine.get_tick();
    uint32_t count = plant.node_count();

    write_val(out, tick);
    write_val(out, count);

    plant.for_each_node([&](const Node& node) {
        write_val(out, node.id);
        uint32_t parent_id = node.parent ? node.parent->id : UINT32_MAX;
        write_val(out, parent_id);
        write_val(out, node.type);
        write_val(out, node.position);
        write_val(out, node.radius);
        write_val(out, node.chemical(ChemicalID::Auxin));
        write_val(out, node.chemical(ChemicalID::Cytokinin));
        write_val(out, node.chemical(ChemicalID::Sugar));
        float ls = node.as_leaf() ? node.as_leaf()->leaf_size : 0.0f;
        write_val(out, ls);
        glm::vec3 fc = node.as_leaf() ? node.as_leaf()->facing : glm::vec3(0.0f, 1.0f, 0.0f);
        write_val(out, fc);
    });
}

TickSnapshot load_tick(std::istream& in) {
    TickSnapshot snap;
    snap.tick_number = read_val<uint32_t>(in);
    uint32_t count = read_val<uint32_t>(in);
    snap.nodes.resize(count);

    for (uint32_t i = 0; i < count; i++) {
        NodeSnapshot& ns = snap.nodes[i];
        ns.id = read_val<uint32_t>(in);
        ns.parent_id = read_val<uint32_t>(in);
        ns.type = read_val<NodeType>(in);
        ns.position = read_val<glm::vec3>(in);
        ns.radius = read_val<float>(in);
        ns.auxin = read_val<float>(in);
        ns.cytokinin = read_val<float>(in);
        ns.sugar = read_val<float>(in);
        ns.leaf_size = read_val<float>(in);
        ns.facing = read_val<glm::vec3>(in);
    }

    return snap;
}

} // namespace botany
