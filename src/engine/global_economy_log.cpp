#include "engine/global_economy_log.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/world_params.h"
#include <iomanip>
#include <filesystem>

namespace botany {

void GlobalEconomyLog::open(const std::string& path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }
    file_.open(path);
    header_written_ = false;
    prev_storage_.clear();
}

void GlobalEconomyLog::close() {
    file_.close();
}

namespace {

// Chemicals we track.  Order matters for CSV output.
const ChemicalID k_tracked_chemicals[] = {
    ChemicalID::Sugar,
    ChemicalID::Water,
    ChemicalID::Auxin,
    ChemicalID::Cytokinin,
    ChemicalID::Gibberellin,
    ChemicalID::Stress,
    ChemicalID::Ethylene,
};

const char* chemical_name(ChemicalID id) {
    switch (id) {
        case ChemicalID::Sugar:       return "Sugar";
        case ChemicalID::Water:       return "Water";
        case ChemicalID::Auxin:       return "Auxin";
        case ChemicalID::Cytokinin:   return "Cytokinin";
        case ChemicalID::Gibberellin: return "Gibberellin";
        case ChemicalID::Stress:      return "Stress";
        case ChemicalID::Ethylene:    return "Ethylene";
        case ChemicalID::Count:       break;
    }
    return "?";
}

} // namespace

void GlobalEconomyLog::log_tick(uint32_t tick, const Plant& plant, const WorldParams& /*world*/) {
    if (!file_.is_open()) return;

    if (!header_written_) {
        file_ << "tick,chemical,storage,produced,consumed,delta,"
              << "storage_local,storage_phloem,storage_xylem,node_count\n";
        header_written_ = true;
    }

    // Accumulators per tracked chemical.
    struct ChemTotals {
        float storage_local  = 0.0f;
        float storage_phloem = 0.0f;
        float storage_xylem  = 0.0f;
        float produced       = 0.0f;
        float consumed       = 0.0f;
    };
    std::unordered_map<ChemicalID, ChemTotals> totals;

    uint32_t node_count = 0;
    plant.for_each_node([&](const Node& n) {
        ++node_count;

        for (ChemicalID id : k_tracked_chemicals) {
            ChemTotals& t = totals[id];
            t.storage_local += n.local().chemical(id);
            if (auto* p = n.phloem()) t.storage_phloem += p->chemical(id);
            if (auto* x = n.xylem())  t.storage_xylem  += x->chemical(id);
        }

        // Known per-chemical production/consumption trackers.
        // Sugar: maintenance is explicit consumption.  tick_sugar_activity
        // encodes net tissue work (positive = photosynthesis production,
        // negative = growth spending), so we split and sum.
        {
            ChemTotals& t = totals[ChemicalID::Sugar];
            t.consumed += n.tick_sugar_maintenance;
            if (n.tick_sugar_activity > 0.0f) t.produced += n.tick_sugar_activity;
            else                              t.consumed += -n.tick_sugar_activity;
        }
        // Auxin and cytokinin have explicit production trackers.
        totals[ChemicalID::Auxin].produced     += n.tick_auxin_produced;
        totals[ChemicalID::Cytokinin].produced += n.tick_cytokinin_produced;
        // Water, Gibberellin, Stress, Ethylene currently have no per-node
        // per-tick production tracker.  Their produced value is left at 0
        // (placeholder).  Consumed for these can still be inferred from
        // mass balance: consumed = produced - delta.
    });

    // Write one row per chemical.
    for (ChemicalID id : k_tracked_chemicals) {
        const ChemTotals& t = totals[id];
        const float storage = t.storage_local + t.storage_phloem + t.storage_xylem;
        const float prev = prev_storage_.count(id) ? prev_storage_[id] : storage;
        const float delta = storage - prev;

        file_ << tick << "," << chemical_name(id) << ","
              << std::fixed << std::setprecision(6)
              << storage << "," << t.produced << "," << t.consumed << ","
              << delta << ","
              << t.storage_local << "," << t.storage_phloem << "," << t.storage_xylem << ","
              << node_count << "\n";

        prev_storage_[id] = storage;
    }
}

} // namespace botany
