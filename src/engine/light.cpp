#include "engine/light.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include <cmath>
#include <unordered_map>
#include <vector>
#include <glm/geometric.hpp>

namespace botany {

struct CellKey {
    int32_t u, v;
    bool operator==(const CellKey& o) const { return u == o.u && v == o.v; }
};

struct CellKeyHash {
    size_t operator()(const CellKey& k) const {
        size_t h = std::hash<int32_t>()(k.u);
        h ^= std::hash<int32_t>()(k.v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct CasterEntry {
    float depth;
    float coverage;
};

void compute_light_exposure(const std::vector<std::unique_ptr<Plant>>& plants,
                            const WorldParams& world,
                            ShadowMapViz* viz_out) {
    const float cell = world.light_cell_size;
    if (cell <= 1e-6f) return;

    const glm::vec3 L = glm::normalize(world.light_direction);

    // Shadow plane basis vectors perpendicular to light direction
    glm::vec3 arbitrary = std::abs(glm::dot(L, glm::vec3(0.0f, 1.0f, 0.0f))) < 0.9f
        ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 plane_u = glm::normalize(glm::cross(L, arbitrary));
    glm::vec3 plane_v = glm::cross(L, plane_u);

    std::unordered_map<CellKey, std::vector<CasterEntry>, CellKeyHash> grid;

    // Phase 1: Stamp all above-ground casters from all plants into the grid
    for (const auto& plant : plants) {
        plant->for_each_node_mut([&](Node& node) {
            if (node.type == NodeType::ROOT || node.is_meristem()) return;

            float shadow_radius = node.radius;
            if (auto* leaf = node.as_leaf()) {
                leaf->light_exposure = 1.0f;
                shadow_radius = leaf->leaf_size;
            }
            if (shadow_radius <= 1e-6f) return;

            float depth = glm::dot(node.position, L);
            float proj_u = glm::dot(node.position, plane_u);
            float proj_v = glm::dot(node.position, plane_v);

            int r_cells = static_cast<int>(std::ceil(shadow_radius / cell));
            int cu = static_cast<int>(std::floor(proj_u / cell));
            int cv = static_cast<int>(std::floor(proj_v / cell));

            for (int du = -r_cells; du <= r_cells; du++) {
                for (int dv = -r_cells; dv <= r_cells; dv++) {
                    float cell_center_u = (cu + du + 0.5f) * cell;
                    float cell_center_v = (cv + dv + 0.5f) * cell;
                    float dist_u = cell_center_u - proj_u;
                    float dist_v = cell_center_v - proj_v;
                    float dist = std::sqrt(dist_u * dist_u + dist_v * dist_v);
                    if (dist < shadow_radius) {
                        grid[{cu + du, cv + dv}].push_back({depth, 1.0f});
                    }
                }
            }
        });
    }

    // Phase 2: Query each leaf's exposure from the grid
    for (const auto& plant : plants) {
        plant->for_each_node_mut([&](Node& node) {
            auto* leaf = node.as_leaf();
            if (!leaf || leaf->leaf_size <= 1e-6f || leaf->senescence_ticks != 0) return;

            float depth = glm::dot(node.position, L);
            float proj_u = glm::dot(node.position, plane_u);
            float proj_v = glm::dot(node.position, plane_v);

            int r_cells = static_cast<int>(std::ceil(leaf->leaf_size / cell));
            int cu = static_cast<int>(std::floor(proj_u / cell));
            int cv = static_cast<int>(std::floor(proj_v / cell));

            float total_exposure = 0.0f;
            int cell_count = 0;

            for (int du = -r_cells; du <= r_cells; du++) {
                for (int dv = -r_cells; dv <= r_cells; dv++) {
                    float cell_center_u = (cu + du + 0.5f) * cell;
                    float cell_center_v = (cv + dv + 0.5f) * cell;
                    float dist_u = cell_center_u - proj_u;
                    float dist_v = cell_center_v - proj_v;
                    float dist = std::sqrt(dist_u * dist_u + dist_v * dist_v);
                    if (dist >= leaf->leaf_size) continue;

                    cell_count++;

                    auto it = grid.find({cu + du, cv + dv});
                    if (it == grid.end()) {
                        total_exposure += 1.0f;
                        continue;
                    }

                    float coverage_above = 0.0f;
                    for (const auto& entry : it->second) {
                        if (entry.depth > depth) {
                            coverage_above += entry.coverage;
                        }
                    }
                    total_exposure += std::max(0.0f, 1.0f - coverage_above);
                }
            }

            leaf->light_exposure = cell_count > 0
                ? total_exposure / static_cast<float>(cell_count)
                : 1.0f;
        });
    }

    // Phase 3: Output visualization data if requested
    if (viz_out) {
        viz_out->cell_size = cell;
        viz_out->cells.clear();
        viz_out->cells.reserve(grid.size());

        for (const auto& [key, entries] : grid) {
            // Convert cell coords back to world XZ (for overhead sun: u=z, v=x)
            float cu = (key.u + 0.5f) * cell;
            float cv = (key.v + 0.5f) * cell;
            // World position = cu * plane_u + cv * plane_v (projected onto ground y=0)
            glm::vec3 world_pos = cu * plane_u + cv * plane_v;

            float total_coverage = 0.0f;
            for (const auto& e : entries) {
                total_coverage += e.coverage;
            }

            viz_out->cells.push_back({world_pos.x, world_pos.z, total_coverage});
        }
    }
}

void compute_light_exposure(Plant& plant, const WorldParams& world) {
    std::vector<std::unique_ptr<Plant>> tmp;
    tmp.emplace_back(&plant);
    compute_light_exposure(tmp, world, nullptr);
    tmp[0].release();
}

} // namespace botany
