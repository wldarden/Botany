// src/app_realtime.cpp
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/geometric.hpp>
#include <cfloat>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "engine/engine.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/node/meristems/helpers.h"
#include "engine/chemical/chemical_registry.h"
#include "engine/ui_helpers.h"
#include "renderer/renderer.h"
#include "format.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>

using namespace botany;

static Genome load_genome_file(const std::string& path) {
    Genome g = default_genome();
    std::ifstream file(path);
    if (!file) return g;

    std::map<std::string, std::string> fields;
    std::string line;
    while (std::getline(file, line)) {
        auto eq = line.find('=');
        if (eq != std::string::npos) {
            fields[line.substr(0, eq)] = line.substr(eq + 1);
        }
    }

    auto get_f = [&](const char* name, float& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = std::stof(it->second);
    };
    auto get_u = [&](const char* name, uint32_t& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = static_cast<uint32_t>(std::stoul(it->second));
    };

    get_f("apical_auxin_baseline", g.apical_auxin_baseline);
    get_f("auxin_diffusion_rate", g.auxin_diffusion_rate);
    get_f("auxin_decay_rate", g.auxin_decay_rate);
    get_f("auxin_threshold", g.auxin_threshold);
    get_f("cytokinin_production_rate", g.cytokinin_production_rate);
    get_f("cytokinin_diffusion_rate", g.cytokinin_diffusion_rate);
    get_f("cytokinin_decay_rate", g.cytokinin_decay_rate);
    get_f("cytokinin_threshold", g.cytokinin_threshold);
    get_f("cytokinin_growth_threshold", g.cytokinin_growth_threshold);
    get_f("growth_rate", g.growth_rate);
    get_u("shoot_plastochron", g.shoot_plastochron);
    get_f("branch_angle", g.branch_angle);
    get_f("cambium_responsiveness", g.cambium_responsiveness);
    get_f("internode_elongation_rate", g.internode_elongation_rate);
    get_f("max_internode_length", g.max_internode_length);
    get_u("internode_maturation_ticks", g.internode_maturation_ticks);
    get_f("root_growth_rate", g.root_growth_rate);
    get_u("root_plastochron", g.root_plastochron);
    get_f("root_branch_angle", g.root_branch_angle);
    get_f("root_internode_elongation_rate", g.root_internode_elongation_rate);
    get_u("root_internode_maturation_ticks", g.root_internode_maturation_ticks);
    get_f("root_gravitropism_strength", g.root_gravitropism_strength);
    get_f("root_gravitropism_depth", g.root_gravitropism_depth);
    get_f("max_leaf_size", g.max_leaf_size);
    get_f("leaf_growth_rate", g.leaf_growth_rate);
    get_f("leaf_bud_size", g.leaf_bud_size);
    get_f("leaf_petiole_length", g.leaf_petiole_length);
    get_f("leaf_opacity", g.leaf_opacity);
    get_f("initial_radius", g.initial_radius);
    get_f("root_initial_radius", g.root_initial_radius);
    get_f("tip_offset", g.tip_offset);
    get_f("growth_noise", g.growth_noise);
    get_f("leaf_phototropism_rate", g.leaf_phototropism_rate);
    get_f("sugar_diffusion_rate", g.sugar_diffusion_rate);
    get_f("seed_sugar", g.seed_sugar);
    get_f("sugar_storage_density_wood", g.sugar_storage_density_wood);
    get_f("sugar_storage_density_leaf", g.sugar_storage_density_leaf);
    get_f("sugar_cap_minimum", g.sugar_cap_minimum);
    get_f("sugar_cap_meristem", g.sugar_cap_meristem);
    get_f("stem_photosynthesis_rate", g.stem_photosynthesis_rate);
    get_f("stem_green_radius_threshold", g.stem_green_radius_threshold);
    get_f("ga_production_rate", g.ga_production_rate);
    get_u("ga_leaf_age_max", g.ga_leaf_age_max);
    get_f("ga_elongation_sensitivity", g.ga_elongation_sensitivity);
    get_f("ga_length_sensitivity", g.ga_length_sensitivity);
    get_f("ga_diffusion_rate", g.ga_diffusion_rate);
    get_f("ga_decay_rate", g.ga_decay_rate);
    get_f("ethylene_starvation_rate", g.ethylene_starvation_rate);
    get_f("ethylene_shade_rate", g.ethylene_shade_rate);
    get_f("ethylene_shade_threshold", g.ethylene_shade_threshold);
    get_f("ethylene_age_rate", g.ethylene_age_rate);
    get_u("ethylene_age_onset", g.ethylene_age_onset);
    get_f("ethylene_crowding_rate", g.ethylene_crowding_rate);
    get_f("ethylene_crowding_radius", g.ethylene_crowding_radius);
    get_f("ethylene_diffusion_radius", g.ethylene_diffusion_radius);
    get_f("ethylene_abscission_threshold", g.ethylene_abscission_threshold);
    get_f("ethylene_elongation_inhibition", g.ethylene_elongation_inhibition);
    get_u("senescence_duration", g.senescence_duration);
    get_f("wood_density", g.wood_density);
    get_f("wood_flexibility", g.wood_flexibility);
    get_f("stress_hormone_threshold", g.stress_hormone_threshold);
    get_f("stress_hormone_production_rate", g.stress_hormone_production_rate);
    get_f("stress_hormone_diffusion_rate", g.stress_hormone_diffusion_rate);
    get_f("stress_hormone_decay_rate", g.stress_hormone_decay_rate);
    get_f("stress_thickening_boost", g.stress_thickening_boost);
    get_f("stress_elongation_inhibition", g.stress_elongation_inhibition);
    get_f("stress_gravitropism_boost", g.stress_gravitropism_boost);
    get_f("vascular_radius_threshold", g.vascular_radius_threshold);

    return g;
}

static std::vector<std::string> scan_genome_files(const std::string& dir) {
    std::vector<std::string> names;
    names.push_back("default");
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() == ".txt") {
            names.push_back(entry.path().stem().string());
        }
    }
    return names;
}

static Renderer* g_renderer = nullptr;
static const Node* g_selected_node = nullptr;
static const Node* g_hovered_node = nullptr;
static bool g_show_node_panel = false;
static bool g_show_economy_modal = false;
static int g_economy_chemical_idx = 3; // default to Sugar

static const char* chemical_name(ChemicalID id) {
    switch (id) {
        case ChemicalID::Auxin:       return "Auxin";
        case ChemicalID::Cytokinin:   return "Cytokinin";
        case ChemicalID::Gibberellin: return "Gibberellin";
        case ChemicalID::Sugar:       return "Sugar";
        case ChemicalID::Ethylene:    return "Ethylene";
        case ChemicalID::Stress:      return "Stress";
        case ChemicalID::Water:       return "Water";
        case ChemicalID::Count:       break;
    }
    return "?";
}

static const char* node_type_label(NodeType t) {
    switch (t) {
        case NodeType::STEM:           return "Stem";
        case NodeType::ROOT:           return "Root";
        case NodeType::LEAF:           return "Leaf";
        case NodeType::APICAL:         return "SAM";
        case NodeType::ROOT_APICAL:    return "RAM";
    }
    return "?";
}

// Find the node closest to a screen-space click via ray casting
static const Node* pick_node(const Plant& plant, const OrbitCamera& camera,
                              int screen_x, int screen_y, int width, int height) {
    float aspect = static_cast<float>(width) / static_cast<float>(height);
    glm::mat4 view = camera.view_matrix();
    glm::mat4 proj = camera.projection_matrix(aspect);
    glm::mat4 vp = proj * view;

    // Convert screen coords to NDC
    float ndc_x = (2.0f * screen_x / width) - 1.0f;
    float ndc_y = 1.0f - (2.0f * screen_y / height);

    // Unproject to get ray direction
    glm::mat4 inv_vp = glm::inverse(vp);
    glm::vec4 near_clip = inv_vp * glm::vec4(ndc_x, ndc_y, -1.0f, 1.0f);
    glm::vec4 far_clip = inv_vp * glm::vec4(ndc_x, ndc_y, 1.0f, 1.0f);
    glm::vec3 near_pt = glm::vec3(near_clip) / near_clip.w;
    glm::vec3 far_pt = glm::vec3(far_clip) / far_clip.w;
    glm::vec3 ray_dir = glm::normalize(far_pt - near_pt);
    glm::vec3 ray_origin = near_pt;

    const Node* closest = nullptr;
    float closest_ray_dist = FLT_MAX;

    // Min distance from ray to segment [seg_a, seg_b].
    // Handles the segment as a whole, not just its endpoints.
    auto ray_seg_dist = [&](glm::vec3 seg_a, glm::vec3 seg_b) -> float {
        glm::vec3 e = seg_b - seg_a;
        float e_len2 = glm::dot(e, e);
        if (e_len2 < 1e-8f) {
            float t = std::max(0.0f, glm::dot(seg_a - ray_origin, ray_dir));
            return glm::length(ray_origin + t * ray_dir - seg_a);
        }
        glm::vec3 w = ray_origin - seg_a;
        float b    = glm::dot(ray_dir, e);
        float ew   = glm::dot(w, e);
        float dw   = glm::dot(ray_dir, w);
        float denom = e_len2 - b * b;
        float s = (std::abs(denom) < 1e-7f) ? 0.0f
                : glm::clamp((ew - b * dw) / denom, 0.0f, 1.0f);
        glm::vec3 seg_pt = seg_a + s * e;
        float t = std::max(0.0f, glm::dot(seg_pt - ray_origin, ray_dir));
        return glm::length(ray_origin + t * ray_dir - seg_pt);
    };

    plant.for_each_node([&](const Node& node) {
        float dist, pick_radius;

        bool is_segment = node.parent &&
                          (node.type == NodeType::STEM || node.type == NodeType::ROOT);
        if (is_segment) {
            // Test the full segment parent→node, not just the tip.
            dist = ray_seg_dist(node.parent->position, node.position);
            pick_radius = std::max(node.radius * 3.0f, 0.05f);
        } else {
            // Point/sphere test for leaves and meristems.
            glm::vec3 to_node = node.position - ray_origin;
            float t = glm::dot(to_node, ray_dir);
            if (t < 0.0f) return;
            dist = glm::length(node.position - (ray_origin + ray_dir * t));
            float eff_r = node.radius;
            if (auto* leaf = node.as_leaf()) {
                if (leaf->leaf_size > 0.0f) eff_r = leaf->leaf_size * 0.5f;
            }
            pick_radius = std::max(eff_r * 2.0f, 0.15f);
        }

        if (dist < pick_radius && dist < closest_ray_dist) {
            closest_ray_dist = dist;
            closest = &node;
        }
    });

    return closest;
}

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main(int argc, char* argv[]) {
    std::string color_chemical;
    std::string world_path;

    bool log_perf = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--color") == 0 && i + 1 < argc) {
            color_chemical = argv[++i];
        } else if (std::strcmp(argv[i], "--world") == 0 && i + 1 < argc) {
            world_path = argv[++i];
        } else if (std::strcmp(argv[i], "--logperf") == 0) {
            log_perf = true;
        }
    }

    Engine engine;
    if (!world_path.empty()) {
        engine.world_params_mut() = load_world_params(world_path);
    }
    const std::string genome_dir = "src/data";
    auto genome_names = scan_genome_files(genome_dir);
    int selected_genome = 0;  // 0 = "default"

    Genome g = default_genome();
    PlantID plant_id = engine.create_plant(g, glm::vec3(0.0f));
    engine.debug_log().open("debug/debug_log.csv");
    if (log_perf) {
        engine.perf_log().open("perf_log.csv");
        std::cout << "Performance logging enabled -> perf_log.csv" << std::endl;
    }

    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;

    // If GPU shadow system initialized, disable the CPU fallback so they don't fight.
    if (renderer.light_system().is_initialized()) {
        engine.world_params_mut().cpu_light_enabled = false;
    }

    if (!color_chemical.empty()) {
        if (color_chemical == "type") {
            renderer.set_color_by_type(true);
        } else {
            ChemicalAccessor accessor;
            if (color_chemical == "auxin") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Auxin); };
            } else if (color_chemical == "cytokinin") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Cytokinin); };
            } else if (color_chemical == "sugar") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Sugar); };
            } else if (color_chemical == "gibberellin") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Gibberellin); };
            } else if (color_chemical == "ethylene") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Ethylene); };
            } else if (color_chemical == "water") {
                accessor = [](const Node& n) { return n.local().chemical(ChemicalID::Water); };
            } else if (color_chemical == "vascular") {
                accessor = [](const Node& n) { return n.get_parent_auxin_flow_bias(); };
            } else {
                std::cerr << "Unknown chemical: " << color_chemical
                          << " (available: auxin, cytokinin, sugar, gibberellin, ethylene, water, vascular, type)" << std::endl;
                renderer.shutdown();
                return 1;
            }
            renderer.set_color_mode(std::move(accessor));
        }
    }

    GLFWwindow* window = renderer.window();
    glfwSetScrollCallback(window, scroll_callback);

    // ImGui setup
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    ImGui::StyleColorsDark();

    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE, STRESS, WATER, GROWTH, VASCULAR };
    Overlay active_overlay = Overlay::NONE;
    bool playing = false;
    int steps_remaining = 0;
    bool space_was_pressed = false;
    bool mouse_was_pressed = false;
    int total_ticks = 0;

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);

        // Arrow keys rotate camera
        const float rotate_speed = 6.0f;
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
            renderer.camera().rotate(-rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
            renderer.camera().rotate(rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
            renderer.camera().rotate(0.0f, rotate_speed);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
            renderer.camera().rotate(0.0f, -rotate_speed);

        // Z/X translate camera target vertically
        const float pan_speed = 0.05f;
        if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() + glm::vec3(0.0f, pan_speed, 0.0f));
        if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() - glm::vec3(0.0f, pan_speed, 0.0f));

        // Space toggles play/pause
        bool space_down = glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS;
        if (space_down && !space_was_pressed) {
            if (playing) {
                playing = false;
            } else {
                playing = true;
                steps_remaining = 0;
            }
        }
        space_was_pressed = space_down;

        // Run sim
        if (playing) {
            engine.tick();
            total_ticks++;
        } else if (steps_remaining > 0) {
            engine.tick();
            total_ticks++;
            steps_remaining--;
        }

        // Sun angle controls — static so they persist across frames.
        static float sun_elevation = 90.0f;  // degrees: 90 = overhead, 5 = grazing
        static float sun_azimuth   =  0.0f;  // degrees: 0..360

        // Apply sun direction every frame so the light update and leaf phototropism
        // always use the latest slider values.
        {
            float el = glm::radians(sun_elevation);
            float az = glm::radians(sun_azimuth);
            glm::vec3 sun_dir(
                std::cos(az) * std::cos(el),
                -std::sin(el),
                std::sin(az) * std::cos(el)
            );
            // sun_dir already normalized by construction (cos²+sin² = 1).
            engine.world_params_mut().sun_direction = sun_dir;
            if (renderer.light_system().is_initialized())
                renderer.light_system().set_sun_direction(sun_dir);
        }

        // GPU light update — every 24 ticks (once per sim-day)
        static int last_light_update = 0;
        if (renderer.light_system().is_initialized() && (total_ticks - last_light_update) >= 24) {
            last_light_update = total_ticks;
            renderer.light_system().update(engine.all_plants());
        }

        // ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Hover + click picking
        double mx, my;
        glfwGetCursorPos(window, &mx, &my);
        int w, h;
        glfwGetWindowSize(window, &w, &h);

        if (!ImGui::GetIO().WantCaptureMouse) {
            g_hovered_node = pick_node(engine.get_plant(plant_id),
                                       renderer.camera(),
                                       static_cast<int>(mx), static_cast<int>(my), w, h);

            bool mouse_down = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
            if (mouse_down && !mouse_was_pressed && g_hovered_node) {
                g_selected_node = g_hovered_node;
                g_show_node_panel = true;
            }
            mouse_was_pressed = mouse_down;
        } else {
            g_hovered_node = nullptr;
        }

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::SetNextWindowSizeConstraints(ImVec2(0, 0), ImVec2(FLT_MAX, FLT_MAX));
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        int hours = total_ticks;
        int days = hours / 24;
        ImGui::Text("Tick: %d (%dd %dh)  Nodes: %d", total_ticks, days, hours % 24,
                     static_cast<int>(engine.get_plant(plant_id).node_count()));

        // Node & sugar stats
        float total_sugar = 0.0f;
        float max_sugar = 0.0f;
        float total_maintenance = 0.0f;
        int stem_count = 0, root_count = 0, leaf_count = 0;
        int shoot_active = 0, shoot_dormant = 0;
        int root_active = 0, root_dormant = 0;
        const Genome& pg = engine.get_plant(plant_id).genome();
        engine.get_plant(plant_id).for_each_node([&](const Node& n) {
            total_sugar += n.local().chemical(ChemicalID::Sugar);
            if (n.local().chemical(ChemicalID::Sugar) > max_sugar) max_sugar = n.local().chemical(ChemicalID::Sugar);
            total_maintenance += n.maintenance_cost(engine.world_params());
            switch (n.type) {
                case NodeType::STEM: stem_count++; break;
                case NodeType::ROOT: root_count++; break;
                case NodeType::LEAF: leaf_count++; break;
                case NodeType::APICAL:
                    if (n.as_apical()->active) shoot_active++; else shoot_dormant++;
                    break;
                case NodeType::ROOT_APICAL:
                    if (n.as_root_apical()->active) root_active++; else root_dormant++;
                    break;
            }
        });

        // Per-tick sugar production (delta of lifetime accumulator)
        float current_produced = engine.get_plant(plant_id).total_sugar_produced();
        static float prev_produced = 0.0f;
        float last_tick_production = current_produced - prev_produced;
        prev_produced = current_produced;

        ImGui::Text("Stem: %d  Root: %d  Leaf: %d", stem_count, root_count, leaf_count);
        ImGui::Text("Shoot: %d active %d dormant  Root: %d active %d dormant",
                     shoot_active, shoot_dormant, root_active, root_dormant);
        ImGui::Text("Sugar: total=%s  max=%s", fmt_mass(total_sugar), fmt_mass(max_sugar));
        ImGui::Text("Production: %s  Maintenance: %s", fmt_mass_rate(last_tick_production), fmt_mass_rate(total_maintenance));
        ImGui::Separator();
        ImGui::SliderFloat("Light Level", &engine.world_params_mut().light_level, 0.0f, 2.0f);
        static bool show_shadow_map = false;
        ImGui::Checkbox("Show Shadow Map", &show_shadow_map);
        static bool show_gpu_shadow_debug = false;
        static int debug_slice_idx = 0;
        ImGui::Checkbox("GPU Shadow Debug", &show_gpu_shadow_debug);
        ImGui::Checkbox("Phloem Debug Log (debug/phloem_log.csv)",
                        &engine.world_params_mut().phloem_debug_log);
        if (engine.world_params().phloem_debug_log)
            ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                               "Logging to debug/phloem_log.csv");
        if (show_gpu_shadow_debug && renderer.light_system().is_initialized()) {
            ImGui::SliderInt("Slice", &debug_slice_idx, 0, LightSystem::NUM_SLICES - 1);
        }
        static bool show_ground_shadow = false;
        if (renderer.light_system().is_initialized()) {
            ImGui::SliderFloat("Sun Elevation", &sun_elevation, 5.0f, 90.0f, "%.0f deg");
            ImGui::SliderFloat("Sun Azimuth",   &sun_azimuth,   0.0f, 360.0f, "%.0f deg");
            ImGui::Checkbox("Show Ground Shadow", &show_ground_shadow);
        }
        ImGui::SeparatorText("Time");
        if (ImGui::Button("Step 1")) {
            playing = false;
            steps_remaining = 1;
        }
        ImGui::SameLine();
        if (ImGui::Button("Step 10")) {
            playing = false;
            steps_remaining = 10;
        }
        ImGui::SameLine();
        if (ImGui::Button("Step 100")) {
            playing = false;
            steps_remaining = 100;
        }
        ImGui::SameLine();
        if (ImGui::Button(playing ? "Pause" : "Play")) {
            if (playing) {
                playing = false;
            } else {
                playing = true;
                steps_remaining = 0;
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            engine.reset();
            plant_id = engine.create_plant(g, glm::vec3(0.0f));
            playing = false;
            steps_remaining = 0;
            total_ticks = 0;
            prev_produced = 0.0f;
            g_selected_node = nullptr;
            g_show_node_panel = false;
        }

        // Genome selector
        if (ImGui::BeginCombo("Genome", genome_names[selected_genome].c_str())) {
            for (int i = 0; i < static_cast<int>(genome_names.size()); i++) {
                bool is_selected = (selected_genome == i);
                if (ImGui::Selectable(genome_names[i].c_str(), is_selected)) {
                    selected_genome = i;
                    if (genome_names[i] == "default") {
                        g = default_genome();
                    } else {
                        g = load_genome_file(genome_dir + "/" + genome_names[i] + ".txt");
                    }
                    engine.reset();
                    plant_id = engine.create_plant(g, glm::vec3(0.0f));
                    playing = false;
                    steps_remaining = 0;
                    total_ticks = 0;
                    prev_produced = 0.0f;
                    g_selected_node = nullptr;
                    g_show_node_panel = false;
                }
                if (is_selected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (ImGui::CollapsingHeader("Overlays")) {
            if (ImGui::Button("None")) {
                renderer.set_color_mode(ChemicalAccessor{});
                renderer.set_color_by_type(false);
                active_overlay = Overlay::NONE;
            }
            ImGui::SameLine();
            if (ImGui::Button("Node Type")) {
                renderer.set_color_mode(ChemicalAccessor{});
                renderer.set_color_by_type(true);
                active_overlay = Overlay::NODE_TYPE;
            }
            ImGui::SameLine();
            if (ImGui::Button("Auxin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Auxin); });
                active_overlay = Overlay::AUXIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Cytokinin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Cytokinin); });
                active_overlay = Overlay::CYTOKININ;
            }
            ImGui::SameLine();
            if (ImGui::Button("Sugar")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Sugar); });
                active_overlay = Overlay::SUGAR;
            }
            ImGui::SameLine();
            if (ImGui::Button("Light")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { auto* l = n.as_leaf(); return l ? l->light_exposure : 0.0f; });
                active_overlay = Overlay::LIGHT;
            }
            if (ImGui::Button("GA")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Gibberellin); });
                active_overlay = Overlay::GIBBERELLIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Ethylene")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Ethylene); });
                active_overlay = Overlay::ETHYLENE;
            }
            ImGui::SameLine();
            if (ImGui::Button("Stress")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.stress; });
                active_overlay = Overlay::STRESS;
            }
            ImGui::SameLine();
            if (ImGui::Button("Water")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.local().chemical(ChemicalID::Water); });
                active_overlay = Overlay::WATER;
            }
            ImGui::SameLine();
            if (ImGui::Button("Growth")) {
                renderer.set_color_by_type(false);
                const Genome& og = engine.get_plant(plant_id).genome();
                const WorldParams& ow = engine.world_params();
                renderer.set_color_mode([og, ow](const Node& n) -> float {
                    using namespace meristem_helpers;
                    if (auto* ap = n.as_apical()) {
                        if (!ap->active) return 0.0f;
                        float max_cost = og.growth_rate * ow.sugar_cost_meristem_growth;
                        float gf = growth_fraction(n.local().chemical(ChemicalID::Sugar), max_cost,
                                                   n.local().chemical(ChemicalID::Cytokinin), og.cytokinin_growth_threshold);
                        float wgf = turgor_fraction(n.local().chemical(ChemicalID::Water), water_cap(n, og));
                        return gf * wgf;
                    } else if (auto* ra = n.as_root_apical()) {
                        if (!ra->active) return 0.0f;
                        float max_cost = og.root_growth_rate * ow.sugar_cost_root_growth;
                        float gf = sugar_growth_fraction(n.local().chemical(ChemicalID::Sugar), max_cost);
                        float wgf = turgor_fraction(n.local().chemical(ChemicalID::Water), water_cap(n, og));
                        return gf * wgf;
                    }
                    return 0.0f;  // non-meristems: dark
                });
                active_overlay = Overlay::GROWTH;
            }
            ImGui::SameLine();
            if (ImGui::Button("Vascular")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.get_parent_auxin_flow_bias(); });
                active_overlay = Overlay::VASCULAR;
            }

            if (active_overlay == Overlay::NODE_TYPE) {
                ImGui::Spacing();
                ImGui::TextColored(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), "Shoot Apical");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.0f, 1.0f), "Shoot Axillary");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.2f, 0.8f, 0.2f, 1.0f), "Stem");
                ImGui::TextColored(ImVec4(0.2f, 0.4f, 1.0f, 1.0f), "Root Apical");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.6f, 0.2f, 0.8f, 1.0f), "Root Axillary");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.8f, 0.6f, 0.2f, 1.0f), "Root");
            } else if (active_overlay != Overlay::NONE) {
                ImGui::Spacing();
                ImGui::TextColored(ImVec4(0.0f, 0.0f, 1.0f, 1.0f), "Low");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.0f, 1.0f, 1.0f, 1.0f), "-");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "-");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "-");
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "High");
            }
        }

        if (ImGui::CollapsingHeader("Info")) {
            if (ImGui::Button("Economy"))
                g_show_economy_modal = true;
        }

        ImGui::End();

        // Economy Modal
        if (g_show_economy_modal) {
            ImVec2 display = ImGui::GetIO().DisplaySize;
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(display);
            ImGui::Begin("Economy", &g_show_economy_modal,
                ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoCollapse);

            // Chemical selector
            const auto& chem_ids = all_chemical_ids;
            if (ImGui::BeginCombo("Chemical", chemical_name(chem_ids[g_economy_chemical_idx]))) {
                for (int i = 0; i < static_cast<int>(chem_ids.size()); i++) {
                    bool selected = (i == g_economy_chemical_idx);
                    if (ImGui::Selectable(chemical_name(chem_ids[i]), selected))
                        g_economy_chemical_idx = i;
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            ChemicalID chem = chem_ids[g_economy_chemical_idx];
            const Genome& eg = engine.get_plant(plant_id).genome();
            const WorldParams& ew = engine.world_params();

            ImGui::Separator();

            // Collect rows
            struct EconRow {
                const Node* node;
                NodeType type; uint32_t id; uint32_t age;
                float level, produced, consumed, net;
            };
            std::vector<EconRow> rows;
            engine.get_plant(plant_id).for_each_node([&](const Node& n) {
                float level = n.local().chemical(chem);
                float produced = 0.0f, consumed = 0.0f;

                if (chem == ChemicalID::Sugar) {
                    consumed = n.maintenance_cost(ew);
                    if (auto* leaf = n.as_leaf()) {
                        if (leaf->leaf_size > 1e-6f && leaf->senescence_ticks == 0) {
                            float angle_eff = 1.0f;
                            float flen = glm::length(leaf->facing);
                            if (flen > 1e-4f) angle_eff = std::max(0.0f, (leaf->facing / flen).y);
                            float area = leaf->leaf_size * leaf->leaf_size;
                            produced = leaf->light_exposure * angle_eff * ew.light_level * area * ew.sugar_production_rate;
                        }
                    }
                } else if (chem == ChemicalID::Auxin) {
                    consumed = level * eg.auxin_decay_rate;
                    if (n.type == NodeType::APICAL) produced = eg.apical_auxin_baseline;
                } else if (chem == ChemicalID::Cytokinin) {
                    consumed = level * eg.cytokinin_decay_rate;
                    if (auto* leaf = n.as_leaf()) {
                        if (leaf->leaf_size > 1e-6f && leaf->senescence_ticks == 0) {
                            float angle_eff = 1.0f;
                            float flen = glm::length(leaf->facing);
                            if (flen > 1e-4f) angle_eff = std::max(0.0f, (leaf->facing / flen).y);
                            float area = leaf->leaf_size * leaf->leaf_size;
                            float sugar_prod = leaf->light_exposure * angle_eff * ew.light_level * area * ew.sugar_production_rate;
                            produced = sugar_prod * eg.cytokinin_production_rate;
                        }
                    }
                } else if (chem == ChemicalID::Gibberellin) {
                    consumed = level * eg.ga_decay_rate;
                    if (auto* leaf = n.as_leaf()) {
                        if (n.age < eg.ga_leaf_age_max && leaf->leaf_size > 1e-6f && leaf->senescence_ticks == 0)
                            produced = leaf->leaf_size * eg.ga_production_rate;
                    }
                } else if (chem == ChemicalID::Ethylene) {
                    produced = level;
                } else if (chem == ChemicalID::Stress) {
                    consumed = level * eg.stress_hormone_decay_rate;
                    float bs = eg.wood_density * ew.break_strength_factor;
                    float sr = (bs > 1e-6f) ? n.stress / bs : 0.0f;
                    if (sr > eg.stress_hormone_threshold) {
                        float excess = (sr - eg.stress_hormone_threshold)
                                     / (1.0f - eg.stress_hormone_threshold);
                        produced = excess * eg.stress_hormone_production_rate;
                    }
                }

                rows.push_back({&n, n.type, n.id, n.age, level, produced, consumed, produced - consumed});
            });

            if (ImGui::BeginTable("economy_table", 7,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_Sortable | ImGuiTableFlags_ScrollY,
                    ImVec2(0, display.y - 100))) {
                ImGui::TableSetupColumn("Type",     ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("ID",       ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("Age",      ImGuiTableColumnFlags_WidthFixed, 50);
                ImGui::TableSetupColumn("Level",    ImGuiTableColumnFlags_WidthFixed, 80);
                ImGui::TableSetupColumn("Produced", ImGuiTableColumnFlags_WidthFixed, 80);
                ImGui::TableSetupColumn("Consumed", ImGuiTableColumnFlags_WidthFixed, 80);
                ImGui::TableSetupColumn("Net",      ImGuiTableColumnFlags_WidthFixed, 80);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                // Sort rows based on clicked column
                if (ImGuiTableSortSpecs* specs = ImGui::TableGetSortSpecs()) {
                    if (specs->SpecsCount > 0) {
                        auto& s = specs->Specs[0];
                        int col = s.ColumnIndex;
                        bool asc = (s.SortDirection == ImGuiSortDirection_Ascending);
                        std::sort(rows.begin(), rows.end(), [col, asc](const EconRow& a, const EconRow& b) {
                            float va = 0, vb = 0;
                            switch (col) {
                                case 0: va = static_cast<float>(a.type); vb = static_cast<float>(b.type); break;
                                case 1: va = static_cast<float>(a.id);   vb = static_cast<float>(b.id); break;
                                case 2: va = static_cast<float>(a.age);  vb = static_cast<float>(b.age); break;
                                case 3: va = a.level;    vb = b.level; break;
                                case 4: va = a.produced; vb = b.produced; break;
                                case 5: va = a.consumed; vb = b.consumed; break;
                                case 6: va = a.net;      vb = b.net; break;
                            }
                            return asc ? va < vb : va > vb;
                        });
                        specs->SpecsDirty = false;
                    }
                }

                // Format helpers — sugar uses metric mass, hormones use raw decimals
                bool is_sugar = (chem == ChemicalID::Sugar);
                auto fmt_lvl = [is_sugar](float v) { return is_sugar ? fmt_mass(v) : fmt_au(v); };
                auto fmt_rate = [is_sugar](float v) { return is_sugar ? fmt_mass_rate(v) : fmt_au(v); };

                for (const auto& r : rows) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);

                    char label[32];
                    std::snprintf(label, sizeof(label), "##node_%u", r.id);
                    if (ImGui::Selectable(label, g_selected_node == r.node,
                            ImGuiSelectableFlags_SpanAllColumns)) {
                        g_selected_node = r.node;
                        g_show_node_panel = true;
                    }
                    ImGui::SameLine();
                    ImGui::Text("%s", node_type_label(r.type));

                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%u", r.id);
                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%u", r.age);
                    ImGui::TableSetColumnIndex(3);
                    ImGui::Text("%s", fmt_lvl(r.level));
                    ImGui::TableSetColumnIndex(4);
                    if (r.produced > 1e-8f) ImGui::TextColored(ImVec4(0.2f, 0.8f, 0.2f, 1.0f), "%s", fmt_rate(r.produced));
                    else ImGui::Text("-");
                    ImGui::TableSetColumnIndex(5);
                    if (r.consumed > 1e-8f) ImGui::TextColored(ImVec4(0.8f, 0.2f, 0.2f, 1.0f), "%s", fmt_rate(r.consumed));
                    else ImGui::Text("-");
                    ImGui::TableSetColumnIndex(6);
                    if (std::abs(r.net) > 1e-8f) {
                        ImVec4 col = r.net > 0 ? ImVec4(0.2f, 0.8f, 0.2f, 1.0f) : ImVec4(0.8f, 0.2f, 0.2f, 1.0f);
                        ImGui::TextColored(col, "%s%s", r.net > 0 ? "+" : "", fmt_rate(r.net));
                    } else {
                        ImGui::Text("-");
                    }
                }

                ImGui::EndTable();
            }

            ImGui::End();
        }

        // Node Inspector Panel (top-right)
        if (g_show_node_panel && g_selected_node) {
            // Verify the node still exists (it may have been pruned)
            bool node_exists = false;
            engine.get_plant(plant_id).for_each_node([&](const Node& n) {
                if (&n == g_selected_node) node_exists = true;
            });
            if (!node_exists) {
                g_selected_node = nullptr;
                g_show_node_panel = false;
            }
        }

        if (g_show_node_panel && g_selected_node) {
            int w, h;
            glfwGetWindowSize(window, &w, &h);
            ImGui::SetNextWindowPos(ImVec2(static_cast<float>(w) - 320, 10), ImGuiCond_Once);
            ImGui::SetNextWindowSize(ImVec2(480, 0), ImGuiCond_Once);
            if (g_show_economy_modal) ImGui::SetNextWindowFocus();

            if (ImGui::Begin("Node Inspector", &g_show_node_panel, ImGuiWindowFlags_AlwaysAutoResize)) {
                const Node& sel = *g_selected_node;

                // --- IDENTITY ---
                {
                    const char* type_str = "?";
                    switch (sel.type) {
                        case NodeType::STEM:        type_str = "STEM"; break;
                        case NodeType::ROOT:        type_str = "ROOT"; break;
                        case NodeType::LEAF:        type_str = "LEAF"; break;
                        case NodeType::APICAL:      type_str = "APICAL"; break;
                        case NodeType::ROOT_APICAL: type_str = "ROOT_APICAL"; break;
                    }
                    ImGui::Text("ID: #%u  Type: %s", sel.id, type_str);
                    ImGui::Text("Age: %u ticks (%.1f days)", sel.age, sel.age / 24.0f);
                    ImGui::Text("Length: %s  Radius: %s", fmt_dist(glm::length(sel.offset)), fmt_dist(sel.radius));
                    float cross_section = 3.14159f * sel.radius * sel.radius;
                    ImGui::Text("Cross-section: %.4f dm\xC2\xB2", cross_section);
                    ImGui::Text("Height (y): %.2f m", sel.position.y / 10.0f);  // dm -> m
                    ImGui::Text("Nodes to seed: %d", nodes_to_seed(sel));

                    const Genome& mg = engine.get_plant(plant_id).genome();
                    if (sel.as_stem() || sel.as_root()) {
                        uint32_t mat_ticks = sel.as_stem()
                            ? mg.internode_maturation_ticks
                            : mg.root_internode_maturation_ticks;
                        bool mature = (sel.parent == nullptr) || (sel.age >= mat_ticks);
                        ImGui::Text("Elongation: %s (age %u/%u)", mature ? "mature" : "growing", sel.age, mat_ticks);
                        float hm = hydraulic_maturity(sel, mg) * 100.0f;
                        ImGui::Text("Hydraulic maturity: %.0f%% closed", hm);
                        ImGui::Text("PIN saturation: %.2f", sel.get_parent_auxin_flow_bias());
                    }
                    ImGui::Text("Starvation: %u ticks", sel.starvation_ticks);
                    if (auto* leaf = sel.as_leaf()) {
                        if (leaf->senescence_ticks > 0)
                            ImGui::Text("Senescence: %u / %u ticks", leaf->senescence_ticks, mg.senescence_duration);
                        else
                            ImGui::Text("Senescence: healthy");
                    }
                }
                ImGui::Separator();
                // --- (existing panel content continues below) ---

                // Node type
                const char* type_str = "?";
                switch (sel.type) {
                    case NodeType::STEM:           type_str = "STEM"; break;
                    case NodeType::ROOT:           type_str = "ROOT"; break;
                    case NodeType::LEAF:           type_str = "LEAF"; break;
                    case NodeType::APICAL:         type_str = "APICAL"; break;
                    case NodeType::ROOT_APICAL:    type_str = "ROOT_APICAL"; break;
                }
                ImGui::Text("Type: %s", type_str);

                // Meristem info with growth factor breakdown
                if (auto* ap = sel.as_apical()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();
                    ImGui::Text("Meristem: %s", ap->active ? "active" : "dormant");
                    if (ap->active) {
                        float max_cost = mg.growth_rate * mw.sugar_cost_meristem_growth;
                        float sugar_gf = (max_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
                        float cyt = sel.local().chemical(ChemicalID::Cytokinin);
                        float cyt_gf = cyt / (cyt + std::max(mg.cytokinin_growth_threshold, 1e-6f));
                        float water_gf = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                        float total = sugar_gf * cyt_gf * water_gf;
                        ImGui::Text("Growth: %.1f%%", total * 100);
                        ImGui::Text("  Sugar: %3.0f%%  Cyt: %3.0f%%  Water: %3.0f%%",
                                    sugar_gf * 100, cyt_gf * 100, water_gf * 100);
                    } else {
                        // Activation conditions for dormant shoot apicals — shown as % toward each threshold
                        float stem_auxin = sel.parent ? sel.parent->local().chemical(ChemicalID::Auxin) : sel.local().chemical(ChemicalID::Auxin);
                        float local_cyt = sel.parent ? sel.parent->local().chemical(ChemicalID::Cytokinin) : sel.local().chemical(ChemicalID::Cytokinin);
                        float sugar = sel.local().chemical(ChemicalID::Sugar);
                        bool auxin_ok = stem_auxin < mg.auxin_threshold;
                        bool cyt_ok = local_cyt >= mg.cytokinin_threshold;
                        bool sugar_ok = sugar >= mw.sugar_cost_activation;
                        // Auxin must be LOW: 100% = at zero (perfect), 0% = at/above threshold
                        float auxin_pct = (mg.auxin_threshold > 1e-6f)
                            ? std::max(0.0f, (mg.auxin_threshold - stem_auxin) / mg.auxin_threshold * 100.0f)
                            : (auxin_ok ? 100.0f : 0.0f);
                        // Cyt must be HIGH: 100% = at/above threshold, lower is worse
                        float cyt_pct = (mg.cytokinin_threshold > 1e-6f)
                            ? std::min(local_cyt / mg.cytokinin_threshold * 100.0f, 100.0f)
                            : (cyt_ok ? 100.0f : 0.0f);
                        // Sugar must cover activation cost: 100% = fully funded
                        float sugar_pct = (mw.sugar_cost_activation > 1e-6f)
                            ? std::min(sugar / mw.sugar_cost_activation * 100.0f, 100.0f)
                            : (sugar_ok ? 100.0f : 0.0f);
                        ImGui::Text("Activation: %s", (auxin_ok && cyt_ok && sugar_ok) ? "READY" : "blocked");
                        ImGui::Text("  Auxin: %3.0f%%  Cyt: %3.0f%%  Sugar: %3.0f%%",
                                    auxin_pct, cyt_pct, sugar_pct);
                    }
                } else if (auto* ra = sel.as_root_apical()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();
                    ImGui::Text("Meristem: %s", ra->active ? "active" : "dormant");
                    if (ra->active) {
                        float max_cost = mg.root_growth_rate * mw.sugar_cost_root_growth;
                        float sugar_gf = (max_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
                        float water_gf = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                        float total = sugar_gf * water_gf;
                        ImGui::Text("Growth: %.1f%%", total * 100);
                        ImGui::Text("  Sugar: %3.0f%%  Water: %3.0f%%",
                                    sugar_gf * 100, water_gf * 100);
                    } else {
                        // Activation conditions for dormant root apicals — shown as % toward each threshold
                        float auxin = sel.local().chemical(ChemicalID::Auxin);
                        float cyt = sel.local().chemical(ChemicalID::Cytokinin);
                        float sugar = sel.local().chemical(ChemicalID::Sugar);
                        bool auxin_ok = auxin >= mg.root_auxin_activation_threshold;
                        bool cyt_ok = cyt <= mg.root_cytokinin_inhibition_threshold;
                        bool sugar_ok = sugar >= mw.sugar_cost_activation;
                        // Auxin must be HIGH: 100% = at/above threshold
                        float auxin_pct = (mg.root_auxin_activation_threshold > 1e-6f)
                            ? std::min(auxin / mg.root_auxin_activation_threshold * 100.0f, 100.0f)
                            : (auxin_ok ? 100.0f : 0.0f);
                        // Cyt must be LOW: 100% = at zero (perfect), 0% = at/above threshold
                        float cyt_pct = (mg.root_cytokinin_inhibition_threshold > 1e-6f)
                            ? std::max(0.0f, (mg.root_cytokinin_inhibition_threshold - cyt) / mg.root_cytokinin_inhibition_threshold * 100.0f)
                            : (cyt_ok ? 100.0f : 0.0f);
                        // Sugar must cover activation cost: 100% = fully funded
                        float sugar_pct = (mw.sugar_cost_activation > 1e-6f)
                            ? std::min(sugar / mw.sugar_cost_activation * 100.0f, 100.0f)
                            : (sugar_ok ? 100.0f : 0.0f);
                        ImGui::Text("Activation: %s", (auxin_ok && cyt_ok && sugar_ok) ? "READY" : "blocked");
                        ImGui::Text("  Auxin: %3.0f%%  Cyt: %3.0f%%  Sugar: %3.0f%%",
                                    auxin_pct, cyt_pct, sugar_pct);
                    }
                } else if (auto* leaf = sel.as_leaf()) {
                    // Leaf expansion growth fraction
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();
                    if (leaf->leaf_size >= mg.max_leaf_size) {
                        ImGui::Text("Growth: 0%% (full size)");
                    } else {
                        float auxin = sel.local().chemical(ChemicalID::Auxin);
                        float auxin_sat = auxin / (auxin + std::max(mg.leaf_auxin_half_saturation, 1e-6f));
                        float auxin_boost = 1.0f + mg.leaf_auxin_max_boost * auxin_sat;
                        float remaining = mg.max_leaf_size - leaf->leaf_size;
                        float max_growth = mg.leaf_growth_rate * auxin_boost;
                        float water_gf = meristem_helpers::turgor_fraction(
                            sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                        float growth_cap = std::min(max_growth, remaining) * water_gf;
                        float sugar_gf = (growth_cap > 1e-6f)
                            ? std::min(sel.local().chemical(ChemicalID::Sugar) / (growth_cap * mw.sugar_cost_leaf_growth), 1.0f)
                            : 1.0f;
                        ImGui::Text("Growth: %.1f%%", water_gf * sugar_gf * 100.0f);
                        ImGui::Text("  Sugar: %3.0f%%  Water: %3.0f%%  Auxin: %3.0f%%",
                                    sugar_gf * 100.0f, water_gf * 100.0f, auxin_sat * 100.0f);
                    }
                } else if (sel.as_stem()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();
                    float density_scale = mg.wood_density / mw.reference_wood_density;

                    // Elongation (rate-limited by sugar + water; modulated by GA/eth/auxin/stress)
                    if (!sel.parent || sel.age >= mg.internode_maturation_ticks) {
                        ImGui::Text("Elongate: 0%% (mature)");
                    } else {
                        float ga_mult = 1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_elongation_sensitivity;
                        float eth_fac = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Ethylene) * mg.ethylene_elongation_inhibition);
                        float str_fac = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Stress) * mg.stress_elongation_inhibition);
                        float axn_boost = meristem_helpers::auxin_growth_factor(
                            sel.local().chemical(ChemicalID::Auxin), mg.stem_auxin_max_boost, mg.stem_auxin_half_saturation);
                        float eff_rate = mg.internode_elongation_rate * ga_mult * eth_fac * str_fac * axn_boost;
                        float eff_cost = eff_rate * mw.sugar_cost_stem_growth * density_scale;
                        float sugar_gf = (eff_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / eff_cost, 1.0f) : 1.0f;
                        float water_gf = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                        ImGui::Text("Elongate: %.1f%%", sugar_gf * water_gf * 100.0f);
                        ImGui::Text("  Sugar: %3.0f%%  Water: %3.0f%%", sugar_gf * 100.0f, water_gf * 100.0f);
                    }

                    // Thickening (bias-driven, sugar-limited)
                    float bias = sel.get_parent_auxin_flow_bias();
                    if (bias < 1e-6f) {
                        ImGui::Text("Thicken:  0%% (no flux)");
                    } else {
                        float thicken_cost = mg.cambium_responsiveness * bias * mw.sugar_cost_stem_growth * density_scale;
                        float sugar_gf = (thicken_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / thicken_cost, 1.0f) : 1.0f;
                        ImGui::Text("Thicken:  %.1f%%", sugar_gf * 100.0f);
                        ImGui::Text("  Sugar: %3.0f%%  Bias: %.4f", sugar_gf * 100.0f, bias);
                    }
                } else if (sel.as_root()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();

                    // Elongation (rate-limited by sugar + water; modulated by GA/eth/auxin)
                    if (!sel.parent || sel.age >= mg.root_internode_maturation_ticks) {
                        ImGui::Text("Elongate: 0%% (mature)");
                    } else {
                        float ga_mult = 1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_elongation_sensitivity;
                        float eth_fac = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Ethylene) * mg.ethylene_elongation_inhibition);
                        float axn_boost = meristem_helpers::auxin_growth_factor(
                            sel.local().chemical(ChemicalID::Auxin), mg.root_auxin_max_boost, mg.root_auxin_half_saturation);
                        float eff_rate = mg.root_internode_elongation_rate * ga_mult * eth_fac * axn_boost;
                        float eff_cost = eff_rate * mw.sugar_cost_stem_growth;
                        float sugar_gf = (eff_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / eff_cost, 1.0f) : 1.0f;
                        float water_gf = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                        ImGui::Text("Elongate: %.1f%%", sugar_gf * water_gf * 100.0f);
                        ImGui::Text("  Sugar: %3.0f%%  Water: %3.0f%%", sugar_gf * 100.0f, water_gf * 100.0f);
                    }

                    // Thickening (bias-driven, sugar-limited; no density scale for roots)
                    float bias = sel.get_parent_auxin_flow_bias();
                    if (bias < 1e-6f) {
                        ImGui::Text("Thicken:  0%% (no flux)");
                    } else {
                        float thicken_cost = mg.cambium_responsiveness * bias * mw.sugar_cost_stem_growth;
                        float sugar_gf = (thicken_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / thicken_cost, 1.0f) : 1.0f;
                        ImGui::Text("Thicken:  %.1f%%", sugar_gf * 100.0f);
                        ImGui::Text("  Sugar: %3.0f%%  Bias: %.4f", sugar_gf * 100.0f, bias);
                    }
                } else if (sel.as_stem()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();
                    float density_scale = mg.wood_density / mw.reference_wood_density;

                    // Thickening: vascular-history driven, proportional to auxin_flow_bias
                    float bias = sel.get_parent_auxin_flow_bias();
                    float thicken_max_cost = mg.cambium_responsiveness * bias * mw.sugar_cost_stem_growth * density_scale;
                    float thicken_sugar_gf = (thicken_max_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / thicken_max_cost, 1.0f) : 1.0f;
                    float stress_boost = 1.0f + sel.local().chemical(ChemicalID::Stress) * mg.stress_thickening_boost;
                    float thicken_rate = (bias >= 1e-6f) ? (mg.cambium_responsiveness * bias * thicken_sugar_gf * stress_boost) : 0.0f;
                    ImGui::Text("Thicken:  bias=%.4f  rate=%s/tick", bias, fmt_dist(thicken_rate));
                    ImGui::Text("  Sugar: %.0f%%  Stress: %.2fx", thicken_sugar_gf * 100.0f, stress_boost);

                    // Elongation: GA/ethylene/auxin/stress modulated, sugar+water gated
                    float current_len = glm::length(sel.offset);
                    float ga_elong = 1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_elongation_sensitivity;
                    float eth_elong = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Ethylene) * mg.ethylene_elongation_inhibition);
                    float stress_elong = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Stress) * mg.stress_elongation_inhibition);
                    float auxin_elong = meristem_helpers::auxin_growth_factor(
                        sel.local().chemical(ChemicalID::Auxin), mg.stem_auxin_max_boost, mg.stem_auxin_half_saturation);
                    float eff_rate = mg.internode_elongation_rate * ga_elong * eth_elong * stress_elong * auxin_elong;
                    float max_len = mg.max_internode_length * (1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_length_sensitivity);
                    float elong_cost = eff_rate * mw.sugar_cost_stem_growth * density_scale;
                    float elong_sug = (elong_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / elong_cost, 1.0f) : 1.0f;
                    float elong_wat = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), water_cap(sel, mg));
                    if (!sel.parent || sel.age >= mg.internode_maturation_ticks) {
                        ImGui::Text("Elongate: mature (age %u / %u)", sel.age, mg.internode_maturation_ticks);
                    } else {
                        ImGui::Text("Elongate: %s / %s  rate=%s/tick",
                                    fmt_dist(current_len), fmt_dist(max_len),
                                    fmt_dist(eff_rate * elong_sug * elong_wat));
                        ImGui::Text("  GA: %.2fx  Eth: %.0f%%  Auxin: %.2fx  Stress: %.0f%%",
                                    ga_elong, eth_elong * 100.0f, auxin_elong, stress_elong * 100.0f);
                        ImGui::Text("  Sugar: %.0f%%  Water: %.0f%%",
                                    elong_sug * 100.0f, elong_wat * 100.0f);
                    }
                } else if (sel.as_root()) {
                    const Genome& mg = engine.get_plant(plant_id).genome();
                    const WorldParams& mw = engine.world_params();

                    // Water absorption: gradient-based, self-limiting when root is full
                    float root_len = std::max(glm::length(sel.offset), 0.01f);
                    float surface_area = 2.0f * 3.14159f * sel.radius * root_len;
                    float rwcap = water_cap(sel, mg);
                    float fill_frac = (rwcap > 1e-6f) ? sel.local().chemical(ChemicalID::Water) / rwcap : 1.0f;
                    float gradient = std::max(0.0f, mw.soil_moisture - fill_frac);
                    float absorbed = mg.water_absorption_rate * surface_area * gradient;
                    ImGui::Text("Absorb: +%s/tick", fmt_vol(absorbed));
                    ImGui::Text("  Area: %.3f dm\xC2\xB2  Gradient: %.0f%%  Fill: %.0f%%",
                                surface_area, gradient * 100.0f, fill_frac * 100.0f);

                    // Thickening: same vascular-driven model as stem (no density scale for roots)
                    float bias = sel.get_parent_auxin_flow_bias();
                    float thicken_max_cost = mg.cambium_responsiveness * bias * mw.sugar_cost_stem_growth;
                    float thicken_sugar_gf = (thicken_max_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / thicken_max_cost, 1.0f) : 1.0f;
                    float thicken_rate = (bias >= 1e-6f) ? (mg.cambium_responsiveness * bias * thicken_sugar_gf) : 0.0f;
                    ImGui::Text("Thicken:  bias=%.4f  rate=%s/tick", bias, fmt_dist(thicken_rate));
                    ImGui::Text("  Sugar: %.0f%%", thicken_sugar_gf * 100.0f);

                    // Elongation: GA/ethylene/auxin modulated, sugar+water gated
                    float current_len = glm::length(sel.offset);
                    float ga_elong = 1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_elongation_sensitivity;
                    float eth_elong = std::max(0.0f, 1.0f - sel.local().chemical(ChemicalID::Ethylene) * mg.ethylene_elongation_inhibition);
                    float auxin_elong = meristem_helpers::auxin_growth_factor(
                        sel.local().chemical(ChemicalID::Auxin), mg.root_auxin_max_boost, mg.root_auxin_half_saturation);
                    float eff_rate = mg.root_internode_elongation_rate * ga_elong * eth_elong * auxin_elong;
                    float max_len = mg.max_internode_length * (1.0f + sel.local().chemical(ChemicalID::Gibberellin) * mg.ga_length_sensitivity);
                    float elong_cost = eff_rate * mw.sugar_cost_stem_growth;
                    float elong_sug = (elong_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / elong_cost, 1.0f) : 1.0f;
                    float elong_wat = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), rwcap);
                    if (!sel.parent || sel.age >= mg.root_internode_maturation_ticks) {
                        ImGui::Text("Elongate: mature (age %u / %u)", sel.age, mg.root_internode_maturation_ticks);
                    } else {
                        ImGui::Text("Elongate: %s / %s  rate=%s/tick",
                                    fmt_dist(current_len), fmt_dist(max_len),
                                    fmt_dist(eff_rate * elong_sug * elong_wat));
                        ImGui::Text("  GA: %.2fx  Eth: %.0f%%  Auxin: %.2fx",
                                    ga_elong, eth_elong * 100.0f, auxin_elong);
                        ImGui::Text("  Sugar: %.0f%%  Water: %.0f%%",
                                    elong_sug * 100.0f, elong_wat * 100.0f);
                    }
                }

                ImGui::Text("ID: %u  Age: %u", sel.id, sel.age);
                ImGui::Text("Radius: %s", fmt_dist(sel.radius));
                ImGui::Text("Length: %s", fmt_dist(glm::length(sel.offset)));
                if (auto* leaf = sel.as_leaf()) {
                    const Genome& lg = engine.get_plant(plant_id).genome();
                    const WorldParams& lw = engine.world_params();
                    float leaf_area = leaf->leaf_size * leaf->leaf_size;
                    float angle_eff = 1.0f;
                    float flen = glm::length(leaf->facing);
                    if (flen > 1e-4f) angle_eff = std::max(0.0f, (leaf->facing / flen).y);

                    // Leaf size + progress
                    ImGui::Text("Leaf Size: %s / %s (%.0f%%)",
                                fmt_dist(leaf->leaf_size), fmt_dist(lg.max_leaf_size),
                                100.0f * leaf->leaf_size / std::max(lg.max_leaf_size, 1e-6f));
                    ImGui::Text("Light: %.1f%%  Angle eff: %.0f%%",
                                leaf->light_exposure * 100.0f, angle_eff * 100.0f);

                    // Photosynthesis with stomatal conductance
                    float lwcap = water_cap(sel, lg);
                    float stomatal = (lwcap > 1e-6f)
                        ? std::min(std::max(sel.local().chemical(ChemicalID::Water) / lwcap, 0.2f), 1.0f)
                        : 1.0f;
                    float production = leaf->light_exposure * angle_eff
                        * lw.light_level * leaf_area
                        * lw.sugar_production_rate * stomatal;
                    ImGui::Text("Photo: %s/tick  Stomatal: %.0f%%",
                                fmt_mass_rate(production), stomatal * 100.0f);

                    // Transpiration
                    float transpiration = lg.transpiration_rate * leaf_area * leaf->light_exposure;
                    ImGui::Text("Transpire: -%s/tick", fmt_vol(transpiration));

                    // Expansion grow choice
                    if (leaf->leaf_size < lg.max_leaf_size) {
                        float auxin_boost = meristem_helpers::auxin_growth_factor(
                            sel.local().chemical(ChemicalID::Auxin), lg.leaf_auxin_max_boost, lg.leaf_auxin_half_saturation);
                        float water_gf = meristem_helpers::turgor_fraction(sel.local().chemical(ChemicalID::Water), lwcap);
                        float max_growth = lg.leaf_growth_rate * auxin_boost;
                        float expand_cost = max_growth * lw.sugar_cost_leaf_growth;
                        float expand_sug = (expand_cost > 1e-6f) ? std::min(sel.local().chemical(ChemicalID::Sugar) / expand_cost, 1.0f) : 1.0f;
                        ImGui::Text("Expand:  Auxin: %.2fx  Water: %.0f%%  Sugar: %.0f%%",
                                    auxin_boost, water_gf * 100.0f, expand_sug * 100.0f);
                    } else {
                        ImGui::Text("Expand: full size");
                    }

                    // Carbon balance and senescence
                    ImGui::Text("Carbon deficit: %u ticks", leaf->deficit_ticks);
                    if (leaf->senescence_ticks > 0) {
                        ImGui::Text("Senescence: %u / %u ticks", leaf->senescence_ticks, lg.senescence_duration);
                    } else {
                        ImGui::Text("Senescence: healthy");
                    }
                }
                ImGui::Text("Starvation: %u ticks", sel.starvation_ticks);
                ImGui::Text("Children: %d", static_cast<int>(sel.children.size()));

                ImGui::Separator();

                // --- Per-edge transport capacity / actual-flow helpers ---
                // Max Trans.: sum the diffusion-formula cap across all edges touching this
                // node (sel->children and parent->sel).  Uses the same radius_factor and
                // base/scale the engine uses inside compute_transport_flow.
                // Transporting: sum of last-tick magnitudes where we track them
                // (auxin via last_auxin_flux, sugar via tick_sugar_transport).
                // Other chemicals currently have no per-tick flux tracking → "-".
                const Genome& gnm = engine.get_plant(plant_id).genome();
                const float ref_r = gnm.initial_radius;
                const float min_conn_r = ref_r * 0.2f;
                auto edge_max = [&](float child_radius, float base, float scale) -> float {
                    float conn_r = std::max(std::min(child_radius, sel.radius), min_conn_r);
                    float radius_factor = conn_r / std::max(ref_r, 1e-6f);
                    return base + radius_factor * scale;
                };
                auto chem_max_trans = [&](float base, float scale) -> float {
                    float total = 0.0f;
                    for (const Node* c : sel.children) total += edge_max(c->radius, base, scale);
                    if (sel.parent) total += edge_max(sel.radius, base, scale);
                    return total;
                };

                float auxin_flux_sum = 0.0f;
                for (const auto& [cptr, flux] : sel.last_auxin_flux) auxin_flux_sum += std::fabs(flux);
                if (sel.parent) {
                    auto it = sel.parent->last_auxin_flux.find(const_cast<Node*>(&sel));
                    if (it != sel.parent->last_auxin_flux.end()) auxin_flux_sum += std::fabs(it->second);
                }

                float auxin_max   = chem_max_trans(gnm.hormone_base_transport, gnm.hormone_transport_scale);
                float cyt_max     = chem_max_trans(gnm.hormone_base_transport, gnm.hormone_transport_scale);
                float ga_max      = chem_max_trans(gnm.hormone_base_transport, gnm.hormone_transport_scale);
                float eth_max     = chem_max_trans(gnm.hormone_base_transport, gnm.hormone_transport_scale);
                float sugar_max   = chem_max_trans(gnm.sugar_base_transport,   gnm.sugar_transport_scale);
                float water_max   = chem_max_trans(gnm.water_base_transport,   gnm.water_transport_scale);

                ImGui::Text("Chemicals:");
                if (ImGui::BeginTable("chemicals", 10, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
                    ImGui::TableSetupColumn("Chem");
                    ImGui::TableSetupColumn("L lvl"); ImGui::TableSetupColumn("L +"); ImGui::TableSetupColumn("L -");
                    ImGui::TableSetupColumn("V lvl"); ImGui::TableSetupColumn("V +"); ImGui::TableSetupColumn("V -");
                    ImGui::TableSetupColumn("T lvl"); ImGui::TableSetupColumn("T +"); ImGui::TableSetupColumn("T -");
                    ImGui::TableHeadersRow();

                    struct Row { const char* name; ChemicalID id; const char* (*fmt)(float); };
                    Row rows[] = {
                        {"Sugar",  ChemicalID::Sugar,       fmt_mass},
                        {"Water",  ChemicalID::Water,       fmt_vol },
                        {"Auxin",  ChemicalID::Auxin,       fmt_au  },
                        {"Cyt",    ChemicalID::Cytokinin,   fmt_au  },
                        {"GA",     ChemicalID::Gibberellin, fmt_au  },
                        {"Eth",    ChemicalID::Ethylene,    fmt_au  },
                        {"Stress", ChemicalID::Stress,      fmt_au  },
                    };
                    for (const Row& r : rows) {
                        const size_t idx = static_cast<size_t>(r.id);
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::Text("%s", r.name);

                        // Local
                        float l_lvl = sel.local().chemical(r.id);
                        float l_p = sel.tick_chem_produced[idx];
                        float l_c = sel.tick_chem_consumed[idx];
                        ImGui::TableSetColumnIndex(1); ImGui::Text("%s", r.fmt(l_lvl));
                        ImGui::TableSetColumnIndex(2); ImGui::Text("%s", r.fmt(l_p));
                        ImGui::TableSetColumnIndex(3); ImGui::Text("%s", r.fmt(l_c));

                        // Vasculature
                        const TransportPool* vp = vascular_scope(sel, r.id);
                        if (vp) {
                            float v_lvl = vp->chemical(r.id);
                            ImGui::TableSetColumnIndex(4); ImGui::Text("%s", r.fmt(v_lvl));
                            ImGui::TableSetColumnIndex(5); ImGui::Text("\xE2\x80\x94");  // em dash — vasc-scope produced/consumed deferred
                            ImGui::TableSetColumnIndex(6); ImGui::Text("\xE2\x80\x94");
                            // Total
                            ImGui::TableSetColumnIndex(7); ImGui::Text("%s", r.fmt(l_lvl + v_lvl));
                            ImGui::TableSetColumnIndex(8); ImGui::Text("%s", r.fmt(l_p));
                            ImGui::TableSetColumnIndex(9); ImGui::Text("%s", r.fmt(l_c));
                        } else {
                            ImGui::TableSetColumnIndex(4); ImGui::Text("\xE2\x80\x94");
                            ImGui::TableSetColumnIndex(5); ImGui::Text("\xE2\x80\x94");
                            ImGui::TableSetColumnIndex(6); ImGui::Text("\xE2\x80\x94");
                            ImGui::TableSetColumnIndex(7); ImGui::Text("%s", r.fmt(l_lvl));
                            ImGui::TableSetColumnIndex(8); ImGui::Text("%s", r.fmt(l_p));
                            ImGui::TableSetColumnIndex(9); ImGui::Text("%s", r.fmt(l_c));
                        }
                    }
                    ImGui::EndTable();
                }
                ImGui::TextDisabled("L=local  V=vascular  T=total  +/- = this-tick produced/consumed");
                ImGui::Separator();

                // --- Sugar budget (last tick) ---
                ImGui::Text("Sugar last tick:");
                ImGui::Text("  Maintenance:  %s", fmt_mass(sel.tick_sugar_maintenance));
                if (sel.tick_sugar_activity >= 0.0f) {
                    ImGui::Text("  Produced:    +%s", fmt_mass(sel.tick_sugar_activity));
                } else {
                    ImGui::Text("  Growth:       %s", fmt_mass(sel.tick_sugar_activity));
                }
                if (sel.tick_sugar_transport < 0.0f) {
                    ImGui::Text("  Exported:     %s", fmt_mass(sel.tick_sugar_transport));
                } else if (sel.tick_sugar_transport > 0.0f) {
                    ImGui::Text("  Imported:    +%s", fmt_mass(sel.tick_sugar_transport));
                } else {
                    ImGui::Text("  Transport:    0");
                }

                ImGui::Separator();

                // Helper: node type as string
                auto node_type_str = [](const Node* n) -> const char* {
                    switch (n->type) {
                        case NodeType::STEM:        return "STEM";
                        case NodeType::ROOT:        return "ROOT";
                        case NodeType::LEAF:        return "LEAF";
                        case NodeType::APICAL:      return "APICAL";
                        case NodeType::ROOT_APICAL: return "ROOT_APICAL";
                    }
                    return "?";
                };

                // Parent navigation
                if (sel.parent) {
                    char label[64];
                    std::snprintf(label, sizeof(label), "Parent: %s", node_type_str(sel.parent));
                    if (ImGui::Button(label)) {
                        g_selected_node = sel.parent;
                    }
                    if (ImGui::IsItemHovered()) {
                        g_hovered_node = sel.parent;
                    }
                }

                // Child navigation
                ImGui::Spacing();
                if (sel.children.empty()) {
                    ImGui::TextDisabled("No children");
                } else {
                    for (int i = 0; i < static_cast<int>(sel.children.size()); ++i) {
                        const Node* child = sel.children[i];
                        char label[64];
                        std::snprintf(label, sizeof(label), "Child %d: %s", i, node_type_str(child));
                        if (ImGui::Button(label)) {
                            g_selected_node = child;
                        }
                        if (ImGui::IsItemHovered()) {
                            g_hovered_node = child;
                        }
                    }
                }
            }
            ImGui::End();
        }

        ImGui::Render();

        renderer.begin_frame();
        if (show_ground_shadow && renderer.light_system().is_initialized())
            renderer.draw_ground_shadow();
        else
            renderer.draw_ground();
        if (show_shadow_map) {
            renderer.draw_shadow_map(engine.shadow_map());
        }
        renderer.draw_plant(engine.get_plant(plant_id));
        if (g_hovered_node) renderer.draw_highlight(*g_hovered_node);
        if (g_selected_node && g_selected_node != g_hovered_node) renderer.draw_highlight(*g_selected_node);
        if (show_gpu_shadow_debug && renderer.light_system().is_initialized()) {
            renderer.light_system().draw_debug_slice(debug_slice_idx);
        }
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        renderer.end_frame();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
