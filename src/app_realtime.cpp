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
#include "engine/node/leaf_node.h"
#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/root_axillary.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/chemical/chemical_registry.h"
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

    get_f("auxin_production_rate", g.auxin_production_rate);
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
    get_f("thickening_rate", g.thickening_rate);
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
    }
    return "?";
}

static const char* node_type_label(NodeType t) {
    switch (t) {
        case NodeType::STEM:           return "Stem";
        case NodeType::ROOT:           return "Root";
        case NodeType::LEAF:           return "Leaf";
        case NodeType::SHOOT_APICAL:   return "SAM";
        case NodeType::SHOOT_AXILLARY: return "SAx";
        case NodeType::ROOT_APICAL:    return "RAM";
        case NodeType::ROOT_AXILLARY:  return "RAx";
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

    plant.for_each_node([&](const Node& node) {
        glm::vec3 test_pos = node.position;
        float effective_radius = node.radius;
        if (auto* leaf = node.as_leaf()) {
            if (leaf->leaf_size > 0.0f)
                effective_radius = leaf->leaf_size * 0.5f;
        }

        glm::vec3 to_node = test_pos - ray_origin;
        float t = glm::dot(to_node, ray_dir);
        if (t < 0.0f) return;

        glm::vec3 closest_on_ray = ray_origin + ray_dir * t;
        float dist = glm::length(test_pos - closest_on_ray);

        // Pick radius: generous for tiny nodes, proportional for larger ones
        float pick_radius = std::max(effective_radius * 2.0f, 0.15f);

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
    engine.debug_log().open("debug_log.csv");
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

    if (!color_chemical.empty()) {
        if (color_chemical == "type") {
            renderer.set_color_by_type(true);
        } else {
            ChemicalAccessor accessor;
            if (color_chemical == "auxin") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Auxin); };
            } else if (color_chemical == "cytokinin") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Cytokinin); };
            } else if (color_chemical == "sugar") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Sugar); };
            } else if (color_chemical == "gibberellin") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Gibberellin); };
            } else if (color_chemical == "ethylene") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Ethylene); };
            } else {
                std::cerr << "Unknown chemical: " << color_chemical
                          << " (available: auxin, cytokinin, sugar, gibberellin, ethylene, type)" << std::endl;
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

    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE, STRESS };
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
        int apical_count = 0, axillary_count = 0;
        int root_apical_count = 0, root_axillary_count = 0;
        const Genome& pg = engine.get_plant(plant_id).genome();
        engine.get_plant(plant_id).for_each_node([&](const Node& n) {
            total_sugar += n.chemical(ChemicalID::Sugar);
            if (n.chemical(ChemicalID::Sugar) > max_sugar) max_sugar = n.chemical(ChemicalID::Sugar);
            total_maintenance += n.maintenance_cost(engine.world_params());
            switch (n.type) {
                case NodeType::STEM: stem_count++; break;
                case NodeType::ROOT: root_count++; break;
                case NodeType::LEAF: leaf_count++; break;
                case NodeType::SHOOT_APICAL:   apical_count++; break;
                case NodeType::SHOOT_AXILLARY: axillary_count++; break;
                case NodeType::ROOT_APICAL:    root_apical_count++; break;
                case NodeType::ROOT_AXILLARY:  root_axillary_count++; break;
            }
        });

        // Per-tick sugar production (delta of lifetime accumulator)
        float current_produced = engine.get_plant(plant_id).total_sugar_produced();
        static float prev_produced = 0.0f;
        float last_tick_production = current_produced - prev_produced;
        prev_produced = current_produced;

        ImGui::Text("Stem: %d  Root: %d  Leaf: %d", stem_count, root_count, leaf_count);
        ImGui::Text("Meristems: SA:%d SX:%d RA:%d RX:%d",
                     apical_count, axillary_count, root_apical_count, root_axillary_count);
        ImGui::Text("Sugar: total=%s  max=%s", fmt_mass(total_sugar), fmt_mass(max_sugar));
        ImGui::Text("Production: %s  Maintenance: %s", fmt_mass_rate(last_tick_production), fmt_mass_rate(total_maintenance));
        ImGui::Separator();
        ImGui::SliderFloat("Light Level", &engine.world_params_mut().light_level, 0.0f, 2.0f);
        static bool show_shadow_map = false;
        ImGui::Checkbox("Show Shadow Map", &show_shadow_map);
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
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Auxin); });
                active_overlay = Overlay::AUXIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Cytokinin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Cytokinin); });
                active_overlay = Overlay::CYTOKININ;
            }
            ImGui::SameLine();
            if (ImGui::Button("Sugar")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Sugar); });
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
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Gibberellin); });
                active_overlay = Overlay::GIBBERELLIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Ethylene")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Ethylene); });
                active_overlay = Overlay::ETHYLENE;
            }
            ImGui::SameLine();
            if (ImGui::Button("Stress")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.stress; });
                active_overlay = Overlay::STRESS;
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
                float level = n.chemical(chem);
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
                    if (n.type == NodeType::SHOOT_APICAL) produced = eg.auxin_production_rate;
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
            ImGui::SetNextWindowSize(ImVec2(310, 0), ImGuiCond_Once);
            if (g_show_economy_modal) ImGui::SetNextWindowFocus();

            if (ImGui::Begin("Node Inspector", &g_show_node_panel, ImGuiWindowFlags_AlwaysAutoResize)) {
                const Node& sel = *g_selected_node;

                // Node type
                const char* type_str = "?";
                switch (sel.type) {
                    case NodeType::STEM:           type_str = "STEM"; break;
                    case NodeType::ROOT:           type_str = "ROOT"; break;
                    case NodeType::LEAF:           type_str = "LEAF"; break;
                    case NodeType::SHOOT_APICAL:   type_str = "SHOOT_APICAL"; break;
                    case NodeType::SHOOT_AXILLARY: type_str = "SHOOT_AXILLARY"; break;
                    case NodeType::ROOT_APICAL:    type_str = "ROOT_APICAL"; break;
                    case NodeType::ROOT_AXILLARY:  type_str = "ROOT_AXILLARY"; break;
                }
                ImGui::Text("Type: %s", type_str);

                // Meristem info
                if (auto* ax = sel.as_shoot_axillary()) {
                    ImGui::Text("Meristem: %s", ax->active ? "active" : "dormant");
                } else if (auto* ax = sel.as_root_axillary()) {
                    ImGui::Text("Meristem: %s", ax->active ? "active" : "dormant");
                } else if (sel.is_meristem()) {
                    ImGui::Text("Meristem: active");
                }

                ImGui::Text("ID: %u  Age: %u", sel.id, sel.age);
                ImGui::Text("Radius: %s", fmt_dist(sel.radius));
                if (auto* leaf = sel.as_leaf()) {
                    ImGui::Text("Leaf Size: %s", fmt_dist(leaf->leaf_size));
                    ImGui::Text("Light Exposure: %.1f%%", leaf->light_exposure * 100.0f);
                    // Sugar production: light * angle * world_light * area * rate
                    const Genome& lg = engine.get_plant(plant_id).genome();
                    float angle_eff = 1.0f;
                    float flen = glm::length(leaf->facing);
                    if (flen > 1e-4f) angle_eff = std::max(0.0f, (leaf->facing / flen).y);
                    float leaf_area = leaf->leaf_size * leaf->leaf_size;
                    float production = leaf->light_exposure * angle_eff
                        * engine.world_params().light_level * leaf_area
                        * engine.world_params().sugar_production_rate;
                    ImGui::Text("Sugar/tick: %s", fmt_mass_rate(production));
                    if (leaf->senescence_ticks > 0) {
                        ImGui::Text("Senescence: %u ticks", leaf->senescence_ticks);
                    }
                }
                ImGui::Text("Starvation: %u ticks", sel.starvation_ticks);
                ImGui::Text("Children: %d", static_cast<int>(sel.children.size()));

                ImGui::Separator();

                // Chemical levels table
                float parent_sugar = sel.parent ? sel.parent->chemical(ChemicalID::Sugar) : 0.0f;
                float parent_auxin = sel.parent ? sel.parent->chemical(ChemicalID::Auxin) : 0.0f;
                float parent_cytokinin = sel.parent ? sel.parent->chemical(ChemicalID::Cytokinin) : 0.0f;

                float child_sugar = 0.0f, child_auxin = 0.0f, child_cytokinin = 0.0f;
                if (!sel.children.empty()) {
                    for (const Node* c : sel.children) {
                        child_sugar += c->chemical(ChemicalID::Sugar);
                        child_auxin += c->chemical(ChemicalID::Auxin);
                        child_cytokinin += c->chemical(ChemicalID::Cytokinin);
                    }
                    float n = static_cast<float>(sel.children.size());
                    child_sugar /= n;
                    child_auxin /= n;
                    child_cytokinin /= n;
                }

                // Compute ratios (upstream/self and downstream/self)
                auto ratio_str = [](float other, float self) -> std::string {
                    if (self < 1e-6f) {
                        return other < 1e-6f ? "-" : "inf";
                    }
                    char buf[16];
                    std::snprintf(buf, sizeof(buf), "%.2fx", other / self);
                    return buf;
                };

                ImGui::Text("Chemical Levels:");
                if (ImGui::BeginTable("chemicals", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("Chemical", ImGuiTableColumnFlags_WidthFixed, 80);
                    ImGui::TableSetupColumn("Seed Side", ImGuiTableColumnFlags_WidthFixed, 65);
                    ImGui::TableSetupColumn("Self", ImGuiTableColumnFlags_WidthFixed, 75);
                    ImGui::TableSetupColumn("Tip Side", ImGuiTableColumnFlags_WidthFixed, 65);
                    ImGui::TableHeadersRow();

                    // Sugar
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Sugar");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_sugar, sel.chemical(ChemicalID::Sugar)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%s", fmt_mass(sel.chemical(ChemicalID::Sugar)));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_sugar, sel.chemical(ChemicalID::Sugar)).c_str());

                    // Auxin (dimensionless signaling units)
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Auxin");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_auxin, sel.chemical(ChemicalID::Auxin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%s", fmt_au(sel.chemical(ChemicalID::Auxin)));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_auxin, sel.chemical(ChemicalID::Auxin)).c_str());

                    // Cytokinin (dimensionless signaling units)
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Cyt");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_cytokinin, sel.chemical(ChemicalID::Cytokinin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%s", fmt_au(sel.chemical(ChemicalID::Cytokinin)));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_cytokinin, sel.chemical(ChemicalID::Cytokinin)).c_str());

                    // Gibberellin
                    float parent_ga = sel.parent ? sel.parent->chemical(ChemicalID::Gibberellin) : 0.0f;
                    float child_ga = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_ga += c->chemical(ChemicalID::Gibberellin);
                        child_ga /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("GA");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_ga, sel.chemical(ChemicalID::Gibberellin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%s", fmt_au(sel.chemical(ChemicalID::Gibberellin)));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_ga, sel.chemical(ChemicalID::Gibberellin)).c_str());

                    // Ethylene
                    float parent_eth = sel.parent ? sel.parent->chemical(ChemicalID::Ethylene) : 0.0f;
                    float child_eth = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_eth += c->chemical(ChemicalID::Ethylene);
                        child_eth /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Eth");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_eth, sel.chemical(ChemicalID::Ethylene)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%s", fmt_au(sel.chemical(ChemicalID::Ethylene)));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_eth, sel.chemical(ChemicalID::Ethylene)).c_str());

                    ImGui::EndTable();
                }
            }
            ImGui::End();
        }

        ImGui::Render();

        renderer.begin_frame();
        renderer.draw_ground();
        if (show_shadow_map) {
            renderer.draw_shadow_map(engine.shadow_map());
        }
        renderer.draw_plant(engine.get_plant(plant_id));
        if (g_hovered_node) renderer.draw_highlight(*g_hovered_node);
        if (g_selected_node && g_selected_node != g_hovered_node) renderer.draw_highlight(*g_selected_node);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        renderer.end_frame();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
