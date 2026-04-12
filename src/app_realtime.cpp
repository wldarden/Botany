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
#include "renderer/renderer.h"

using namespace botany;

static Renderer* g_renderer = nullptr;
static const Node* g_selected_node = nullptr;
static const Node* g_hovered_node = nullptr;
static bool g_show_node_panel = false;

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
            if (leaf->leaf_size > 0.0f && node.parent) {
                glm::vec3 outward = glm::normalize(node.position - node.parent->position);
                test_pos = node.position + outward * (leaf->leaf_size * 0.5f);
                effective_radius = leaf->leaf_size * 0.5f;
            }
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

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--color") == 0 && i + 1 < argc) {
            color_chemical = argv[++i];
        } else if (std::strcmp(argv[i], "--world") == 0 && i + 1 < argc) {
            world_path = argv[++i];
        }
    }

    Engine engine;
    if (!world_path.empty()) {
        engine.world_params_mut() = load_world_params(world_path);
    }
    Genome g = default_genome();
    PlantID plant_id = engine.create_plant(g, glm::vec3(0.0f));
    engine.debug_log().open("debug_log.csv");

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

    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE };
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
            total_maintenance += n.maintenance_cost(pg);
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
        ImGui::Text("Sugar: total=%.2fg  max=%.4fg", total_sugar, max_sugar);
        ImGui::Text("Production: %.4fg/tick  Maintenance: %.4fg/tick", last_tick_production, total_maintenance);
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
            g_selected_node = nullptr;
            g_show_node_panel = false;
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

        ImGui::End();

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
                ImGui::Text("Radius: %.4f", sel.radius);
                if (auto* leaf = sel.as_leaf()) {
                    ImGui::Text("Leaf Size: %.3f", leaf->leaf_size);
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
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Sugar (g)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_sugar, sel.chemical(ChemicalID::Sugar)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.chemical(ChemicalID::Sugar));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_sugar, sel.chemical(ChemicalID::Sugar)).c_str());

                    // Auxin (dimensionless signaling units)
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Auxin (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_auxin, sel.chemical(ChemicalID::Auxin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.chemical(ChemicalID::Auxin));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_auxin, sel.chemical(ChemicalID::Auxin)).c_str());

                    // Cytokinin (dimensionless signaling units)
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Cyt (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_cytokinin, sel.chemical(ChemicalID::Cytokinin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.chemical(ChemicalID::Cytokinin));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_cytokinin, sel.chemical(ChemicalID::Cytokinin)).c_str());

                    // Gibberellin
                    float parent_ga = sel.parent ? sel.parent->chemical(ChemicalID::Gibberellin) : 0.0f;
                    float child_ga = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_ga += c->chemical(ChemicalID::Gibberellin);
                        child_ga /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("GA (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_ga, sel.chemical(ChemicalID::Gibberellin)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.chemical(ChemicalID::Gibberellin));
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_ga, sel.chemical(ChemicalID::Gibberellin)).c_str());

                    // Ethylene
                    float parent_eth = sel.parent ? sel.parent->chemical(ChemicalID::Ethylene) : 0.0f;
                    float child_eth = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_eth += c->chemical(ChemicalID::Ethylene);
                        child_eth /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Eth (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_eth, sel.chemical(ChemicalID::Ethylene)).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.chemical(ChemicalID::Ethylene));
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
