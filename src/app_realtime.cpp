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
#include <cstdio>
#include <string>
#include "engine/engine.h"
#include "engine/node.h"
#include "engine/world_params.h"
#include "renderer/renderer.h"

using namespace botany;

static Renderer* g_renderer = nullptr;
static const Node* g_selected_node = nullptr;
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
    float closest_dist = FLT_MAX;
    float max_pick_dist = 0.5f; // max distance from ray to count as a hit

    plant.for_each_node([&](const Node& node) {
        // Point-to-ray distance
        glm::vec3 to_node = node.position - ray_origin;
        float t = glm::dot(to_node, ray_dir);
        if (t < 0.0f) return; // behind camera

        glm::vec3 closest_on_ray = ray_origin + ray_dir * t;
        float dist = glm::length(node.position - closest_on_ray);

        // Use a pick radius that scales with node radius (but has a minimum)
        float pick_radius = std::max(node.radius * 3.0f, max_pick_dist);

        if (dist < pick_radius && t < closest_dist) {
            closest_dist = t;
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
                accessor = [](const Node& n) { return n.auxin; };
            } else if (color_chemical == "cytokinin") {
                accessor = [](const Node& n) { return n.cytokinin; };
            } else if (color_chemical == "sugar") {
                accessor = [](const Node& n) { return n.sugar; };
            } else {
                std::cerr << "Unknown chemical: " << color_chemical
                          << " (available: auxin, cytokinin, sugar, type)" << std::endl;
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

    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR };
    Overlay active_overlay = Overlay::NONE;
    bool playing = false;
    int steps_remaining = 0;
    bool space_was_pressed = false;
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

        // Mouse picking (left click when ImGui doesn't want the mouse)
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            if (!ImGui::GetIO().WantCaptureMouse) {
                double mx, my;
                glfwGetCursorPos(window, &mx, &my);
                int w, h;
                glfwGetFramebufferSize(window, &w, &h);
                const Node* picked = pick_node(engine.get_plant(plant_id),
                                                renderer.camera(),
                                                static_cast<int>(mx), static_cast<int>(my), w, h);
                if (picked) {
                    g_selected_node = picked;
                    g_show_node_panel = true;
                }
            }
        }

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::SetNextWindowSizeConstraints(ImVec2(0, 0), ImVec2(FLT_MAX, FLT_MAX));
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        ImGui::Text("Tick: %d  Nodes: %d", total_ticks,
                     static_cast<int>(engine.get_plant(plant_id).node_count()));

        // Sugar stats
        float total_sugar = 0.0f;
        float max_sugar = 0.0f;
        int leaf_count = 0;
        engine.get_plant(plant_id).for_each_node([&](const Node& n) {
            total_sugar += n.sugar;
            if (n.sugar > max_sugar) max_sugar = n.sugar;
            if (n.type == NodeType::LEAF) leaf_count++;
        });
        ImGui::Text("Leaves: %d", leaf_count);
        ImGui::Text("Sugar: total=%.1f max=%.3f", total_sugar, max_sugar);
        ImGui::Separator();
        ImGui::SliderFloat("Light", &engine.world_params_mut().light_level, 0.0f, 2.0f);
        ImGui::SliderInt("Diffusion Iters",
            &engine.world_params_mut().sugar_diffusion_iterations, 1, 20);

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
                renderer.set_color_mode([](const Node& n) { return n.auxin; });
                active_overlay = Overlay::AUXIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Cytokinin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.cytokinin; });
                active_overlay = Overlay::CYTOKININ;
            }
            ImGui::SameLine();
            if (ImGui::Button("Sugar")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.sugar; });
                active_overlay = Overlay::SUGAR;
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
            glfwGetFramebufferSize(window, &w, &h);
            ImGui::SetNextWindowPos(ImVec2(static_cast<float>(w) - 320, 10), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(310, 0), ImGuiCond_Always);

            if (ImGui::Begin("Node Inspector", &g_show_node_panel, ImGuiWindowFlags_AlwaysAutoResize)) {
                const Node& sel = *g_selected_node;

                // Node type
                const char* type_str = "?";
                switch (sel.type) {
                    case NodeType::STEM: type_str = "STEM"; break;
                    case NodeType::ROOT: type_str = "ROOT"; break;
                    case NodeType::LEAF: type_str = "LEAF"; break;
                }
                ImGui::Text("Type: %s", type_str);

                // Meristem info
                if (sel.meristem) {
                    const char* mer_str = "?";
                    switch (sel.meristem->type()) {
                        case MeristemType::APICAL:        mer_str = "Shoot Apical"; break;
                        case MeristemType::AXILLARY:      mer_str = "Shoot Axillary"; break;
                        case MeristemType::ROOT_APICAL:   mer_str = "Root Apical"; break;
                        case MeristemType::ROOT_AXILLARY: mer_str = "Root Axillary"; break;
                    }
                    ImGui::Text("Meristem: %s (%s)", mer_str, sel.meristem->active ? "active" : "dormant");
                }

                ImGui::Text("ID: %u  Age: %u", sel.id, sel.age);
                ImGui::Text("Radius: %.4f", sel.radius);
                if (sel.type == NodeType::LEAF) {
                    ImGui::Text("Leaf Size: %.3f", sel.leaf_size);
                }
                ImGui::Text("Starvation: %u ticks", sel.starvation_ticks);
                ImGui::Text("Children: %d", static_cast<int>(sel.children.size()));

                ImGui::Separator();

                // Chemical levels table
                float parent_sugar = sel.parent ? sel.parent->sugar : 0.0f;
                float parent_auxin = sel.parent ? sel.parent->auxin : 0.0f;
                float parent_cytokinin = sel.parent ? sel.parent->cytokinin : 0.0f;

                float child_sugar = 0.0f, child_auxin = 0.0f, child_cytokinin = 0.0f;
                if (!sel.children.empty()) {
                    for (const Node* c : sel.children) {
                        child_sugar += c->sugar;
                        child_auxin += c->auxin;
                        child_cytokinin += c->cytokinin;
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
                    ImGui::TableSetupColumn("Upstream", ImGuiTableColumnFlags_WidthFixed, 65);
                    ImGui::TableSetupColumn("Self", ImGuiTableColumnFlags_WidthFixed, 75);
                    ImGui::TableSetupColumn("Downstream", ImGuiTableColumnFlags_WidthFixed, 65);
                    ImGui::TableHeadersRow();

                    // Sugar
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Sugar");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_sugar, sel.sugar).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.3f", sel.sugar);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_sugar, sel.sugar).c_str());

                    // Auxin
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Auxin");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_auxin, sel.auxin).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.auxin);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_auxin, sel.auxin).c_str());

                    // Cytokinin
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Cytokinin");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_cytokinin, sel.cytokinin).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.cytokinin);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_cytokinin, sel.cytokinin).c_str());

                    ImGui::EndTable();
                }
            }
            ImGui::End();
        }

        ImGui::Render();

        renderer.begin_frame();
        renderer.draw_ground();
        renderer.draw_plant(engine.get_plant(plant_id));
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        renderer.end_frame();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
