// src/app_realtime.cpp
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "engine/engine.h"
#include "engine/node.h"
#include "renderer/renderer.h"

using namespace botany;

static Renderer* g_renderer = nullptr;

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main(int argc, char* argv[]) {
    std::string color_chemical;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--color") == 0 && i + 1 < argc) {
            color_chemical = argv[++i];
        }
    }

    Engine engine;
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
            }
            ImGui::SameLine();
            if (ImGui::Button("Node Type")) {
                renderer.set_color_mode(ChemicalAccessor{});
                renderer.set_color_by_type(true);
            }
            ImGui::SameLine();
            if (ImGui::Button("Auxin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.auxin; });
            }
            ImGui::SameLine();
            if (ImGui::Button("Cytokinin")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.cytokinin; });
            }
        }

        ImGui::End();
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
