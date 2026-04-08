// src/app_playback.cpp
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "serialization/serializer.h"
#include "renderer/renderer.h"

using namespace botany;

static bool mouse_pressed = false;
static double last_mouse_x = 0, last_mouse_y = 0;
static Renderer* g_renderer = nullptr;
static bool imgui_wants_mouse = false;

static void mouse_button_callback(GLFWwindow*, int button, int action, int) {
    if (imgui_wants_mouse) return;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mouse_pressed = (action == GLFW_PRESS);
    }
}

static void cursor_pos_callback(GLFWwindow*, double x, double y) {
    if (imgui_wants_mouse) return;
    if (mouse_pressed && g_renderer) {
        float dx = static_cast<float>(x - last_mouse_x);
        float dy = static_cast<float>(y - last_mouse_y);
        g_renderer->camera().rotate(dx, dy);
    }
    last_mouse_x = x;
    last_mouse_y = y;
}

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (imgui_wants_mouse) return;
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: botany_playback <recording.bin>" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];

    std::ifstream file(input_path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open " << input_path << std::endl;
        return 1;
    }

    RecordingHeader header = load_recording_header(file);
    std::vector<TickSnapshot> ticks;
    ticks.reserve(header.num_ticks);
    for (uint32_t i = 0; i < header.num_ticks; i++) {
        ticks.push_back(load_tick(file));
    }
    std::cout << "Loaded " << ticks.size() << " ticks" << std::endl;

    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;

    GLFWwindow* window = renderer.window();
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetScrollCallback(window, scroll_callback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    ImGui::StyleColorsDark();

    int current_tick = 0;
    bool playing = false;
    int speed = 1;
    float timer = 0.0f;
    float last_time = static_cast<float>(glfwGetTime());

    while (!glfwWindowShouldClose(window)) {
        float now = static_cast<float>(glfwGetTime());
        float dt = now - last_time;
        last_time = now;

        if (playing) {
            timer += dt * speed * 30.0f;
            while (timer >= 1.0f && current_tick < static_cast<int>(ticks.size()) - 1) {
                current_tick++;
                timer -= 1.0f;
            }
        }

        renderer.begin_frame();
        renderer.draw_ground();
        if (!ticks.empty()) {
            renderer.draw_snapshot(ticks[current_tick]);
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        imgui_wants_mouse = ImGui::GetIO().WantCaptureMouse;

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::Begin("Playback Controls");
        ImGui::Text("Tick: %d / %d", current_tick + 1, static_cast<int>(ticks.size()));
        ImGui::Text("Nodes: %d", static_cast<int>(ticks[current_tick].nodes.size()));
        ImGui::SliderInt("##tick", &current_tick, 0, static_cast<int>(ticks.size()) - 1);
        if (ImGui::Button(playing ? "Pause" : "Play")) {
            playing = !playing;
        }
        ImGui::SameLine();
        ImGui::SliderInt("Speed", &speed, 1, 10);
        if (ImGui::Button("Reset")) {
            current_tick = 0;
            playing = false;
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        renderer.end_frame();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
