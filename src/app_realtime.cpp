// src/app_realtime.cpp
#include <GLFW/glfw3.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "engine/engine.h"
#include "engine/node.h"
#include "renderer/renderer.h"

using namespace botany;

static bool mouse_pressed = false;
static double last_mouse_x = 0, last_mouse_y = 0;
static Renderer* g_renderer = nullptr;

static void mouse_button_callback(GLFWwindow*, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mouse_pressed = (action == GLFW_PRESS);
    }
}

static void cursor_pos_callback(GLFWwindow* window, double x, double y) {
    if (mouse_pressed && g_renderer) {
        float dx = static_cast<float>(x - last_mouse_x);
        float dy = static_cast<float>(y - last_mouse_y);
        g_renderer->camera().rotate(dx, dy);
    }
    last_mouse_x = x;
    last_mouse_y = y;
}

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main(int argc, char* argv[]) {
    int stop_at = -1;
    std::string color_chemical;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--color") == 0 && i + 1 < argc) {
            color_chemical = argv[++i];
        } else {
            stop_at = std::atoi(argv[i]);
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
            std::cout << "Color mode: type (green=shoot, orange=root)" << std::endl;
        } else {
            ChemicalAccessor accessor;
            if (color_chemical == "auxin") {
                accessor = [](const Node& n) { return n.auxin; };
            } else if (color_chemical == "cytokinin") {
                accessor = [](const Node& n) { return n.cytokinin; };
            } else {
                std::cerr << "Unknown chemical: " << color_chemical
                          << " (available: auxin, cytokinin, type)" << std::endl;
                renderer.shutdown();
                return 1;
            }
            renderer.set_color_mode(std::move(accessor));
            std::cout << "Color mode: " << color_chemical << std::endl;
        }
    }

    GLFWwindow* window = renderer.window();
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetScrollCallback(window, scroll_callback);

    int ticks_per_frame = 1;
    bool paused = false;
    bool space_was_pressed = false;

    // If stop_at given, run those ticks immediately before rendering
    if (stop_at > 0) {
        for (int i = 0; i < stop_at; i++) {
            engine.tick();
        }
        paused = true;
        std::cout << "Ran " << stop_at << " ticks (" << engine.get_plant(plant_id).node_count() << " nodes). Paused — press Space to continue." << std::endl;
    }

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);

        // Debounce space — toggle on press, not hold
        bool space_down = glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS;
        if (space_down && !space_was_pressed) {
            paused = !paused;
        }
        space_was_pressed = space_down;

        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
            ticks_per_frame = std::min(ticks_per_frame + 1, 100);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
            ticks_per_frame = std::max(ticks_per_frame - 1, 1);

        if (!paused) {
            for (int i = 0; i < ticks_per_frame; i++) {
                engine.tick();
            }
        }

        renderer.begin_frame();
        renderer.draw_ground();
        renderer.draw_plant(engine.get_plant(plant_id));
        renderer.end_frame();
    }

    renderer.shutdown();
    return 0;
}
