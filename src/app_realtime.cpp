// src/app_realtime.cpp
#include <GLFW/glfw3.h>
#include <iostream>
#include "engine/engine.h"
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

int main() {
    Engine engine;
    Genome g = default_genome();
    PlantID plant_id = engine.create_plant(g, glm::vec3(0.0f));

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

    int ticks_per_frame = 1;
    bool paused = false;

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            paused = !paused;
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
