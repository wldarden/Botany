#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <glm/vec3.hpp>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include "engine/engine.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "evolution/evolution_runner.h"
#include "evolution/genome_bridge.h"
#include "renderer/renderer.h"

using namespace botany;

static Renderer* g_renderer = nullptr;

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) g_renderer->camera().zoom(static_cast<float>(yoffset));
}

static void save_genome(const Genome& g, const std::string& path) {
    std::ofstream out(path);
    if (!out) return;
    out << "auxin_production_rate=" << g.auxin_production_rate << "\n";
    out << "auxin_transport_rate=" << g.auxin_transport_rate << "\n";
    out << "auxin_directional_bias=" << g.auxin_directional_bias << "\n";
    out << "auxin_decay_rate=" << g.auxin_decay_rate << "\n";
    out << "auxin_threshold=" << g.auxin_threshold << "\n";
    out << "cytokinin_production_rate=" << g.cytokinin_production_rate << "\n";
    out << "cytokinin_transport_rate=" << g.cytokinin_transport_rate << "\n";
    out << "cytokinin_directional_bias=" << g.cytokinin_directional_bias << "\n";
    out << "cytokinin_decay_rate=" << g.cytokinin_decay_rate << "\n";
    out << "cytokinin_threshold=" << g.cytokinin_threshold << "\n";
    out << "growth_rate=" << g.growth_rate << "\n";
    out << "max_internode_length=" << g.max_internode_length << "\n";
    out << "min_internode_length=" << g.min_internode_length << "\n";
    out << "branch_angle=" << g.branch_angle << "\n";
    out << "thickening_rate=" << g.thickening_rate << "\n";
    out << "internode_elongation_rate=" << g.internode_elongation_rate << "\n";
    out << "internode_maturation_ticks=" << g.internode_maturation_ticks << "\n";
    out << "root_growth_rate=" << g.root_growth_rate << "\n";
    out << "root_max_internode_length=" << g.root_max_internode_length << "\n";
    out << "root_min_internode_length=" << g.root_min_internode_length << "\n";
    out << "root_branch_angle=" << g.root_branch_angle << "\n";
    out << "root_internode_elongation_rate=" << g.root_internode_elongation_rate << "\n";
    out << "root_internode_maturation_ticks=" << g.root_internode_maturation_ticks << "\n";
    out << "root_gravitropism_strength=" << g.root_gravitropism_strength << "\n";
    out << "root_gravitropism_depth=" << g.root_gravitropism_depth << "\n";
    out << "max_leaf_size=" << g.max_leaf_size << "\n";
    out << "leaf_growth_rate=" << g.leaf_growth_rate << "\n";
    out << "leaf_bud_size=" << g.leaf_bud_size << "\n";
    out << "initial_radius=" << g.initial_radius << "\n";
    out << "root_initial_radius=" << g.root_initial_radius << "\n";
    out << "tip_offset=" << g.tip_offset << "\n";
    out << "growth_noise=" << g.growth_noise << "\n";
    out << "leaf_phototropism_rate=" << g.leaf_phototropism_rate << "\n";
    out << "sugar_production_rate=" << g.sugar_production_rate << "\n";
    out << "sugar_transport_conductance=" << g.sugar_transport_conductance << "\n";
    out << "sugar_maintenance_leaf=" << g.sugar_maintenance_leaf << "\n";
    out << "sugar_maintenance_stem=" << g.sugar_maintenance_stem << "\n";
    out << "sugar_maintenance_root=" << g.sugar_maintenance_root << "\n";
    out << "sugar_maintenance_meristem=" << g.sugar_maintenance_meristem << "\n";
    out << "seed_sugar=" << g.seed_sugar << "\n";
    out << "sugar_storage_density_wood=" << g.sugar_storage_density_wood << "\n";
    out << "sugar_storage_density_leaf=" << g.sugar_storage_density_leaf << "\n";
    out << "sugar_cap_minimum=" << g.sugar_cap_minimum << "\n";
    out << "sugar_cap_meristem=" << g.sugar_cap_meristem << "\n";
    out << "sugar_activation_shoot=" << g.sugar_activation_shoot << "\n";
    out << "sugar_activation_root=" << g.sugar_activation_root << "\n";
    out << "ga_production_rate=" << g.ga_production_rate << "\n";
    out << "ga_leaf_age_max=" << g.ga_leaf_age_max << "\n";
    out << "ga_elongation_sensitivity=" << g.ga_elongation_sensitivity << "\n";
    out << "ga_length_sensitivity=" << g.ga_length_sensitivity << "\n";
    out << "ga_transport_rate=" << g.ga_transport_rate << "\n";
    out << "ga_directional_bias=" << g.ga_directional_bias << "\n";
    out << "ga_decay_rate=" << g.ga_decay_rate << "\n";
    out << "ethylene_starvation_rate=" << g.ethylene_starvation_rate << "\n";
    out << "ethylene_shade_rate=" << g.ethylene_shade_rate << "\n";
    out << "ethylene_shade_threshold=" << g.ethylene_shade_threshold << "\n";
    out << "ethylene_age_rate=" << g.ethylene_age_rate << "\n";
    out << "ethylene_age_onset=" << g.ethylene_age_onset << "\n";
    out << "ethylene_crowding_rate=" << g.ethylene_crowding_rate << "\n";
    out << "ethylene_crowding_radius=" << g.ethylene_crowding_radius << "\n";
    out << "ethylene_diffusion_radius=" << g.ethylene_diffusion_radius << "\n";
    out << "ethylene_abscission_threshold=" << g.ethylene_abscission_threshold << "\n";
    out << "ethylene_elongation_inhibition=" << g.ethylene_elongation_inhibition << "\n";
    out << "senescence_duration=" << g.senescence_duration << "\n";
}

int main() {
    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;
    GLFWwindow* window = renderer.window();
    glfwSetScrollCallback(window, scroll_callback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    ImGui::StyleColorsDark();

    EvolutionConfig evo_config;
    evo_config.population_size = 50;
    evo_config.max_ticks = 5000;
    evo_config.num_threads = std::max(1u, std::thread::hardware_concurrency());
    EvolutionRunner runner(evo_config);

    // Display engine: re-simulates the best plant for rendering
    Engine display_engine;
    PlantID display_plant = display_engine.create_plant(default_genome(), glm::vec3(0.0f));
    bool display_needs_update = false;

    bool running = false;
    std::atomic<bool> gen_in_progress{false};
    std::thread evo_thread;

    // Config state (editable while paused)
    int pop_size = static_cast<int>(evo_config.population_size);
    int max_ticks = static_cast<int>(evo_config.max_ticks);
    int num_threads = static_cast<int>(evo_config.num_threads);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);

        // Camera controls
        const float rotate_speed = 6.0f;
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) renderer.camera().rotate(-rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) renderer.camera().rotate(rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) renderer.camera().rotate(0.0f, rotate_speed);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) renderer.camera().rotate(0.0f, -rotate_speed);
        const float pan_speed = 0.05f;
        if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() + glm::vec3(0.0f, pan_speed, 0.0f));
        if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() - glm::vec3(0.0f, pan_speed, 0.0f));

        // Check if background generation finished
        if (!gen_in_progress && evo_thread.joinable()) {
            evo_thread.join();
            display_needs_update = true;

            // Autosave on improvement
            if (runner.fitness_improved()) {
                save_genome(runner.best_as_botany_genome(), "best_genome.txt");
                std::cout << "Gen " << runner.generation()
                          << ": new best fitness " << runner.best_fitness()
                          << " -> saved best_genome.txt" << std::endl;
            }
        }

        // Auto-launch next generation if running
        if (running && !gen_in_progress && !evo_thread.joinable()) {
            gen_in_progress = true;
            evo_thread = std::thread([&runner, &gen_in_progress]() {
                runner.run_generation();
                gen_in_progress = false;
            });
        }

        // Re-simulate best plant for display
        if (display_needs_update && runner.generation() > 0) {
            display_engine.reset();
            Genome best_g = runner.best_as_botany_genome();
            display_plant = display_engine.create_plant(best_g, glm::vec3(0.0f));
            uint32_t replay_ticks = runner.best_stats().survival_ticks;
            for (uint32_t i = 0; i < replay_ticks; i++) {
                display_engine.tick();
            }
            display_needs_update = false;
        }

        // ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::Begin("Evolution", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

        // --- Config section ---
        if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen)) {
            bool paused = !running;
            if (paused) {
                ImGui::SliderInt("Population", &pop_size, 10, 500);
                ImGui::SliderInt("Max Ticks", &max_ticks, 100, 20000);
                ImGui::SliderInt("Threads", &num_threads, 1, 16);
            } else {
                ImGui::Text("Population: %d  Max Ticks: %d  Threads: %d", pop_size, max_ticks, num_threads);
            }

            ImGui::SeparatorText("Fitness Weights");
            FitnessWeights& w = runner.config_mut().weights;
            ImGui::SliderFloat("Survival", &w.survival, 0.0f, 5.0f);
            ImGui::SliderFloat("Biomass", &w.biomass, 0.0f, 5.0f);
            ImGui::SliderFloat("Sugar", &w.sugar, 0.0f, 5.0f);
            ImGui::SliderFloat("Leaves", &w.leaves, 0.0f, 5.0f);
            ImGui::SliderFloat("Height", &w.height, 0.0f, 5.0f);
            ImGui::SliderFloat("Crown Ratio", &w.crown_ratio, 0.0f, 5.0f);
            ImGui::SliderFloat("Branch Depth", &w.branch_depth, 0.0f, 5.0f);
            ImGui::SliderFloat("Leaf Spread", &w.leaf_spread, 0.0f, 5.0f);
        }

        // --- Controls ---
        ImGui::SeparatorText("Controls");
        if (!running) {
            if (ImGui::Button("Start")) {
                runner.config_mut().population_size = static_cast<uint32_t>(pop_size);
                runner.config_mut().max_ticks = static_cast<uint32_t>(max_ticks);
                runner.config_mut().num_threads = static_cast<uint32_t>(num_threads);
                if (runner.generation() == 0) runner.reset();
                running = true;
            }
        } else {
            if (ImGui::Button("Pause")) {
                running = false;
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            running = false;
            if (evo_thread.joinable()) evo_thread.join();
            gen_in_progress = false;
            runner.config_mut().population_size = static_cast<uint32_t>(pop_size);
            runner.config_mut().max_ticks = static_cast<uint32_t>(max_ticks);
            runner.config_mut().num_threads = static_cast<uint32_t>(num_threads);
            runner.reset();
            display_engine.reset();
            display_plant = display_engine.create_plant(default_genome(), glm::vec3(0.0f));
        }

        // --- Stats ---
        if (runner.generation() > 0) {
            ImGui::SeparatorText("Stats");
            ImGui::Text("Generation: %u", runner.generation());
            ImGui::Text("Best Fitness: %.3f", runner.best_fitness());
            if (gen_in_progress) ImGui::Text("(evaluating...)");

            const PlantStats& s = runner.best_stats();
            if (ImGui::BeginTable("stats", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("Metric", ImGuiTableColumnFlags_WidthFixed, 120);
                ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 100);
                ImGui::TableHeadersRow();

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Survival");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%u ticks", s.survival_ticks);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Nodes");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%u", s.node_count);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Leaves");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%u", s.leaf_count);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Sugar");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.1f g", s.total_sugar_produced);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Height");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.2f dm", s.height);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Crown Ratio");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.2f", s.crown_ratio);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Branch Depth");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%u", s.branch_depth);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Leaf Spread");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.2f dm", s.leaf_height_spread);

                ImGui::EndTable();
            }

            // Fitness history plot
            const auto& hist = runner.fitness_history();
            if (hist.size() > 1) {
                ImGui::PlotLines("Fitness", hist.data(), static_cast<int>(hist.size()),
                                 0, nullptr, FLT_MAX, FLT_MAX, ImVec2(0, 80));
            }

            // Export Best button
            if (ImGui::Button("Export Best")) {
                save_genome(runner.best_as_botany_genome(), "best_genome.txt");
                std::cout << "Exported best genome to best_genome.txt" << std::endl;
            }
        }

        ImGui::End();
        ImGui::Render();

        // Render
        renderer.begin_frame();
        renderer.draw_ground();
        renderer.draw_plant(display_engine.get_plant(display_plant));
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        renderer.end_frame();
    }

    // Cleanup
    running = false;
    if (evo_thread.joinable()) evo_thread.join();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
