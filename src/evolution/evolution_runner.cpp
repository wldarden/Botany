#include "evolution/evolution_runner.h"
#include "evolution/genome_bridge.h"
#include <algorithm>
#include <cmath>
#include <thread>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

namespace botany {

EvolutionRunner::EvolutionRunner(const EvolutionConfig& config, uint32_t seed)
    : config_(config), rng_(seed)
{
    genome_template_ = build_genome_template(default_genome(), config_.mutation_strength_pct);
    init_population();
}

void EvolutionRunner::init_population() {
    population_.clear();
    population_.resize(config_.population_size);

    for (auto& ind : population_) {
        ind.genome = genome_template_;
        evolve::mutate(ind.genome, rng_);
    }
    // Keep one individual as the unperturbed default
    population_[0].genome = genome_template_;
}

void EvolutionRunner::reset() {
    generation_ = 0;
    best_fitness_ = 0.0f;
    best_stats_ = {};
    fitness_history_.clear();
    init_population();
}

WorldParams EvolutionRunner::randomize_world() {
    WorldParams world = default_world_params();

    std::uniform_real_distribution<float> level_dist(config_.light_level_min, config_.light_level_max);
    world.light_level = level_dist(rng_);

    std::uniform_real_distribution<float> angle_dist(0.0f, config_.light_tilt_max);
    std::uniform_real_distribution<float> azimuth_dist(0.0f, 6.2831853f);
    float tilt = angle_dist(rng_);
    float azimuth = azimuth_dist(rng_);
    float st = std::sin(tilt);
    world.light_direction = glm::normalize(glm::vec3(st * std::cos(azimuth), std::cos(tilt), st * std::sin(azimuth)));

    return world;
}

void EvolutionRunner::evaluate_all() {
    WorldParams world = randomize_world();
    uint32_t pop_size = static_cast<uint32_t>(population_.size());
    uint32_t num_threads = std::min(config_.num_threads, pop_size);

    if (num_threads <= 1) {
        for (auto& ind : population_) {
            Genome g = from_structured(ind.genome);
            ind.stats = evaluate_plant(g, world, config_.max_ticks);
        }
        return;
    }

    std::vector<std::thread> threads;
    uint32_t chunk = pop_size / num_threads;
    uint32_t remainder = pop_size % num_threads;
    uint32_t start = 0;

    for (uint32_t t = 0; t < num_threads; t++) {
        uint32_t end = start + chunk + (t < remainder ? 1 : 0);
        threads.emplace_back([this, &world, start, end]() {
            for (uint32_t i = start; i < end; i++) {
                Genome g = from_structured(population_[i].genome);
                population_[i].stats = evaluate_plant(g, world, config_.max_ticks);
            }
        });
        start = end;
    }

    for (auto& t : threads) t.join();
}

void EvolutionRunner::score_all() {
    PlantStats gen_max;
    for (const auto& ind : population_) {
        gen_max.survival_ticks = std::max(gen_max.survival_ticks, ind.stats.survival_ticks);
        gen_max.node_count = std::max(gen_max.node_count, ind.stats.node_count);
        gen_max.leaf_count = std::max(gen_max.leaf_count, ind.stats.leaf_count);
        gen_max.total_sugar_produced = std::max(gen_max.total_sugar_produced, ind.stats.total_sugar_produced);
        gen_max.height = std::max(gen_max.height, ind.stats.height);
        gen_max.crown_ratio = std::max(gen_max.crown_ratio, ind.stats.crown_ratio);
        gen_max.branch_depth = std::max(gen_max.branch_depth, ind.stats.branch_depth);
        gen_max.leaf_height_spread = std::max(gen_max.leaf_height_spread, ind.stats.leaf_height_spread);
    }

    for (auto& ind : population_) {
        ind.fitness = compute_fitness(ind.stats, gen_max, config_.weights);
    }

    std::sort(population_.begin(), population_.end(),
              [](const Individual& a, const Individual& b) { return a.fitness > b.fitness; });

    float prev_best = best_fitness_;
    best_fitness_ = population_[0].fitness;
    best_stats_ = population_[0].stats;
    best_genome_ = population_[0].genome;
    fitness_history_.push_back(best_fitness_);
    fitness_improved_ = (best_fitness_ > prev_best + 1e-6f);
}

const Individual& EvolutionRunner::tournament_select() {
    std::uniform_int_distribution<uint32_t> dist(0, static_cast<uint32_t>(population_.size()) - 1);
    uint32_t best_idx = dist(rng_);

    for (uint32_t i = 1; i < config_.tournament_size; i++) {
        uint32_t idx = dist(rng_);
        if (population_[idx].fitness > population_[best_idx].fitness) {
            best_idx = idx;
        }
    }
    return population_[best_idx];
}

void EvolutionRunner::evolve_population() {
    std::vector<Individual> next_gen;
    next_gen.reserve(config_.population_size);

    uint32_t elite = std::min(config_.elitism_count, static_cast<uint32_t>(population_.size()));
    for (uint32_t i = 0; i < elite; i++) {
        next_gen.push_back(population_[i]);
    }

    while (next_gen.size() < config_.population_size) {
        const Individual& parent_a = tournament_select();
        const Individual& parent_b = tournament_select();

        Individual child;
        child.genome = evolve::crossover(parent_a.genome, parent_b.genome, rng_);
        evolve::mutate(child.genome, rng_);
        next_gen.push_back(std::move(child));
    }

    population_ = std::move(next_gen);
}

void EvolutionRunner::run_generation() {
    evaluate_all();
    score_all();
    evolve_population();
    generation_++;
}

Genome EvolutionRunner::best_as_botany_genome() const {
    return from_structured(best_genome_);
}

} // namespace botany
