#pragma once

#include <atomic>
#include <cstdint>
#include <random>
#include <vector>
#include "evolution/fitness.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include <evolve/structured_genome.h>

namespace botany {

struct EvolutionConfig {
    uint32_t population_size = 100;
    uint32_t max_ticks = 17520;
    uint32_t num_threads = 4;
    uint32_t elitism_count = 2;
    uint32_t tournament_size = 5;
    float light_level_min = 0.5f;
    float light_level_max = 1.0f;
    uint32_t competitors = 1;      // plants per sim (1 = solo, >1 = shared sim with competition)
    float plant_spacing = 2.0f;    // dm between plants when competitors > 1
    float light_tilt_max = 0.52f;  // ~30 degrees max tilt
    float mutation_strength_pct = 0.03f;  // mutation stddev as fraction of gene's valid range (3%)
    FitnessWeights weights;
};

struct Individual {
    evolve::StructuredGenome genome;
    PlantStats stats;
    float fitness = 0.0f;
    uint32_t group_id = 0;
};

class EvolutionRunner {
public:
    explicit EvolutionRunner(const EvolutionConfig& config, uint32_t seed = 42);

    void run_generation();
    void reset();

    uint32_t generation() const { return generation_; }
    float best_fitness() const { return best_fitness_; }
    bool fitness_improved() const { return fitness_improved_; }
    const PlantStats& best_stats() const { return best_stats_; }
    const evolve::StructuredGenome& best_genome() const { return best_genome_; }
    const std::vector<float>& fitness_history() const { return fitness_history_; }
    const EvolutionConfig& config() const { return config_; }
    EvolutionConfig& config_mut() { return config_; }

    Genome best_as_botany_genome() const;

    // Genomes of competitors that were in the same sim group as the best plant.
    const std::vector<evolve::StructuredGenome>& best_competitor_genomes() const { return best_competitor_genomes_; }

    // Progress: how many individuals evaluated so far this generation (thread-safe read)
    uint32_t evaluated_count() const { return evaluated_count_.load(std::memory_order_relaxed); }

private:
    void init_population();
    void evaluate_all();
    void score_all();
    void evolve_population();
    WorldParams randomize_world();
    const Individual& tournament_select();

    EvolutionConfig config_;
    std::mt19937 rng_;
    uint32_t generation_ = 0;

    std::vector<Individual> population_;
    evolve::StructuredGenome genome_template_;

    float best_fitness_ = 0.0f;
    bool fitness_improved_ = false;
    std::atomic<uint32_t> evaluated_count_{0};
    PlantStats best_stats_;
    evolve::StructuredGenome best_genome_;
    std::vector<float> fitness_history_;
    std::vector<evolve::StructuredGenome> best_competitor_genomes_;
};

} // namespace botany
