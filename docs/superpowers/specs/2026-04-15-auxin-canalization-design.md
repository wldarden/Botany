# Auxin Canalization (Auxin Enhancement #2)

Connections that carry lots of auxin become bigger pipes for all chemicals.

## Context

Currently transport is stateless — every tick, `compute_transport_flow()` calculates desired flow purely from current chemical levels, radii, and diffusion params. A parent distributes to children proportionally by these computed weights. There is no memory of past flow.

In real plants, auxin flow is self-reinforcing. Cells that experience high auxin flux upregulate their transport proteins (PIN proteins), making that path easier to flow through. Over longer timescales, sustained auxin flow triggers actual vascular tissue development (xylem/phloem), permanently increasing the transport capacity of that route. This is how trunks, dominant branches, and stable architecture emerge — not from explicit rules, but from flow history.

This is the second of several planned auxin enhancements (see `docs/auxin_enhancement.txt`).

## The core idea

Each parent-child connection remembers how much auxin has flowed through it. Two kinds of memory:

1. **Transient bias** — builds up quickly when auxin flows, fades when flow stops. Like PIN protein polarization. Responds within hours/days.
2. **Structural bias** — builds up slowly when auxin flow is sustained, never fades. Like xylem/phloem tissue being built. Accumulates over days/weeks.

These biases don't increase total chemical flow — they **redistribute** it. When a parent splits chemicals among its children, biased connections get a bigger share of the same pie. This applies to ALL chemicals (sugar, auxin, cytokinin, GA, stress), not just auxin.

## Data model

### Per-connection state on the parent node

Two maps on `Node`, keyed by child pointer:

```cpp
std::unordered_map<Node*, float> auxin_flow_bias;      // transient — fast, decays
std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
```

Both initialize to 0 for new connections. The effective weight multiplier for any connection:

```
bias_multiplier = 1.0 + canalization_weight * (auxin_flow_bias[child] + structural_flow_bias[child])
```

At zero bias the connection behaves identically to today (multiplier = 1.0).

### Lifecycle management

- **`add_child(child)`** — biases start at 0, 0 (natural default from map).
- **`replace_child(old, new)`** — transfer both bias entries from old key to new key. This is critical for chain growth: when a meristem inserts an internode, the parent's accumulated biases for that branch path transfer to the new internode entry.
- **`die()`** — remove bias entries for dying children from the parent's maps. No special handling needed beyond normal cleanup.
- **Chain growth specifics:** When an apical spawns an internode between itself and its parent, the parent's biases transfer to the internode (via `replace_child`). The internode's biases for its new children (apical, axillary, leaf) all start at 0, 0 — new tissue, fresh connections.

## Auxin flux capture

During `transport_with_children()`, when processing the Auxin chemical specifically, we record the absolute amount of auxin that moved through each parent-child connection. This includes both directions — auxin flowing from child to parent (Phase 1) and from parent to child (Phase 2) both count as flux through the pipe.

The flux is stored temporarily (e.g., a local map or inline computation) and used immediately after transport to update the biases. It does not need to persist between ticks.

## Bias update mechanics

Runs after all chemical transport completes, before chemical decay. Uses the captured auxin flux per child.

### Transient bias update

Exponential chase — the bias tracks a target derived from current flux:

```
target = auxin_flux * transient_gain
auxin_flow_bias[child] += (target - auxin_flow_bias[child]) * transient_rate
```

- High flux: target is high, bias climbs toward it.
- Flux drops to zero: target is 0, bias decays back toward 0.
- `transient_rate` of 0.2 means ~87% of the gap is closed in 8 ticks (8 hours). Responsive enough to track day-to-day changes.

### Structural bias update

Slow ratchet — only grows, never shrinks:

```
if (auxin_flux > structural_threshold) {
    structural_flow_bias[child] += structural_growth_rate
    structural_flow_bias[child] = min(structural_flow_bias[child], structural_max)
}
// no decay — built vasculature stays permanently
```

- Only accumulates when flux exceeds a minimum threshold (no vascular development from a trickle).
- Grows at a fixed slow rate, capped at a maximum.
- At default rate (0.005/tick), reaching 1.0 takes ~200 ticks (~8 days of sustained flow).

## How biases affect transport

Applied to **all chemicals** in `transport_with_children()`. The raw `desired` flow from `compute_transport_flow()` is unchanged — concentration gradients, radii, and diffusion params still drive the physics. The bias only modifies the **proportional weight** used when splitting flow among siblings.

### Phase 1: Children giving to parent (limited parent headroom)

When multiple children want to give chemicals to the parent but the parent can't accept everything, each child's contribution is scaled proportionally by its effective weight:

```
effective_weight[child] = |desired[child]| * bias_multiplier[child]
child_share = total_headroom * (effective_weight[child] / sum_effective_weights)
```

High-bias children get priority to push their chemicals through.

### Phase 2: Parent giving to children

When the parent distributes chemicals to children, each child's share is proportional to its bias-adjusted weight:

```
effective_weight[child] = desired[child] * bias_multiplier[child]
child_share = to_distribute * (effective_weight[child] / sum_effective_weights)
```

High-bias children receive a larger share of the same total distribution.

### Key property: redistribution, not amplification

The total amount of chemicals moved is determined by the raw `desired` values and available supply/headroom — exactly as today. Biases only change **who gets what fraction** of that total. A connection with zero bias (multiplier = 1.0) works identically to the current system — it just competes against biased siblings for the same pool.

## Tick ordering

1. Transport all chemicals (biases applied as distribution weight multipliers)
2. During auxin transport specifically, record per-child flux amounts
3. After all transport completes, update both biases from recorded auxin flux
4. Decay chemicals

Biases from tick N affect transport in tick N. Flux from tick N's transport updates biases for tick N+1. Clean one-tick feedback loop.

## New genome parameters

| Field | Default | Evolvable range | Description |
|---|---|---|---|
| `transient_gain` | 2.0 | [0, 20] | Target bias per unit of auxin flux. Higher = connections respond more strongly to flow. |
| `transient_rate` | 0.2 | [0, 1] | How fast transient bias chases its target. 0.2 = ~87% in 8 hours. |
| `structural_threshold` | 0.05 | [0, 1] | Minimum auxin flux required to grow structural bias. Filters noise. |
| `structural_growth_rate` | 0.005 | [0, 0.1] | Structural bias increment per tick when flux exceeds threshold. ~8 days to reach 1.0. |
| `structural_max` | 2.0 | [0, 10] | Maximum structural bias. At max, this connection alone gets 1 + 2.0 = 3.0x weight. |
| `canalization_weight` | 1.0 | [0, 5] | Global scaling on the combined bias effect. 0 = canalization disabled. |

All 6 parameters go in a new "canalization" linkage group for evolution crossover.

## Files changed

### `src/engine/node/node.h`
- Add `auxin_flow_bias` and `structural_flow_bias` maps (both `std::unordered_map<Node*, float>`)

### `src/engine/node/node.cpp`
- `add_child()` — no code change needed (map default is 0)
- `replace_child()` — transfer bias entries from old child key to new child key
- `die()` — clean up bias entries for dying children from parent's maps
- `transport_with_children()` — apply bias as weight multiplier in Phase 1 and Phase 2 distribution; capture auxin flux per child during auxin iteration
- New method `update_canalization(const Genome& g)` — called after transport, updates both biases from captured auxin flux

### `src/engine/genome.h`
- Add 6 new fields to `Genome` struct
- Update `default_genome()` with defaults

### `src/evolution/genome_bridge.cpp`
- Register 6 new genes
- Add "canalization" linkage group
- Update `from_structured()` to read new genes

### `CLAUDE.md`
- Add canalization section describing the two-layer bias system
- Add new genome parameters to tuning section

### Tests
- Transient bias builds with auxin flux and decays without it
- Structural bias builds with sustained flux and does NOT decay
- Structural bias does NOT build below threshold flux
- Biased connection gets larger share of distributed sugar (redistribution test)
- Biases transfer on `replace_child` (chain growth inheritance)
- Zero `canalization_weight` disables all bias effects
- Full-size structural bias is capped at `structural_max`

## Expected emergent behavior

- **Dominant branches emerge:** A branch producing more auxin (more growing tips/leaves) builds transport preference, receives more sugar, grows more, produces more auxin — positive feedback loop.
- **Trunk formation:** The main stem accumulates high structural bias from sustained auxin flow through it. Even after seasonal leaf drop, the trunk remains the preferred transport path.
- **Graceful competition:** When a branch is shaded and its auxin drops, transient bias fades within days. But structural bias persists — the branch retains some resource access from its past success, declining gradually only because other branches' biases continue growing.
- **Crown shifting:** If light conditions change, new growth on the sunny side starts building fresh transport preference. Over time it overtakes the old dominant branch through the transient bias system, while structural bias ensures the old pathway doesn't collapse instantly.
