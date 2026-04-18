# Münch Pressure Flow: Design Document

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current global-allocation phloem model with a local, pressure-driven Münch flow system where sugar transport emerges from osmotic pressure differentials and physical phloem velocity, without a central scheduler, reserved-sugar bookkeeping, or explicit source/sink classification.

**Architecture:** Nodes do all local work (produce, consume, grow) first during the DFS tick, with no transport. After the DFS, a phloem resolve pass reads the resulting post-consumption sugar landscape, computes osmotic pressure at every vascular node, and runs a distance-limited flow from every high-pressure node outward — with unloading at each traversed sink. Xylem runs as a separate resolve pass after phloem.

**Scope:** Phloem only. The xylem pass (water + cytokinin) keeps its Phase 1/Phase 2 structure.

---

## 1. The Problem with the Current Model

The current `vascular_transport` is a two-phase global scheduler that runs *before* nodes have done their local work for the tick.

### 1.1 The tick order is backwards

The current order in `tick_tree`:

```
pre_transport_growth()     // nodes pre-declare how much sugar they'll need for growth
vascular_transport()       // global allocation — deducts from sources, delivers to sinks
pin_transport()
tick_recursive()           // DFS: nodes actually grow and produce
```

Sugar is allocated before it is produced. To work around this, every node must pre-declare a `sugar_reserved_for_growth` estimate and the vascular pass must exclude that reservation from the leaf's available supply. This is a bookkeeping hack — an approximation of what nodes will actually need, computed before the tick's physics have run.

The correct biological order is: leaves produce first, meristems consume first, then the phloem equilibrates the resulting pressure landscape. Loading is a consequence of high local sugar concentration; unloading is a consequence of low local concentration after consumption. Neither requires pre-declaration.

### 1.2 Sugar teleports across the plant in one tick

Both a canopy leaf and a deep root meristem contribute to the seed's global supply/demand totals in Phase 1, then draw from the same pool in Phase 2 — in a single tick. The seed mediates all exchange globally. Real phloem moves at 0.3–1.5 m/hr. A large plant has meters of phloem between source and sink. The transport model should reflect that distance costs time.

### 1.3 Meristems need artificial demand caps

`meristem_sink_fraction` caps how much sugar an active meristem can demand from phloem per tick. This exists because in the allocation model, a meristem can demand its full deficit and the planner tries to fill it, outcompeting leaves and stem conduits. Under Münch flow, a meristem's pull is bounded by the concentration gradient from its parent conduit and the pipe conductance — no explicit cap needed.

### 1.4 The leaf reserve floor is a static hack

`phloem_reserve_fraction` prevents leaves from contributing all their sugar to Phase 1 supply. Under the new model, leaves produce sugar during the DFS tick, consume what they need, and whatever remains is their actual available concentration for phloem. No separate floor is needed — the leaf naturally retains whatever sugar it didn't consume.

### 1.5 Flow direction is fixed by classification

Every node is classified as source, sink, or conduit at the start of each tick. Stems cannot be sources. In Münch flow, a stem node with high sugar (from mobilized starch, for example) simply has high osmotic pressure and exports naturally. Direction follows pressure, not classification.

### 1.6 Water doesn't affect phloem flow

Phloem requires water to build osmotic pressure in sieve tubes. The current model moves sugar regardless of water status. Münch coupling is free: `phloem_pressure = (sugar/volume) × osmotic_coeff × water_fraction`. No extra code — drought automatically suppresses phloem flow.

---

## 2. How Münch Pressure Flow Works (Biology)

### 2.1 The physical mechanism

Phloem sieve tubes run from leaves to sinks as continuous pipes. At source leaves, companion cells pump sucrose from mesophyll into sieve tubes — raising sucrose concentration there, drawing water in from xylem by osmosis, and building turgor pressure (typically 1–3 MPa at a productive leaf vs. 0.1–0.5 MPa at a growing sink).

At sinks, sucrose is consumed. Sieve tube concentration drops. Water exits back to xylem. Turgor drops. The pressure gradient between source and sink drives bulk flow through the sieve tube. Sucrose dissolved in the water moves with the stream at 0.3–1.5 m/hr in real trees.

Key properties:

- **Local:** Each sieve element only "knows" its pressure and the pressure of its immediate neighbors. No scheduler.
- **Bidirectional:** Flow goes from high to low pressure, regardless of plant orientation.
- **Distance-sensitive:** Pressure drops along the path. A nearby sink gets a steeper gradient; a distant sink gets a diluted stream that has already unloaded at intermediate sinks.
- **Water-coupled:** A dehydrated source cannot build osmotic pressure, so it cannot export even if it has sugar.

### 2.2 Loading

Loading is **implicit** — there is no separate loading step. A leaf with high post-photosynthesis sugar concentration simply has high osmotic pressure in its phloem sieve elements. That elevated pressure IS the loading signal. The companion cell pump is approximated by the osmotic equation directly; there is no additional loading boost, no loading cost, and no loading rate parameter. The leaf's sugar concentration after consuming for its own growth is what drives export.

### 2.3 Unloading

Unloading is **passive** and happens during phloem resolution at each node the pressure wave passes through. The driving force is the concentration difference between the phloem stream and the node's local cytoplasmic sugar. Fast-growing nodes (meristems, expanding leaves) have depleted their sugar during the DFS tick → steep gradient → fast unloading. Mature stems with stable sugar are near equilibrium with the stream → slow unloading → act as conduits.

Unloading permeabilities differ by tissue type. Meristems have high permeability (active sieve tube unloaders); stem conduits have very low permeability (sealed phloem passing through with minimal leakage).

---

## 3. Finalized Tick Order

```
1. pin_transport()
   — PIN-mediated polar auxin transport (global pass)

2. tick_recursive()   [DFS: LOCAL WORK ONLY]
   — produce: photosynthesis adds sugar to self
   — consume: maintenance, growth, leaf expansion all deduct from local sugar
   — NO transport during this step (sugar stays local)
   — auxin/GA/stress local diffusion still happens here (transport_with_children
     skips Sugar and Water; those are handled by resolve passes)

3. phloem_resolve()   [three sub-steps]

   3a. Leaf loading pass  (pre-BFS, endpoint sources)
       — For each leaf: load = max(0, leaf_conc − parent_ring_conc) × perm_leaf × pipe_cap
       — Capped at parent ring room (ring_vol × max_conc − parent.sugar)
       — leaf.sugar reduced; parent stem.sugar increased
       — Leaves are NOT in the BFS graph

   3b. BFS on stem/root/seed network
       — pressure = (sugar/phloem_vol) × osmotic_coeff × water_fraction
       — phloem_vol: ring for stems/roots; TOTAL cylinder volume for seed (storage organ)
       — distance-limited BFS from high-pressure conduits outward
       — source capacity guard: flow_vol ≤ available[source] + running_delta[source]
       — destination cap guard: flow_vol ≤ room in destination ring
       — transit dead-end credit: when stream cannot continue, credit transit back to conduit
       — multi-source: shared running delta so later BFS walks see earlier fills
       — unloading at each traversed conduit via permeability formula

   3c. Meristem unloading pass  (post-BFS, terminal sinks)
       — For each meristem: unload = max(0, parent_conc − mer_conc) × perm_mer × pipe_cap(bottleneck)
       — Capped at min(parent.sugar, meristem_room)  ← meristem room cap essential
       — parent.sugar reduced; meristem.sugar increased
       — Meristems are NOT in the BFS graph

4. xylem_resolve()
   — Phase 1/Phase 2 for Water and Cytokinin (unchanged)
   — runs after phloem so it sees the final sugar state for the tick
```

**Why phloem before xylem:** phloem pressure uses last-tick's water values (the water level after the previous xylem resolve). This one-tick lag is correct: roots absorb water, xylem transports it over the course of an hour, and it becomes available for phloem sieve tube osmosis in the subsequent hour. The lag is not an approximation — it is biologically accurate discrete-time physics.

**Why "LOCAL WORK ONLY" in the DFS:** Grow-before-transport is now inherent. Each node consumes whatever sugar it needs during the DFS tick. When phloem_resolve runs, it sees the real post-consumption concentration — not a pre-declared reservation. No `sugar_reserved_for_growth` field. No `pre_transport_growth()` pass. No `phloem_reserve_fraction` floor. The local consumption step produces the correct source/sink landscape automatically.

---

## 4. Proposed Algorithm

### 4.1 Pressure computation

After the DFS tick, compute phloem osmotic pressure at every vascular node using **bulk sugar concentration** (g/dm³) — sugar divided by total node volume. Using total volume rather than ring-only volume prevents the 85× concentration asymmetry between leaves and stems that inverts flow direction (see trial calculation and worked examples below):

```
for each vascular node n:
    node_vol    = node_volume(n)           // dm³ — total node volume (full cylinder for stems/roots)
    sugar_conc  = n.sugar / node_vol       // g/dm³ — bulk sugar concentration
    water_frac  = clamp(n.water / water_cap(n, g), 0.0, 1.0)
    phloem_pressure[n] = sugar_conc × g.phloem_osmotic_coefficient × water_frac
```

**Two volume concepts — used for different purposes:**

- **`node_volume(n)`** — total node volume, used for **concentration** (pressure) and **sugar cap**. Full cylinder for stems/roots (including heartwood — the dilution is intentional; it prevents ring-volume concentrations from making stems appear as higher-pressure than leaves at realistic sugar levels).
- **`phloem_ring_area(r, t)`** — phloem ring cross-section, used exclusively for **pipe capacity** in the flow formula (`flow_vol = velocity × ring_area`). Heartwood excluded. Scales linearly with r.

| Node type | `node_volume` formula | `phloem_ring_area` | Role in phloem_resolve |
|-----------|---------|------------------------|-------|
| STEM, ROOT | `π × r² × \|offset\|` (full cylinder) | `π × (2r×t - t²)` | **BFS conduit** — pressure from full cylinder volume; pipe flow uses ring area |
| SEED | `π × r² × \|offset\|` | `π × (2r×t - t²)` | **BFS source** — total cylinder (storage organ; no heartwood to exclude) |
| LEAF | `leaf_size² × world.leaf_thickness` | N/A — uses parent stem's ring for pipe_cap | **Endpoint source** — not in BFS graph; loads into parent stem pre-BFS |
| APICAL, ROOT_APICAL | `(4/3)π × radius³` | N/A — uses parent stem's ring for pipe_cap | **Terminal sinks** — not in BFS graph; unload from parent stem post-BFS |

**Ring area helper:**  `phloem_ring_area(r, t) = π × (2r×t - t²)`. For r >> t, this ≈ 2π×r×t.

**Worked examples** — showing why full node volume fixes flow direction:
- Young stem (r=0.015 dm, L=0.3 dm): ring area ≈ 0.000176 dm². node_vol = π × 0.015² × 0.3 ≈ 2.12×10⁻⁴ dm³. With 0.001g sugar: **bulk conc ≈ 4.72 g/dm³** (vs. 19 g/dm³ with ring vol — 4× lower).
- Mature stem (r=0.1 dm, L=1.0 dm): ring area ≈ 0.00124 dm². node_vol = π × 0.1² × 1.0 ≈ 0.0314 dm³. With 0.1g sugar: **bulk conc ≈ 3.18 g/dm³** (vs. 80 g/dm³ with ring vol — 25× lower).
- Leaf (leaf_size=0.5 dm): node_vol = 0.25 × 0.003 = 7.5×10⁻⁴ dm³. With 0.03g post-photosynthesis sugar: **conc ≈ 40 g/dm³**. Leaf pressure (~40) is ~9× higher than adjacent young stem (~4.7) — flow goes leaf→stem. ✓

`node_volume` and `phloem_ring_area` are derived from live geometry fields — no stored parameters, no per-tissue density lookups. If volume is zero or degenerate (new node with zero offset), return 0 pressure.

- Leaves with high post-photosynthesis sugar → high bulk g/dm³ concentration → high pressure → they are sources
- Meristems after consuming for growth → low g/dm³ → low pressure → they are sinks
- Stems equilibrate between their upstream source and downstream sink
- Drought (low water_frac) suppresses pressure even with abundant sugar — phloem stalls

No loading boost. No classification. No reserve floor. No per-tissue cap lookup. Pressure is computed purely from node geometry and the actual post-consumption chemical state of the node.

### 4.2 Distance-limited BFS from sources

Each high-pressure node is a source. From each source, walk outward toward lower-pressure neighbors, spending a time budget of 1.0 tick. Each edge consumes time based on its physical length and the local phloem velocity.

**Phloem velocity at an edge** (velocity-capped model):

Pressure drives velocity, but sieve plate resistance and sap viscosity impose a physical upper limit:

```
pressure_diff = phloem_pressure[node] - phloem_pressure[neighbor]
velocity = pressure_diff × g.conductance_per_pressure    // dm/tick per (g/dm³) gradient
velocity = min(velocity, world.max_phloem_velocity)      // physical cap: ~10 dm/tick (= 1 m/hr)
```

`conductance_per_pressure` maps osmotic pressure gradient [g/dm³] to sap velocity [dm/tick]. `max_phloem_velocity` is the hard physical ceiling — prevents astronomical flow rates at high-gradient edges (e.g., a leaf at 40 g/dm³ adjacent to an almost-empty stem at 0.5 g/dm³). Wider pipes carry more sugar per tick at the same velocity via the ring area in the flow formula, not via higher velocity.

**Time cost per edge:**

```
edge_length = |child.position - parent.position|     // dm
time_cost   = edge_length / velocity                 // ticks — same velocity drives BFS traversal and flow volume
```

**BFS walk from one source node:**

```
queue: [(source_node, time_remaining=1.0, stream_conc=pressure[source]/g.phloem_osmotic_coefficient)]

while queue not empty:
    (node, budget, stream_conc) = dequeue

    for each neighbor of node with phloem_pressure[neighbor] < phloem_pressure[node]:
        edge = edge between node and neighbor
        time_cost = edge_length / speed(edge)

        pressure_diff = phloem_pressure[node] - phloem_pressure[neighbor]
        if time_cost >= budget:
            time_fraction = budget / time_cost
            // partial edge — flow arrives but budget exhausted
            apply_flow(edge, source=node, dest=neighbor,
                       pressure_diff=pressure_diff, stream_conc=stream_conc, fraction=time_fraction)
            // do NOT enqueue neighbor (budget used up)
        else:
            time_fraction = 1.0
            stream_conc_after = apply_flow(edge, source=node, dest=neighbor,
                                           pressure_diff=pressure_diff, stream_conc=stream_conc, fraction=1.0)
            enqueue(neighbor, budget - time_cost, stream_conc_after)
```

All deltas are accumulated in a per-node buffer. After all BFS walks complete, apply deltas atomically:

```
for each node:
    n.sugar = clamp(n.sugar + delta[n], 0.0, sugar_cap(n, g))
```

### 4.3 Flow formula: two limits, one transfer

Per-edge sugar transfer is governed by two independent limits. **The lesser applies.**

```
function apply_flow(edge, source, dest, pressure_diff, stream_conc, fraction):

    // ── Limit 1: Velocity (pipe physical capacity) ─────────────────────────────────────────
    r_eff     = min(source.radius, dest.radius)
    velocity  = min(pressure_diff × g.conductance_per_pressure, world.max_phloem_velocity)
    ring_area = phloem_ring_area(r_eff, t)       // pipe cross-section (ring only, not π×r²)
    flow_vol  = velocity × ring_area × fraction  // dm³ — sap volume moving this partial tick
    max_by_velocity = flow_vol × stream_conc     // g — sugar dissolved in that sap

    // ── Limit 2: Equalization (thermodynamic endpoint) ─────────────────────────────────────
    // If both nodes pooled their sugar across their combined volume, each would settle to:
    //     equil_conc = (source.sugar + dest.sugar) / (source.volume + dest.volume)
    // The sender CANNOT drop below equil_conc — the gradient has vanished, flow stops.
    // This is the same equalization clamp used in Node::compute_transport_flow() for diffusion.
    src_vol    = node_volume(source)
    dest_vol   = node_volume(dest)
    equil_conc = (source.sugar + dest.sugar) / max(src_vol + dest_vol, 1e-8)
    max_by_equalization = max(0, source.sugar - equil_conc × src_vol)

    // ── Transfer: lesser of pipe capacity and thermodynamic limit ──────────────────────────
    // Velocity limits HOW FAST per tick. Equalization limits HOW MUCH in total.
    // Adjacent nodes at similar concentrations: equalization is binding.
    // Long-distance flow through many hops: velocity is binding.
    desired_sugar = min(max_by_velocity, max_by_equalization)
    desired_sugar = min(desired_sugar, source.sugar)  // floating-point safety: should rarely fire

    // ── Unloading: fraction of sugar that leaves sieve tube into dest tissue ───────────────
    // gradient × perm is the unloaded fraction (perm has implicit units 1/(g/dm³)).
    dest_local_conc = (dest_vol > 0) ? dest.sugar / dest_vol : 0.0
    gradient = max(0, stream_conc - dest_local_conc)
    unload   = min(gradient × unloading_permeability(dest.type) × desired_sugar, desired_sugar)

    delta[dest]   += unload
    delta[source] -= desired_sugar   // source pays the full transfer;
                                     // transit remainder credited back via dead-end credit below

    // Stream dilutes by what unloaded — less sugar per volume for next hop
    stream_conc_after = (desired_sugar > 1e-6f) ? stream_conc × (1.0 - unload / desired_sugar) : 0.0f
    return stream_conc_after
}
```

**Why equalization is the key concept:**

- A tiny meristem with high concentration next to a large stem: equalization allows the meristem to push a small amount of sugar (small volume × small gradient = small max_equil). The stem's concentration barely changes.
- A large stem with stored sugar next to a small empty leaf: equalization allows substantial flow — but the leaf fills quickly (small volume → conc rises fast) and max_equil shrinks toward zero as concentrations converge. Flow stops naturally without any explicit cap.
- Two stems of equal volume, one with 1g and the other 0g: max_equil = 0.5g. They settle to equal concentration.
- Flow stops when concentrations equalize — the gradient disappears, `max_by_equalization` hits zero, `desired_sugar` = 0. This is the self-regulating property: no node gets over-drained or over-filled, even without explicit per-node caps.

This is the same equalization clamp in `Node::compute_transport_flow()` — the principle is proven in the existing diffusion code. Here it's applied to bulk phloem flow.

**Transit dead-end credit (conservation fix):** When BFS cannot continue from `dest`
(no downhill neighbours or budget exhausted), the transit fraction `desired_sugar − unload`
has nowhere to go and disappears from the sugar balance. Credit it back to `dest`:

```
bool will_continue = remaining > 1e-4f && stream_after > 1e-8f && has_downhill_neighbours(dest);
if (!will_continue) {
    delta[dest] += (desired_sugar - unload);   // transit stays in dest's sieve tube
}
```

This ensures `delta[source] + delta[dest] = 0` for every hop, giving exact conservation. The `min(desired_sugar, source.sugar)` safety clamp should almost never fire — equalization naturally prevents over-draft before the clamp is needed.

**Destination cap guard:** Cap `desired_sugar` at destination room before recording the delta, to prevent over-saturation when multiple BFS walks target the same conduit:

```
float dest_cap  = node_volume(dest) * world.max_sugar_concentration;
float dest_room = dest_cap - (dest.sugar + accumulated_delta[dest]);
desired_sugar = std::min(desired_sugar, std::max(0.0f, dest_room));
```

The stream starts at source concentration. At each sink it passes through, some sugar is unloaded. The first hungry sink (meristem, young leaf) gets the richest stream. A distant sink receives what's left. Priority emerges from physics — no explicit priority list needed.

**Unloading permeabilities by node type** (these are genome params — see Section 6b for calibration):

| Node type | Permeability | Biological basis |
|---|---|---|
| APICAL, ROOT_APICAL (active) | HIGH — `g.phloem_unloading_meristem` = 0.08 | Actively growing, fast metabolic consumption |
| LEAF (young, still expanding) | MODERATE — `g.phloem_unloading_leaf` = 0.10 | Transitioning sink→source; still net consumer |
| ROOT (mature conduit root) | LOW — `g.phloem_unloading_root` = 0.008 | Mainly conduit; some maintenance leakage |
| STEM (mature conduit) | VERY LOW — `g.phloem_unloading_stem` = 0.002 | Sealed phloem passing through; minimal loss |

Inactive (dormant) meristems use `phloem_unloading_stem` — they are not consuming and have equilibrated with surrounding sugar, so the gradient is shallow and almost nothing unloads. They don't need to be explicitly excluded.

### 4.4 Multi-source interaction

Multiple source nodes each run sequential BFS walks sharing a single accumulated delta array.
A conduit stem may receive contributions from both a leaf-loaded stem above it and the seed
below; both contribute to the stem's running delta. **Important:** each BFS walk must use the
running delta to compute destination room:

```
nxt_room = sugar_cap(nxt) - (nxt.sugar + accumulated_delta[nxt])
```

Without this, a second source walk would not know the first already filled a destination and
could add sugar beyond the ring cap. The atomic design (all flows computed from start-of-tick
state, applied together at end) fails for multiple sources competing for the same ring —
a shared running delta is required.

### 4.5 Last-mile delivery to non-vascular nodes

Leaves and meristems are not in the vascular network (no `has_vasculature` membership). They are connected to their parent vascular conduit via local diffusion during the DFS tick. This is unchanged: `transport_with_children` handles the final hop from nearest vascular stem into the leaf mesophyll or meristem tissue. The sugar that arrived at the conduit stem via phloem_resolve is available to diffuse into non-vascular children during the next DFS tick.

---

### 4.6 Phase 2 Revision: Jacobi Simultaneous Resolution

> **Status: planned replacement.** The BFS-per-source approach in sections 4.2–4.4 is the current implementation target. This section describes the planned next iteration that fixes the ordering-bias problem. The swap is isolated to `phloem_resolve()` — all data structures, pressure formulas, and endpoint passes (leaf loading pre-BFS, meristem unloading post-BFS) remain unchanged.

#### The problem with BFS-per-source

Sources are processed sequentially. The first source gets priority — later sources find depleted capacity or already-filled nodes. Processing order determines who "wins," which is arbitrary. In real phloem, all sources push simultaneously.

#### The Jacobi approach

Each tick of `phloem_resolve` runs `phloem_iterations` inner iterations (WorldParam, default 3–5). Each iteration resolves **all edges simultaneously** from the current sugar state:

```cpp
for (int iter = 0; iter < world.phloem_iterations; ++iter) {

    // 1. Recompute pressure at every vascular node from current sugar
    for (int i = 0; i < N; ++i)
        pressure[i] = has_vasculature(*flat[i].node, g)
            ? compute_phloem_pressure(*flat[i].node, g, world) : 0.0f;

    // 2. Compute desired flow for EVERY edge simultaneously.
    //    All reads from start-of-iteration sugar — no intra-iteration updates.
    std::vector<float> iter_delta(N, 0.0f);

    for (int i = 0; i < N; ++i) {
        if (!has_vasculature(*flat[i].node, g)) continue;
        Node& cur = *flat[i].node;

        auto process_edge = [&](int j) {
            if (!has_vasculature(*flat[j].node, g)) return;
            float dp = pressure[i] - pressure[j];
            if (dp <= 0.0f) return;   // cur is higher pressure: cur is sender

            Node& next = *flat[j].node;
            float r_eff      = std::min(cur.radius, next.radius);
            float ring_area  = phloem_ring_area(r_eff, world.phloem_ring_thickness);
            float velocity   = std::min(dp * g.conductance_per_pressure,
                                        world.max_phloem_velocity);
            float flow_vol   = velocity * ring_area;  // dm³/tick — velocity × cross-section

            float cur_vol    = node_volume(cur, world);
            float cur_conc   = (cur_vol > 0)
                               ? cur.chemical(ChemicalID::Sugar) / cur_vol : 0.0f;
            float max_velocity = flow_vol * cur_conc;  // g — velocity-limited sugar

            // Equalization: cur cannot drop below the equilibrium concentration.
            // Same principle as Node::compute_transport_flow() diffusion equalization clamp.
            float next_vol   = node_volume(next, world);
            float equil_conc = (cur.chemical(ChemicalID::Sugar) + next.chemical(ChemicalID::Sugar))
                               / std::max(cur_vol + next_vol, 1e-8f);
            float max_equil  = std::max(0.0f, cur.chemical(ChemicalID::Sugar) - equil_conc * cur_vol);

            // Lesser of pipe capacity and thermodynamic limit
            float desired_sugar = std::min(max_velocity, max_equil);
            desired_sugar = std::min(desired_sugar, cur.chemical(ChemicalID::Sugar));  // safety net

            // Apply unloading: split desired_sugar into conduit transit + tissue delivery.
            // Same permeability formula as Section 4.3.
            float next_local_conc = (next_vol > 0)
                                   ? next.chemical(ChemicalID::Sugar) / next_vol : 0.0f;
            float gradient = std::max(0.0f, cur_conc - next_local_conc);
            float perm     = unloading_permeability(next.type, g);
            float unload   = std::min(gradient * perm * desired_sugar, desired_sugar);

            iter_delta[i] -= desired_sugar;   // cur pays for full flow
            iter_delta[j] += unload;          // next receives unloaded fraction
            // transit remainder (desired_sugar - unload) stays in the sieve tube and is
            // credited back to cur at the cap-apply step (it never actually left cur's
            // compartment — only the unloaded fraction crosses the membrane)
        };

        if (flat[i].parent_idx >= 0) process_edge(flat[i].parent_idx);
        for (int ci : flat[i].child_idxs) process_edge(ci);
    }

    // 3. Fairness scaling: if any node's total outflow exceeds its sugar, scale down.
    //    Ensures conservation — no node goes negative.
    for (int i = 0; i < N; ++i) {
        if (iter_delta[i] >= 0.0f) continue;   // net receiver — no scaling needed
        float available = flat[i].node->chemical(ChemicalID::Sugar);
        float total_out = -iter_delta[i];
        if (total_out > available + 1e-8f) {
            float scale = available / total_out;
            iter_delta[i] = -available;         // node sends exactly what it has
            // Scale back all inflows that originated from this node.
            // Because we iterate edges from the sender's perspective, receivers
            // credited from this node need their share scaled proportionally.
            // Implementation: track sender per receiver in a separate pass, or
            // accept slight over-credit (bounded by ring cap) and clamp at apply.
            // The clamp at step 4 is the safety net; per-edge scaling is the refinement.
        }
    }

    // 4. Apply iteration delta — clamp to [0, cap]
    for (int i = 0; i < N; ++i) {
        if (std::abs(iter_delta[i]) < 1e-10f) continue;
        Node& n = *flat[i].node;
        float cap = node_volume(n, world) * world.max_sugar_concentration;
        n.chemical(ChemicalID::Sugar) =
            std::clamp(n.chemical(ChemicalID::Sugar) + iter_delta[i], 0.0f, cap);
    }
}
```

**`phloem_iterations` — WorldParam, default 3–5.** Each iteration propagates sugar one pipe section and runs one round of equalization on every edge. 4 iterations = sugar moves up to 4 edges per tick and the concentration landscape converges toward equilibrium over 4 steps per tick. The iteration count replaces the BFS time budget as the distance control. Unlike the BFS budget, iterations are discrete hops — simpler to reason about, though less physically precise about long vs. short edges. Equalization per iteration means the system self-corrects: nodes that overshoot in one iteration get corrected in the next.

#### Why this is better than BFS-per-source

| Property | BFS-per-source | Jacobi |
|---|---|---|
| Ordering bias | Yes — first source wins | None — all sources push simultaneously |
| Multiple sinks from same conduit | Sequential fills; early walkers deplete capacity | Proportional by gradient × conductance_per_pressure |
| Conservation | Transit dead-end credit + shared running delta required | Fairness scaling; cap clamp is safety net |
| Implementation complexity | Queue, source sorting, running delta bookkeeping | Two O(edges) passes per iteration; no queue |
| Distance control | Continuous (time budget × r² speed) | Discrete (iteration count × hops) |
| Per-tick cost | O(nodes) amortized over BFS walks | O(edges × iterations) — same order |

#### Swap footprint

The Jacobi revision is contained entirely within `phloem_resolve()`. The BFS implementation (Task 5) and the Jacobi revision share the same function signature, the same pressure formula (`compute_phloem_pressure`), the same velocity-capped flow formula (`velocity × phloem_ring_area` with `max_phloem_velocity` cap), the same volume helpers (`node_volume`, `phloem_ring_area`), and the same endpoint passes (leaf loading pre-BFS, meristem unloading post-BFS). Swapping is a targeted rewrite of the resolution loop, not an architectural change.

#### New WorldParam

Add `phloem_iterations` alongside the other Münch WorldParams in Task 3:

```cpp
uint32_t phloem_iterations = 4;   // inner Jacobi iterations per tick — controls how many
                                   // pipe sections sugar can travel per tick.
                                   // 3 = conservative (short-range flow); 5 = aggressive
                                   // (sugar travels further per tick, approaches BFS reach).
```

---

## 5. What Gets Deleted from the Codebase

These are removed entirely — not simplified or refactored, deleted:

| Deleted | Location | Why |
|---|---|---|
| `sugar_reserved_for_growth` field | `src/engine/node/node.h` | Grow-before-transport is now inherent to tick order |
| `pre_transport_growth()` function | `src/engine/plant.h/cpp` | Same — no longer needed |
| `compute_growth_reserve()` virtual | `src/engine/node/node.h` + all subclasses | Same |
| `phloem_reserve_fraction` genome param | `src/engine/genome.h` | Leaf retains its real post-consumption sugar naturally |
| `meristem_sink_fraction` genome param | `src/engine/genome.h` | Meristem pull bounded by gradient × conductance |
| Sugar Phase 1/Phase 2 logic in `run_vascular` | `src/engine/vascular.cpp` | Entire sugar pass replaced by `phloem_resolve` |
| Sugar source/sink classification (leaf supply, meristem demand, starvation sinks) | `src/engine/vascular.cpp` | No classification needed — pressure landscape drives flow |
| Sugar from `transport_with_children` | `src/engine/node/node.cpp` | Sugar no longer diffuses locally; phloem handles it |
| `phloem_loading_rate`, `phloem_loading_cost`, `phloem_water_coupling` | (were in v1 of this design) | Superseded — loading is implicit, water coupling is built into pressure formula |

Xylem Phase 1/Phase 2 (Water and Cytokinin) is **kept unchanged**.

---

## 6. New Parameters

### WorldParams additions

| Parameter | Default | Unit | Notes |
|---|---|---|---|
| `max_phloem_velocity` | 10.0 | dm/tick (= 1 m/hr) | Physical upper bound on phloem sap velocity — set by sieve plate resistance and sap viscosity. Real phloem: 0.3–1.5 m/hr; 10 dm/tick is the high end. Caps `velocity = pressure_diff × conductance_per_pressure` per edge. At max velocity through a young stem ring (~0.000176 dm²): max throughput = 10 × 0.000176 × 300 ≈ 0.53 g/tick — bounded but generous. |
| `base_phloem_speed` | 5.0 | dm/tick | _Superseded by velocity-capped model (Section 4.2). Not used in phloem_resolve. Retained for reference — remove in a future cleanup._ |
| `phloem_reference_radius` | 0.05 | dm | _Superseded by velocity-capped model. Not used. Retained for reference._ |
| `phloem_ring_thickness` | 0.002 | dm (= 0.2 mm) | Thickness of the **full active sugar-handling zone**: sieve tubes + companion cells + phloem parenchyma. Real dicots: 0.1–0.5 mm wide, constant regardless of stem diameter. Used in `phloem_ring_area(r, t) = π × (2r×t − t²)` for pipe capacity only — not for concentration. |
| `max_sugar_concentration` | 300.0 | g/dm³ | Upper limit of bulk sugar concentration per node — roughly a 30% sucrose solution, matching real phloem sap. Used to cap the atomic delta apply: `sugar_cap = node_volume(n) × max_sugar_concentration`. Replaces per-tissue density params with a single universal physical constant. |
| `leaf_thickness` | 0.003 | dm | Leaf mesophyll thickness used in `node_volume` for LEAF nodes. `volume = leaf_size² × leaf_thickness`. |
| `phloem_iterations` | 4 | iterations/tick | Inner Jacobi iterations per phloem_resolve call (Section 4.6 — planned replacement for BFS). Each iteration propagates sugar one edge. Default 4 = sugar travels up to 4 pipe sections per tick. |

These are WorldParams (physical constants), not genome params — they represent phloem sieve tube biology and plant cell physics, not plant strategy.

### Genome additions

| Parameter | Default | Notes |
|---|---|---|
| `phloem_osmotic_coefficient` | 1.0 | Maps absolute sugar concentration [g/dm³] to osmotic pressure. At 300 g/dm³ (max saturation) with full water: pressure = 300. The primary flow-rate calibration knob — evolution tunes how strongly concentration gradients drive flow. |
| `phloem_unloading_meristem` | 0.08 | Permeability for APICAL and ROOT_APICAL nodes. Active sink — growing tissue consumes continuously, keeping local sugar low. See calibration table below. |
| `phloem_unloading_leaf` | 0.10 | Permeability for LEAF nodes. Used for both loading (mature leaf → phloem) and unloading (young leaf ← phloem). Direction determined by gradient. See calibration table below. |
| `phloem_unloading_root` | 0.008 | Permeability for mature ROOT conduit nodes. Low — mainly sealed pipe with some maintenance leakage. See calibration table below. |
| `phloem_unloading_stem` | 0.002 | Permeability for STEM conduit nodes. Very low — sealed phloem, minimal leakage for surface maintenance only. See calibration table below. |

**`conductance_per_pressure`** — replaces `phloem_conductance` (old default 8.0; recalibrated to ~0.5). Maps osmotic pressure gradient [g/dm³] to phloem sap velocity [dm/tick]. Together with `max_phloem_velocity`, determines how quickly the phloem stream accelerates to its physical speed limit. At `conductance_per_pressure = 0.5`, a gradient of 20 g/dm³ gives velocity = 10 dm/tick (at cap); a gradient of 5 g/dm³ gives 2.5 dm/tick. Acts as the primary phloem throughput calibration knob — evolution tunes how strongly concentration gradients drive flow.

### Parameters removed (were in genome, now deleted)

`meristem_sink_fraction`, `phloem_reserve_fraction` — see Section 5.

`sugar_storage_density_wood`, `sugar_storage_density_leaf`, `sugar_cap_meristem`, `sugar_cap_minimum` — replaced by `node_volume(n) × world.max_sugar_concentration`. Sugar capacity is now derived from total node volume and a single universal concentration cap, not from per-tissue density lookups.

**Broader cleanup note:** The existing `sugar_cap()` helper and per-tissue density params exist throughout the sim outside the Münch code. The Münch implementation is the first adopter of `node_volume`-based caps. The rest of the sim (maintenance costs, local diffusion headroom checks, etc.) still uses `sugar_cap()`. These should eventually be unified — replace `sugar_cap()` everywhere with `node_volume(n) × max_sugar_concentration`. This is out of scope for the Münch implementation but the design here is intentionally forward-compatible: `max_sugar_concentration` as a WorldParam is the consolidation point.

---

## 6b. Permeability Calibration

### Anchor point

Leaf photosynthesis rate: `0.02 g/(dm²·hr)` at full sun (WorldParams). A full-size leaf (1.5 dm²) produces `0.03 g/tick` gross. After ~`0.01 g` maintenance + growth consumption during the DFS tick, it has roughly `0.02 g` surplus available to load into phloem. Leaf phloem volume ≈ `1.5² × 0.003 ≈ 0.007 dm³` (all living tissue), so surplus concentration ≈ `0.02 / 0.007 ≈ 3 g/dm³` above phloem stream. Permeability values are set so this surplus clears in approximately one tick.

For stems in the calibration table, pressure uses **full node volume** (total cylinder), not ring volume. A mature stem conduit (r=0.1, L=1.0 dm) has node_vol = π × 0.1² × 1.0 ≈ 0.0314 dm³. If it holds 0.1g sugar: bulk concentration ≈ 3.18 g/dm³. A young stem (r=0.015, L=0.3 dm) node_vol ≈ 2.12×10⁻⁴ dm³. With 0.001g sugar: bulk concentration ≈ 4.72 g/dm³. A leaf (leaf_size=0.5 dm) with 0.03g surplus: concentration ≈ 40 g/dm³. Leaf pressure (~40) exceeds stem pressure (~5) — flow goes leaf→stem as expected. Ring area is used only for pipe capacity in the flow formula, not for concentration.

### The formula

```
exchange = (phloem_concentration - local_concentration) × permeability
```

Both concentrations are in **g/dm³** (physical). The formula is **bidirectional**: when `local > phloem`, sugar enters the pipe (loading). When `phloem > local`, sugar exits (unloading). The same permeability param handles both directions. Direction is determined entirely by the concentration gradient — no separate loading vs. unloading classification.

### Calibration table

| Tissue | Permeability | At typical gradient (g/dm³) | Exchange rate | Rationale |
|--------|-------------|-------------------|---------------|-----------|
| Mature leaf | 0.10 | ~3 g/dm³ (leaf ~5 g/dm³, phloem ~2 g/dm³) | ~0.03 g/tick per dm³ flow | Efficient loader with companion cells; clears full surplus in ~1 tick |
| Young leaf | ~0.06 | ~2 g/dm³ (growing, depleted local) | ~0.012 g/tick | Developing companion cells; still transitioning from sink to source |
| Meristem | 0.08 | ~2 g/dm³ (consumed, phloem moderate) | ~0.016 g/tick | Active sink; covers growth cost (~0.015g/tick) plus maintenance comfortably |
| Root tip | 0.08 | ~2 g/dm³ | ~0.016 g/tick | Same growth demand as shoot meristem |
| Mature stem | 0.002 | ~1 g/dm³ (conduit near equilibrium) | ~0.002 g/tick | Sealed pipe; covers surface-area maintenance, passes most flow through |
| Mature root | 0.008 | ~1 g/dm³ | ~0.008 g/tick | More living tissue than stem (ray parenchyma, root hairs); slightly more permeable |

**Key ratio:** Leaf loading should be ~20–100× stem leakage. At these values the ratio is 50× (0.10 / 0.002). This ensures most sugar flows *through* stem conduits to distant sinks rather than being absorbed locally by the pipe itself.

### Young leaf distinction

The table shows mature leaf (0.10) and young leaf (~0.06) as separate rows. In code, a single `phloem_unloading_leaf = 0.10` is used for all LEAF nodes. Young leaves have depleted local sugar from active expansion — this naturally steepens the gradient in the unloading direction, so they receive more sugar per tick than the permeability alone suggests. Mature leaves have local sugar near photosynthetic saturation, which flattens the gradient in the loading direction, so the effective loading rate is gradient-limited rather than permeability-limited. The single param handles both cases through gradient dynamics.

If this distinction becomes important to tune independently, add `phloem_unloading_leaf_young` and gate by `leaf_size < g.max_leaf_size × 0.9`.

### On directionality and starch mobilization

The bidirectional formula is how starch mobilization in stems creates local sources. When a stem mobilizes starch to sugar, local sugar concentration rises above the phloem stream concentration. The gradient reverses: `local > phloem`. The same permeability formula now loads sugar *into* the pipe from that stem node. This creates a local pressure peak that drives outward flow toward adjacent sinks — without any source classification, without any special case. The Münch model handles it identically to a productive leaf.

### Saturation (future)

Real leaf loading uses membrane transporters (SUC/SUT family) with Michaelis-Menten kinetics. When phloem sucrose concentration is high (saturated phloem), the effective permeability decreases — less sugar loads even if the leaf has surplus. The constant-permeability model is a valid first pass. Add saturation if leaves over-load when phloem is full:

```
effective_permeability = phloem_unloading_leaf × (1.0 - phloem_conc / world.max_sugar_concentration)
```

Leave this for after initial calibration confirms the behavior is actually a problem.

---

## 7. Integration with Existing Systems

### 7.1 Xylem (water + cytokinin)

Xylem runs as a separate `xylem_resolve()` pass after phloem, using the same Phase 1/Phase 2 structure. The ordering is:

```
phloem_resolve()  →  xylem_resolve()
```

Phloem uses the water values from last tick's xylem pass. This one-tick lag is biologically accurate — xylem-transported water becomes available for phloem osmosis one hour later.

The coupling is naturally unidirectional: phloem flow does not return water to xylem in this implementation (a future refinement). In real plants, phloem sieve tubes are surrounded by xylem; water cycles between them locally. This can be added by having `phloem_resolve` record total sugar flow per node and returning `flow × water_per_sugar_ratio` to xylem at each sink. Leave this for when the hydraulic integration becomes a priority.

### 7.2 PIN transport and canalization

The Münch model uses `velocity × phloem_ring_area(r, t)` for pipe flow volume, where velocity is bounded by pressure gradient and `max_phloem_velocity`. The ring area `π × (2r×t − t²)` grows linearly with radius (for r >> t) rather than quadratically, but the canalization feedback is structurally identical:

```
PIN auxin flux → auxin_flow_bias → cambium activation → radius grows →
ring_area increases → more sap volume per tick at the same velocity →
sugar delivery improves → growth improves → more auxin → more flux
```

Thicker pipes carry more sugar volume per tick because ring area grows linearly with radius — at the same pressure-driven velocity, a wider pipe moves proportionally more sap. The linear-vs-quadratic difference between ring area and full cross-section means the capacity advantage of wide pipes is somewhat moderated — a 10× radius increase gives ~10× ring capacity vs. 100× full-circle capacity. This is biologically accurate: phloem capacity scales with girth perimeter, not girth area.

### 7.3 Starch storage (storage/plan.md)

Starch mobilization integrates for free under Münch:

**Starch synthesis:** A node converts surplus sugar to starch. Local sugar drops. Phloem pressure drops. The node imports less from adjacent high-pressure sources and exports less to adjacent sinks. It naturally "steps back" from the phloem stream while banking reserves.

**Starch mobilization:** A node mobilizes starch to sugar (triggered by low local sugar + GA). Local sugar rises. Phloem pressure rises. The node becomes a local source, pushing sugar outward to adjacent lower-pressure sinks. A mid-trunk node mobilizing winter starch reserves becomes a source automatically — no allocation classification change required. Münch routes sugar from wherever pressure is highest to wherever it is lowest.

**Seed bootstrap:** The seed's starch reserve mobilizes as GA from the shoot apical triggers amylase. Seed sugar rises. Phloem pressure at the seed rises above adjacent nodes. Sugar flows outward to early growing regions. This replaces the current seed-as-special-source logic in Phase 1/Phase 2 entirely.

### 7.4 `sugar_reserved_for_growth` and `pre_transport_growth()` — deleted

Both are removed. The new tick order makes them unnecessary:

**Old:** `pre_transport_growth()` → `vascular_transport()` → `tick_recursive(grow)`  
**New:** `tick_recursive(grow + produce)` → `phloem_resolve()`

Nodes grow and consume during the DFS tick. Phloem sees the resulting landscape. No reservation needed.

### 7.5 Anti-teleportation correctness

The distance-based BFS naturally limits how far sugar propagates per tick — it can only reach as far as the time budget allows. The atomic delta application (all flows computed from start-of-phloem-resolve state, applied together at the end) prevents order-dependent mid-pass state corruption. This is stricter than the current Phase 2, which modifies source sugar as it walks down.

---

## 8. Implementation Steps

### Task 1: Reorder tick phases in `plant.cpp` — no behavior change

**Files:**
- Modify: `src/engine/plant.cpp`
- Modify: `src/engine/plant.h`

The DFS tick currently runs `transport_chemicals()` which calls `transport_with_children()` for sugar. Sugar must stop diffusing locally during the DFS — it now moves only via phloem_resolve. All other chemicals (auxin, GA, cytokinin, stress) keep their existing local diffusion.

- [ ] In `node.cpp`, find the `transport_with_children` call. Confirm it handles Sugar via local diffusion. Add a guard: skip Sugar in `transport_with_children` (it will be handled by phloem_resolve). Check that Water is also skipped (handled by xylem). The guard: `if (chem_id == ChemicalID::Sugar || chem_id == ChemicalID::Water) continue;`
- [ ] Stub out `phloem_resolve(Plant&, const Genome&, const WorldParams&)` and `xylem_resolve(Plant&, const Genome&, const WorldParams&)` as empty functions in `vascular.cpp` / `vascular.h`
- [ ] In `plant.cpp::tick_tree()`, replace current call order:
  ```cpp
  // OLD:
  pre_transport_growth(world);
  vascular_transport(*this, genome_, world);
  pin_transport(*this, genome_);
  tick_recursive(*nodes_[0], *this, world);

  // NEW:
  pin_transport(*this, genome_);
  tick_recursive(*nodes_[0], *this, world);
  phloem_resolve(*this, genome_, world);
  xylem_resolve(*this, genome_, world);
  ```
- [ ] Delete `pre_transport_growth()` from `plant.h` and `plant.cpp`
- [ ] Delete `compute_growth_reserve()` virtual from `node.h` and all subclass implementations
- [ ] Delete `sugar_reserved_for_growth` field from `node.h` and all references to it
- [ ] Run build: `/usr/local/bin/cmake --build build` — fix all compilation errors
- [ ] Run tests: `./build/botany_tests` — many sugar tests will break (expected — phloem is stubbed to no-op)
- [ ] Note which tests break and which pass, for comparison in Task 5
- [ ] Commit: `refactor: reorder tick phases — DFS first, phloem/xylem resolve after`

### Task 2: Move xylem into `xylem_resolve()`

**Files:**
- Modify: `src/engine/vascular.cpp`
- Modify: `src/engine/vascular.h`

Move the existing xylem Phase 1/Phase 2 from `vascular_transport` into the new `xylem_resolve` stub.

- [ ] In `vascular.cpp`, rename the xylem portion of `vascular_transport` into `xylem_resolve`. Keep `build_flat` and all `VascNodeInfo` logic. Remove all sugar-related Phase 1/Phase 2 code from `run_vascular` (the sugar path calls will be gone — only Water and Cytokinin remain in `run_vascular`)
- [ ] Delete `run_vascular` calls for `ChemicalID::Sugar` entirely from the old `vascular_transport` body
- [ ] Populate `xylem_resolve` to call `build_flat`, then run `run_vascular` for Water, then for Cytokinin, with supply/demand resets between passes — same as before
- [ ] Delete `vascular_transport()` from `vascular.h/cpp` and its call in `plant.cpp` (already replaced by `xylem_resolve` in Task 1)
- [ ] Run build and tests
- [ ] Commit: `refactor: extract xylem Phase1/Phase2 into xylem_resolve()`

### Task 3: Add WorldParams and Genome parameters

**Files:**
- Modify: `src/engine/world_params.h`
- Modify: `src/engine/genome.h`
- Modify: `src/evolution/genome_bridge.cpp`

- [ ] Add to `WorldParams` struct and `default_world_params()`:
  ```cpp
  float max_phloem_velocity     = 10.0f;   // dm/tick (= 1 m/hr) — physical upper bound on phloem velocity.
                                            //   Caps: velocity = min(dp × conductance_per_pressure, max_phloem_velocity)
                                            //   Real phloem: 0.3–1.5 m/hr; 10 dm/tick is the high end.
  float base_phloem_speed       = 5.0f;    // superseded by velocity-capped model — kept for reference, not used
  float phloem_reference_radius = 0.05f;   // superseded by velocity-capped model — kept for reference, not used
  float phloem_ring_thickness   = 0.002f;  // dm (= 0.2 mm) — thickness of the living phloem+cambium ring.
                                            //   Used in phloem_ring_area(r, t) = π × (2r×t − t²).
                                            //   Used ONLY for pipe capacity (flow_vol = velocity × ring_area).
                                            //   Real dicots: ~0.1–0.5 mm. Constant regardless of stem diameter.
  float max_sugar_concentration = 300.0f;  // g/dm³ — upper limit of bulk sugar concentration (~30% solution)
                                            //   Sugar cap = node_volume(n) × max_sugar_concentration.
                                            //   Replaces per-tissue density params.
  float leaf_thickness          = 0.003f;  // dm — leaf mesophyll depth for node_volume calculation
                                            //   volume = leaf_size² × leaf_thickness
  uint32_t phloem_iterations    = 4;       // inner Jacobi iterations per tick (Section 4.6, planned replacement)
  ```
- [ ] Add to `Genome` struct and `default_genome()`:
  ```cpp
  float phloem_osmotic_coefficient = 1.0f;  // maps sugar concentration [g/dm³] to osmotic pressure
  float conductance_per_pressure   = 0.5f;  // dm/tick per (g/dm³) — velocity per pressure gradient unit.
                                             //   velocity = min(dp × conductance_per_pressure, max_phloem_velocity)
                                             //   At gradient=20 g/dm³: vel=10 dm/tick (at cap).
                                             //   Replaces phloem_conductance (old default 8.0, different units).
  float phloem_unloading_meristem  = 0.08f; // APICAL, ROOT_APICAL — calibrated active sink
  float phloem_unloading_leaf      = 0.10f; // LEAF — bidirectional; loading & unloading
  float phloem_unloading_root      = 0.008f;// mature ROOT conduit — low leakage
  float phloem_unloading_stem      = 0.002f;// STEM conduit — very low, mostly sealed pipe
  ```
- [ ] Remove from `Genome` and `default_genome()`:
  - `meristem_sink_fraction`
  - `phloem_reserve_fraction`
  - `sugar_storage_density_wood`
  - `sugar_storage_density_leaf`
  - `sugar_cap_meristem`
  - `sugar_cap_minimum`
  - Grep to find every use site: `grep -r "meristem_sink_fraction\|phloem_reserve_fraction\|sugar_storage_density\|sugar_cap_meristem\|sugar_cap_minimum" src/ tests/`
  - Fix all compilation errors (existing `sugar_cap()` calls outside Münch still compile — they get updated separately per the broader cleanup note in Section 6)
- [ ] Add four unloading params to `genome_bridge.cpp` under `sugar_economy` linkage group. Remove the two deleted params from the bridge.
- [ ] Run build: `/usr/local/bin/cmake --build build`
- [ ] Run tests: `./build/botany_tests`
- [ ] Commit: `feat: add Münch genome/world params, remove allocation params`

### Task 4: Implement `phloem_resolve()` — pressure computation

**Files:**
- Modify: `src/engine/vascular.cpp`
- Modify: `src/engine/vascular.h`
- Create: `tests/test_munch_phloem.cpp` (or add to `test_vascularization.cpp`)

- [ ] Add static helpers `phloem_ring_area` and `node_volume` to `vascular.cpp`. Two separate volume concepts — see Section 4.1:
  ```cpp
  // Phloem ring cross-section area — used ONLY for pipe capacity (flow_vol = velocity × ring_area).
  // For stems/roots: living phloem+cambium layer; heartwood excluded. Scales linearly with r.
  static float phloem_ring_area(float r, float t) {
      constexpr float pi = 3.14159265f;
      return pi * (2.0f * r * t - t * t);
  }

  // Total node volume — used for concentration (pressure) and sugar cap.
  // Full cylinder for stems/roots: using total volume rather than ring volume prevents the
  // 85× asymmetry that inverts leaf→stem flow direction (see Section 4.1 worked examples).
  static float node_volume(const Node& n, const WorldParams& world) {
      constexpr float pi = 3.14159265f;
      switch (n.type) {
          case NodeType::STEM:
          case NodeType::ROOT: {
              float length = glm::length(n.offset);
              return pi * n.radius * n.radius * length;   // full cylinder
          }
          case NodeType::LEAF: {
              float ls = n.as_leaf()->leaf_size;
              return ls * ls * world.leaf_thickness;
          }
          case NodeType::APICAL:
          case NodeType::ROOT_APICAL: {
              return (4.0f / 3.0f) * pi * n.radius * n.radius * n.radius;
          }
          default:
              return 0.0f;
      }
  }
  ```
- [ ] Add static helper `compute_phloem_pressure` to `vascular.cpp`:
  ```cpp
  static float compute_phloem_pressure(const Node& n, const Genome& g, const WorldParams& world) {
      float vol = node_volume(n, world);      // full node volume — not ring volume
      if (vol <= 0.0f) return 0.0f;
      float sugar_conc = n.chemical(ChemicalID::Sugar) / vol;  // g/dm³ — bulk concentration
      float water_cap_val = water_cap(n, g);
      float water_frac = (water_cap_val > 0)
          ? std::clamp(n.chemical(ChemicalID::Water) / water_cap_val, 0.0f, 1.0f)
          : 0.0f;
      return sugar_conc * g.phloem_osmotic_coefficient * water_frac;
  }
  ```
- [ ] Write test `test_phloem_pressure_high_sugar_high_pressure`: create two STEM nodes with identical geometry (same radius, same offset length). Give one sugar = `0.9 × node_vol × max_conc`, the other `0.1 × node_vol × max_conc`. Both with full water. Verify the first has ~9× the pressure of the second. Note: node_vol = π × r² × length (full cylinder).
- [ ] Write test `test_phloem_pressure_zero_water_zero_pressure`: STEM node with high sugar but zero water. Verify `compute_phloem_pressure` returns 0.0f.
- [ ] Write test `test_phloem_pressure_meristem_after_consumption`: set up a meristem (APICAL, r=0.01dm) with very low sugar. Verify pressure is low relative to a leaf node with full post-photosynthesis sugar. Note: meristem uses sphere volume `(4/3)π×r³ ≈ 4.2×10⁻⁶ dm³`, leaf uses `leaf_size²×leaf_thickness`.
- [ ] Run tests: `./build/botany_tests` — new tests pass
- [ ] Commit: `feat: add compute_phloem_pressure helper`

### Task 5: Implement `phloem_resolve()` — distance-based BFS flow

**Files:**
- Modify: `src/engine/vascular.cpp`

This is the core of the implementation. `phloem_resolve` populates `delta[]` via BFS from each high-pressure source.

- [ ] Add `phloem_edge_velocity` helper (replaces `phloem_local_speed` — velocity now from pressure, not r²):
  ```cpp
  static float phloem_edge_velocity(float pressure_diff, const Genome& g, const WorldParams& w) {
      return std::min(pressure_diff * g.conductance_per_pressure, w.max_phloem_velocity);
  }
  ```
- [ ] Add `unloading_permeability` helper:
  ```cpp
  static float unloading_permeability(NodeType t, const Genome& g) {
      switch (t) {
          case NodeType::APICAL:      return g.phloem_unloading_meristem;
          case NodeType::ROOT_APICAL: return g.phloem_unloading_meristem;
          case NodeType::LEAF:        return g.phloem_unloading_leaf;
          case NodeType::ROOT:        return g.phloem_unloading_root;
          case NodeType::STEM:        return g.phloem_unloading_stem;
          default:                    return g.phloem_unloading_stem;
      }
  }
  ```
- [ ] Implement `phloem_resolve`:
  ```cpp
  void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {
      Node* seed = plant.seed_mut();
      if (!seed) return;

      // Build flat pre-order array (same helper as xylem_resolve)
      std::vector<VascNodeInfo> flat;
      flat.reserve(plant.node_count());
      build_flat(seed, -1, flat);
      int N = (int)flat.size();

      // Step 1: compute phloem pressure at every vascular node
      // pressure = (sugar / node_volume) × osmotic_coeff × water_frac   [g/dm³ × coeff]
      // node_volume is total node volume (full cylinder for stems/roots, sphere for meristems, slab for leaves)
      std::vector<float> pressure(N, 0.0f);
      for (int i = 0; i < N; ++i) {
          if (has_vasculature(*flat[i].node, g))
              pressure[i] = compute_phloem_pressure(*flat[i].node, g, world);
      }

      // Step 2: BFS from each vascular source node outward toward lower-pressure neighbors
      // Accumulate sugar deltas. Apply atomically after all walks complete.
      std::vector<float> delta(N, 0.0f);

      // A node is a source for BFS if it has higher pressure than at least one neighbor.
      // stream_conc is phloem sugar concentration [g/dm³] of the moving phloem sap.
      struct BFSEntry { int idx; float budget; float stream_conc; };
      std::vector<BFSEntry> queue;
      queue.reserve(N);

      for (int i = 0; i < N; ++i) {
          if (!has_vasculature(*flat[i].node, g)) continue;
          bool is_source = false;
          if (flat[i].parent_idx >= 0 && has_vasculature(*flat[flat[i].parent_idx].node, g))
              if (pressure[i] > pressure[flat[i].parent_idx]) is_source = true;
          for (int ci : flat[i].child_idxs)
              if (has_vasculature(*flat[ci].node, g) && pressure[i] > pressure[ci])
                  is_source = true;
          if (is_source) {
              float vol = node_volume(*flat[i].node, world);
              float src_conc = (vol > 0) ? flat[i].node->chemical(ChemicalID::Sugar) / vol : 0.0f;
              queue.push_back({i, 1.0f, src_conc});
          }
      }

      for (auto& entry : queue) {
          struct WalkState { int from_idx; int to_idx; float budget; float stream_conc; };
          std::vector<WalkState> walk = {{-1, entry.idx, entry.budget, entry.stream_conc}};

          while (!walk.empty()) {
              auto [from_i, cur_i, budget, stream_conc] = walk.back();
              walk.pop_back();

              Node& cur = *flat[cur_i].node;

              std::vector<int> nexts;
              if (flat[cur_i].parent_idx >= 0 && flat[cur_i].parent_idx != from_i
                  && has_vasculature(*flat[flat[cur_i].parent_idx].node, g)
                  && pressure[flat[cur_i].parent_idx] < pressure[cur_i])
                  nexts.push_back(flat[cur_i].parent_idx);
              for (int ci : flat[cur_i].child_idxs)
                  if (ci != from_i && has_vasculature(*flat[ci].node, g)
                      && pressure[ci] < pressure[cur_i])
                      nexts.push_back(ci);

              for (int next_i : nexts) {
                  Node& next = *flat[next_i].node;

                  float edge_len = glm::length(next.offset);
                  float pressure_diff = pressure[cur_i] - pressure[next_i];
                  float velocity = phloem_edge_velocity(pressure_diff, g, world);
                  float time_cost = (velocity > 1e-8f) ? edge_len / velocity : 1e30f;

                  float time_fraction = std::min(1.0f, budget / time_cost);
                  float remaining = budget - time_cost * time_fraction;

                  // Velocity-limited sugar: pipe physically moves at most this much per tick
                  float r_eff = std::min(cur.radius, next.radius);
                  float ring_area = phloem_ring_area(r_eff, world.phloem_ring_thickness);
                  float flow_vol = velocity * ring_area * time_fraction;  // dm³
                  float max_velocity = flow_vol * stream_conc;            // g (sugar in sap)

                  // Equalization limit: cur cannot drop below the equilibrium concentration.
                  // Same principle as the diffusion equalization clamp in compute_transport_flow().
                  float next_vol = node_volume(next, world);
                  float cur_vol  = node_volume(cur, world);
                  float equil_conc = (cur.chemical(ChemicalID::Sugar) + next.chemical(ChemicalID::Sugar))
                                     / std::max(cur_vol + next_vol, 1e-8f);
                  float max_equil = std::max(0.0f,
                      cur.chemical(ChemicalID::Sugar) - equil_conc * cur_vol);

                  // Actual transfer: lesser of pipe capacity and thermodynamic limit
                  float desired_sugar = std::min(max_velocity, max_equil);
                  // Running-delta guard: don't drain past what's accumulated so far
                  float src_available = cur.chemical(ChemicalID::Sugar) + delta[cur_i];
                  desired_sugar = std::min(desired_sugar, std::max(0.0f, src_available));

                  // Passive unloading: gradient in g/dm³ between stream and destination
                  float next_local_conc = (next_vol > 0)
                      ? next.chemical(ChemicalID::Sugar) / next_vol : 0.0f;
                  float gradient = std::max(0.0f, stream_conc - next_local_conc);
                  float perm = unloading_permeability(next.type, g);
                  float unload = std::min(gradient * perm * desired_sugar, desired_sugar);

                  delta[next_i] += unload;
                  delta[cur_i]  -= desired_sugar;

                  float stream_after = (desired_sugar > 1e-6f)
                      ? stream_conc * (1.0f - unload / desired_sugar) : 0.0f;

                  if (remaining > 1e-4f && stream_after > 1e-6f) {
                      walk.push_back({cur_i, next_i, remaining, stream_after});
                  }
              }
          }
      }

      // Step 3: apply deltas atomically
      // Sugar cap = node_volume × max_sugar_concentration (total node volume for all types)
      for (int i = 0; i < N; ++i) {
          if (delta[i] == 0.0f) continue;
          Node& n = *flat[i].node;
          float cap = node_volume(n, world) * world.max_sugar_concentration;
          float new_val = std::clamp(n.chemical(ChemicalID::Sugar) + delta[i], 0.0f, cap);
          n.chemical(ChemicalID::Sugar) = new_val;
      }
  }
  ```
- [ ] Run build: `/usr/local/bin/cmake --build build`
- [ ] Run tests: `./build/botany_tests`
- [ ] Commit: `feat: implement phloem_resolve() — Münch distance-based BFS`

### Task 6: Write Münch-specific unit tests

**Files:**
- Modify: `tests/test_munch_phloem.cpp` (new file) or `tests/test_vascularization.cpp`

- [ ] Write test `test_munch_flow_leaf_to_stem`: 2-node chain (leaf → stem). Give leaf high sugar + full water. Give stem low sugar. Call `phloem_resolve` once. Verify stem gained sugar, leaf lost sugar.
- [ ] Write test `test_munch_flow_bidirectional`: same chain reversed — stem has high sugar (starch mobilization scenario), leaf has low sugar. Verify flow goes stem→leaf (reversed direction).
- [ ] Write test `test_munch_no_flow_without_water`: leaf has high sugar, zero water. Verify stem receives no sugar (pressure = 0 because water_frac = 0).
- [ ] Write test `test_munch_distance_limits_reach`: 6-node chain with slow (thin) internodes. Source leaf at tip. Run one tick of phloem_resolve with calibrated speed so sugar cannot reach 5 hops away. Verify middle nodes have gained sugar, far node has not.
- [ ] Write test `test_munch_unloading_priority_near_over_far`: source leaf at one end, two sinks — one adjacent, one 3 hops away. Verify adjacent sink received more sugar and richer stream concentration than the far sink.
- [ ] Write test `test_munch_conduit_stem_low_unloading`: place a mature STEM between source and sink. Verify the stem's sugar change is < 5% of the amount transferred through it (low unloading permeability).
- [ ] Run: `./build/botany_tests` — all pass
- [ ] Commit: `test: add Münch phloem unit tests`

### Task 7: Tune WorldParams defaults and run visual validation

**Files:**
- Modify: `src/engine/world_params.h` (adjust defaults)

This step is calibration. Goals:
1. A leaf exporting at peak production can supply 1–2 active meristems
2. Distance penalty is visible — meristems adjacent to leaves are supplied before distant ones
3. Drought (soil_moisture = 0.1) visibly reduces canopy-to-root sugar flow
4. Young seedling (thin stems) transports sugar noticeably slower than an established plant

- [ ] Run `./build/botany_realtime --color sugar` and observe gradient from canopy to roots
- [ ] If sugar equilibrates too quickly between nodes: lower `base_phloem_speed` or `phloem_osmotic_coefficient`
- [ ] If meristems starve despite healthy leaves: increase `base_phloem_speed` or `phloem_unloading_meristem`
- [ ] If trunk vs. thin branch shows no speed difference: verify the r² term in `phloem_local_speed`
- [ ] Run `./build/botany_sugar_test --tree seedling --ticks 100` — production/maintenance ratio should remain ~4x
- [ ] Commit: `tune: calibrate Münch phloem speed defaults`

### Task 8: Update documentation

**Files:**
- Modify: `src/engine/vascular.cpp` (block comment)
- Modify: `src/engine/vascular.h`
- Modify: `CLAUDE.md` (Chemical Transport Model section, Tick Control Flow section)

- [ ] Update top block comment in `vascular.cpp` to describe the two separate resolve functions, explain Münch vs. xylem, note tick ordering
- [ ] Update `CLAUDE.md` "Tick Control Flow" section to reflect the new 4-phase order
- [ ] Update `CLAUDE.md` "Chemical Transport Model" table: Sugar row should read `Phloem (Münch pressure resolve)` not `Vascular (phloem)`
- [ ] Update `CLAUDE.md` "Key Design Decisions" to note `sugar_reserved_for_growth` is gone and why
- [ ] Commit: `docs: update Chemical Transport Model and tick flow for Münch phloem`

---

## 9. Expected Behavior Changes

### What improves

**Source-proximity advantage with physical speed.** A meristem two nodes from a leaf receives sugar in 2 ticks. A meristem ten nodes away may need many ticks to be supplied. Thin bottleneck stems slow delivery proportionally. This matches competitive dynamics in real plant branching.

**Bidirectional flow, no classification.** Any node with higher post-consumption sugar than its neighbors becomes a source. Spring starch mobilization, green stem photosynthesis, or a well-supplied side branch can all push sugar outward without being "classified" as sources.

**Drought stalls phloem automatically.** Water stress reduces osmotic pressure at leaves → slower phloem → root meristems receive less sugar → growth slows. No additional code. The hydraulic feedback loop is complete.

**Elimination of bookkeeping artifacts.** `sugar_reserved_for_growth` was a one-tick-stale estimate that sometimes caused leaves to under-supply phloem (over-reserved) or over-supply (under-reserved). The new model uses the actual post-consumption state — no estimation.

**Unloading priority by proximity.** The first hungry sink on a path from source gets the richest stream. Distant sinks get diluted leftovers. This creates natural priority without an explicit priority system.

### What may need monitoring

**Early seedling transport.** Young stems (r ≈ 0.015 dm) at reference_radius = 0.05 dm run at ~9% of base phloem speed. A 5-node seedling may require multiple ticks for sugar to transit from cotyledon to root tip. If the seedling cannot survive this delay, increase `base_phloem_speed` or reduce `phloem_reference_radius` to match initial radius more closely.

**Starvation from slow propagation.** Long chains of thin nodes act as phloem bottlenecks. A 30-node root with all internodes at initial radius could starve the tip unless phloem speed is calibrated realistically. This is the price of physical accuracy. If root tips stall, the calibration lever is `base_phloem_speed`.

**Tests with hardcoded sugar delivery amounts.** Any test that asserts `node.sugar == X` after transport will need recalibration. The physical delivery is now bounded by speed and time budget, not by supply/demand matching.

---

## 10. What Is Not Changing

- Xylem (water + cytokinin) keeps Phase 1/Phase 2 allocation. Transpiration pull is a different physical mechanism.
- PIN transport pass is unchanged.
- Canalization model (`auxin_flow_bias`, radius-driven thickening) is unchanged.
- `has_vasculature()` admission criteria unchanged.
- Local diffusion for Auxin, GA, Stress, Cytokinin during DFS — unchanged.
- Sugar local diffusion during DFS is removed (phloem handles sugar now).
- Gibberellin, ethylene, stress, growth, photosynthesis models — all unchanged.

---

## 11. Relationship to Other Milestone 2 Work

**Starch storage (storage/plan.md):** Münch is a prerequisite. Starch mobilization at a stem node should push sugar outward — this only works if the stem can be a source, which requires pressure-driven flow. Implement Münch first, then starch.

**World physics / soil model:** When the soil model adds spatial moisture variation, root water uptake becomes spatially heterogeneous. Phloem's water_frac term automatically responds: roots in dry zones have lower water, which propagates via xylem to the connected leaves, which then have lower phloem pressure. No additional coupling code needed.

**Production throttle:** When stem conduits are full (high sugar after starch storage is saturated), the pressure gradient from leaf to conduit collapses. Leaf phloem pressure barely exceeds conduit pressure → almost no flow → leaf sugar accumulates → photosynthesis naturally downregulates via the sugar signaling pathway (future). The Münch model provides the physical mechanism for this saturation signal implicitly.

**Evolution:** Four unloading permeability params (`phloem_unloading_meristem`, `_leaf`, `_root`, `_stem`) + `phloem_osmotic_coefficient` = 5 new evolvable parameters in the `sugar_economy` linkage group. Net change from deleted params: `meristem_sink_fraction` and `phloem_reserve_fraction` removed. Net: +3 parameters. All five new params have direct biological interpretations and clear effect directions for evolution to work with.
