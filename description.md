# botany

an evolution sim where plants grow and evolve. we get to watch evolution occur over huge amounts of time. we review their lineages and dna and see how they changed over time.

# engine

a world is started and seeded by random plants. the engine can run the world headless. 

# plants

plants have these genes:
- height
- leaf size
- root depth
- root width
- color
- branch length


- germination time
- growth speed
  - used to derive water requirements and metabolic speed

## How plants grow
- they start as a seed that falls on the ground somewhere. this is their starting location
- there are "growth points" aka meristems. for the initial plants theres 2: up and down. for growth, each growth point is simulated per tick, and runs the plants genome to know what to do (grow/extend, branch, differentiate into leaf/branch/seed/whatever).
  - apical meristem: the growth point at the tip of the plant that allows it to grow taller. 
  - lateral meristem: the growth point along the sides of the plant that allows it to grow wider. it can also branch by creating new growth points.
  - root meristem: the growth point at the bottom of the plant that allows it to grow roots. it can also branch by creating new growth points.
  - axilary meristem: the growth point at the junction of a leaf and stem that allows it to grow branches. it can also branch by creating new growth points.
- every tick, every active meristem in a plant "takes a turn". 
- the plant grows by extending its growth points. it can also branch by creating new growth points. the plant can also differentiate its growth points into different types (leaf, branch, seed, etc) which have different properties and functions.
  - so a seed on tick one may grow its root growth point, and grow its stem growth point. on tick two, it may branch its stem growth point into two new stem growth points, and grow its root point more. on tick three, it may grow one of the stem growth points, and differentiate the other into a leaf. etc.

## plant parts
- trunk: supports the plant, transports water and nutrients. can grow and branch.
- branch: extends from the trunk or other branches, supports leaves and seeds. can grow and branch.
- leaf: photosynthesizes and produces energy for the plant. can grow but not branch.
- root: anchors the plant and absorbs water and nutrients. can grow but not branch.
- seed: reproductive unit that can grow into a new plant. can grow but not branch.

the difference between a trunk part and a branch is chemically depenent. 

# Genome

- apical meristem: radially distal
  - growth()
  - spawn axilary meristem node if axin level high enough
  - produce axin()
- axilary meristem
  - stay dormant or activate
  - suppress axin()
  - growth
  - spawn axilary meristem
- Leaf
- root meristem
  - produce cytokinin()
  - extend()

# apical meristem tick cycle

1. extend
2. produce auxin
3. form axilary meristem if auxin low enough

# axilary meristem tick cycle




# ui

a 3D world space 
