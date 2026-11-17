# AGENTS.md — palabos-skill

Open-source counterpart of the private `mstar-skill` repo: same
discussion → case → run → report workflow, same Couette example,
Palabos solver instead of licensed M-Star. No license, no GUI.

- Skill: `skills/palabos-simulation/SKILL.md`, mirrored byte-for-byte
  at `.github/skills/palabos-simulation/SKILL.md` (the copy Copilot
  CLI discovers). Keep the two copies in sync; edit both per commit.
- No `allowed-tools` in the skill frontmatter, deliberately — Copilot
  must ask before shell/solver runs.
- `case.yaml` is the contract. Defaults are labelled; confirm before
  running. Never invent a number.
- Lattice units are the agent's job: `lattice.py` converts physical
  → lattice and enforces u_lb < 0.1. Never run with omega outside
  roughly [0.8, 1.9] without saying so.
- Validation gate: Newtonian torque vs analytic Couette torque within
  ~10% before any non-Newtonian curve. `couette_post.py` enforces it.
- v0.1 status: `lattice.py` and `couette_post.py` are tested here;
  `couette.cpp` is a template, not yet compiled against Palabos —
  see its `TODO (verify)` markers.
- Palabos itself is AGPLv3; this repo's own code is MIT.
