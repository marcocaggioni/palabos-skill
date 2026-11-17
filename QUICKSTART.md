# Quickstart — palabos-skill

```bash
git clone https://github.com/marcocaggioni/palabos-skill.git
cd palabos-skill
uv sync            # Python side: numpy, matplotlib, pyyaml (see note below)
bash scripts/check-setup.sh
```

Palabos itself is C++, not Python — there is no PyPI package for
it, so `uv sync` cannot install it. Fetch and build it separately:

```bash
bash scripts/install-palabos.sh   # clones Palabos to ~/palabos
export PALABOS_ROOT=~/palabos
```

Then build the Couette case and run it:

```bash
cd examples/couette && mkdir -p build && cd build
cmake .. -DPALABOS_ROOT="$PALABOS_ROOT" && make
cd .. && ./../../run-couette.sh --run
```

Then start Copilot CLI in the repo root and invoke the skill:

```bash
copilot
```

```
/palabos-simulation run the Couette example in examples/couette/case.yaml
```

Or in one go:

```bash
copilot -i "Use the /palabos-simulation skill to run the Couette example in examples/couette/case.yaml"
```

Copilot must be started from the repo root for skill discovery.
If the skill was added during a session, run `/skills reload` first.
The skill has no `allowed-tools` frontmatter on purpose — Copilot
must ask before running shell commands or a solver run.
