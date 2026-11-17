# Quickstart — palabos-skill

```bash
git clone https://github.com/marcocaggioni/palabos-skill.git
cd palabos-skill
bash scripts/check-setup.sh
```

If Palabos is missing:

```bash
git clone https://gitlab.com/unigespc/palabos.git ~/palabos
export PALABOS_ROOT=~/palabos
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
