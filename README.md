# donana-hydroperiod-evidence

[![CI](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/ci.yml/badge.svg)](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/ci.yml)
[![Jupyter Book](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/jupyter-book.yml/badge.svg)](https://annefou.github.io/donana-hydroperiod-evidence/)
[![Docker](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/docker.yml/badge.svg)](https://github.com/annefou/donana-hydroperiod-evidence/pkgs/container/donana-hydroperiod-evidence)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/{{ZENODO_DOI}}.svg)]({{ZENODO_DOI}})
[![FAIR4RS](https://img.shields.io/badge/FAIR4RS-conformant-brightgreen)](docs/fair4rs-checklist.md)
[![FORRT](https://img.shields.io/badge/FORRT-replication-blue)](https://forrt.org/)
[![Science Live](https://img.shields.io/badge/Science%20Live-nanopub%20chain-purple)](nanopubs/PUBLISHED.md)
[![RO-Crate](https://img.shields.io/badge/RO--Crate-1.2-orange)](ro-crate-metadata.json)
[![Software Heritage](https://archive.softwareheritage.org/badge/origin/https://github.com/annefou/donana-hydroperiod-evidence/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/annefou/donana-hydroperiod-evidence)

> **Fast and Automatic Data-Driven Thresholding for Inundation Mapping with Sentinel-2 Data** — replication study.
> Reference paper: [10.3390/rs10060910](https://doi.org/10.3390/rs10060910)

**DRAFT, not reviewed.** LifeWatch ERIC evidence for Spain, and the technical side of the FAIR2Adapt Doñana transfer case. Two analyses for Natura 2000 site ES0000024 (Doñana):

1. **Flooding regime** (`pipeline/flooding/`). Every Sentinel-2 scene (2017/18 onwards) is resampled onto HEALPix cells (healpix-geo, WGS84 ellipsoid, NESTED, depth 19, about 12 m) before classification: reflectance conservatively, scene classification by majority. Flooding is mapped with a re-implementation of Kordelas et al. (2018), Sentinel-1 is tested for cloudy dates, and results are validated against the LAST-EBD Landsat flood masks and hydroperiod maps (Díaz-Delgado et al., Zenodo, CC BY 4.0). Outputs: hydroperiod per cell and hydrological year, anomaly against 1984/85–2016/17, lagoon hydroperiod against groundwater (CHG piezometers) and rainfall, and inputs for the [Doñana Digital Twin](https://gitlab.lifewatch.dev/lfw024-donana-digital-twin).
2. **Climate DT fitness test** (`pipeline/climate_dt/`). Whether the DestinE Climate Change Adaptation Digital Twin (IFS-NEMO, IFS-FESOM; historical then SSP3-7.0) reproduces Doñana rainfall and temperature statistics (E-OBS, ERA5), before and after bias correction. Verdict: usable after bias correction, as a complement to the EURO-CORDEX ensemble (`pipeline/climate_dt/REPORT.txt`).

**Status (23 October 2026).** The scripts ran end to end in a work folder on 22 October 2026; summary results are in `results/` and figures in `figures/`. A full rerun from a clean clone has not been done yet. The Copernicus Land Water Cover Duration validation (`pixi run flood-wcd`) has not run: it needs a Copernicus Data Space token.

---

## Quick start

```bash
git clone https://github.com/annefou/donana-hydroperiod-evidence.git
cd donana-hydroperiod-evidence
pixi install
pixi run flooding     # open data; ERA5 needs a CDS account (~/.cdsapirc). Long: every Sentinel-2 date 2017-2026
pixi run climate-dt   # DestinE upgraded access (~/.polytopeapirc) and CDS account
```

This repository publishes code, plots and summary statistics only; input data stay in `work/` (gitignored) and are downloaded from their sources. Figure sources and caveats: `figures/README.md`. CHG piezometer data are reproduced citing the source under the CHG legal notice and are not official data. Raw Climate DT extracts are not redistributed (DestinE terms); only derived statistics are reported.

**Inputs and licences.** Natura 2000 boundary (EEA). LAST-EBD flood masks and hydroperiod maps, ground-truth points (Zenodo, CC BY 4.0). ICTS-Doñana marsh water levels (CC BY 4.0). CHG piezometers (no explicit licence; CHG legal notice authorises reproduction citing the source; data "not official"). OpenStreetMap water polygons (ODbL). Sentinel-2 and Sentinel-1 (Copernicus terms; Earth Search and Microsoft Planetary Computer hosting). ERA5 (Copernicus licence). E-OBS (ECA&D data policy: non-commercial research and education only; used here only as a research cross-check, cite Cornes et al. 2018). DestinE Climate DT (DestinE terms).

## Built from a template

This repository was created from [`sciencelivehub/forrt-replication-template`](https://github.com/sciencelivehub/forrt-replication-template). The template ships an operating manual for AI assistants ([`CLAUDE.md`](CLAUDE.md), [`AGENTS.md`](AGENTS.md)), domain conventions ([`DOMAIN.md`](DOMAIN.md)), and reference docs (`docs/`) so that an AI working only inside this repository can guide a researcher from "paper PDF + GitHub repo" to "published FORRT chain + Zenodo DOI" with no other context.

If you are reading this in a fresh fork, run [`/init-template`](.claude/skills/init-template/SKILL.md) inside Claude Code to substitute the placeholder tokens with your details. (For other AI tools, see [`docs/ai-portability.md`](docs/ai-portability.md).)

After `/init-template`, do these one-time setup steps to enable the full CI/CD path:

- **Enable GitHub Pages** at *Settings → Pages → Source: GitHub Actions*. Until enabled, the Jupyter Book build runs but the deploy step is skipped (CI stays green).
- All three workflows share one **readiness guard** (`.github/actions/check-ready`). Before `/init-template` runs, the `.template-uninitialised` sentinel makes them skip with an informative `::notice::` (badges stay green); `/init-template` deletes the sentinel, which activates them. They also skip while `notebooks/*.py` are still scaffolds (Phase 2). **Once you've published a nanopub chain** (real URIs in `nanopubs/PUBLISHED.md`), a skip is treated as a bug and **fails the run loudly** — so a finished replication can't sit on silently-green-but-empty CI.

## Repository structure

```
.
├── CLAUDE.md / AGENTS.md       # operating manual for AI assistants
├── DOMAIN.md                   # domain flavour (current: biodiversity + earth observation)
├── USER_PREFERENCES.md         # per-user style (edit on first clone)
├── README.md                   # this file
├── LICENSE                     # MIT
├── CITATION.cff                # how to cite
├── codemeta.json               # software metadata (CodeMeta-2.0)
├── ro-crate-metadata.json      # research object packaging (RO-Crate 1.2)
├── pixi.toml + pixi.lock       # pinned dependencies (single source of truth; lockfile is per-platform)
├── Dockerfile                  # container build
├── Snakefile                   # pipeline orchestration
├── myst.yml + index.md         # Jupyter Book scaffold
├── paper/                      # the source paper PDF
├── data/                       # downloaded artefacts (gitignored)
├── notebooks/                  # jupytext .py pipeline (01–04)
├── nanopubs/                   # FORRT chain drafts + published-URI registry
├── docs/                       # reference material
├── figures/                    # curated figures used in the Jupyter Book
├── .github/workflows/          # CI, Jupyter Book, Docker
└── .claude/                    # Claude Code agents, skills, sandbox config
```

## What you get

This template bakes in conventions that took multiple replications to discover. By using it, you inherit:

- **FAIR4RS conformance** — see [`docs/fair4rs-checklist.md`](docs/fair4rs-checklist.md) for the principle-by-principle mapping.
- **Self-contained data downloads** — the first notebook fetches everything; no manual data prep.
- **`pixi.toml` + `pixi.lock` as single source of truth** — local dev, Docker, and CI all install the same per-platform-pinned env.
- **`prefix-dev/setup-pixi`-based CI** — caches the env, runs the pipeline with `pixi run`, executes notebooks via a glob, fails fast on a stale lockfile.
- **Jupyter Book deployment** — auto-deploys to GitHub Pages with `BASE_URL` set correctly. (Don't put `base_url` in `myst.yml` — MyST silently ignores it.)
- **Docker + GHCR + Zenodo image archival** — `release` trigger pushes to GHCR and (optionally) archives to Zenodo for long-term preservation.
- **RO-Crate packaging** — the entire repo is a navigable Research Object via `ro-crate-metadata.json` (Process Run Crate + Workflow RO-Crate profiles).
- **Six-step FORRT chain workspace** — `nanopubs/drafts/` has a field-by-field skeleton for each step. `nanopubs/PUBLISHED.md` is the URI registry.
- **Layered AI guidance** — `CLAUDE.md` (universal) + `DOMAIN.md` (swappable per field) + `USER_PREFERENCES.md` (per-user). See [`docs/ai-portability.md`](docs/ai-portability.md) for non-Claude AI tools.
- **Sandbox by default** — `.claude/settings.json` denies file ops outside the repo, so a fresh AI session can't accidentally read `~/.ssh/` or write to `/etc/`.

## The six FORRT chain steps

A complete FORRT chain has six steps published on [platform.sciencelive4all.org](https://platform.sciencelive4all.org):

```
Quote-with-comment  →  AIDA  →  FORRT Claim  →  Replication Study  →  Replication Outcome  →  CiTO Citation
```

(For question-rooted chains with no upstream paper, replace step 1 with PICO or PCC. See [`docs/chain-decision-tree.md`](docs/chain-decision-tree.md).)

Drafts live in [`nanopubs/drafts/`](nanopubs/drafts/) field-by-field. Published URIs go into [`nanopubs/PUBLISHED.md`](nanopubs/PUBLISHED.md).

Optional further layers:

- **Research Software nanopub** — for reusable upstream tools (not demo repos). See [`docs/forrt-form-fields.md`](docs/forrt-form-fields.md) § Research Software.
- **Research Synthesis nanopub** — when this chain is part of a multi-chain story. See [`docs/forrt-form-fields.md`](docs/forrt-form-fields.md) § Research Synthesis.

## After publishing

When the chain is live and the FAIR4RS checklist is green, drafting an announcement post is the next step. See [`docs/announcement-template.md`](docs/announcement-template.md) for the structural template (vision-piece-first; the worked replication is the payoff, not the lead).

For lower-level nanopub work — retraction, superseding, batch publishing — see [`docs/programmatic-nanopubs.md`](docs/programmatic-nanopubs.md).

## Citation

If you use this work, please cite both:

- This software: [`CITATION.cff`](CITATION.cff) → DOI [{{ZENODO_DOI}}]({{ZENODO_DOI}})
- The original paper: [10.3390/rs10060910](https://doi.org/10.3390/rs10060910)

## Acknowledgements

This repository was built from [`sciencelivehub/forrt-replication-template`](https://github.com/sciencelivehub/forrt-replication-template), part of the [Science Live platform](https://platform.sciencelive4all.org). The template is licensed MIT and contributions (especially new domain flavours under [`docs/domain-flavours/`](docs/domain-flavours/)) are welcome.
