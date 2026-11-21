# donana-hydroperiod-evidence

[![CI](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/ci.yml/badge.svg)](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/ci.yml)
[![Jupyter Book](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/jupyter-book.yml/badge.svg)](https://annefou.github.io/donana-hydroperiod-evidence/)
[![Docker](https://github.com/annefou/donana-hydroperiod-evidence/actions/workflows/docker.yml/badge.svg)](https://github.com/annefou/donana-hydroperiod-evidence/pkgs/container/donana-hydroperiod-evidence)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23945812.svg)](https://doi.org/10.5281/zenodo.23945812)
[![FAIR4RS](https://img.shields.io/badge/FAIR4RS-conformant-brightgreen)](docs/fair4rs-checklist.md)
[![FORRT](https://img.shields.io/badge/FORRT-replication-blue)](https://forrt.org/)
[![Science Live](https://img.shields.io/badge/Science%20Live-nanopub%20chain-purple)](nanopubs/PUBLISHED.md)
[![RO-Crate](https://img.shields.io/badge/RO--Crate-1.2-orange)](ro-crate-metadata.json)
[![Software Heritage](https://archive.softwareheritage.org/badge/origin/https://github.com/annefou/donana-hydroperiod-evidence/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/annefou/donana-hydroperiod-evidence)

> **Fast and Automatic Data-Driven Thresholding for Inundation Mapping with Sentinel-2 Data** — replication study.
> Reference paper: [10.3390/rs10060910](https://doi.org/10.3390/rs10060910)

**DRAFT, not reviewed.** LifeWatch ERIC evidence for Spain, and the technical side of the FAIR2Adapt Doñana transfer case. Two analyses for Natura 2000 site ES0000024 (Doñana):

1. **Flooding regime** (`pipeline/flooding/`). Every Sentinel-2 scene (2017/18 onwards) is resampled onto HEALPix cells (healpix-geo, WGS84 ellipsoid, NESTED, depth 19, about 12 m) before classification: reflectance conservatively, scene classification by majority. Flooding is mapped with a re-implementation of Kordelas et al. (2018), Sentinel-1 is tested for cloudy dates, and results are validated against the LAST-EBD Landsat flood masks and hydroperiod maps (Díaz-Delgado et al., Zenodo, CC BY 4.0). Outputs: hydroperiod per cell and hydrological year, anomaly against 1984/85–2016/17, lagoon hydroperiod against groundwater (CHG piezometers) and rainfall, and inputs for the [Doñana Digital Twin](https://gitlab.lifewatch.dev/lfw024-donana-digital-twin). Against the LAST-EBD Landsat flood masks (Díaz-Delgado et al., doi:10.5281/zenodo.3518182), on 37 matched dates, the seasonal-marsh flood maps reach Cohen's kappa 0.78 with all cells and 0.89 when cells within about 37 m of a flood edge in the reference mask are left out (`results/flooding/results_validation.json`).
2. **Climate DT fitness test** (`pipeline/climate_dt/`). Whether the DestinE Climate Change Adaptation Digital Twin (IFS-NEMO, IFS-FESOM; historical then SSP3-7.0) reproduces Doñana rainfall and temperature statistics (E-OBS, ERA5), before and after bias correction. Verdict: usable after bias correction, as a complement to the EURO-CORDEX ensemble (`pipeline/climate_dt/REPORT.txt`).

**Status (23 October 2026).** The scripts ran end to end in a work folder on 22 October 2026; summary results are in `results/` and figures in `figures/`. A full rerun from a clean clone has not been done yet. The Copernicus Land Water Cover Duration validation (`pixi run flood-wcd`) has not run: it needs a Copernicus Data Space token.

---

## Install and run

You can read the results without installing anything: the summary tables are in `results/`, the figures in
`figures/`, and the [Jupyter Book](https://annefou.github.io/donana-hydroperiod-evidence/) shows both. Rerunning the analysis needs the
steps below. Status: prototype (TRL 3). The scripts ran end to end on 22 October 2026 in a work folder; a full rerun
from a clean clone has not been done yet.

### Install

1. Install pixi (<https://pixi.sh>): `curl -fsSL https://pixi.sh/install.sh | bash`.
2. Get the code and the pinned environment:

   ```bash
   git clone https://github.com/annefou/donana-hydroperiod-evidence.git
   cd donana-hydroperiod-evidence
   pixi install
   ```

   Pixi installs everything from `pixi.toml` and `pixi.lock` into `.pixi/`; `pixi run <task>` runs a task in that
   environment without activating it.

### Accounts each step needs

| What you run | Data | Account needed |
| --- | --- | --- |
| `pixi run snakemake --cores 1` (flooding regime) | Sentinel-2 L2A and Sentinel-1 (Earth Search, Microsoft Planetary Computer), LAST-EBD flood masks (Zenodo), ICTS-Doñana levels, CHG piezometers, Natura 2000 boundary, OpenStreetMap | none |
| | ERA5 (rainfall for the hydroperiod drivers) | **Copernicus Climate Data Store** (free) |
| `pixi run snakemake --cores 1 climate_dt` (Climate DT fitness test) | Climate DT generation 2, IFS-NEMO and IFS-FESOM | **DestinE upgraded access** |
| | ERA5; E-OBS v31.0e (ECA&D, non-commercial research only) | **Copernicus Climate Data Store**; none for E-OBS |
| `pixi run snakemake --cores 1 flood_wcd` (not run yet) | Copernicus Land Water Cover Duration | **Copernicus Data Space** (free) |

Accounts are read from your home directory or the environment, never from the repository, and nothing in `work/`
(the downloads) is committed. Set up only the ones for the steps you run.

### Copernicus Climate Data Store (ERA5)

1. Create a free account at <https://cds.climate.copernicus.eu> and log in.
2. Accept the licence of each dataset used, on the dataset page under *Download* → *Terms of use*: `reanalysis-era5-single-levels-monthly-means` (flooding), `reanalysis-era5-single-levels` and `reanalysis-era5-single-levels-timeseries` (Climate DT test).
3. Copy your personal access token from your profile page into `~/.cdsapirc`:

   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your personal access token>
   ```

   `cdsapi` also reads `CDSAPI_URL` and `CDSAPI_KEY` from the environment. Requests wait in the CDS queue (minutes to
   hours); the downloads are resumable.

### Destination Earth: Climate DT (DestinE upgraded access)

1. Create an account on the DestinE Platform, <https://platform.destine.eu/>, and request **upgraded access**, which
   includes the Climate Change Adaptation Digital Twin. Eligibility and terms are set by the platform (at the time of
   writing: users in EU Member States and in countries associated with the Digital Europe Programme; check the current
   terms of use before you apply). Upgraded access is not automatic: allow time for the request to be reviewed.
2. Get the Polytope key (an offline token) with ECMWF's script; it asks for your DestinE user name and password and
   writes `{"user_key": "<offline token>"}` to `~/.polytopeapirc`:

   ```bash
   curl -sO https://raw.githubusercontent.com/destination-earth-digital-twins/polytope-examples/main/desp-authentication.py
   pixi run python desp-authentication.py
   ```

   `polytope-client` also reads `POLYTOPE_USER_KEY` from the environment. A *DestinE Data Lake* (DEDL) API key is a
   different credential and is refused by Polytope (HTTP 401).
3. Authentication is not authorisation: a valid key can still be refused the data if upgraded access has not been
   granted yet; the first request tells you.
4. Terms: raw Climate DT data may not be redistributed. The pipeline keeps downloads in `work/` (not in git) and
   commits only derived statistics; keep it that way in any fork.

### Copernicus Data Space (only for `flood_wcd`)

The Copernicus Land Water Cover Duration check (`pipeline/flooding/13_wcd_validate.py`) needs a free Copernicus Data
Space account (<https://dataspace.copernicus.eu>) and its access token in the environment variable `CDSE_TOKEN`. This
check has not run yet.

### Run

```bash
pixi run snakemake --cores 1              # flooding regime: every Sentinel-2 date 2017-2026 (long); CDS account for ERA5
pixi run snakemake --cores 1 climate_dt   # Climate DT fitness test: DestinE upgraded access and CDS account
pixi run snakemake --cores 1 flood_wcd    # Copernicus Land Water Cover Duration check: CDSE_TOKEN
```

The same steps are pixi tasks (`pixi run flooding`, `pixi run climate-dt`, `pixi run flood-wcd`).

Every step is a script in `pipeline/` with its usage in its docstring; `Snakefile` gives the order, inputs and
outputs, and copies the summary results from `work/` into `results/` and `figures/`.
`pixi run snakemake --cores 1 --dry-run` lists what would run without downloading anything.

### With Docker

The image holds the environment and the code, never an account. Put the accounts the steps need in a file that
stays on your machine (for example `~/credentials.env`, outside the repository) and pass it at run time:

```
CDSAPI_URL=https://cds.climate.copernicus.eu/api
CDSAPI_KEY=<CDS personal access token>
COPERNICUSMARINE_SERVICE_USERNAME=<Copernicus Marine user name>
COPERNICUSMARINE_SERVICE_PASSWORD=<Copernicus Marine password>
POLYTOPE_USER_KEY=<DestinE offline token>
```

```bash
docker run --rm --env-file ~/credentials.env -v "$PWD/work:/app/work" -v "$PWD/results:/app/results" \
    ghcr.io/annefou/donana-hydroperiod-evidence:latest pixi run snakemake --cores 1 climate_dt
```

Leave out the lines a step does not need. The image is published on GitHub's container registry with each release.

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

- This software: [`CITATION.cff`](CITATION.cff) → DOI [10.5281/zenodo.23945812](https://doi.org/10.5281/zenodo.23945812)
- The original paper: [10.3390/rs10060910](https://doi.org/10.3390/rs10060910)

## Acknowledgements

This repository was built from [`sciencelivehub/forrt-replication-template`](https://github.com/sciencelivehub/forrt-replication-template), part of the [Science Live platform](https://platform.sciencelive4all.org). The template is licensed MIT and contributions (especially new domain flavours under [`docs/domain-flavours/`](docs/domain-flavours/)) are welcome.
