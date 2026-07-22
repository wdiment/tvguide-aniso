# TVGuide-Aniso

TVGuide-Aniso is an interactive web app for visualizing seismic velocities derived from
elastic tensors. It ships with a built-in database of 96 published rock tensors and lets
you plot velocity surfaces, average tensors together, run a fold model, decompose tensors
into symmetry components, and add your own tensors. The app is
built with Python, Plotly Dash, and gunicorn.

## Running it

```bash
./run.sh
```

On first run this calls `build.sh`, which creates a self-contained conda environment in
`./tvguide` (downloading Miniforge if no conda/mamba is on the system). It then starts
gunicorn and serves the app at `http://0.0.0.0:8000`.

Optional settings can go in a gitignored `.env` file in the repo root:

- `TVGUIDE_BIND` — address:port to bind (default `0.0.0.0:8000`)
- `MAPBOX_TOKEN` — Mapbox access token, if map features are used

Re-running `run.sh` stops the previous instance and starts a fresh one.

## Using the tool

When the page loads, the **Vp** velocity plot for tensor 1 is displayed. Large plots take
a few seconds to compute — the spinner shows while a plot is being generated.

**Selecting a tensor.** Pick a tensor number (1–96) from the dropdown and click
**Submit Tensor** to plot it. The **Database** button in the top navbar shows the full
tensor table (filterable, and exportable to CSV) so you can look up which number is which.

**Tools** (buttons under the title):

- **Visualize Velocities** — the main tool. Tabs switch between plot types: Vp, Vs1, Vs2,
  Vp/Vs1, Vs1 polarization & splitting time, the 3D versions of each, a back-azimuthal
  plot, and radial plots.
- **Average Tensors** — average 2–5 tensors. Choose how many, select each tensor, give
  each a weight (default 1), and submit to get the averaged tensor components.
- **Calculate Fold Model** — runs the fold model for the selected tensor, with its own
  set of 2D and 3D plot tabs.
- **Decompose Tensors** — breaks the selected tensor into isotropic, hexagonal, and
  orthorhombic components, with a percentage breakdown table and exportable components.

**Enter Your Own Tensor** (top navbar) — paste a tensor as 22 comma-separated values:
density (g/cm³) followed by the 21 independent stiffness components in Mbar, in the order
C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46,
C55, C56, C66. Click **Submit Input Tensor** and it is added to the session's database as
tensor 97 (then 98, and so on). User tensors work in every tool, same as database tensors.

## Repo layout

- `tvguide.py` — Dash app: layout and callbacks
- `vector_functions.py` — computation and plotting functions
- `databases/` — tensor database CSV and precomputed sphere-grid files
- `assets/` — CSS and logos served by Dash
- `run.sh` / `build.sh` — launcher and environment setup

For a code-level walkthrough, see [code_overview.md](code_overview.md).
