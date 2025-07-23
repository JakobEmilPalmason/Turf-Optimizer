
# Turf-Optimizer

> **Optimal roll layout for synthetic turf, in seconds.**  <br>
> The turf-optimizer lays fixed-width strips on a given area, scans through 180 rotation angles, accounts for offset, and returns all Pareto optimal solutions (installing efficiency vs. least waste)

<p align="center">
  <img src="docs/demo.gif" width="400" alt="Demo: draw polygon → JSON → run script → PDF">
</p>


Built in collaboration with Reykjavík University and Metatron, Iceland's leading artificial turf installer. The project aimed to minimize waste and increase efficiency in turf installation. Shown to reduce waste for up to 15%, the program also managed to shorten quote-making from 1-2 hours to a couple of minutes. 

---

## Features

- Minimizes turf waste on polygons, concave or “unorthodox” shapes. Up to 15% waste reduction, significantly smoother than manual planning, and provides detailed instructions for the workforce.
- Outputs a shareable PDF summary (tables + layout figure).
- Handles the `{ "body": "<json string>" }` wrapper format used by the web tool.
- Comes with two example inputs (rectangle + irregular polygon).
- MIT licensed.

---
## System Architecture

> **Browser → JSON → Python optimizer → PDF**.  
> Full architecture is documented in [`docs/workflow.md`](docs/workflow.md).
---
## Quick start

```bash
# 1. Install deps
pip install -r requirements.txt   # or: py -m pip install -r requirements.txt

# 2. Run on the examples
python turf_optimizer.py examples/event_sample_rectangle.json
python turf_optimizer.py examples/event_sample_unorthodox.json
```

---
## Draw your area

**Draw your area and copy the JSON at:**  
**https://pjk.cse.mybluehost.me/turf-optimizer/**

1. Draw the polygon  
2. Pick turf + extras  
3. Copy the JSON shown on the page  
4. Run:

```bash
# 1. Install deps
pip install -r requirements.txt   # or: py -m pip install -r requirements.txt

# 2. Run your area
python turf_optimizer.py my_area.json

