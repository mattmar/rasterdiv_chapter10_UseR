# Chapter 10: The rasterdiv package for measuring diversity from space

Author: Matteo Marcantonio  
Date: 2025-06-18  
Contact: marcantoniomatteo@gmail.com

---

## Overview

This repository contains the R script and data organisation needed to reproduce all analyses and figures (03–08) in Chapter 10. We apply the rasterdiv package to:

- Stebbins Cold Canyon (CA, USA) — NAIP CIR imagery  
- Macchiarvana (IT) — PlanetScope NDVI time series  

Key outputs:
1. NDVI maps for years 2016, 2020, 2022  
2. Rényi entropy (α = 0, 1, 2)  
3. Rao’s quadratic entropy (global & area‐based)  
4. Helical plots of seasonal NDVI  

---

## Repository Structure

├── cold_canyon_CIR/         # NAIP CIR & NDVI subset TIFFs  (from OSF link see below)
├── macchiarvana_NDVI/       # PlanetScope NDVI time series  (from OSF link see below)
├── results/                 # Precomputed .RDS files (optional) (from OSF link see below)
├── figures/                 # Output PNGs (figure03–figure08)  (from OSF link see below)
└── Chapter10_rasterdiv.R    # Main analysis script  

---

## Data & Licensing

- Download materials:
  https://osf.io/pgbmq/?view_only=73fed91d69b6409d95a72b39d6f24e5f  
- PlanetScope imagery is under a research licence and not included here. Contact the author for infos.

---

## License

This project is released under the MIT License. See LICENSE for details.
