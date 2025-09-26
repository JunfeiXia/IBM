# Individual-Based Model for Plant Invasion Dynamics

Code repository for the individual-based model (IBM) used in the following publications:

* **Lu, Y., DeAngelis, D. L., Xia, J., & Jiang, J. (2022).**
  *Modeling the impact of invasive species litter on conditions affecting its spread and potential regime shift.*
  *Ecological Modelling, 468, 109962.*
  [https://doi.org/10.1016/j.ecolmodel.2022.109962](https://doi.org/10.1016/j.ecolmodel.2022.109962)

* **Lu, Y., Xia, J., Magee, L. J., & DeAngelis, D. L. (2023).**
  *Seed dispersal and tree legacies influence spatial patterns of plant invasion dynamics.*
  *Frontiers in Applied Mathematics and Statistics, 9, 1086781.*
  [https://doi.org/10.3389/fams.2023.1086781](https://doi.org/10.3389/fams.2023.1086781)

* **Lu, Y., Xia, J., Holt, R. D., & DeAngelis, D. L. (2024).**
  *Modeling the Effects of Spatial Distribution on Dynamics of an Invading Melaleuca quinquenervia (Cav.) Blake Population.*
  *Forests, 15(8), 1308.*
  [https://doi.org/10.3390/f15081308](https://doi.org/10.3390/f15081308)

---

## 📖 Overview

This repository contains MATLAB code for an **individual-based modeling framework** developed to study invasive plant dynamics, specifically focusing on:

* Litter feedback effects and potential regime shifts in invaded ecosystems.
* Seed dispersal mechanisms and the influence of legacy trees on invasion spread.
* Spatial distribution and population-level consequences for *Melaleuca quinquenervia*.

The IBM explicitly represents individuals, their growth, mortality, dispersal, and interactions with the environment, enabling flexible exploration of invasion dynamics across ecological contexts.

---

## ⚙️ Requirements

* MATLAB R2020a or later (earlier versions may also work).
* No specialized toolboxes required beyond base MATLAB functions.
* Scripts are self-contained; datasets are either generated within the code or included as `.mat` files.

---

## 🚀 Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/<your-username>/IndividualBasedModel-Invasion.git
   cd IndividualBasedModel-Invasion
   ```

2. Open MATLAB and run the main simulation scripts:

   ```matlab
   % Example: run invasion simulation with default parameters
   main_simulation
   ```

3. Adjust ecological and simulation parameters in the configuration section of the scripts (e.g., initial density, dispersal kernel, litter effects, simulation duration).

4. Outputs include:

   * Time series of population density and composition.
   * Spatial maps of individuals across the landscape.
   * Figures matching those in the related publications.

---

## 📌 Citation

If you use this code, please cite the corresponding papers:

```
Lu, Y., DeAngelis, D. L., Xia, J., & Jiang, J. (2022). 
Modeling the impact of invasive species litter on conditions affecting its spread and potential regime shift. 
Ecological Modelling, 468, 109962.

Lu, Y., Xia, J., Magee, L. J., & DeAngelis, D. L. (2023).
Seed dispersal and tree legacies influence spatial patterns of plant invasion dynamics. 
Frontiers in Applied Mathematics and Statistics, 9, 1086781.

Lu, Y., Xia, J., Holt, R. D., & DeAngelis, D. L. (2024).
Modeling the Effects of Spatial Distribution on Dynamics of an Invading Melaleuca quinquenervia (Cav.) Blake Population. 
Forests, 15(8), 1308.
```

---

## 📬 Contact

For questions or collaborations, please contact:
**Junfei Xia** – [junfei.xia@outlook.com](mailto:junfei.xia@outlook.com) • [www.junfeixia.com](https://www.junfeixia.com)
---

