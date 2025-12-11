# Nuclear-Membrane-Binding-Analysis

This repository contains the image-analysis pipeline and code used to quantify nuclear-membrane invaginations/folds, nuclear volume, and protein nuclear‑membrane adsorption for the manuscript:

**Endoplasmic reticulum disruption stimulates nuclear membrane mechanotransduction.**
<br>**DOI: https://www.nature.com/articles/s41556-025-01820-9**
---

The pipeline processes 3D fluorescence microscopy datasets, extracts morphology and intensity metrics, and exports publication‑ready measurements and visualizations.

---

## 📌 Overview

Key capabilities

- Robust 3D nuclear segmentation from z‑stacks/time‑lapse movies
- Morphometrics: volume and nuclear‑membrane folds/invaginations
- Protein rim↔center binding quantification
- Protein rim↔LMNB1 (median across nuclear mid‑section) binding quantification


## Visualization of Workflow and Segmentation Results

  ### Nuclear Segmentation Workflow ###
![Nuclear Segmentation Workflow](Pipeline%20Image/Nucleus%20Analysis%20Pipeline%201.png)

  
 ### Calculation of Nuclear Volume and Protein Adsorption ###
![Segmentation to Metrics Workflow ](Pipeline%20Image/Nucleus%20Analysis%20Pipeline%202.png)


 ### Sample Segmentation Masks Overlay Movie ###
![Sample Movies with Segmentation Masks Overlay](Pipeline%20Image/Lamin%20Folds%20Analysis%20Visualization.png)


 ### Quantification of Nuclear Membrane Invaginations ### 
![Nuclear Membrane Invagination Analysis Pipeline](Pipeline%20Image/workflow.gif)



## Data flow

1. `Image_Import.py` converts raw **.nd2** movies into a single **.hdf5** (HDF5) file that stores the image stack as a NumPy array.
2. `Nucleus_3D_Analysis_Program.py` loads that `.hdf5`, performs segmentation/quantification (via `Nucleus_Analysis_Module.py`), and writes:
   - Segmentation results (**.hdf5**)
   - Measurements (**.csv**)
3. `segmentation_result_visualization.py` (optional) generates segmentation‑mask and raw images overlay in Napari Viewer for the users to assess segmentation accuracy.

---

## ✨ Organization

This repository contains **three analysis programs** with the same layout:

| Directory                              | Description                                                                   |
| -------------------------------------- | ----------------------------------------------------------------------------------- |
|`Nucleus Analysis Lamin Analysis`|**Invagination/Fold analysis +Rim vs nucleoplasm binding ratio+ Volume Analysis**|
|`Nucleus Analysis Binding Over Nucleoplasm` | **Rim vs nucleoplasm binding ratio + Volume Analysis**|
|`Nucleus Analysis Binding Over Lamin Control`|**Rim vs LMNB1 (LMNB1 intensity across mid‑section) binding ratio + Volume Analysis**| 


> Each program directory includes the files below.

| File                                   | Purpose                                                                             |
| -------------------------------------- | ----------------------------------------------------------------------------------- |
| `environment.yaml`                     | Dependency specification (Conda/Mamba).                                             |
| `Image_Import.py`                      | Import raw `.nd2` movies → export a consolidated `.hdf5` (NumPy array inside HDF5). |
| `Nucleus_3D_Analysis_Program.py`       | Load the `.hdf5`, run segmentation & metrics via the module, save results.          |
| `Nucleus_Analysis_Module.py`           | Reusable functions for I/O, preprocessing, segmentation, feature extraction.        |
| `segmentation_result_visualization.py` | Present segmentation alongside intermediate results and raw data in the Napari viewer, enabling users to compare and assess segmentation accuracy.                                 |
|`aicssegmentation`                      | Contains accessory functions used by the program| 

---

## ⚙ Installation

We recommend **micromamba** for environment management. (Conda works as well.) We have provided an environmenmt.yaml file that contains all the packages (with specific versions) for the analysis pipeline.

```bash
# Micromamba
micromamba create -n nucleus_analysis --file environment.yaml
micromamba activate nucleus_analysis
```

```bash
conda env create -f environment.yaml
conda activate nucleus_analysis
```

---

## System Requirement

The scripts were tested with Python 3.11 on Mac OS Sequoia Version 15.4.1

## 🚀 Quick Start (GUI file selection)

1. **Open the folder** in **Visual Studio Code** or **PyCharm**.

2. **Run** **Image_Import.py**\
   • You will be prompted to select the folder that contains all raw **.nd2** movie(s).\
   • For each ND2 file, the script writes a single **.hdf5** file that contains the image stack as a NumPy array.\
3. **Run** **Nucleus_3D_Analysis_Program.py**\
   • Select one `.hdf5` file from step 2.\
   • The program produces:\
   – **Segmentation results** (.hdf5)\
   – **Measurement results** (.csv)
   
4. *(Optional)* **Run** **segmentation_result_visualization.py** to inspect overlays and accuracy of segmentation mask.


---

## 📤 Outputs

- **HDF5 (segmentation)**: labeled/processed image data for downstream checks
- **CSV (metrics)**: per‑object and/or per‑frame measurements (e.g., volume and protein adsorption ratio )


---

## 🧪 Reproducibility

- The **.hdf5** exported by `Image_Import.py` is the canonical intermediate used by the analysis program.
- Keep acquisition metadata (voxel size, z‑step, channels) with each dataset.
- Pin exact versions in `environment.yaml`; consider exporting with `--from-history` to keep it minimal.

---

## Demo and Walkthrough

- ### We have provided the demo data and screenshots of the key steps to demonstrate the analysis workflow.
  - The demo data is stored in the **Demo** directory. It is cropped from the original timelapse movies and only contains 3 timepoints.

- ### Below are the walkthrough of this analysis program.
   
   - **Step 1: Run Nucleus_3D_Anlaysis_Program.py**
  ![Step1](Demo%20Data%20Step/step1.png)
    
   - **Step 2: Choose Segmentation channel based on nuclear marker**
  ![Step 2](Demo%20Data%20Step/Step2.png)

   - **Step 3: Choose whether to analyze only on subset of total timelapse**
  ![Step 3](Demo%20Data%20Step/Step3.png)

   - **Step 4: Program is running**
  ![Step 4](Demo%20Data%20Step/Step4.png)

   - **Step 5: Use Napari viewer to visualize all intermediate steps and generate quantitative plots**
  ![Step 5](Demo%20Data%20Step/Step5.png)
  
   - **Step 5.5: If **Lamin Analysis** is selected, the pipeline generates additional visualizations and quantitative plots characterizing nuclear membrane folds**
  ![Step 5.5](Demo%20Data%20Step/Lamin%20Step.png)
   
   - **Step 6: Choose to delete specific tracks (e.g. broken) before saving the intermediate visualizations as hdf5 and analysis resluts as csv**
   ![Step 6](Demo%20Data%20Step/Step6.png)

## ❓ Troubleshooting

If you have any questions, feel free to email joeshenz123@gmail.com

---


