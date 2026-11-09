# Progress in Neurobiology 2026

Image-analysis scripts from the Radulovic Lab for quantifying hippocampal cFos expression, rabies-labeled inputs to CA2, and CA2 projections.

## Requirements

[Fiji / ImageJ](https://fiji.sc/), with Python (Jython) support for the `.py` scripts. These scripts use the ImageJ API and run within Fiji rather than standalone Python.

## Contents

| Folder | Description |
| --- | --- |
| `cFos Quantification Code` | Region-specific macros for cropping and merging images, plus a batch macro for quantifying nuclear cFos expression. |
| `Rabies Quantification Code` | Scripts for counting the CA2 rabies starting population and calculating regional ipsilateral and contralateral connection-strength indices (CSI). |
| `CA2 Projection Quantification` | Script for measuring green-channel fluorescence within manually drawn ipsilateral and contralateral anatomical regions of interest (ROIs), normalized to local background. |

## Usage

Download the files and extract any ZIP archives before use.

### cFos quantification

Open the `.ijm` files in Fiji's script editor, select **IJ1 Macro**, and click **Run**. The region-specific macros prompt for separate blue- and red-channel image folders and an output folder, then allow manual positioning of the crop ROI. Input files are paired by sorted filename order.

Run `Quantification_annotated_v9 (CFOS).ijm` on multichannel TIFF images with **channel 1 = cFos** and **channel 2 = Hoechst**. It segments Hoechst-positive nuclei and measures background-subtracted cFos signal, exporting per-cell results, batch summaries, ROI sets, and quality-control images.

### Rabies quantification

Open the `.py` files in Fiji's script editor, select **Python (Jython)**, and click **Run**.

- `Rabies Starting Population Quantification (CA2).py`: Select a parent folder containing per-mouse subfolders of three-channel, single-plane TIFF images. Follow the prompts to draw CA2 ROIs and establish detection settings. Exports cell counts, per-mouse summaries, ROIs, and quality-control images.
- `Rabies CSI Quantification.py`: Analyzes one anatomical region per run and expects folders `M1`, `M2`, `M3`, and `M4`, each containing exactly three TIFF images. Draw ipsilateral and contralateral ROIs when prompted. Exports counts, CSI summaries, and quality-control images. CSI is the summed input-cell count for each side divided by that mouse's CA2 starting population.

### CA2 projection quantification

Run `CA2_Projection_ROI_Quantification.py` in Fiji using **Python (Jython)**. Select multichannel TIFF images with one 2D plane per channel, enter the region and channel settings, and draw anatomical and background ROIs when prompted. Exports fluorescence measurements and ROI sets. Normalized fluorescence is mean green fluorescence divided by local background mean.

## Analysis settings

Settings are specific to the original dataset. Review channel assignments, image calibration, cropping parameters, and detection thresholds before analyzing other data. The CSI script contains hard-coded starting populations (`M1=99`, `M2=227`, `M3=9`, `M4=167`) and fixed segmentation settings; update these for a different dataset.
