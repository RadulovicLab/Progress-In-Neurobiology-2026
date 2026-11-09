"""
CA2 Rabies Starting-Population Counter (Fiji / ImageJ Jython)

Purpose
-------
For each mouse folder inside a selected parent directory, process a fixed number
of 3-channel, single-plane TIFF images. The user manually draws a CA2 ROI on a
3-color composite (red CA2 marker + green virus + blue Hoechst). Only the GREEN
channel is used for cell detection/counting. Red and blue are visualization only.

Each image is rotated 90 degrees clockwise for display and analysis so that the
manually drawn ROI remains aligned with the green detection channel. The source
TIFF on disk is not modified.

The first image is used to set an absolute green-intensity threshold interactively.
The user also chooses a soma-cleanup strength on the first image. Cleanup consists
of N binary erosions followed by N binary dilations, which suppresses thin neurites
while approximately restoring thicker soma-sized objects. Threshold and cleanup
strength are then LOCKED for the entire run. Particle-size/circularity filters are
also fixed for the run. A detected object is counted only when its centroid lies
inside the manually drawn CA2 ROI.

Outputs
-------
1) User-named detailed CSV: one row per image.
2) Auto-generated *_Summary.csv: one row per mouse, with summed count.
3) QC folder containing, for every image (saved automatically):
   - CA2 ROI (.roi)
   - detected-cell ROIs (.zip)
   - overlay preview (.png)

Recommended Fiji usage
----------------------
File > New > Script, choose language "Python", open/paste this file, then Run.
This script is intended for Fiji's Jython/Python interpreter (ImageJ 1.x API).
"""

from ij import IJ, ImagePlus, Prefs
from ij.gui import GenericDialog, WaitForUserDialog, YesNoCancelDialog, Overlay, TextRoi
from ij.io import DirectoryChooser, SaveDialog, FileSaver, RoiEncoder
from ij.plugin import Duplicator, RGBStackMerge
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import Measurements, ResultsTable
from ij.process import ImageProcessor
from java.awt import Color, Font
from java.io import File
import os
import re
import csv


# ---------------------------- Helper functions ---------------------------- #

def natural_key(text):
    parts = re.split(r'(\d+)', text)
    out = []
    for p in parts:
        if p.isdigit():
            out.append(int(p))
        else:
            out.append(p.lower())
    return out


def safe_name(name):
    """Filesystem-friendly basename."""
    base = os.path.splitext(os.path.basename(name))[0]
    return re.sub(r'[^A-Za-z0-9._-]+', '_', base)


def list_tiffs(folder):
    names = []
    for name in os.listdir(folder):
        low = name.lower()
        if low.endswith('.tif') or low.endswith('.tiff'):
            if os.path.isfile(os.path.join(folder, name)):
                names.append(name)
    names.sort(key=natural_key)
    return names


def get_mouse_folders(parent):
    """Return immediate mouse folders named M1, M2, M3, ..."""
    folders = []
    for name in os.listdir(parent):
        p = os.path.join(parent, name)
        if os.path.isdir(p) and re.match(r'^M\d+$', name, re.IGNORECASE):
            folders.append((name, p))
    folders.sort(key=lambda x: natural_key(x[0]))
    return folders


def close_imp(imp):
    if imp is not None:
        try:
            imp.changes = False
            imp.close()
        except:
            pass


def channel_duplicate(imp, channel):
    """Duplicate one channel from a single-plane or hyperstack image."""
    nC = max(1, imp.getNChannels())
    nZ = max(1, imp.getNSlices())
    nT = max(1, imp.getNFrames())
    if channel < 1 or channel > nC:
        raise ValueError("Requested channel %d, but image has %d channels." % (channel, nC))
    # User specified single-plane images. If Z/T > 1, use current first plane only.
    return Duplicator().run(imp, channel, channel, 1, 1, 1, 1)


def rotate_clockwise_90(imp):
    """Rotate an ImagePlus 90 degrees clockwise in place."""
    IJ.run(imp, "Rotate 90 Degrees Right", "")
    return imp


def make_display_composite(imp, red_ch, green_ch, blue_ch, title):
    """Build a display-only RGB-style composite in R,G,B order, rotated clockwise."""
    r = channel_duplicate(imp, red_ch)
    g = channel_duplicate(imp, green_ch)
    b = channel_duplicate(imp, blue_ch)
    rotate_clockwise_90(r)
    rotate_clockwise_90(g)
    rotate_clockwise_90(b)
    comp = RGBStackMerge.mergeChannels([r, g, b], True)
    comp.setTitle(title)
    close_imp(r)
    close_imp(g)
    close_imp(b)
    return comp


def prompt_for_ca2_roi(comp, image_name, mouse_name):
    """Wait until user draws a valid area ROI on the displayed composite."""
    while True:
        comp.show()
        IJ.setTool("freehand")
        msg = (
            "Image: %s\nMouse: %s\n\n"
            "Draw the CA2 anatomical ROI on the composite.\n"
            "Use the red CA2 marker and/or Hoechst for anatomy.\n"
            "Red intensity is NOT quantified.\n\n"
            "Zoom and pan as needed. When the ROI is complete, click OK here."
        ) % (image_name, mouse_name)
        WaitForUserDialog("Draw CA2 ROI", msg).show()
        roi = comp.getRoi()
        if roi is not None and roi.isArea():
            return roi.clone()
        IJ.showMessage("ROI required", "Please draw a closed area ROI around CA2 before continuing.")


def prepare_detection_image(green_imp, rolling_ball, median_radius):
    """Duplicate the green channel and apply optional fixed preprocessing."""
    work = green_imp.duplicate()
    work.setTitle(green_imp.getTitle() + " - detection working copy")
    if rolling_ball > 0:
        IJ.run(work, "Subtract Background...", "rolling=%s" % str(rolling_ball))
    if median_radius > 0:
        IJ.run(work, "Median...", "radius=%s" % str(median_radius))
    return work


def calibrate_threshold(work, ca2_roi, auto_method, image_name):
    """
    Interactively set a threshold on the first image. The returned low/high
    values are absolute intensity thresholds and should be reused for all images.
    """
    work.show()
    work.setRoi(ca2_roi)

    # Give the user a reproducible automated starting point, then allow manual tuning.
    if auto_method != "None":
        try:
            IJ.setAutoThreshold(work, auto_method + " dark")
        except:
            IJ.setAutoThreshold(work, "Default dark")

    while True:
        IJ.run(work, "Threshold...", "")
        WaitForUserDialog(
            "Set global green threshold",
            "FIRST IMAGE CALIBRATION: %s\n\n"
            "Adjust the threshold so true green virus-positive SOMATA are selected.\n"
            "Do not click Apply in the Threshold window.\n"
            "When satisfied, click OK here.\n\n"
            "This absolute threshold will be locked for ALL images in this run."
            % image_name
        ).show()

        ip = work.getProcessor()
        low = ip.getMinThreshold()
        high = ip.getMaxThreshold()
        if low == ImageProcessor.NO_THRESHOLD or high == ImageProcessor.NO_THRESHOLD:
            IJ.showMessage("No threshold", "No threshold is active. Please set a threshold and try again.")
            continue
        return float(low), float(high)


def make_binary_mask(work, low, high, cleanup_passes, watershed):
    """Create a binary mask, suppress thin neurites, and optionally watershed."""
    mask = work.duplicate()
    # Detection is performed on the full image; CA2 inclusion is imposed later
    # by the centroid-inside-ROI rule. Do not clip particles at the ROI boundary.
    mask.deleteRoi()
    ip = mask.getProcessor()
    ip.setThreshold(low, high, ImageProcessor.NO_LUT_UPDATE)

    # Ensure bright selected objects become white foreground on black background.
    old_black = Prefs.blackBackground
    Prefs.blackBackground = True
    IJ.run(mask, "Convert to Mask", "")

    # Soma cleanup: N erosions followed by N dilations. This is intentionally
    # applied after thresholding so thin neurites are preferentially removed, while
    # thicker soma-sized objects are approximately restored. 0 = no cleanup.
    for i in range(int(cleanup_passes)):
        IJ.run(mask, "Erode", "")
    for i in range(int(cleanup_passes)):
        IJ.run(mask, "Dilate", "")

    if watershed:
        IJ.run(mask, "Watershed", "")
    Prefs.blackBackground = old_black
    return mask


def prompt_cleanup_passes(current_value, image_name):
    """Choose soma-cleanup strength on the first calibration image."""
    gd = GenericDialog("Soma cleanup calibration")
    gd.addMessage(
        "FIRST IMAGE CALIBRATION: %s\n\n"
        "Cleanup removes thin green processes before cell counting.\n"
        "0 = none; 1 = mild; 2 = stronger.\n"
        "The script performs N erosions followed by N dilations.\n\n"
        "Start with 1. After reviewing the numbered detections, choose NO if you want to try another value."
        % image_name
    )
    gd.addNumericField("Soma_cleanup_passes_(0-5):", float(current_value), 0)
    gd.showDialog()
    if gd.wasCanceled():
        return None
    value = int(round(gd.getNextNumber()))
    if value < 0:
        value = 0
    if value > 5:
        value = 5
    return value


def detect_cells(mask, ca2_roi, min_area, max_area, min_circ, max_circ):
    """
    Detect particles on the full image, then keep only objects whose geometric
    centroid falls inside the CA2 ROI. Returns a list of ROI clones.
    """
    rt = ResultsTable()
    rm = RoiManager(True)  # private/hidden manager
    ParticleAnalyzer.setRoiManager(rm)

    options = ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.CLEAR_WORKSHEET | ParticleAnalyzer.SHOW_NONE
    measurements = Measurements.AREA | Measurements.CENTROID
    pa = ParticleAnalyzer(options, measurements, rt, min_area, max_area, min_circ, max_circ)
    ok = pa.analyze(mask)
    if not ok:
        try:
            rm.reset()
        except:
            pass
        try:
            rm.close()
        except:
            pass
        return []

    all_rois = rm.getRoisAsArray()
    cal = mask.getCalibration()
    kept = []

    n = min(len(all_rois), rt.getCounter())
    for i in range(n):
        try:
            x_cal = rt.getValue("X", i)
            y_cal = rt.getValue("Y", i)
            x_raw = int(round(cal.getRawX(x_cal)))
            y_raw = int(round(cal.getRawY(y_cal)))
        except:
            # Fallback: center of particle bounding box.
            b = all_rois[i].getBounds()
            x_raw = int(round(b.x + b.width / 2.0))
            y_raw = int(round(b.y + b.height / 2.0))

        if ca2_roi.contains(x_raw, y_raw):
            kept.append(all_rois[i].clone())

    # Clear the temporary ROI Manager before closing so ImageJ does not
    # ask whether to save displayed ROIs. The kept ROIs have already been cloned.
    try:
        rm.reset()
    except:
        pass
    try:
        rm.close()
    except:
        pass
    return kept


def apply_review_overlay(comp, ca2_roi, cell_rois, count):
    overlay = Overlay()

    ca = ca2_roi.clone()
    ca.setStrokeColor(Color.yellow)
    ca.setStrokeWidth(2.0)
    ca.setName("CA2")
    overlay.add(ca)

    for i, roi in enumerate(cell_rois):
        r = roi.clone()
        r.setStrokeColor(Color.white)
        r.setStrokeWidth(1.0)
        r.setName("Cell_%04d" % (i + 1))
        overlay.add(r)

    label = TextRoi(10, 10, "Detected green+ cells in CA2: %d" % count,
                    Font("SansSerif", Font.BOLD, 16))
    label.setStrokeColor(Color.white)
    overlay.add(label)

    comp.setOverlay(overlay)
    comp.updateAndDraw()
    return overlay


def review_count(comp, image_name, mouse_name, count, is_first=False):
    if is_first:
        text = (
            "Image: %s\nMouse: %s\n\n"
            "Detected green-positive cells in CA2: %d\n\n"
            "Inspect the white cell outlines and yellow CA2 boundary.\n"
            "YES = accept this threshold + soma cleanup and lock them for the entire run.\n"
            "NO = return to threshold/cleanup calibration for this first image.\n"
            "CANCEL = stop the script."
        ) % (image_name, mouse_name, count)
    else:
        text = (
            "Image: %s\nMouse: %s\n\n"
            "Detected green-positive cells in CA2: %d\n\n"
            "Inspect the white cell outlines and yellow CA2 boundary.\n"
            "YES = accept this image.\n"
            "NO = stop the run so you can restart with revised GLOBAL settings.\n"
            "CANCEL = stop the script."
        ) % (image_name, mouse_name, count)

    dlg = YesNoCancelDialog(comp.getWindow(), "Review detected cells", text)
    if dlg.cancelPressed():
        return "cancel"
    if dlg.yesPressed():
        return "accept"
    return "retry" if is_first else "stop"


def save_roi_file(roi, path):
    enc = RoiEncoder(path)
    enc.write(roi)


def save_cell_rois(cell_rois, zip_path):
    rm = RoiManager(True)
    for i, roi in enumerate(cell_rois):
        r = roi.clone()
        r.setName("Cell_%04d" % (i + 1))
        rm.addRoi(r)
    # Save as ZIP if there is at least one detected cell.
    if len(cell_rois) > 0:
        rm.runCommand("Save", zip_path)
    else:
        # Make the absence explicit; no fake/empty ZIP is created.
        note_path = os.path.splitext(zip_path)[0] + "_NO_CELLS.txt"
        f = open(note_path, 'w')
        f.write("No particles passed the detection criteria within the CA2 ROI.\n")
        f.close()
    # Detected-cell ROIs are already saved above. Clear the temporary manager
    # before closing to suppress ImageJ's "save displayed ROIs?" prompt.
    try:
        rm.reset()
    except:
        pass
    try:
        rm.close()
    except:
        pass


def save_qc_png(comp, path):
    flat = comp.flatten()
    flat.setTitle("QC")
    FileSaver(flat).saveAsPng(path)
    close_imp(flat)


def csv_write(path, rows, header):
    f = open(path, 'wb')
    try:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    finally:
        f.close()


def validate_input(parent, expected_per_mouse):
    mice = get_mouse_folders(parent)
    if len(mice) == 0:
        raise ValueError("No immediate subfolders were found in the selected parent directory.")

    problems = []
    for mouse_name, mouse_path in mice:
        n = len(list_tiffs(mouse_path))
        if n != expected_per_mouse:
            problems.append("%s: found %d TIFF(s), expected %d" % (mouse_name, n, expected_per_mouse))
    if len(problems) > 0:
        raise ValueError(
            "Input validation failed. Each mouse folder must contain exactly %d TIFF images.\n\n%s"
            % (expected_per_mouse, "\n".join(problems))
        )
    return mice


# ------------------------------- Main script ------------------------------ #

def main():
    IJ.showStatus("Rabies starting-population counter: setup")

    dc = DirectoryChooser("Choose parent directory containing mouse subfolders (e.g., M1, M2, M3, M4)")
    parent = dc.getDirectory()
    if parent is None:
        return

    # Read first TIFF only to discover channel count and calibration for setup messaging.
    mouse_folders = get_mouse_folders(parent)
    if len(mouse_folders) == 0:
        IJ.showMessage("No folders", "No mouse subfolders were found.")
        return

    first_tiffs = list_tiffs(mouse_folders[0][1])
    if len(first_tiffs) == 0:
        IJ.showMessage("No TIFFs", "The first mouse folder contains no TIFF images.")
        return

    first_path = os.path.join(mouse_folders[0][1], first_tiffs[0])
    probe = IJ.openImage(first_path)
    if probe is None:
        IJ.showMessage("Open error", "Could not open the first TIFF:\n" + first_path)
        return

    n_channels = probe.getNChannels()
    cal = probe.getCalibration()
    area_unit = cal.getUnit() + "^2" if cal.scaled() else "pixel^2"
    close_imp(probe)

    if n_channels < 3:
        IJ.showMessage(
            "Channel error",
            "The first image reports only %d channel(s). This script expects at least 3 channels." % n_channels
        )
        return

    methods = ["Default", "Otsu", "Triangle", "Yen", "Li", "Huang", "Moments", "RenyiEntropy", "None"]

    gd = GenericDialog("CA2 rabies starting-population setup")
    gd.addMessage("Red and blue channels are DISPLAY ONLY. Only the green channel is quantified.")
    gd.addStringField("Analysis_label:", "CA2_Starting_Population", 28)
    gd.addNumericField("Green_virus_channel_number:", 2, 0)
    gd.addNumericField("Red_CA2_marker_channel_number:", 1, 0)
    gd.addNumericField("Blue_Hoechst_channel_number:", 3, 0)
    gd.addNumericField("Expected_images_per_mouse:", 3, 0)
    gd.addMessage("Particle size units for this dataset: " + area_unit)
    gd.addNumericField("Minimum_cell_area_(%s):" % area_unit, 25.0, 2)
    gd.addNumericField("Maximum_cell_area_(%s):" % area_unit, 5000.0, 2)
    gd.addNumericField("Minimum_circularity_(0-1):", 0.00, 2)
    gd.addNumericField("Maximum_circularity_(0-1):", 1.00, 2)
    gd.addChoice("Auto-threshold_starting_method:", methods, "Otsu")
    gd.addNumericField("Rolling_ball_radius_pixels_(0=off):", 0.0, 1)
    gd.addNumericField("Median_filter_radius_pixels_(0=off):", 0.0, 1)
    gd.addNumericField("Starting_soma_cleanup_passes_(0-5):", 1, 0)
    gd.addCheckbox("Watershed_touching_objects", True)
    gd.addMessage(
        "Starting defaults were estimated from the representative cropped dataset.\n"
        "Verify the detected-cell overlay on the first full image before proceeding."
    )
    gd.showDialog()
    if gd.wasCanceled():
        return

    analysis_label = gd.getNextString()
    green_ch = int(gd.getNextNumber())
    red_ch = int(gd.getNextNumber())
    blue_ch = int(gd.getNextNumber())
    expected_per_mouse = int(gd.getNextNumber())
    min_area = float(gd.getNextNumber())
    max_area = float(gd.getNextNumber())
    min_circ = float(gd.getNextNumber())
    max_circ = float(gd.getNextNumber())
    auto_method = gd.getNextChoice()
    rolling_ball = float(gd.getNextNumber())
    median_radius = float(gd.getNextNumber())
    cleanup_passes = int(round(gd.getNextNumber()))
    watershed = gd.getNextBoolean()

    if expected_per_mouse < 1:
        IJ.showMessage("Invalid setting", "Expected images per mouse must be at least 1.")
        return
    if min_area < 0 or max_area <= min_area:
        IJ.showMessage("Invalid setting", "Cell-area limits are invalid.")
        return
    if min_circ < 0 or max_circ > 1 or max_circ < min_circ:
        IJ.showMessage("Invalid setting", "Circularity must satisfy 0 <= minimum <= maximum <= 1.")
        return
    if cleanup_passes < 0 or cleanup_passes > 5:
        IJ.showMessage("Invalid setting", "Starting soma-cleanup passes must be between 0 and 5.")
        return
    for c, label in [(green_ch, "green"), (red_ch, "red"), (blue_ch, "blue")]:
        if c < 1 or c > n_channels:
            IJ.showMessage("Channel error", "%s channel %d is outside 1-%d." % (label.capitalize(), c, n_channels))
            return
    if len(set([green_ch, red_ch, blue_ch])) != 3:
        IJ.showMessage("Channel error", "Green, red, and blue channel numbers must be different.")
        return

    try:
        mouse_folders = validate_input(parent, expected_per_mouse)
    except Exception as e:
        IJ.showMessage("Input validation failed", str(e))
        return

    sd = SaveDialog("Name the detailed output CSV", "Rabies_Starting_Population", ".csv")
    out_dir = sd.getDirectory()
    out_name = sd.getFileName()
    if out_dir is None or out_name is None:
        return
    if not out_name.lower().endswith('.csv'):
        out_name += '.csv'

    detailed_csv = os.path.join(out_dir, out_name)
    base = os.path.splitext(out_name)[0]
    summary_csv = os.path.join(out_dir, base + "_Summary.csv")
    qc_dir = os.path.join(out_dir, base + "_QC")
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)

    detailed_rows = []
    summary_rows = []
    threshold_low = None
    threshold_high = None
    threshold_locked = False
    cleanup_locked = False
    total_images = len(mouse_folders) * expected_per_mouse
    image_counter = 0
    aborted = False

    detail_header = [
        "Mouse", "Image", "Source_Path", "Analysis_Label",
        "CA2_ROI_Area", "Area_Unit", "Green_Positive_Cell_Count",
        "Threshold_Low", "Threshold_High",
        "Min_Cell_Area", "Max_Cell_Area",
        "Min_Circularity", "Max_Circularity",
        "Rolling_Ball_Pixels", "Median_Radius_Pixels", "Soma_Cleanup_Passes", "Watershed",
        "CA2_ROI_File", "Detected_Cell_ROIs_File", "QC_Overlay_File"
    ]
    summary_header = [
        "Mouse", "Analysis_Label", "Image_1_Count", "Image_2_Count", "Image_3_Count",
        "Summed_Starting_Population", "N_Images"
    ]

    for mouse_name, mouse_path in mouse_folders:
        mouse_counts = []
        tiffs = list_tiffs(mouse_path)

        mouse_qc_dir = os.path.join(qc_dir, safe_name(mouse_name))
        if not os.path.exists(mouse_qc_dir):
            os.makedirs(mouse_qc_dir)

        for image_name in tiffs:
            image_counter += 1
            path = os.path.join(mouse_path, image_name)
            IJ.showStatus("Processing %d/%d: %s / %s" % (image_counter, total_images, mouse_name, image_name))

            src = None
            comp = None
            green = None
            work = None
            mask = None

            try:
                src = IJ.openImage(path)
                if src is None:
                    raise RuntimeError("Could not open image: " + path)
                if src.getNChannels() < 3:
                    raise RuntimeError("Image has fewer than 3 channels: " + image_name)
                if src.getNSlices() != 1 or src.getNFrames() != 1:
                    raise RuntimeError("Image is not a single-plane single-timepoint image: " + image_name)

                comp = make_display_composite(
                    src, red_ch, green_ch, blue_ch,
                    "%s | %s | draw/review" % (mouse_name, image_name)
                )
                comp.show()

                ca2_roi = prompt_for_ca2_roi(comp, image_name, mouse_name)
                comp.setRoi(ca2_roi)
                roi_stats = comp.getStatistics(Measurements.AREA)
                roi_area = roi_stats.area

                green = channel_duplicate(src, green_ch)
                rotate_clockwise_90(green)
                green.setTitle("%s | %s | GREEN virus (rotated 90 deg CW)" % (mouse_name, image_name))
                work = prepare_detection_image(green, rolling_ball, median_radius)

                is_first_global = not threshold_locked

                while True:
                    if not threshold_locked:
                        threshold_low, threshold_high = calibrate_threshold(
                            work, ca2_roi, auto_method, image_name
                        )

                    if not cleanup_locked:
                        chosen_cleanup = prompt_cleanup_passes(cleanup_passes, image_name)
                        if chosen_cleanup is None:
                            aborted = True
                            break
                        cleanup_passes = chosen_cleanup

                    mask = make_binary_mask(work, threshold_low, threshold_high, cleanup_passes, watershed)
                    cell_rois = detect_cells(mask, ca2_roi, min_area, max_area, min_circ, max_circ)
                    count = len(cell_rois)

                    apply_review_overlay(comp, ca2_roi, cell_rois, count)
                    comp.show()
                    if comp.getWindow() is not None:
                        comp.getWindow().toFront()

                    decision = review_count(comp, image_name, mouse_name, count, is_first_global)

                    if decision == "accept":
                        if is_first_global:
                            threshold_locked = True
                            cleanup_locked = True
                        break
                    elif decision == "retry" and is_first_global:
                        comp.setOverlay(None)
                        close_imp(mask)
                        mask = None
                        # Return to threshold and cleanup calibration on the same first image.
                        continue
                    else:
                        aborted = True
                        break

                if aborted:
                    break

                stem = safe_name(image_name)
                roi_path = os.path.join(mouse_qc_dir, stem + "__CA2.roi")
                cells_path = os.path.join(mouse_qc_dir, stem + "__DetectedCells.zip")
                qc_png = os.path.join(mouse_qc_dir, stem + "__QC.png")

                save_roi_file(ca2_roi, roi_path)
                save_cell_rois(cell_rois, cells_path)
                save_qc_png(comp, qc_png)

                detailed_rows.append([
                    mouse_name, image_name, path, analysis_label,
                    "%.6f" % roi_area, area_unit, count,
                    "%.6f" % threshold_low, "%.6f" % threshold_high,
                    min_area, max_area, min_circ, max_circ,
                    rolling_ball, median_radius, cleanup_passes, str(bool(watershed)),
                    roi_path, cells_path if count > 0 else os.path.splitext(cells_path)[0] + "_NO_CELLS.txt",
                    qc_png
                ])
                mouse_counts.append(count)

                # Write partial detailed results after every accepted image for crash protection.
                csv_write(detailed_csv, detailed_rows, detail_header)

            except Exception as e:
                IJ.showMessage("Error", "Problem while processing:\n%s\n\n%s" % (path, str(e)))
                aborted = True

            finally:
                close_imp(mask)
                close_imp(work)
                close_imp(green)
                close_imp(comp)
                close_imp(src)

            if aborted:
                break

        if aborted:
            break

        if len(mouse_counts) == expected_per_mouse:
            # The current study design uses 3 images per mouse. Keep the CSV columns
            # explicit for 3 counts, while still recording N_Images.
            c1 = mouse_counts[0] if len(mouse_counts) > 0 else ""
            c2 = mouse_counts[1] if len(mouse_counts) > 1 else ""
            c3 = mouse_counts[2] if len(mouse_counts) > 2 else ""
            summary_rows.append([
                mouse_name, analysis_label, c1, c2, c3,
                sum(mouse_counts), len(mouse_counts)
            ])
            csv_write(summary_csv, summary_rows, summary_header)

    # Final writes.
    if len(detailed_rows) > 0:
        csv_write(detailed_csv, detailed_rows, detail_header)
    if len(summary_rows) > 0:
        csv_write(summary_csv, summary_rows, summary_header)

    if aborted:
        IJ.showMessage(
            "Run stopped",
            "The run was stopped before completion.\n\n"
            "Accepted image-level results were preserved in:\n%s\n\n"
            "Completed mouse summaries, if any, were preserved in:\n%s"
            % (detailed_csv, summary_csv)
        )
    else:
        IJ.showMessage(
            "Finished",
            "Rabies starting-population quantification is complete.\n\n"
            "Detailed results:\n%s\n\n"
            "Mouse-level summed results:\n%s\n\n"
            "Saved ROIs and QC overlays:\n%s\n\n"
            "Global green threshold used: %.3f to %.3f\n"
            "Soma-cleanup passes used: %d"
            % (detailed_csv, summary_csv, qc_dir, threshold_low, threshold_high, cleanup_passes)
        )
    IJ.showStatus("Rabies starting-population counter finished")


main()
