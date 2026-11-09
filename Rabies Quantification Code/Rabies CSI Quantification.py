"""
Rabies Regional Input Counter + Ipsilateral/Contralateral Connection Strength Index (CSI)
Fiji / ImageJ Jython

Purpose
-------
Analyze one anatomical region per run. The selected parent directory must contain
mouse folders M1, M2, M3, and M4, with exactly 3 TIFF images in each folder.

For each image:
1) Build a red/green/blue composite and rotate it 90 degrees clockwise.
2) User draws the IPSILATERAL ROI for the region.
3) User draws the CONTRALATERAL ROI for the same region.
4) Count green rabies-positive objects with the SAME fixed segmentation settings
   used for the finalized starting-population analysis.
5) Show both sides together for visual QC.
6) After acceptance, save ONE combined ipsilateral/contralateral QC JPEG,
   downsampled to 50% width and 50% height, then record the numerical results.
   All detection/counting is still performed at full resolution.

For each mouse, the 3 image counts are summed separately for each side:

    Ipsi CSI   = summed ipsilateral input cells / starting population
    Contra CSI = summed contralateral input cells / starting population

Hard-coded starting populations:
    M1 = 99
    M2 = 227
    M3 = 9
    M4 = 167

Fixed detection settings copied from the finalized starting-population workflow:
    Green channel: 2
    Red channel: 1
    Blue channel: 3
    Green threshold: 14135 to 65535
    Minimum particle area: 25
    Maximum particle area: 5000
    Circularity: 0.00 to 1.00
    Rolling-ball background subtraction: OFF
    Median filter: OFF
    Soma cleanup: 1 erosion followed by 1 dilation
    Watershed: ON
    Rotation: 90 degrees clockwise

Outputs
-------
Saved automatically in the selected parent directory:
1) Rabies_CSI_<region>_IpsiContra_Detailed.csv
2) Rabies_CSI_<region>_IpsiContra_Summary.csv
3) Rabies_CSI_<region>_IpsiContra_QC/ containing ONE combined QC JPEG per accepted image:
       *__IpsiContra_QC.jpg

QC images are saved at 50% linear resolution and JPEG quality 90 to speed writing.
No ROI files or detected-cell ROI ZIPs are saved.

Recommended Fiji usage
----------------------
File > New > Script, choose language "Python", open/paste this file, then Run.
This script is intended for Fiji's Jython/Python interpreter (ImageJ 1.x API).
"""

from ij import IJ, Prefs, ImagePlus
from ij.gui import GenericDialog, WaitForUserDialog, YesNoCancelDialog, Overlay, TextRoi
from ij.io import DirectoryChooser, FileSaver
from ij.plugin import Duplicator, RGBStackMerge, JpegWriter
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import Measurements, ResultsTable
from ij.process import ImageProcessor
from java.awt import Color, Font
import os
import re
import csv


# --------------------------- Fixed study settings -------------------------- #

STARTING_POPULATIONS = {
    "M1": 99,
    "M2": 227,
    "M3": 9,
    "M4": 167,
}

EXPECTED_MICE = ["M1", "M2", "M3", "M4"]
EXPECTED_IMAGES_PER_MOUSE = 3

RED_CHANNEL = 1
GREEN_CHANNEL = 2
BLUE_CHANNEL = 3

THRESHOLD_LOW = 14135.0
THRESHOLD_HIGH = 65535.0

MIN_AREA = 25.0
MAX_AREA = 5000.0
MIN_CIRCULARITY = 0.00
MAX_CIRCULARITY = 1.00

ROLLING_BALL = 0.0
MEDIAN_RADIUS = 0.0
SOMA_CLEANUP_PASSES = 1
WATERSHED = True

# QC export settings only. These do NOT affect detection or counts.
QC_SCALE = 0.50
QC_JPEG_QUALITY = 90


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
    base = os.path.splitext(os.path.basename(name))[0]
    return re.sub(r'[^A-Za-z0-9._-]+', '_', base)


def list_tiffs(folder):
    names = []
    for name in os.listdir(folder):
        low = name.lower()
        if low.endswith('.tif') or low.endswith('.tiff'):
            path = os.path.join(folder, name)
            if os.path.isfile(path):
                names.append(name)
    names.sort(key=natural_key)
    return names


def get_mouse_folders(parent):
    folders = []
    for name in os.listdir(parent):
        path = os.path.join(parent, name)
        if os.path.isdir(path) and re.match(r'^M\d+$', name, re.IGNORECASE):
            folders.append((name.upper(), path))
    folders.sort(key=lambda x: natural_key(x[0]))
    return folders


def validate_input(parent):
    mouse_folders = get_mouse_folders(parent)
    found_names = [x[0] for x in mouse_folders]

    missing = [m for m in EXPECTED_MICE if m not in found_names]
    unsupported = [m for m in found_names if m not in STARTING_POPULATIONS]

    problems = []
    if len(missing) > 0:
        problems.append("Missing mouse folder(s): %s" % ", ".join(missing))
    if len(unsupported) > 0:
        problems.append(
            "Mouse folder(s) without a hard-coded starting population: %s"
            % ", ".join(unsupported)
        )

    for mouse_name, mouse_path in mouse_folders:
        if mouse_name in STARTING_POPULATIONS:
            n = len(list_tiffs(mouse_path))
            if n != EXPECTED_IMAGES_PER_MOUSE:
                problems.append(
                    "%s: found %d TIFF(s), expected exactly %d"
                    % (mouse_name, n, EXPECTED_IMAGES_PER_MOUSE)
                )

    if len(problems) > 0:
        raise ValueError("Input validation failed:\n\n" + "\n".join(problems))

    lookup = dict(mouse_folders)
    return [(m, lookup[m]) for m in EXPECTED_MICE]


def close_imp(imp):
    if imp is not None:
        try:
            imp.changes = False
            imp.close()
        except:
            pass


def channel_duplicate(imp, channel):
    n_channels = max(1, imp.getNChannels())
    if channel < 1 or channel > n_channels:
        raise ValueError(
            "Requested channel %d, but image has %d channel(s)."
            % (channel, n_channels)
        )
    return Duplicator().run(imp, channel, channel, 1, 1, 1, 1)


def rotate_clockwise_90(imp):
    IJ.run(imp, "Rotate 90 Degrees Right", "")
    return imp


def make_display_composite(imp, title):
    r = channel_duplicate(imp, RED_CHANNEL)
    g = channel_duplicate(imp, GREEN_CHANNEL)
    b = channel_duplicate(imp, BLUE_CHANNEL)

    rotate_clockwise_90(r)
    rotate_clockwise_90(g)
    rotate_clockwise_90(b)

    comp = RGBStackMerge.mergeChannels([r, g, b], True)
    comp.setTitle(title)

    close_imp(r)
    close_imp(g)
    close_imp(b)
    return comp


def make_single_roi_overlay(roi, side_label):
    overlay = Overlay()
    r = roi.clone()
    if side_label == "Ipsilateral":
        r.setStrokeColor(Color.yellow)
    else:
        r.setStrokeColor(Color.cyan)
    r.setStrokeWidth(2.0)
    r.setName(side_label)
    overlay.add(r)
    return overlay


def prompt_for_side_roi(comp, image_name, mouse_name, region_label, side_label, show_existing_roi=None):
    while True:
        comp.show()
        if show_existing_roi is not None:
            comp.setOverlay(make_single_roi_overlay(show_existing_roi, "Ipsilateral"))
        else:
            comp.setOverlay(None)

        try:
            comp.deleteRoi()
        except:
            comp.setRoi(None)

        IJ.setTool("freehand")
        msg = (
            "Region: %s\n"
            "Mouse: %s\n"
            "Image: %s\n"
            "Side: %s\n\n"
            "Draw the %s ROI for %s on the composite.\n"
            "Only the GREEN channel will be counted.\n\n"
            "Zoom and pan as needed. When the ROI is complete, click OK here."
        ) % (
            region_label,
            mouse_name,
            image_name,
            side_label,
            side_label.lower(),
            region_label
        )
        WaitForUserDialog("Draw %s ROI" % side_label, msg).show()

        roi = comp.getRoi()
        if roi is not None and roi.isArea():
            return roi.clone()

        IJ.showMessage(
            "ROI required",
            "Please draw a closed area ROI for the %s side before continuing."
            % side_label.lower()
        )


def prepare_detection_image(green_imp):
    work = green_imp.duplicate()
    work.setTitle(green_imp.getTitle() + " - detection working copy")

    if ROLLING_BALL > 0:
        IJ.run(work, "Subtract Background...", "rolling=%s" % str(ROLLING_BALL))
    if MEDIAN_RADIUS > 0:
        IJ.run(work, "Median...", "radius=%s" % str(MEDIAN_RADIUS))

    return work


def make_binary_mask(work):
    mask = work.duplicate()
    mask.deleteRoi()
    ip = mask.getProcessor()
    ip.setThreshold(THRESHOLD_LOW, THRESHOLD_HIGH, ImageProcessor.NO_LUT_UPDATE)

    old_black = Prefs.blackBackground
    Prefs.blackBackground = True

    try:
        IJ.run(mask, "Convert to Mask", "")

        for i in range(int(SOMA_CLEANUP_PASSES)):
            IJ.run(mask, "Erode", "")
        for i in range(int(SOMA_CLEANUP_PASSES)):
            IJ.run(mask, "Dilate", "")

        if WATERSHED:
            IJ.run(mask, "Watershed", "")
    finally:
        Prefs.blackBackground = old_black

    return mask


def detect_all_particles(mask):
    """Run particle detection once and retain each particle plus its raw-pixel centroid."""
    rt = ResultsTable()
    rm = RoiManager(True)
    ParticleAnalyzer.setRoiManager(rm)

    options = (
        ParticleAnalyzer.ADD_TO_MANAGER |
        ParticleAnalyzer.CLEAR_WORKSHEET |
        ParticleAnalyzer.SHOW_NONE
    )
    measurements = Measurements.AREA | Measurements.CENTROID

    pa = ParticleAnalyzer(
        options,
        measurements,
        rt,
        MIN_AREA,
        MAX_AREA,
        MIN_CIRCULARITY,
        MAX_CIRCULARITY
    )

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
    particles = []

    n = min(len(all_rois), rt.getCounter())
    for i in range(n):
        try:
            x_cal = rt.getValue("X", i)
            y_cal = rt.getValue("Y", i)
            x_raw = int(round(cal.getRawX(x_cal)))
            y_raw = int(round(cal.getRawY(y_cal)))
        except:
            b = all_rois[i].getBounds()
            x_raw = int(round(b.x + b.width / 2.0))
            y_raw = int(round(b.y + b.height / 2.0))

        particles.append((all_rois[i].clone(), x_raw, y_raw))

    try:
        rm.reset()
    except:
        pass
    try:
        rm.close()
    except:
        pass

    return particles


def particles_in_roi(particles, roi):
    kept = []
    for particle_roi, x_raw, y_raw in particles:
        if roi.contains(x_raw, y_raw):
            kept.append(particle_roi.clone())
    return kept


def make_review_overlay(ipsi_roi, contra_roi, ipsi_cells, contra_cells,
                        ipsi_count, contra_count, region_label, mouse_name):
    overlay = Overlay()

    ipsi = ipsi_roi.clone()
    ipsi.setStrokeColor(Color.yellow)
    ipsi.setStrokeWidth(2.0)
    ipsi.setName("Ipsilateral_%s" % region_label)
    overlay.add(ipsi)

    contra = contra_roi.clone()
    contra.setStrokeColor(Color.cyan)
    contra.setStrokeWidth(2.0)
    contra.setName("Contralateral_%s" % region_label)
    overlay.add(contra)

    for i, roi in enumerate(ipsi_cells):
        r = roi.clone()
        r.setStrokeColor(Color.white)
        r.setStrokeWidth(1.0)
        r.setName("Ipsi_Cell_%04d" % (i + 1))
        overlay.add(r)

    for i, roi in enumerate(contra_cells):
        r = roi.clone()
        r.setStrokeColor(Color.magenta)
        r.setStrokeWidth(1.0)
        r.setName("Contra_Cell_%04d" % (i + 1))
        overlay.add(r)

    label_text = "%s | %s | Ipsi: %d | Contra: %d" % (
        mouse_name,
        region_label,
        ipsi_count,
        contra_count
    )
    label = TextRoi(10, 10, label_text, Font("SansSerif", Font.BOLD, 16))
    label.setStrokeColor(Color.white)
    overlay.add(label)

    return overlay


def review_count(comp, image_name, mouse_name, region_label, ipsi_count, contra_count):
    text = (
        "Region: %s\n"
        "Mouse: %s\n"
        "Image: %s\n\n"
        "Ipsilateral green-positive input cells: %d\n"
        "Contralateral green-positive input cells: %d\n\n"
        "QC colors:\n"
        "  Yellow boundary = ipsilateral ROI\n"
        "  Cyan boundary = contralateral ROI\n"
        "  White outlines = ipsilateral detected cells\n"
        "  Magenta outlines = contralateral detected cells\n\n"
        "YES = accept, save one combined downsampled QC JPEG, and continue.\n"
        "NO = stop the run.\n"
        "CANCEL = stop the run."
    ) % (
        region_label,
        mouse_name,
        image_name,
        ipsi_count,
        contra_count
    )

    dlg = YesNoCancelDialog(comp.getWindow(), "Review ipsilateral/contralateral detections", text)
    if dlg.yesPressed():
        return True
    return False


def save_combined_qc_jpeg(comp, path):
    """
    Save the CURRENT combined ipsi/contra overlay as a downsampled JPEG.

    Important: this function is called only after all full-resolution analysis
    and counting are complete. Downsampling changes the QC file only and has
    no effect on ROI geometry, particle detection, or CSI.
    """
    flat = comp.flatten()
    flat.setTitle("QC full resolution")

    try:
        scale = float(QC_SCALE)
        if scale <= 0 or scale > 1:
            scale = 1.0

        new_w = max(1, int(round(flat.getWidth() * scale)))
        new_h = max(1, int(round(flat.getHeight() * scale)))

        if new_w != flat.getWidth() or new_h != flat.getHeight():
            # useAveraging=True gives cleaner downsampling for QC review.
            resized_ip = flat.getProcessor().resize(new_w, new_h, True)
            qc = ImagePlus("QC", resized_ip)
        else:
            qc = flat.duplicate()
            qc.setTitle("QC")

        try:
            JpegWriter.setQuality(int(QC_JPEG_QUALITY))
        except:
            # If an older ImageJ build lacks the setter, FileSaver still saves
            # a JPEG using ImageJ's current/default JPEG quality.
            pass

        ok = FileSaver(qc).saveAsJpeg(path)
        close_imp(qc)
        if not ok:
            raise RuntimeError("Could not save QC JPEG: " + path)
    finally:
        close_imp(flat)


def csv_write(path, rows, header):
    f = open(path, 'wb')
    try:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    finally:
        f.close()


def unique_output_paths(parent, region_label):
    region_safe = safe_name(region_label)
    if len(region_safe) == 0:
        region_safe = "Region"

    base = "Rabies_CSI_%s_IpsiContra" % region_safe
    candidate = base
    index = 2

    while True:
        detailed = os.path.join(parent, candidate + "_Detailed.csv")
        summary = os.path.join(parent, candidate + "_Summary.csv")
        qc_dir = os.path.join(parent, candidate + "_QC")

        if (
            not os.path.exists(detailed) and
            not os.path.exists(summary) and
            not os.path.exists(qc_dir)
        ):
            return detailed, summary, qc_dir

        candidate = "%s_%d" % (base, index)
        index += 1


# -------------------------------- Main ------------------------------------ #

def main():
    IJ.showStatus("Rabies ipsi/contra CSI counter: setup")

    dc = DirectoryChooser("Choose region folder containing M1, M2, M3, and M4")
    parent = dc.getDirectory()
    if parent is None:
        return

    gd = GenericDialog("Rabies regional ipsilateral/contralateral CSI analysis")
    gd.addStringField("Region_label:", "Region", 24)
    gd.addMessage(
        "For every image you will draw TWO ROIs:\n"
        "1) Ipsilateral ROI\n"
        "2) Contralateral ROI\n\n"
        "Fixed detection settings:\n"
        "Green threshold: 14135-65535\n"
        "Area: 25-5000\n"
        "Circularity: 0.00-1.00\n"
        "Soma cleanup: 1 pass\n"
        "Watershed: ON\n"
        "Channels R/G/B: 1/2/3\n\n"
        "Starting populations: M1=99, M2=227, M3=9, M4=167\n"
        "One combined QC JPEG per accepted image will be saved automatically.\n"
        "QC export: 50% width/height, JPEG quality 90. Analysis remains full resolution."
    )
    gd.showDialog()
    if gd.wasCanceled():
        return

    region_label = gd.getNextString().strip()
    if len(region_label) == 0:
        IJ.showMessage("Region label required", "Please enter a region label.")
        return

    try:
        mouse_folders = validate_input(parent)
    except Exception as e:
        IJ.showMessage("Input validation failed", str(e))
        return

    first_mouse_name, first_mouse_path = mouse_folders[0]
    first_tiffs = list_tiffs(first_mouse_path)
    probe_path = os.path.join(first_mouse_path, first_tiffs[0])
    probe = IJ.openImage(probe_path)
    if probe is None:
        IJ.showMessage("Open error", "Could not open first TIFF:\n" + probe_path)
        return

    n_channels = probe.getNChannels()
    bit_depth = probe.getBitDepth()
    cal = probe.getCalibration()
    area_unit = cal.getUnit() + "^2" if cal.scaled() else "pixel^2"
    close_imp(probe)

    if n_channels < 3:
        IJ.showMessage(
            "Channel error",
            "The first image has %d channel(s); at least 3 are required." % n_channels
        )
        return

    if bit_depth != 16:
        IJ.showMessage(
            "Bit-depth mismatch",
            "The fixed threshold 14135-65535 was established on 16-bit images, "
            "but the first image reports %d-bit data. The run has been stopped."
            % bit_depth
        )
        return

    detailed_csv, summary_csv, qc_dir = unique_output_paths(parent, region_label)
    os.makedirs(qc_dir)

    detail_header = [
        "Mouse",
        "Region",
        "Image",
        "Source_Path",
        "Ipsi_ROI_Area",
        "Contra_ROI_Area",
        "Area_Unit",
        "Ipsi_Green_Positive_Input_Count",
        "Contra_Green_Positive_Input_Count",
        "Starting_Population",
        "Threshold_Low",
        "Threshold_High",
        "Min_Cell_Area",
        "Max_Cell_Area",
        "Min_Circularity",
        "Max_Circularity",
        "Soma_Cleanup_Passes",
        "Watershed",
        "QC_Overlay_File"
    ]

    summary_header = [
        "Mouse",
        "Region",
        "Starting_Population",
        "Ipsi_Image_1_Count",
        "Ipsi_Image_2_Count",
        "Ipsi_Image_3_Count",
        "Ipsi_Summed_Input_Cells",
        "Ipsi_CSI",
        "Contra_Image_1_Count",
        "Contra_Image_2_Count",
        "Contra_Image_3_Count",
        "Contra_Summed_Input_Cells",
        "Contra_CSI",
        "N_Images"
    ]

    detailed_rows = []
    summary_rows = []
    total_images = len(mouse_folders) * EXPECTED_IMAGES_PER_MOUSE
    image_counter = 0
    aborted = False

    for mouse_name, mouse_path in mouse_folders:
        ipsi_counts = []
        contra_counts = []
        starting_population = STARTING_POPULATIONS[mouse_name]
        tiffs = list_tiffs(mouse_path)

        mouse_qc_dir = os.path.join(qc_dir, mouse_name)
        if not os.path.exists(mouse_qc_dir):
            os.makedirs(mouse_qc_dir)

        for image_name in tiffs:
            image_counter += 1
            source_path = os.path.join(mouse_path, image_name)
            IJ.showStatus(
                "Processing %d/%d: %s / %s / %s"
                % (image_counter, total_images, region_label, mouse_name, image_name)
            )

            src = None
            comp = None
            green = None
            work = None
            mask = None

            try:
                src = IJ.openImage(source_path)
                if src is None:
                    raise RuntimeError("Could not open image: " + source_path)

                if src.getNChannels() < 3:
                    raise RuntimeError("Image has fewer than 3 channels: " + image_name)

                if src.getNSlices() != 1 or src.getNFrames() != 1:
                    raise RuntimeError(
                        "Image is not a single-plane, single-timepoint image: " + image_name
                    )

                if src.getBitDepth() != 16:
                    raise RuntimeError(
                        "Image is not 16-bit, so the fixed threshold is not valid: " + image_name
                    )

                comp = make_display_composite(
                    src,
                    "%s | %s | %s" % (region_label, mouse_name, image_name)
                )
                comp.show()

                ipsi_roi = prompt_for_side_roi(
                    comp,
                    image_name,
                    mouse_name,
                    region_label,
                    "Ipsilateral"
                )

                contra_roi = prompt_for_side_roi(
                    comp,
                    image_name,
                    mouse_name,
                    region_label,
                    "Contralateral",
                    show_existing_roi=ipsi_roi
                )

                comp.setOverlay(None)

                comp.setRoi(ipsi_roi)
                ipsi_stats = comp.getStatistics(Measurements.AREA)
                ipsi_area = ipsi_stats.area

                comp.setRoi(contra_roi)
                contra_stats = comp.getStatistics(Measurements.AREA)
                contra_area = contra_stats.area

                try:
                    comp.deleteRoi()
                except:
                    comp.setRoi(None)

                green = channel_duplicate(src, GREEN_CHANNEL)
                rotate_clockwise_90(green)
                green.setTitle(
                    "%s | %s | %s | GREEN virus"
                    % (region_label, mouse_name, image_name)
                )

                work = prepare_detection_image(green)
                mask = make_binary_mask(work)
                particles = detect_all_particles(mask)

                ipsi_cells = particles_in_roi(particles, ipsi_roi)
                contra_cells = particles_in_roi(particles, contra_roi)
                ipsi_count = len(ipsi_cells)
                contra_count = len(contra_cells)

                review_overlay = make_review_overlay(
                    ipsi_roi,
                    contra_roi,
                    ipsi_cells,
                    contra_cells,
                    ipsi_count,
                    contra_count,
                    region_label,
                    mouse_name
                )
                comp.setOverlay(review_overlay)
                comp.updateAndDraw()
                comp.show()
                if comp.getWindow() is not None:
                    comp.getWindow().toFront()

                if not review_count(
                    comp,
                    image_name,
                    mouse_name,
                    region_label,
                    ipsi_count,
                    contra_count
                ):
                    aborted = True
                    break

                stem = safe_name(image_name)
                qc_jpg = os.path.join(mouse_qc_dir, stem + "__IpsiContra_QC.jpg")

                # The review overlay already contains BOTH ipsilateral and
                # contralateral ROIs/detections. Save that single combined view.
                # Export downsampling happens only here, after full-resolution counting.
                comp.setOverlay(review_overlay)
                comp.updateAndDraw()
                save_combined_qc_jpeg(comp, qc_jpg)

                detailed_rows.append([
                    mouse_name,
                    region_label,
                    image_name,
                    source_path,
                    "%.6f" % ipsi_area,
                    "%.6f" % contra_area,
                    area_unit,
                    ipsi_count,
                    contra_count,
                    starting_population,
                    "%.0f" % THRESHOLD_LOW,
                    "%.0f" % THRESHOLD_HIGH,
                    MIN_AREA,
                    MAX_AREA,
                    MIN_CIRCULARITY,
                    MAX_CIRCULARITY,
                    SOMA_CLEANUP_PASSES,
                    str(bool(WATERSHED)),
                    qc_jpg
                ])

                ipsi_counts.append(ipsi_count)
                contra_counts.append(contra_count)

                # For speed, do not rewrite the CSV after every image.
                # Results are written after each completed mouse and again at the end.

            except Exception as e:
                IJ.showMessage(
                    "Error",
                    "Problem while processing:\n%s\n\n%s"
                    % (source_path, str(e))
                )
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

        if (
            len(ipsi_counts) == EXPECTED_IMAGES_PER_MOUSE and
            len(contra_counts) == EXPECTED_IMAGES_PER_MOUSE
        ):
            ipsi_sum = sum(ipsi_counts)
            contra_sum = sum(contra_counts)
            ipsi_csi = float(ipsi_sum) / float(starting_population)
            contra_csi = float(contra_sum) / float(starting_population)

            summary_rows.append([
                mouse_name,
                region_label,
                starting_population,
                ipsi_counts[0],
                ipsi_counts[1],
                ipsi_counts[2],
                ipsi_sum,
                "%.6f" % ipsi_csi,
                contra_counts[0],
                contra_counts[1],
                contra_counts[2],
                contra_sum,
                "%.6f" % contra_csi,
                len(ipsi_counts)
            ])

            csv_write(summary_csv, summary_rows, summary_header)
            csv_write(detailed_csv, detailed_rows, detail_header)

    if len(detailed_rows) > 0:
        csv_write(detailed_csv, detailed_rows, detail_header)
    if len(summary_rows) > 0:
        csv_write(summary_csv, summary_rows, summary_header)

    if aborted:
        IJ.showMessage(
            "Run stopped",
            "The run stopped before completion.\n\n"
            "Accepted image-level results were preserved in:\n%s\n\n"
            "Completed mouse CSI summaries, if any, were preserved in:\n%s\n\n"
            "Saved combined QC JPEGs are in:\n%s"
            % (detailed_csv, summary_csv, qc_dir)
        )
    else:
        IJ.showMessage(
            "Finished",
            "Rabies regional ipsilateral/contralateral CSI analysis is complete for: %s\n\n"
            "Detailed image counts:\n%s\n\n"
            "Mouse-level CSI summary:\n%s\n\n"
            "Combined QC JPEGs (50% size, quality 90):\n%s\n\n"
            "CSI = summed side-specific input cells / mouse starting population"
            % (region_label, detailed_csv, summary_csv, qc_dir)
        )

    IJ.showStatus("Rabies ipsi/contra CSI counter finished")


main()
