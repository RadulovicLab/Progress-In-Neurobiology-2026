#@ ImageJ ij
# Fiji / ImageJ Jython script for semi-automated CA2 projection quantification
#
# Workflow per run:
#   - Select multiple multichannel TIFF images.
#   - Enter one anatomical region label for the whole run (e.g. Dorsal_CA1).
#   - For each image, draw ipsilateral and contralateral anatomical ROIs.
#   - Place fixed-size ipsilateral and contralateral background ROIs.
#   - Measure the GREEN viral-expression channel only.
#   - Save raw mean intensity, background mean, normalized intensity, ROI area,
#     and an ROI .zip for each image.
#
# Normalization:
#   normalized fluorescence = mean green fluorescence / local background mean
#
# IMPORTANT:
#   This script expects a true multichannel TIFF/hyperstack with one 2D plane
#   per channel (not a flattened RGB image and not a Z-stack/time series).

from ij import IJ, CompositeImage
from ij.gui import Roi, WaitForUserDialog
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import LUT
from java.awt import Color
from java.io import File, BufferedWriter, FileWriter
from javax.swing import JFileChooser, JOptionPane
import os


def safe_name(text):
    """Return a filename-safe version of a label."""
    bad = '<>:"/\\|?*'
    out = str(text)
    for ch in bad:
        out = out.replace(ch, '_')
    return out.strip() or 'unnamed'


def csv_escape(value):
    """Minimal CSV escaping that is reliable in Jython/Java file writing."""
    if value is None:
        s = ''
    else:
        s = str(value)
    if ',' in s or '"' in s or '\n' in s or '\r' in s:
        s = '"' + s.replace('"', '""') + '"'
    return s


def write_csv_row(writer, values):
    writer.write(','.join([csv_escape(v) for v in values]))
    writer.newLine()
    writer.flush()


def choose_images():
    chooser = JFileChooser()
    chooser.setDialogTitle('Select TIFF images for this anatomical region')
    chooser.setMultiSelectionEnabled(True)
    result = chooser.showOpenDialog(None)
    if result != JFileChooser.APPROVE_OPTION:
        return []
    selected = list(chooser.getSelectedFiles())
    selected.sort(key=lambda f: f.getName().lower())
    return selected


def choose_output_csv(default_name):
    chooser = JFileChooser()
    chooser.setDialogTitle('Name the quantification output CSV')
    chooser.setSelectedFile(File(default_name))
    result = chooser.showSaveDialog(None)
    if result != JFileChooser.APPROVE_OPTION:
        return None

    path = chooser.getSelectedFile().getAbsolutePath()
    if not path.lower().endswith('.csv'):
        path += '.csv'

    if File(path).exists():
        answer = JOptionPane.showConfirmDialog(
            None,
            'The file already exists:\n\n%s\n\nOverwrite it?' % path,
            'Confirm overwrite',
            JOptionPane.YES_NO_OPTION
        )
        if answer != JOptionPane.YES_OPTION:
            return None
    return path


def get_run_settings():
    from ij.gui import GenericDialog

    gd = GenericDialog('CA2 projection quantification settings')
    gd.addStringField('Anatomical_region_label:', 'Dorsal_CA1', 25)
    gd.addNumericField('Green_projection_channel_number:', 1, 0)
    gd.addNumericField('Blue_Hoechst_channel_number:', 2, 0)
    gd.addNumericField('Background_ROI_width_pixels:', 300, 0)
    gd.addNumericField('Background_ROI_height_pixels:', 300, 0)
    gd.addMessage(
        'One run should contain only one anatomical category.\n'
        'Example: all images in this run are Dorsal_CA1.\n\n'
        'The blue+green composite is used for ROI drawing.\n'
        'Only the green channel is quantified.'
    )
    gd.showDialog()
    if gd.wasCanceled():
        return None

    region = gd.getNextString().strip()
    green_ch = int(gd.getNextNumber())
    blue_ch = int(gd.getNextNumber())
    bg_w = int(gd.getNextNumber())
    bg_h = int(gd.getNextNumber())

    if region == '':
        region = 'Unlabeled_region'
    if green_ch < 1 or blue_ch < 1:
        raise ValueError('Channel numbers must be >= 1.')
    if green_ch == blue_ch:
        raise ValueError('Green and blue channel numbers must be different.')
    if bg_w < 1 or bg_h < 1:
        raise ValueError('Background ROI dimensions must be positive.')

    return region, green_ch, blue_ch, bg_w, bg_h


def prepare_composite(imp, green_ch, blue_ch):
    n_channels = imp.getNChannels()

    if imp.getBitDepth() == 24:
        raise ValueError(
            'This image is a flattened RGB image. The script requires a true '
            'multichannel TIFF/hyperstack so the green channel can be measured '
            'independently.'
        )

    if n_channels < 2:
        raise ValueError(
            'Fiji reports only %d channel(s). The script requires a true '
            'multichannel TIFF/hyperstack.' % n_channels
        )

    if green_ch > n_channels or blue_ch > n_channels:
        raise ValueError(
            'This image has %d channels, but the selected channel numbers are '
            'green=%d and blue=%d.' % (n_channels, green_ch, blue_ch)
        )

    if imp.getNSlices() != 1 or imp.getNFrames() != 1:
        raise ValueError(
            'This script expects one 2D optical plane per file. Fiji reports '
            '%d Z-slices and %d time frames.' % (imp.getNSlices(), imp.getNFrames())
        )

    if isinstance(imp, CompositeImage):
        comp = imp
        comp.setMode(CompositeImage.COMPOSITE)
    else:
        comp = CompositeImage(imp, CompositeImage.COMPOSITE)

    # Assign display LUTs only. This does NOT change pixel values.
    try:
        comp.setChannelLut(LUT.createLutFromColor(Color.GREEN), green_ch)
        comp.setChannelLut(LUT.createLutFromColor(Color.BLUE), blue_ch)
    except:
        # The measurements are still valid if LUT assignment is unsupported.
        pass

    # If there are extra channels, display only the chosen green and blue channels.
    try:
        active = ['0'] * n_channels
        active[green_ch - 1] = '1'
        active[blue_ch - 1] = '1'
        comp.setActiveChannels(''.join(active))
    except:
        pass

    comp.setTitle(imp.getTitle() + ' - ROI drawing')
    comp.show()
    comp.updateAndDraw()
    return comp


def get_freehand_area_roi(comp, side_label):
    while True:
        comp.killRoi()
        try:
            IJ.setTool('freehand')
        except:
            pass

        dialog = WaitForUserDialog(
            'Draw %s ROI' % side_label,
            'Draw the %s anatomical ROI on the blue+green composite.\n\n'
            'Use the Hoechst signal to follow anatomical boundaries.\n'
            'The ROI may be any size; it should follow the anatomical region.\n\n'
            'When the ROI is complete, click OK.' % side_label
        )
        dialog.show()

        roi = comp.getRoi()
        if roi is not None:
            try:
                is_area = roi.isArea()
            except:
                is_area = True
            if is_area:
                saved = roi.clone()
                saved.setName('%s_Anatomical' % side_label)
                return saved

        JOptionPane.showMessageDialog(
            None,
            'No valid area ROI was detected. Please draw an area selection and try again.',
            'ROI required',
            JOptionPane.WARNING_MESSAGE
        )


def clamp(value, low, high):
    if value < low:
        return low
    if value > high:
        return high
    return value


def place_fixed_background_roi(comp, side_label, bg_w, bg_h):
    img_w = comp.getWidth()
    img_h = comp.getHeight()

    if bg_w > img_w or bg_h > img_h:
        raise ValueError(
            'Background ROI (%d x %d px) is larger than the image (%d x %d px).'
            % (bg_w, bg_h, img_w, img_h)
        )

    # Start near image center. User moves the box to the desired tissue location.
    start_x = int((img_w - bg_w) / 2)
    start_y = int((img_h - bg_h) / 2)
    comp.setRoi(Roi(start_x, start_y, bg_w, bg_h))
    try:
        IJ.setTool('rectangle')
    except:
        pass

    dialog = WaitForUserDialog(
        'Place %s background ROI' % side_label,
        'A fixed-size %d x %d pixel background box has been created.\n\n'
        'DRAG the box into your predetermined %s background tissue region.\n'
        'Avoid empty space outside the tissue, tears, large vessels, folds, and '
        'obvious fluorescent artifacts.\n\n'
        'Do not worry if the box is accidentally resized: the script will restore '
        'the specified fixed dimensions around its final center.\n\n'
        'Click OK when it is positioned.' % (bg_w, bg_h, side_label)
    )
    dialog.show()

    current = comp.getRoi()
    if current is None:
        raise ValueError('The background ROI was removed before confirmation.')

    bounds = current.getBounds()
    cx = bounds.x + bounds.width / 2.0
    cy = bounds.y + bounds.height / 2.0

    x = int(round(cx - bg_w / 2.0))
    y = int(round(cy - bg_h / 2.0))
    x = clamp(x, 0, img_w - bg_w)
    y = clamp(y, 0, img_h - bg_h)

    fixed = Roi(x, y, bg_w, bg_h)
    fixed.setName('%s_Background' % side_label)
    comp.setRoi(fixed)
    comp.updateAndDraw()
    return fixed.clone()


def measure_roi_on_green(comp, green_ch, roi):
    old_c = comp.getC()
    comp.setC(green_ch)
    comp.setRoi(roi)
    stats = comp.getStatistics(Measurements.AREA | Measurements.MEAN)
    comp.setC(old_c)
    return stats


def save_roi_zip(roi_dir, index, image_name, region, rois):
    base = os.path.splitext(image_name)[0]
    zip_name = '%03d_%s_%s_ROIs.zip' % (
        index,
        safe_name(base),
        safe_name(region)
    )
    zip_path = os.path.join(roi_dir, zip_name)

    rm = RoiManager(False)
    try:
        for roi in rois:
            rm.addRoi(roi)
        rm.runCommand('Save', zip_path)
    finally:
        try:
            rm.reset()
            rm.close()
        except:
            pass

    return zip_path


def process_image(image_file, index, region, green_ch, blue_ch, bg_w, bg_h,
                  roi_dir, writer):
    path = image_file.getAbsolutePath()
    image_name = image_file.getName()
    imp = IJ.openImage(path)
    if imp is None:
        raise ValueError('Fiji could not open the image.')

    # Immediately rotate the complete multichannel image 90 degrees clockwise.
    # ImageJ's built-in command applies the transform before the composite is
    # displayed or any ROIs are drawn/measured.
    IJ.run(imp, 'Rotate 90 Degrees Right', '')

    comp = None
    try:
        comp = prepare_composite(imp, green_ch, blue_ch)

        # Both anatomical sides are drawn before moving on to backgrounds/image.
        ipsi_roi = get_freehand_area_roi(comp, 'Ipsilateral')
        contra_roi = get_freehand_area_roi(comp, 'Contralateral')

        ipsi_bg = place_fixed_background_roi(comp, 'Ipsilateral', bg_w, bg_h)
        contra_bg = place_fixed_background_roi(comp, 'Contralateral', bg_w, bg_h)

        ipsi_stats = measure_roi_on_green(comp, green_ch, ipsi_roi)
        contra_stats = measure_roi_on_green(comp, green_ch, contra_roi)
        ipsi_bg_stats = measure_roi_on_green(comp, green_ch, ipsi_bg)
        contra_bg_stats = measure_roi_on_green(comp, green_ch, contra_bg)

        ipsi_norm = float('nan') if ipsi_bg_stats.mean == 0 else ipsi_stats.mean / ipsi_bg_stats.mean
        contra_norm = float('nan') if contra_bg_stats.mean == 0 else contra_stats.mean / contra_bg_stats.mean

        roi_zip = save_roi_zip(
            roi_dir, index, image_name, region,
            [ipsi_roi, contra_roi, ipsi_bg, contra_bg]
        )

        cal = comp.getCalibration()
        unit = cal.getUnit()
        if unit is None or unit == '':
            unit = 'pixel'
        area_unit = unit + '^2'

        write_csv_row(writer, [
            image_name, path, region, 'Ipsilateral',
            int(ipsi_stats.pixelCount), '%.6f' % ipsi_stats.area, area_unit,
            '%.6f' % ipsi_stats.mean,
            '%.6f' % ipsi_bg_stats.mean,
            '%.6f' % ipsi_norm,
            bg_w, bg_h,
            ipsi_bg.getBounds().x, ipsi_bg.getBounds().y,
            roi_zip
        ])

        write_csv_row(writer, [
            image_name, path, region, 'Contralateral',
            int(contra_stats.pixelCount), '%.6f' % contra_stats.area, area_unit,
            '%.6f' % contra_stats.mean,
            '%.6f' % contra_bg_stats.mean,
            '%.6f' % contra_norm,
            bg_w, bg_h,
            contra_bg.getBounds().x, contra_bg.getBounds().y,
            roi_zip
        ])

    finally:
        if comp is not None:
            try:
                comp.changes = False
                comp.close()
            except:
                pass
        elif imp is not None:
            try:
                imp.changes = False
                imp.close()
            except:
                pass


def main():
    try:
        settings = get_run_settings()
        if settings is None:
            return
        region, green_ch, blue_ch, bg_w, bg_h = settings

        images = choose_images()
        if len(images) == 0:
            return

        output_path = choose_output_csv(safe_name(region) + '_projection_quantification.csv')
        if output_path is None:
            return

        output_dir = os.path.dirname(output_path)
        output_base = os.path.splitext(os.path.basename(output_path))[0]
        roi_dir = os.path.join(output_dir, safe_name(output_base) + '_ROIs')
        if not os.path.exists(roi_dir):
            os.makedirs(roi_dir)

        writer = BufferedWriter(FileWriter(output_path, False))
        try:
            write_csv_row(writer, [
                'Image',
                'Source_Path',
                'Region',
                'Side',
                'ROI_Area_Pixels',
                'ROI_Area_Calibrated',
                'ROI_Area_Unit',
                'Mean_Green',
                'Background_Mean_Green',
                'Normalized_Green_ROI_div_Background',
                'Background_Width_Pixels',
                'Background_Height_Pixels',
                'Background_X',
                'Background_Y',
                'ROI_Zip'
            ])

            errors = []
            total = len(images)
            for i, image_file in enumerate(images):
                IJ.showStatus('Processing image %d of %d: %s' % (i + 1, total, image_file.getName()))
                try:
                    process_image(
                        image_file, i + 1, region, green_ch, blue_ch,
                        bg_w, bg_h, roi_dir, writer
                    )
                except Exception as e:
                    errors.append('%s: %s' % (image_file.getName(), str(e)))
                    JOptionPane.showMessageDialog(
                        None,
                        'This image could not be completed:\n\n%s\n\nReason:\n%s\n\n'
                        'The script will continue to the next selected image.'
                        % (image_file.getName(), str(e)),
                        'Image skipped',
                        JOptionPane.WARNING_MESSAGE
                    )
        finally:
            writer.close()

        message = (
            'Quantification complete.\n\n'
            'Results CSV:\n%s\n\n'
            'Saved ROI sets:\n%s\n\n'
            'Normalization used:\nMean green ROI / local mean green background'
            % (output_path, roi_dir)
        )
        if len(errors) > 0:
            message += '\n\n%d image(s) were skipped because of errors.' % len(errors)

        JOptionPane.showMessageDialog(
            None,
            message,
            'CA2 projection quantification complete',
            JOptionPane.INFORMATION_MESSAGE
        )
        IJ.showStatus('CA2 projection quantification complete')

    except Exception as e:
        JOptionPane.showMessageDialog(
            None,
            'The script stopped because of an error:\n\n%s' % str(e),
            'Quantification error',
            JOptionPane.ERROR_MESSAGE
        )


main()
