// ============================================================
// Fiji Macro: Robust batch quantification of nuclear cFos
//
// Expected input:
//   Channel 1 = cFos
//   Channel 2 = Hoechst
//
// For each input image a directory is created:
//   outputDir/imageName/
//
// All per-image outputs go there.
// The batch summary CSV stays in outputDir.
// ============================================================

// ---------------- USER SETTINGS ----------------
minSize = 20;
maxSize = 2000;
circularityMin = 0.10;
circularityMax = 1.00;

rollingBallHoechst = 10;
gaussSigma = 0.5;

useCLAHE = true;
claheBlocksize = 63;
claheHistogram = 256;
claheMaxSlope = 3.0;

thresholdMethod = "Default";
doWatershed = true;

ringInner = 2;
ringOuter = 5;

savePreprocessedHoechst = true;
saveAnnotatedPositiveImage = true;
positiveOutlineWidth = 2;
labelPositiveNuclei = false;
saveSegmentationQCImage = true;
segmentationOutlineWidth = 2;
saveMergedPositiveImage = true;
mergedPositiveOutlineWidth = 2;

cfosPositiveThreshold = 1125;

// ---------------- SELECT FOLDERS ----------------
inputDir = getDirectory("Choose folder containing TIFF images");
outputDir = getDirectory("Choose folder for results");

// ---------------- SUMMARY FILE ----------------
outputDirName = getLastFolderName(outputDir);
summaryFile = outputDir + outputDirName + "_Batch_cFos_summary.csv";

File.delete(summaryFile);
File.append("Image,NucleiCount,cFosPosCount,Percent_cFosPos,Mean_cFos_Raw,Mean_cFos_BgSub,Total_cFos_RawIntDen,Total_cFos_BgSubIntDen\n", summaryFile);

// ---------------- FILE LIST ----------------
list = getFileList(inputDir);

// ============================================================
// BATCH LOOP
// ============================================================

for (f = 0; f < list.length; f++) {

    name = list[f];
    lower = toLowerCase(name);

    if (!(endsWith(lower, ".tif") || endsWith(lower, ".tiff")))
        continue;

    baseName = stripExt(name);

    // Create per-image output folder
    imageOutputDir = outputDir + baseName + File.separator;

    File.makeDirectory(imageOutputDir);

    open(inputDir + name);
    origID = getImageID();

    getDimensions(width, height, channels, slices, frames);

    cfosID = -1;
    hoechstID = -1;
    preID = -1;
    maskID = -1;
    annotID = -1;
    qcID = -1;
    mergedID = -1;

    // ---------------- CHANNEL EXTRACTION ----------------
    if (channels == 2) {

        selectImage(origID);
        run("Duplicate...", "title=[cFos] duplicate channels=1");
        cfosID = getImageID();

        selectImage(origID);
        run("Duplicate...", "title=[Hoechst] duplicate channels=2");
        hoechstID = getImageID();

        safeClose(origID);

    } else if (channels == 1 && slices == 2) {

        selectImage(origID);
        setSlice(1);
        run("Duplicate...", "title=[cFos]");
        cfosID = getImageID();

        selectImage(origID);
        setSlice(2);
        run("Duplicate...", "title=[Hoechst]");
        hoechstID = getImageID();

        safeClose(origID);

    } else if (channels == 1 && frames == 2) {

        selectImage(origID);
        run("Duplicate...", "title=[cFos] duplicate frames=1");
        cfosID = getImageID();

        selectImage(origID);
        run("Duplicate...", "title=[Hoechst] duplicate frames=2");
        hoechstID = getImageID();

        safeClose(origID);

    } else {

        print("Skipping " + name);
        safeClose(origID);
        continue;

    }

    // ---------------- HOECHST PREPROCESSING ----------------
    selectImage(hoechstID);

    run("Duplicate...", "title=[Hoechst_pre]");
    preID = getImageID();

    selectImage(preID);

    run("Subtract Background...", "rolling=" + rollingBallHoechst);

    if (useCLAHE) {

        run("Enhance Local Contrast (CLAHE)",
            "blocksize=" + claheBlocksize +
            " histogram=" + claheHistogram +
            " maximum=" + claheMaxSlope +
            " mask=*None* fast_(less_accurate)");

    }

    run("Gaussian Blur...", "sigma=" + gaussSigma);

    if (savePreprocessedHoechst)

        saveAs("Tiff", imageOutputDir + baseName + "_Hoechst_preprocessed.tif");

    // ---------------- SEGMENTATION ----------------
    run("8-bit");
    setAutoThreshold(thresholdMethod);
    run("Convert to Mask");

    run("Invert");

    run("Fill Holes");

    if (doWatershed)
        run("Watershed");

    maskID = getImageID();

    saveAs("Tiff", imageOutputDir + baseName + "_Hoechst_mask.tif");

    // ---------------- NUCLEI ROIS ----------------
    roiManager("Reset");
    run("Clear Results");

    selectImage(maskID);

    run("Set Measurements...",
        "area mean min centroid center shape integrated decimal=3");

    run("Analyze Particles...",
        "size=" + minSize + "-" + maxSize +
        " circularity=" + circularityMin + "-" + circularityMax +
        " show=Nothing exclude clear add");

    nucleiCount = roiManager("count");

    if (nucleiCount == 0) {

        File.append(name + ",0,0,0,0,0,0,0\n", summaryFile);
        cleanupByID(cfosID, hoechstID, preID, maskID, annotID, qcID, mergedID);
        continue;

    }

    roiManager("Save",
        imageOutputDir + baseName + "_Nuclei_ROIs.zip");

    // ---------------- SEGMENTATION QC IMAGE ----------------
    if (saveSegmentationQCImage) {

        selectImage(hoechstID);
        run("Select None");
        run("Duplicate...", "title=[Hoechst_segmentation_QC]");
        qcID = getImageID();

        selectImage(qcID);
        run("RGB Color");
        run("Remove Overlay");
        run("Select None");
        setForegroundColor(0, 255, 0);
        setColor(0, 255, 0);
        setLineWidth(segmentationOutlineWidth);

        for (i = 0; i < nucleiCount; i++) {
            selectImage(qcID);
            roiManager("Select", i);
            run("Add Selection...");
            run("Select None");
        }

        selectImage(qcID);
        run("Select None");
        run("Flatten");
        flatID = getImageID();
        safeClose(qcID);
        qcID = flatID;

        selectImage(qcID);
        run("Select None");
        saveAs("Png", imageOutputDir + baseName + "_Hoechst_segmentation_QC.png");
    }

    // ---------------- CFOS MEASUREMENTS ----------------
    selectImage(cfosID);

    run("16-bit");

    run("Clear Results");

    run("Set Measurements...",
        "area mean min max std integrated centroid center decimal=3");

    roiManager("Measure");

    n = nResults;

    rawArea = newArray(n);
    rawMean = newArray(n);
    rawIntDen = newArray(n);

    for (i = 0; i < n; i++) {

        rawArea[i] = getResult("Area", i);
        rawMean[i] = getResult("Mean", i);
        rawIntDen[i] = getResult("IntDen", i);

    }

    bgMean = newArray(n);
    bgSubMean = newArray(n);
    bgSubIntDen = newArray(n);

    isPositive = newArray(n);

    // ---------------- BACKGROUND RING ----------------
    for (i = 0; i < n; i++) {

        roiManager("Select", i);
        run("Enlarge...", "enlarge=" + ringOuter);
        run("Measure");

        outerMean = getResult("Mean", nResults - 1);
        outerArea = getResult("Area", nResults - 1);

        roiManager("Select", i);
        run("Enlarge...", "enlarge=" + ringInner);
        run("Measure");

        innerMean = getResult("Mean", nResults - 1);
        innerArea = getResult("Area", nResults - 1);

        ringArea = outerArea - innerArea;

        if (ringArea > 0) {

            ringSum = outerMean * outerArea - innerMean * innerArea;
            bgMean[i] = ringSum / ringArea;

        } else {

            bgMean[i] = 0;

        }

        bgSubMean[i] = rawMean[i] - bgMean[i];

        if (bgSubMean[i] < 0)
            bgSubMean[i] = 0;

        bgSubIntDen[i] = bgSubMean[i] * rawArea[i];

        if (bgSubMean[i] >= cfosPositiveThreshold)
            isPositive[i] = 1;
        else
            isPositive[i] = 0;

    }

    // ---------------- ANNOTATED POSITIVE NUCLEI IMAGE ----------------
    if (saveAnnotatedPositiveImage) {

        // Important: clear any active ROI on the cFos image first.
        // Otherwise Duplicate can crop to the selected nucleus instead of the full field.
        selectImage(cfosID);
        run("Select None");
        run("Duplicate...", "title=[cFos_positive_overlay]");
        annotID = getImageID();

        selectImage(annotID);
        run("RGB Color");
        run("Remove Overlay");
        run("Select None");
        setForegroundColor(255, 0, 0);
        setColor(255, 0, 0);
        setLineWidth(positiveOutlineWidth);
        setFont("SansSerif", 12, "antialiased");

        for (i = 0; i < n; i++) {
            if (isPositive[i] == 1) {

                // Make the annotation image active before applying the ROI from ROI Manager.
                // That attaches the selection/overlay to the cFos duplicate, not the Hoechst mask.
                selectImage(annotID);
                roiManager("Select", i);
                getSelectionBounds(x, y, w, h);
                run("Add Selection...");

                if (labelPositiveNuclei) {
                    drawString("" + (i+1), x, y-2);
                }

                run("Select None");
            }
        }

        selectImage(annotID);
        run("Select None");
        run("Flatten");
        flatID = getImageID();
        saveAs("Png", imageOutputDir + baseName + "_cFos_positive_annotated.png");
        if (imageExists(annotID)) safeClose(annotID);
        annotID = flatID;
    }

    // ---------------- MERGED CFOS + HOECHST POSITIVE OVERLAY ----------------
    if (saveMergedPositiveImage) {

        // Capture the current source display ranges so the merged PNG more closely
        // matches what Fiji shows for the extracted source channels.
        selectImage(cfosID);
        getMinAndMax(cfosDispMin, cfosDispMax);
        run("Select None");
        run("Duplicate...", "title=[cFos_merge_tmp]");
        cfosMergeID = getImageID();

        selectImage(hoechstID);
        getMinAndMax(hoechstDispMin, hoechstDispMax);
        run("Select None");
        run("Duplicate...", "title=[Hoechst_merge_tmp]");
        hoechstMergeID = getImageID();

        // Reapply the original display ranges to the duplicated channels instead of
        // autoscaling them with Enhance Contrast. This makes the merged rendering
        // visually more comparable to the source image.
        selectImage(cfosMergeID);
        setMinAndMax(cfosDispMin, cfosDispMax);
        run("Grays");

        selectImage(hoechstMergeID);
        setMinAndMax(hoechstDispMin, hoechstDispMax);
        run("Grays");

        run("Merge Channels...", "c1=[cFos_merge_tmp] c3=[Hoechst_merge_tmp] create");
        mergedID = getImageID();

        if (imageExists(cfosMergeID)) safeClose(cfosMergeID);
        if (imageExists(hoechstMergeID)) safeClose(hoechstMergeID);

        selectImage(mergedID);
        run("Remove Overlay");
        run("Select None");
        setForegroundColor(255, 255, 0);
        setColor(255, 255, 0);
        setLineWidth(mergedPositiveOutlineWidth);
        setFont("SansSerif", 12, "antialiased");

        for (i = 0; i < n; i++) {
            if (isPositive[i] == 1) {
                selectImage(mergedID);
                roiManager("Select", i);
                getSelectionBounds(x, y, w, h);
                run("Add Selection...");

                if (labelPositiveNuclei) {
                    drawString("" + (i+1), x, y-2);
                }

                run("Select None");
            }
        }

        selectImage(mergedID);
        run("Select None");
        run("Flatten");
        flatID = getImageID();
        saveAs("Png", imageOutputDir + baseName + "_cFos_positive_on_merged.png");
        if (imageExists(mergedID)) safeClose(mergedID);
        mergedID = flatID;
    }

    // ---------------- PER-CELL CSV ----------------
    cellFile = imageOutputDir + baseName + "_cFos_cells.csv";

    File.delete(cellFile);

    File.append("Cell,Area,RawMean,RawIntDen,LocalBgMean,BgSubMean,BgSubIntDen,cFosPositive\n", cellFile);

    sumRawMean = 0;
    sumBgSubMean = 0;
    sumRawIntDen = 0;
    sumBgSubIntDen = 0;

    cfosPosCount = 0;

    for (i = 0; i < n; i++) {

        File.append(
            (i+1) + "," +
            rawArea[i] + "," +
            rawMean[i] + "," +
            rawIntDen[i] + "," +
            bgMean[i] + "," +
            bgSubMean[i] + "," +
            bgSubIntDen[i] + "," +
            isPositive[i] + "\n",
            cellFile);

        sumRawMean += rawMean[i];
        sumBgSubMean += bgSubMean[i];
        sumRawIntDen += rawIntDen[i];
        sumBgSubIntDen += bgSubIntDen[i];

        cfosPosCount += isPositive[i];

    }

    meanRawMean = sumRawMean / n;
    meanBgSubMean = sumBgSubMean / n;

    percentPos = 100 * cfosPosCount / n;

    // ---------------- SUMMARY CSV ----------------
    File.append(
        name + "," +
        nucleiCount + "," +
        cfosPosCount + "," +
        percentPos + "," +
        meanRawMean + "," +
        meanBgSubMean + "," +
        sumRawIntDen + "," +
        sumBgSubIntDen + "\n",
        summaryFile);

    cleanupByID(cfosID, hoechstID, preID, maskID, annotID, qcID, mergedID);

}

print("Batch analysis complete.");

// ============================================================
// HELPERS
// ============================================================

function stripExt(nm) {

    lowername = toLowerCase(nm);

    if (endsWith(lowername, ".tiff"))
        return substring(nm, 0, lengthOf(nm)-5);

    if (endsWith(lowername, ".tif"))
        return substring(nm, 0, lengthOf(nm)-4);

    return nm;

}

function imageExists(id) {

    n = nImages;

    for (k = 1; k <= n; k++) {

        selectImage(k);

        if (getImageID() == id)
            return true;

    }

    return false;

}

function safeClose(id) {

    if (id == -1)
        return;

    if (imageExists(id)) {

        selectImage(id);
        close();

    }

}

function cleanupByID(cfosID, hoechstID, preID, maskID, annotID, qcID, mergedID) {

    run("Clear Results");
    roiManager("Reset");

    if (annotID != -1) safeClose(annotID);

    if (qcID != -1 && qcID != annotID) safeClose(qcID);

    if (mergedID != -1 && mergedID != annotID && mergedID != qcID) safeClose(mergedID);

    if (maskID != -1 && maskID != annotID && maskID != qcID && maskID != mergedID) safeClose(maskID);

    if (preID != -1 && preID != maskID && preID != annotID && preID != qcID && preID != mergedID)
        safeClose(preID);

    if (hoechstID != -1 && hoechstID != preID && hoechstID != maskID && hoechstID != annotID && hoechstID != qcID && hoechstID != mergedID)
        safeClose(hoechstID);

    if (cfosID != -1 && cfosID != hoechstID && cfosID != preID && cfosID != maskID && cfosID != annotID && cfosID != qcID && cfosID != mergedID)
        safeClose(cfosID);

}

function getLastFolderName(path) {

    p = path;

    if (endsWith(p, "\\") || endsWith(p, "/"))
        p = substring(p, 0, lengthOf(p)-1);

    lastSlash = lastIndexOf(p, "\\");
    if (lastSlash == -1)
        lastSlash = lastIndexOf(p, "/");

    if (lastSlash == -1)
        return p;

    return substring(p, lastSlash + 1, lengthOf(p));

}