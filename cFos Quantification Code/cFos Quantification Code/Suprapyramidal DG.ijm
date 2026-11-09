dir1 = getDirectory("Choose Source Directory For Blue");
dir2 = getDirectory("Choose Source Directory For Red");
list1 = getFileList(dir1);
list2 = getFileList(dir2);
Array.sort(list1);
Array.sort(list2);
dir5 = getDirectory("Choose Destination Directory");

run("Set Measurements...", "area mean min centroid redirect=None decimal=3");

for (i = 0; i < list1.length; i++) {
    showProgress(i + 1, list1.length);

    // Open blue channel
    open(dir1 + list1[i]);
    blueTitle = getTitle();
    run("Rotate 90 Degrees Right");
    run("Set Scale...", "distance=1 known=0.650 unit=μm");

    // Open red channel
    open(dir2 + list2[i]);
    redTitle = getTitle();
    run("Rotate 90 Degrees Right");
    run("Set Scale...", "distance=1 known=0.650 unit=μm");

    // Merge channels so ROI is chosen on both channels together
    run("Merge Channels...", "c1=[" + redTitle + "] c3=[" + blueTitle + "] create");
    run("Rotate...", "angle=-20 interpolation=Bilinear");
    mergedTitle = getTitle();

    // Draw ROI on merged image
    makeRectangle(1300, 3500, 850, 200);
    waitForUser("Move ROI on the merged image, then hit OK");
    Roi.getBounds(x, y, width, height);

    // Crop merged image and save as a single merged TIFF
    makeRectangle(x, y, width, height);
    run("Crop");
    saveAs("TIFF", dir5 + list1[i]);

    // Close all open images
    close(); // cropped merged image

}