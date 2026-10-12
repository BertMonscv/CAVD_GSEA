// make_bodipy_figure.ijm
// Panel A: BODIPY 493/503 representative fields (PM, PM+Lep, PM+Lep+Rapa)
// Row 1: BODIPY (green); Row 2: merge of BODIPY (green) + DAPI (blue). No grayscale rows.
// All three fields use the SAME fixed display range (no per-image Auto).
// Output: figure_images.tif in the figure folder.

bodipyDir = getDirectory("Folder containing BODIPY TIFFs");
dapiDir = getDirectory("Folder containing matching DAPI TIFFs");
outDir = getDirectory("Output folder for the representative panel");

// Representative fields retained from the original image-assembly source.
// Field-to-group assignments require the associated laboratory records.
// Each panel is a uniform 512 x 512 px (166.4 x 166.4 um) crop of the quantified field.
// Crops were placed so that the quantified cells inside (ROI bbox >60% inside the crop)
// have median droplet number / area close to the group mean of field medians:
//   PM          s19 crop (288,160): 4 cells, 33.5 droplets, 168 um2  (group 34.0 / 168)
//   PM+Lep      s15 crop (512,192): 5 cells, 40.0 droplets, 487 um2  (group 48.5 / 442)
//                s15 is the PM+Lep field with the highest median droplet area (701 um2);
//                chosen for illustration (s17 has a dividing cell, s35/s37 high DAPI background).
//   PM+Lep+Rapa s09 crop (288,512): 6 cells, 15.0 droplets,  26 um2  (group 13.5 / 31)
names  = newArray("Project001_s19_Z01_BODIPY", "Project001_s15_Z01_BODIPY", "Project001_s09_Z01_BODIPY");
cropX  = newArray(288, 512, 288);
cropY  = newArray(160, 192, 512);
cropS  = 512;
labels = newArray("PM", "PM+Lep", "PM+Lep+Rapa");

// Fixed display ranges (16-bit raw values), identical for all three panels (no Auto).
bMin = 490;  bMax = 2900;   // BODIPY
dMin = 455;  dMax = 3650;   // DAPI

pixelUm  = 0.3251;
barUm    = 20;

// Layout (px)
W = cropS; H = cropS;
gap = 8; left = 60; top = 50;
fontSize = 32;
nRows = 2;

setBatchMode(true);

canvasW = left + 3*W + 2*gap;
canvasH = top + nRows*H + (nRows-1)*gap;
newImage("figure_images", "RGB white", canvasW, canvasH, 1);
canvas = getImageID();

for (i = 0; i < 3; i++) {
    bPath = bodipyDir + names[i] + ".tif";
    dName = replace(names[i], "_BODIPY", "_DAPI");
    dPath = dapiDir + dName + ".tif";
    if (!File.exists(dPath)) dPath = bodipyDir + dName + ".tif";
    if (!File.exists(bPath)) exit("Missing: " + bPath);
    if (!File.exists(dPath)) exit("Missing: " + dPath);
    x0 = left + i*(W+gap);

    open(bPath); b = getTitle();
    makeRectangle(cropX[i], cropY[i], cropS, cropS); run("Crop");
    setMinAndMax(bMin, bMax); run("8-bit");
    open(dPath); d = getTitle();
    makeRectangle(cropX[i], cropY[i], cropS, cropS); run("Crop");
    setMinAndMax(dMin, dMax); run("8-bit");

    // Row 1: BODIPY in green (same green LUT as in the merge)
    selectImage(b); run("Duplicate...", "title=tmp"); run("Green"); run("RGB Color");
    run("Select All"); run("Copy"); close("tmp");
    selectImage(canvas); makeRectangle(x0, top, W, H); run("Paste");

    // Row 2: merge, green = BODIPY, blue = DAPI
    run("Merge Channels...", "c2=[" + b + "] c3=[" + d + "]");
    merged = getImageID();

    // Scale bar only on the bottom-right panel
    if (i == 2) {
        barPx = round(barUm / pixelUm);
        barH = 7; margin = 20;
        setColor(255, 255, 255);
        fillRect(W - margin - barPx, H - margin - barH, barPx, barH);
        setFont("Arial", fontSize, "antialiased");
        setJustification("right");
        drawString(barUm + " " + fromCharCode(181) + "m", W - margin, H - margin - barH - 6);
    }

    run("Select All"); run("Copy");
    selectImage(canvas); makeRectangle(x0, top + (H+gap), W, H); run("Paste");
    selectImage(merged); close();
}

// Column labels
selectImage(canvas); run("Select None");
setColor(0, 0, 0);
setFont("Arial", fontSize, "antialiased");
setJustification("center");
for (i = 0; i < 3; i++)
    drawString(labels[i], left + i*(W+gap) + W/2, top - 10);

// Row labels (rotated)
rowLab = newArray("BODIPY", "Merge");
for (k = 0; k < nRows; k++) {
    newImage("lab", "RGB white", H, left, 1);
    setColor(0, 0, 0); setFont("Arial", fontSize, "antialiased"); setJustification("center");
    drawString(rowLab[k], H/2, left - 15);
    run("Rotate 90 Degrees Left");
    run("Select All"); run("Copy"); close("lab");
    selectImage(canvas);
    makeRectangle(0, top + k*(H+gap), left, H); run("Paste");
}
run("Select None");

// Calibrate output so image size is meaningful
run("Set Scale...", "distance=1 known=" + pixelUm + " unit=micron");
saveAs("Tiff", outDir + "figure_images.tif");
close("*");
setBatchMode(false);
print("Saved " + outDir + "figure_images.tif");
