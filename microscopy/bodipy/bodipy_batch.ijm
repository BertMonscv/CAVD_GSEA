// bodipy_batch.ijm
// One click for the whole folder. Do not open TIFFs one by one.
//
// A  Batch prepare
//      LIF or multi-plane TIFF -> one single-plane BODIPY TIFF per field
//      (and a matching DAPI TIFF). Z = the plane with the highest BODIPY mean.
//      No autoscale, no 8-bit stretch.
//
// C  Draw cells, script saves
//      You only draw and press T. The macro names and saves the zip.
//      Project001_s01_Z03_BODIPY.tif -> Project001_s01_Z03_BODIPY.zip

//
// Install: Plugins > Macros > Install...
// Mac folder dialog: click Open to choose the folder. Do not double-click a file.

var MIN_AREA_UM2 = 0.1;
var MAX_AREA_UM2 = 1000000;
var THR_PREF = "bodipy.rapa.thr";

macro "A Batch prepare single-plane TIFFs" {
    batchPrepare();
}

macro "B Batch quantify all fields" {
    batchQuantify();
}

macro "C Draw cells, script saves" {
    drawCells();
}

function stemOf(name) {
    low = toLowerCase(name);
    if (endsWith(low, ".tiff")) return substring(name, 0, lengthOf(name) - 5);
    if (endsWith(low, ".tif")) return substring(name, 0, lengthOf(name) - 4);
    return name;
}

function twoDigits(n) {
    if (n < 10) return "0" + n;
    return "" + n;
}

function drawCells() {
    dir = getDirectory("Folder of *_BODIPY.tif. Click Open, do not open a file");
    if (dir == "") exit("Cancelled");
    list = getFileList(dir);
    names = "";
    nFiles = 0;
    for (i = 0; i < list.length; i++) {
        low = toLowerCase(list[i]);
        if (endsWith(low, "_bodipy.tif") || endsWith(low, "_bodipy.tiff")) {
            names = names + list[i] + "|";
            nFiles++;
        }
    }
    if (nFiles == 0) exit("No *_BODIPY.tif in that folder. Run A first. BODIPY channel = 1.");
    files = split(names, "|");
    saved = 0;
    already = 0;
    for (i = 0; i < lengthOf(files); i++) {
        image = files[i];
        zip = dir + stemOf(image) + ".zip";
        if (File.exists(zip)) {
            already++;
            continue;
        }
        dapi = replace(image, "_BODIPY", "_DAPI");
        open(dir + image);
        bodTitle = getTitle();
        if (File.exists(dir + dapi)) {
            open(dir + dapi);
            dapiTitle = getTitle();
            run("Merge Channels...", "c2=[" + bodTitle + "] c3=[" + dapiTitle + "] create");
        }
        run("ROI Manager...");
        roiManager("reset");
        roiManager("Show All with labels");
        setTool("polygon");
        waitForUser("Draw cells  " + (i + 1) + "/" + nFiles, image + "\n\nPolygon is selected.\nOutline one cell, press T. Repeat.\nThen click OK. The zip name is written for you.\n\nCancel stops. Fields already saved are kept.");
        nRoi = roiManager("count");
        if (nRoi > 0) {
            roiManager("Save", zip);
            print("Saved " + nRoi + " -> " + zip);
            saved++;
        } else {
            print("No cells, not saved: " + image);
        }
        close();
    }
    showMessage("Drawing pass done", "Saved: " + saved + "\nAlready had a zip: " + already + "\n\nNext: B Batch quantify all fields");
}

function batchPrepare() {
    dir = getDirectory("Input folder: LIF or TIFF. Click Open, do not open a file");
    if (dir == "") exit("Cancelled");
    out = getDirectory("Empty output folder. Click Open");
    if (out == "") exit("Cancelled");
    if (dir == out) exit("Input and output folders must be different.");

    Dialog.create("Channels");
    Dialog.addNumber("BODIPY channel (1-based)", 1);
    Dialog.addNumber("DAPI channel (0 = do not save)", 2);
    Dialog.show();
    bodCh = Dialog.getNumber();
    dapiCh = Dialog.getNumber();
    if (bodCh < 1) exit("BODIPY channel must be >= 1");

    list = getFileList(dir);
    run("Bio-Formats Macro Extensions");
    logPath = out + "prepare_log.csv";
    File.append("source,series,channels,z_total,z_kept,bodipy_mean,pixel_um,unit,bit_depth,bodipy_file,dapi_file\n", logPath);

    nOk = 0;
    nSkip = 0;
    setBatchMode(true);
    for (i = 0; i < list.length; i++) {
        name = list[i];
        if (endsWith(name, "/")) continue;
        low = toLowerCase(name);
        if (!(endsWith(low, ".lif") || endsWith(low, ".tif") || endsWith(low, ".tiff"))) continue;
        path = dir + name;
        Ext.setId(path);
        Ext.getSeriesCount(nSeries);
        if (nSeries < 1) nSeries = 1;
        for (s = 1; s <= nSeries; s++) {
            opts = "open=[" + path + "] color_mode=Grayscale view=Hyperstack stack_order=XYCZT series_" + s + " windowless=true";
            run("Bio-Formats Importer", opts);
            if (nImages == 0) {
                print("SKIP open failed: " + name + " series " + s);
                nSkip++;
                continue;
            }
            kept = savePlanes(out, name, s, bodCh, dapiCh, logPath);
            close();
            if (kept) nOk++;
            else nSkip++;
        }
    }
    setBatchMode(false);
    showMessage("Prepare done", "Fields written: " + nOk + "\nSkipped: " + nSkip + "\n\n" + out + "\n\nLog: prepare_log.csv");
}

function savePlanes(out, name, series, bodCh, dapiCh, logPath) {
    getDimensions(w, h, c, z, f);
    getVoxelSize(px, py, pz, unit);
    bd = bitDepth();
    if (bodCh > c) {
        print("SKIP " + name + " series " + series + ": BODIPY channel " + bodCh + " but image has " + c + " channels");
        File.append(name + "," + series + "," + c + "," + z + ",,,," + unit + "," + bd + ",SKIP,\n", logPath);
        return 0;
    }
    best = 1;
    bestMean = -1;
    for (zi = 1; zi <= z; zi++) {
        if (c > 1) Stack.setPosition(bodCh, zi, 1);
        else setSlice(zi);
        getStatistics(area, mean);
        if (mean > bestMean) {
            bestMean = mean;
            best = zi;
        }
    }
    base = replace(replace(replace(name, ".lif", ""), ".LIF", ""), ".tif", "");
    base = replace(replace(base, ".TIF", ""), ".tiff", "");
    base = replace(base, ".TIFF", "");
    stem = base + "_s" + twoDigits(series) + "_Z" + twoDigits(best);
    bodName = stem + "_BODIPY.tif";
    dapName = "";
    src = getImageID();

    selectImage(src);
    run("Duplicate...", "title=bodipy duplicate channels=" + bodCh + "-" + bodCh + " slices=" + best + "-" + best + " frames=1-1");
    saveAs("Tiff", out + bodName);
    close();

    if (dapiCh >= 1 && dapiCh <= c && dapiCh != bodCh) {
        selectImage(src);
        dapName = stem + "_DAPI.tif";
        run("Duplicate...", "title=dapi duplicate channels=" + dapiCh + "-" + dapiCh + " slices=" + best + "-" + best + " frames=1-1");
        saveAs("Tiff", out + dapName);
        close();
    }
    selectImage(src);
    File.append(name + "," + series + "," + c + "," + z + "," + best + "," + bestMean + "," + px + "," + unit + "," + bd + "," + bodName + "," + dapName + "\n", logPath);
    print(bodName + "  Z " + best + "/" + z);
    return 1;
}

function batchQuantify() {
    dir = getDirectory("Folder of *_BODIPY.tif and matching .zip ROIs");
    if (dir == "") exit("Cancelled");
    sheet = dir + "sample_sheet.csv";
    if (!File.exists(sheet)) {
        writeSheet(dir, sheet);
        showMessage("Fill the sheet", "Wrote\n" + sheet + "\n\nFill experiment (1, 2 or 3) and group (PM, PM+Lep or PM+Lep+Rapa).\nSave it, then run this command again.\n\nROI file for each row:\nsame name as the TIFF, extension .zip");
        return;
    }

    lastThr = parseFloat(call("ij.Prefs.get", THR_PREF, "100"));
    Dialog.create("Frozen threshold");
    Dialog.addNumber("Raw gray threshold", lastThr);
    Dialog.show();
    thr = Dialog.getNumber();
    if (thr < 1) exit("Threshold must be >= 1");
    call("ij.Prefs.set", THR_PREF, "" + thr);

    csv = dir + "bodipy_rapamycin_cells.csv";
    if (!File.exists(csv))
        File.append("experiment,group,field,cell,n_droplets,total_area_um2,mean_droplet_area_um2,pixel_um,unit,threshold,bit_depth,image\n", csv);

    text = File.openAsString(sheet);
    text = replace(text, "\r", "");
    lines = split(text, "\n");
    nOk = 0;
    nSkip = 0;
    setBatchMode(true);
    for (i = 0; i < lines.length; i++) {
        line = lines[i];
        if (line == "" || startsWith(line, "image")) continue;
        cols = split(line, ",");
        if (lengthOf(cols) < 3) { nSkip++; continue; }
        image = cols[0];
        expt = cols[1];
        grp = cols[2];
        if (image == "" || expt == "" || grp == "") { nSkip++; continue; }
        if (grp != "PM" && grp != "PM+Lep" && grp != "PM+Lep+Rapa") {
            print("SKIP bad group: " + line);
            nSkip++;
            continue;
        }
        tif = dir + image;
        zip = dir + replace(replace(image, ".tif", ".zip"), ".tiff", ".zip");
        if (!File.exists(tif)) { print("SKIP missing TIFF: " + image); nSkip++; continue; }
        if (!File.exists(zip)) { print("SKIP missing ROI zip: " + zip); nSkip++; continue; }
        open(tif);
        if (!calibratedMicrons()) {
            print("SKIP not calibrated in microns: " + image);
            close();
            nSkip++;
            continue;
        }
        roiManager("reset");
        roiManager("Open", zip);
        nRoi = roiManager("count");
        if (nRoi < 1) {
            print("SKIP no ROIs: " + zip);
            close();
            nSkip++;
            continue;
        }
        field = replace(image, ".tif", "");
        measureOpenImage(csv, expt, grp, field, thr);
        close();
        nOk++;
    }
    setBatchMode(false);
    showMessage("Quantify done", "Fields measured: " + nOk + "\nSkipped: " + nSkip + "\n\n" + csv);
}

function writeSheet(dir, sheet) {
    list = getFileList(dir);
    File.append("image,experiment,group\n", sheet);
    for (i = 0; i < list.length; i++) {
        low = toLowerCase(list[i]);
        if (endsWith(low, "_bodipy.tif") || endsWith(low, "_bodipy.tiff"))
            File.append(list[i] + ",,\n", sheet);
    }
}

function calibratedMicrons() {
    getVoxelSize(px, py, pz, unit);
    if (px == 0 || py == 0) return 0;
    unitLow = toLowerCase(unit);
    if (indexOf(unitLow, "micron") < 0 && indexOf(unit, "µm") < 0 && indexOf(unitLow, "um") < 0) return 0;
    return 1;
}

function measureOpenImage(csv, expt, grp, field, thr) {
    getDimensions(w, h, c, z, f);
    if (c != 1 || z != 1) {
        print("SKIP not a single plane: " + field);
        return;
    }
    getVoxelSize(px, py, pz, unit);
    if ((MIN_AREA_UM2 / (px * py)) < 1)
        print("Warning: 0.1 um^2 is below 1 pixel in " + field);
    maxThr = 255;
    bd = bitDepth();
    if (bd == 16) maxThr = 65535;
    if (bd == 32) maxThr = 1e30;
    src = getImageID();
    nRoi = roiManager("count");
    for (i = 0; i < nRoi; i++) {
        selectImage(src);
        roiManager("select", i);
        run("Duplicate...", "title=celltmp");
        dup = getImageID();
        run("Clear Outside");
        run("Select None");
        setThreshold(thr, maxThr, "raw");
        if (isOpen("Results")) run("Clear Results");
        run("Analyze Particles...", "size=" + MIN_AREA_UM2 + "-" + MAX_AREA_UM2 + " display clear");
        count = nResults;
        totalUm = 0;
        for (r = 0; r < count; r++) totalUm += getResult("Area", r);
        meanUm = 0;
        if (count > 0) meanUm = totalUm / count;
        img = replace(field, ",", "_");
        row = expt + "," + grp + "," + img + "," + (i + 1) + "," + count + "," + d2s(totalUm, 6) + "," + d2s(meanUm, 6) + "," + px + "," + unit + "," + d2s(thr, 2) + "," + bd + "," + img;
        File.append(row + "\n", csv);
        selectImage(dup);
        close();
    }
    selectImage(src);
    resetThreshold();
    print(field + "  cells=" + nRoi);
}
