// Publication adaptation (2026-10-12): input/output directory selectors only.
// Quantification and QC logic unchanged; see README and provenance manifest.
// ============================================================================
//  IF_identity_percell.ijm      v2.4 (2026-08-08)
//
//  Per-cell quantification of the two IDENTITY markers (vimentin, alpha-SMA)
//  and their secondary-only controls.
//
//  ---------------------------------------------------------------------------
//  EXCLUSION POLICY  (changed in v2 - read this before running)
//  ---------------------------------------------------------------------------
//  This macro makes NO judgement-based exclusion. Every raw series is attempted
//  and every raw series gets exactly one row in perfield.csv.
//
//  A series yields no cells ONLY on a hard technical failure, each of which is
//  recorded in the `status` column:
//      open_fail   Bio-Formats returned no image
//      not16bit    not 16-bit (pipeline assumes 16-bit histograms)
//      nC<2        fewer than 2 channels
//      no_nuclei   segmentation returned 0 particles after the size/circ filter
//      svcc        series name marks it as computationally cleared / LIGHTNING
//                  (a DATA-SELECTION rule, not a QC rule: only unprocessed raw
//                   series are analysed - this must be stated in Methods)
//
//  Everything that used to be a skip is now only a FLAG in the `flags` column.
//  The exclusion decision is therefore made at the CSV stage, prospectively,
//  and can be reported both ways:
//      dapi_sat    C0_max = 65535 (nuclei may be clipped/merged)
//      dapi_weak   C0_max < C0_WEAK
//      fitc_sat    C1_max = 65535 (binary positivity unaffected; intensity is)
//      low_nuclei  fewer than LOW_NUC cells segmented
//
//  WHY: for a positivity readout, exclusion changes both numerator and
//  denominator, and any rule that removes fields unequally from the stained
//  and the control coverslip biases the threshold itself. Never define an
//  exclusion rule on the FITC (readout) channel - that is the error that
//  invalidated the focusC1 metric in the p-S6 dataset. Focus is reported here
//  from the DAPI counterstain only (nucEdgeRatio) and is NOT used to exclude.
//
//  Per-cell attrition is also an exclusion criterion and is made auditable:
//      nObjAll   all nuclear objects, no filters at all
//      nObjEdge  after edge-touching removal only
//      nCells    after 30-600 um^2 + circularity >= 0.35   <- analysed
//  ---------------------------------------------------------------------------
//
//  Expected layout (DIR is scanned recursively):
//      <identity_raw_root>/
//          Vimentin/  Vimentin.lif      -> VIM / stain
//                     Vimentin Sac.lif  -> VIM / control
//          SMA/       SMA.lif           -> SMA / stain
//                     SMA-Sac.lif       -> SMA / control
//      .lifext companions must stay beside their .lif.
//
//  marker  from path      : "vimentin" -> VIM ; "sma" -> SMA
//  role    from file name : "sac" -> control ; otherwise stain
//
//  OUTPUT -> ~/if_identity_out/  : percell.csv, perfield.csv, run_log.txt,
//                                  roi/, overlay/
//
//  RUN: set DIR, leave PROBE = true, run, read the Log. If sane, PROBE = false.
//  BATCH must stay false (Bio-Formats + setBatchMode(true) silently yields
//  nImages = 0 and every series gets skipped).
// ============================================================================

// ---------------------------- USER SETTINGS ---------------------------------
var DIR          = getDirectory("Select VIC identity LIF root with Vimentin and SMA folders");
var PROBE        = true;    // true = first PROBE_SERIES series of EVERY file
var PROBE_SERIES = 3;
var BATCH        = false;   // KEEP FALSE

var MIN_AREA = 30;          // um^2   |  unified pipeline, identical to p-S6
var MAX_AREA = 600;         // um^2   |
var MIN_CIRC = 0.35;        //        |
var RING_UM  = 5;           // perinuclear ring width (calibrated units)

var SATVAL   = 65535;       // 16-bit ceiling
var C0_WEAK  = 6000;        // DAPI flag threshold only - does NOT skip
var LOW_NUC  = 5;           // low_nuclei flag threshold only - does NOT skip

var SAVE_ROI         = true;
var SAVE_OVERLAY     = true;
var OVERLAY_PER_FILE = 2;   // -1 = every field

// ---------------------------- GLOBALS ---------------------------------------
var OUT = "", pathPer = "", pathFld = "", pathLog = "";
var nRowPer = 0, nRowFld = 0;

var fMarker, fRole, fFile, fSeries, fName, fStatus, fFlags;
var fW, fH, fNC, fNZ, fBD, fPx;
var fExp0, fExp1, fLed, fPosX, fPosY, fPosZ, fTime;
var fObjAll, fObjEdge, fSmall, fLarge, fLowCirc, fCells;
var f0min, f0med, f0max, f0sat;
var f1min, f1med, f1mode, f1p99, f1max, f1sat, f1zero;
var fRingMed, fNetMed, fEdgeMed;

// ============================================================================
//  MAIN
// ============================================================================
requires("1.53");
run("Set Measurements...", "area mean standard modal min centroid median shape redirect=None decimal=4");
setOption("BlackBackground", true);
setBatchMode(BATCH);

OUT = getDirectory("Select empty VIC identity quantification output folder");
File.makeDirectory(OUT);
File.makeDirectory(OUT + "roi");
File.makeDirectory(OUT + "overlay");

pathPer = OUT + "percell.csv";
pathFld = OUT + "perfield.csv";
pathLog = OUT + "run_log.txt";

File.saveString("marker,role,file,series,seriesName,cell,nucArea_um2,nucX_um,nucY_um," +
                "ringArea_um2,ringMean,ringMedian,ringMin,ringMax,fieldMode,net,nucEdgeRatio\n", pathPer);
File.saveString("marker,role,file,series,seriesName,status,flags,w,h,nC,nZ,bitDepth,px_um," +
                "expC0_s,expC1_s,ledC1,posX,posY,posZ,timestamp,nObjAll,nObjEdge,nSmall,nLarge,nLowCirc,nCells," +
                "C0_min,C0_med,C0_max,C0_satPct,C1_min,C1_med,C1_mode,C1_p99,C1_max,C1_satPct,C1_zeroPct," +
                "ringMedian_field,netMedian_field,dapiNucEdge\n", pathFld);
File.saveString("IF_identity_percell.ijm v2.4 - " + getTimeStamp() + "\n", pathLog);

// File.append self-test: an earlier macro silently produced header-only CSVs
len0 = lengthOf(File.openAsString(pathLog));
File.append("__SELFTEST__", pathLog);
len1 = lengthOf(File.openAsString(pathLog));
if (len1 <= len0) {
    print("!! FATAL: File.append does not work on " + pathLog + " - choose another OUT");
    exit();
}
File.saveString("IF_identity_percell.ijm v2.4 - " + getTimeStamp() + "\n", pathLog);
logLine("append self-test OK (" + len0 + " -> " + len1 + " bytes)");

if (!File.exists(DIR)) exit("DIR does not exist: " + DIR);

lifs = newArray(0);
lifs = listLif(DIR, lifs);
logLine("found " + lifs.length + " .lif under " + DIR);
for (i = 0; i < lifs.length; i++)
    logLine("   [" + i + "] " + lifs[i] + "   marker=" + markerOf(lifs[i]) + "  role=" + roleOf(lifs[i]));

for (f = 0; f < lifs.length; f++) processLif(lifs[f]);

logLine("");
logLine("DONE.  percell rows = " + nRowPer + "  (" + File.length(pathPer) + " bytes)");
logLine("       perfield rows = " + nRowFld + "  (" + File.length(pathFld) + " bytes)");
logLine("       output folder = " + OUT);
setBatchMode(false);


// ============================================================================
//  PER-FILE
// ============================================================================
function processLif(path) {
    stem   = stripExt(baseName(path));
    marker = markerOf(path);
    role   = roleOf(path);

    logLine("");
    logLine("=== " + baseName(path) + "   marker=" + marker + "  role=" + role);

    run("Bio-Formats Macro Extensions");
    Ext.setId(path);
    Ext.getSeriesCount(nSer);
    logLine("    series count = " + nSer);

    nDo = nSer;
    if (PROBE && nDo > PROBE_SERIES) nDo = PROBE_SERIES;

    ovCount = 0;

    for (s = 0; s < nDo; s++) {
        Ext.setSeries(s);
        Ext.getSeriesName(sName);

        resetField();
        fMarker = marker;  fRole = role;  fFile = stem;
        fSeries = s;       fName = csvSafe(sName);

        if (isProcessed(sName)) {
            fStatus = "svcc";
            emitField();
            logLine("    s" + s + " '" + sName + "' -> not raw (computationally cleared)");
            continue;
        }

        close("*");
        roiManager("reset");
        run("Bio-Formats Importer",
            "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (s + 1));

        if (nImages == 0) {
            fStatus = "open_fail";
            emitField();
            logLine("    s" + s + " -> OPEN FAILED");
            continue;
        }

        idOrig = getImageID();
        getDimensions(w, h, nC, nZ, nT);
        getPixelSize(unit, pw, ph);
        info = getMetadata("Info");

        fW = w;  fH = h;  fNC = nC;  fNZ = nZ;  fBD = bitDepth();  fPx = d4(pw);
        fExp0 = "" + grabNth(info, "Exposure", 1);
        fExp1 = "" + grabNth(info, "Exposure", 2);
        fLed  = "" + grab(info, "IntensityDev");
        if (fLed == "NA") fLed = "" + grab(info, "LightSourceIntensity");
        if (fLed == "NA") fLed = "" + grab(info, "Intensity");
        fPosX = "" + grab(info, "PosX");
        fPosY = "" + grab(info, "PosY");
        fPosZ = "" + grab(info, "PosZ");
        fTime = "" + grab(info, "Timestamp");
        if (fTime == "NA") fTime = "" + grab(info, "AcquiredDate");
        if (fTime == "NA") fTime = "" + grab(info, "StartTime");

        if (fBD != 16) { fStatus = "not16bit"; emitField();
                         logLine("    s" + s + " -> " + fBD + "-bit, skipped"); continue; }
        if (nC < 2)    { fStatus = "nC<2";     emitField();
                         logLine("    s" + s + " -> " + nC + " channel(s), skipped"); continue; }

        // ------------- whole-field channel statistics ------------------------
        selectImage(idOrig);
        Stack.setChannel(1);
        run("Select None");
        getRawStatistics(nPix, mean0, min0, max0, sd0, h0);
        f0min = min0;
        f0med = pctFromHist(h0, nPix, 50);
        f0max = max0;
        f0sat = 0;
        if (h0.length > SATVAL) f0sat = d4(100.0 * h0[SATVAL] / nPix);

        Stack.setChannel(2);
        run("Select None");
        getRawStatistics(nPix, mean1, min1, max1, sd1, h1);
        mode1  = argmaxHist(h1);
        f1min  = min1;
        f1med  = pctFromHist(h1, nPix, 50);
        f1mode = mode1;
        f1p99  = pctFromHist(h1, nPix, 99);
        f1max  = max1;
        f1sat  = 0;
        if (h1.length > SATVAL) f1sat = d4(100.0 * h1[SATVAL] / nPix);
        f1zero = d4(100.0 * h1[0] / nPix);

        // ------------- FLAGS ONLY - no field is dropped here ------------------
        if (max0 >= SATVAL) addFlag("dapi_sat");
        if (max0 <  C0_WEAK) addFlag("dapi_weak");
        if (max1 >= SATVAL) addFlag("fitc_sat");

        // ------------- DAPI segmentation --------------------------------------
        selectImage(idOrig);
        Stack.setChannel(1);
        run("Duplicate...", "title=dapiRaw duplicate channels=1");
        idDapi = getImageID();

        run("Duplicate...", "title=dapiEdge");
        idEdge = getImageID();
        run("32-bit");
        run("Find Edges");

        selectImage(idDapi);
        run("Duplicate...", "title=dapiMask");
        idMask = getImageID();
        run("Gaussian Blur...", "sigma=2");
        setAutoThreshold("Otsu dark");
        run("Convert to Mask");
        run("Fill Holes");
        run("Watershed");

        // attrition audit: no filter -> edge only -> full filter
        run("Clear Results");
        run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 display clear");
        fObjAll = nResults;
        run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 exclude display clear");
        fObjEdge = nResults;
        nSmall = 0; nLarge = 0; nLowCirc = 0;
        for (r = 0; r < nResults; r++) {
            ar = getResult("Area", r);
            ci = getResult("Circ.", r);
            if (ar < MIN_AREA) nSmall++;
            if (ar > MAX_AREA) nLarge++;
            if (ci < MIN_CIRC) nLowCirc++;
        }
        fSmall = nSmall; fLarge = nLarge; fLowCirc = nLowCirc;   // categories overlap
        run("Clear Results");

        roiManager("reset");
        run("Analyze Particles...",
            "size=" + MIN_AREA + "-" + MAX_AREA + " circularity=" + MIN_CIRC + "-1.00 exclude add");
        nNuc = roiManager("count");
        fCells = nNuc;

        if (nNuc == 0) {
            fStatus = "no_nuclei";
            emitField();
            logLine("    s" + s + " -> 0 nuclei after filter (objects: all=" + fObjAll +
                    " edge=" + fObjEdge + ")");
            close("*");
            continue;
        }
        if (nNuc < LOW_NUC) addFlag("low_nuclei");

        // ------------- rings ---------------------------------------------------
        selectImage(idOrig);
        for (k = 0; k < nNuc; k++) {
            roiManager("Select", k);
            run("Enlarge...", "enlarge=" + RING_UM);
            roiManager("Add");
        }

        ringVals = newArray(nNuc);
        netVals  = newArray(nNuc);
        edgeVals = newArray(nNuc);

        for (k = 0; k < nNuc; k++) {
            selectImage(idOrig);
            Stack.setChannel(1);
            roiManager("Select", k);
            nucArea = getValue("Area");
            nucX    = getValue("X");
            nucY    = getValue("Y");

            selectImage(idDapi);
            run("Restore Selection");
            dapiMean = getValue("Mean");

            selectImage(idEdge);
            run("Restore Selection");
            edgeMean = getValue("Mean");
            run("Select None");

            edgeRatio = 0;
            if (dapiMean > 0) edgeRatio = edgeMean / dapiMean;
            edgeVals[k] = edgeRatio;

            selectImage(idOrig);
            roiManager("Select", newArray(k, nNuc + k));
            roiManager("XOR");
            Stack.setChannel(2);
            ringArea = getValue("Area");
            ringMean = getValue("Mean");
            ringMed  = getValue("Median");
            ringMin  = getValue("Min");
            ringMax  = getValue("Max");

            net = ringMean - mode1;
            ringVals[k] = ringMean;
            netVals[k]  = net;

            row = marker + "," + role + "," + stem + "," + s + "," + fName + "," + (k + 1) + "," +
                  d4(nucArea) + "," + d4(nucX) + "," + d4(nucY) + "," +
                  d4(ringArea) + "," + d4(ringMean) + "," + d4(ringMed) + "," +
                  d4(ringMin) + "," + d4(ringMax) + "," + mode1 + "," + d4(net) + "," + d4(edgeRatio);
            File.append(row, pathPer);
            nRowPer++;
        }

        fRingMed = d4(medianOf(ringVals));
        fNetMed  = d4(medianOf(netVals));
        fEdgeMed = d4(medianOf(edgeVals));
        fStatus  = "ok";
        emitField();

        logLine("    s" + s + " '" + sName + "'  ok  cells=" + nNuc +
                " (objects all=" + fObjAll + " edge=" + fObjEdge + ")" +
                "  C0max=" + f0max + "  C1 mode=" + f1mode + " med=" + f1med + " max=" + f1max +
                "  ring=" + fRingMed + "  net=" + fNetMed + "  focus=" + fEdgeMed +
                "  flags=" + fFlags);

        if (SAVE_ROI) {
            roiManager("Deselect");
            roiManager("Save", OUT + "roi" + File.separator + stem + "_s" + pad3(s) + "_RoiSet.zip");
        }

        okOverlay = false;
        if (SAVE_OVERLAY) {
            if (OVERLAY_PER_FILE < 0) okOverlay = true;
            if (ovCount < OVERLAY_PER_FILE) okOverlay = true;
        }
        if (okOverlay) {
            // CRITICAL: the per-cell loop leaves a selection active on idDapi, and
            // run("Duplicate...") CROPS TO THE ACTIVE SELECTION. Without this
            // Select None the overlay is a thumbnail of the last nucleus, not the field.
            selectImage(idDapi);
            run("Select None");
            run("Duplicate...", "title=ov");
            idOv = getImageID();
            run("Enhance Contrast", "saturated=0.35");
            run("Grays");                 // Bio-Formats hands back a GREEN LUT
            run("RGB Color");
            setForegroundColor(0, 255, 255);              // cyan  = nucleus (analysed)
            for (k = 0; k < nNuc; k++) { roiManager("Select", k); run("Draw", "slice"); }
            setForegroundColor(255, 255, 0);              // yellow = ring outer edge
            for (k = 0; k < nNuc; k++) { roiManager("Select", nNuc + k); run("Draw", "slice"); }

            // red = objects that survived the edge filter but were REJECTED by the
            // size / circularity filter - lets you see whether the rejects are
            // watershed fragments or real nuclei. ROI manager is safe to reuse here:
            // rings are already measured and the RoiSet has already been saved.
            selectImage(idMask);
            roiManager("reset");
            run("Clear Results");
            run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 exclude display add");
            selectImage(idOv);
            setForegroundColor(255, 0, 0);
            for (r = 0; r < roiManager("count"); r++) {
                ar = getResult("Area", r);
                ci = getResult("Circ.", r);
                rej = false;
                if (ar < MIN_AREA) rej = true;
                if (ar > MAX_AREA) rej = true;
                if (ci < MIN_CIRC) rej = true;
                if (rej) { roiManager("Select", r); run("Draw", "slice"); }
            }
            run("Clear Results");
            selectImage(idOv);
            run("Select None");
            saveAs("PNG", OUT + "overlay" + File.separator + stem + "_s" + pad3(s) + "_seg.png");
            close();
            ovCount++;
        }

        close("*");
    }

    Ext.close();
}


// ============================================================================
//  FIELD ROW - single assembly point, so column count can never drift
// ============================================================================
function resetField() {
    fMarker = "NA"; fRole = "NA"; fFile = "NA"; fSeries = "NA"; fName = "NA";
    fStatus = "NA"; fFlags = "";
    fW = "NA"; fH = "NA"; fNC = "NA"; fNZ = "NA"; fBD = "NA"; fPx = "NA";
    fExp0 = "NA"; fExp1 = "NA"; fLed = "NA"; fPosX = "NA"; fPosY = "NA"; fPosZ = "NA"; fTime = "NA";
    fObjAll = "NA"; fObjEdge = "NA"; fSmall = "NA"; fLarge = "NA"; fLowCirc = "NA"; fCells = 0;
    f0min = "NA"; f0med = "NA"; f0max = "NA"; f0sat = "NA";
    f1min = "NA"; f1med = "NA"; f1mode = "NA"; f1p99 = "NA"; f1max = "NA";
    f1sat = "NA"; f1zero = "NA";
    fRingMed = "NA"; fNetMed = "NA"; fEdgeMed = "NA";
}

function addFlag(t) {
    if (lengthOf(fFlags) == 0) fFlags = t;
    else fFlags = fFlags + ";" + t;
}

function emitField() {
    fl = fFlags;
    if (lengthOf(fl) == 0) fl = "none";
    row = "" + fMarker + "," + fRole + "," + fFile + "," + fSeries + "," + fName + "," +
          fStatus + "," + fl + "," + fW + "," + fH + "," + fNC + "," + fNZ + "," + fBD + "," + fPx + "," +
          fExp0 + "," + fExp1 + "," + fLed + "," + fPosX + "," + fPosY + "," + fPosZ + "," + fTime + "," +
          fObjAll + "," + fObjEdge + "," + fSmall + "," + fLarge + "," + fLowCirc + "," + fCells + "," +
          f0min + "," + f0med + "," + f0max + "," + f0sat + "," +
          f1min + "," + f1med + "," + f1mode + "," + f1p99 + "," + f1max + "," + f1sat + "," + f1zero + "," +
          fRingMed + "," + fNetMed + "," + fEdgeMed;
    File.append(row, pathFld);
    nRowFld++;
}

function logLine(s) {
    print(s);
    File.append(s, pathLog);
}


// ============================================================================
//  UTILITIES
// ============================================================================
function listLif(dir, acc) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        p = dir + list[i];
        if (endsWith(list[i], "/")) {
            acc = listLif(p, acc);
        } else {
            lw = toLowerCase(list[i]);
            if (endsWith(lw, ".lif")) {
                if (!startsWith(list[i], ".")) acc = Array.concat(acc, p);
            }
        }
    }
    return acc;
}

function baseName(p) {
    i = lastIndexOf(p, "/");
    if (i < 0) i = lastIndexOf(p, "\\");
    return substring(p, i + 1);
}

function stripExt(n) {
    i = lastIndexOf(n, ".");
    if (i < 0) return n;
    return substring(n, 0, i);
}

function markerOf(p) {
    lw = toLowerCase(p);
    if (indexOf(lw, "viment") >= 0) return "VIM";
    if (indexOf(lw, "sma") >= 0) return "SMA";
    return "UNK";
}

function roleOf(p) {
    lw = toLowerCase(baseName(p));
    if (indexOf(lw, "sac") >= 0) return "control";
    return "stain";
}

function isProcessed(nm) {
    lw = toLowerCase(nm);
    if (indexOf(lw, "_lng") >= 0) return true;
    if (indexOf(lw, "svcc") >= 0) return true;
    if (indexOf(lw, "lightning") >= 0) return true;
    if (indexOf(lw, "deconv") >= 0) return true;
    return false;
}

function csvSafe(s) {
    s = replace(s, ",", ";");
    s = replace(s, "\n", " ");
    s = replace(s, "\r", " ");
    return s;
}

function grab(info, key) {
    v = grabNth(info, key, 1);      // must go via a variable: the interpreter
    return v;                       // reads `return userFunc(...)` as numeric
}

function grabNth(info, key, n) {
    from  = 0;
    count = 0;
    while (true) {
        rest = substring(info, from);
        idx  = indexOf(rest, key);
        if (idx < 0) return "NA";
        hitPos = from + idx;
        sub = substring(info, hitPos);
        e = indexOf(sub, "\n");
        if (e < 0) e = lengthOf(sub);
        line = substring(sub, 0, e);
        c = indexOf(line, "=");
        count++;
        if (count >= n) {
            if (c < 0) return "NA";
            v = substring(line, c + 1);   // one step per line, then return the
            v = trimStr(v);               // variable - never return a call to a
            v = csvSafe(v);               // user-defined string function directly
            return v;
        }
        from = hitPos + lengthOf(key);
    }
    return "NA";
}

function trimStr(s) {
    while (lengthOf(s) > 0) {
        if (startsWith(s, " ") || startsWith(s, "\t")) s = substring(s, 1);
        else break;
    }
    while (lengthOf(s) > 0) {
        if (endsWith(s, " ") || endsWith(s, "\t") || endsWith(s, "\r")) s = substring(s, 0, lengthOf(s) - 1);
        else break;
    }
    return s;
}

function argmaxHist(h) {
    best = -1; bi = 0;
    for (i = 0; i < h.length; i++) { if (h[i] > best) { best = h[i]; bi = i; } }
    return bi;
}

// 16-bit getRawStatistics returns a histogram of length max+1, NOT 65536 -
// never index it with a hard-coded 65535 without a length check.
function pctFromHist(h, nPix, pct) {
    target = nPix * pct / 100.0;
    cum = 0;
    for (i = 0; i < h.length; i++) {
        cum += h[i];
        if (cum >= target) return i;
    }
    return h.length - 1;
}

function medianOf(a) {
    b = Array.copy(a);
    Array.sort(b);
    n = b.length;
    if (n == 0) return 0;
    if (n % 2 == 1) return b[floor(n / 2)];
    return (b[n / 2 - 1] + b[n / 2]) / 2.0;
}

function pad3(n) {
    s = "" + n;
    while (lengthOf(s) < 3) s = "0" + s;
    return s;
}

function d4(x) {
    v = d2s(x, 4);
    return v;
}

function getTimeStamp() {
    getDateAndTime(y, mo, dw, dm, hh, mm, ss, ms);
    return "" + y + "-" + (mo + 1) + "-" + dm + " " + hh + ":" + mm + ":" + ss;
}
