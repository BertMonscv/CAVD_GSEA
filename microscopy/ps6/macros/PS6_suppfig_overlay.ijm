// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  PS6_suppfig_overlay.ijm  —— 補充圖 S10：分割稽核疊圖
//
//  對「Figure 4C 用的那 5 個視野」重跑一次與定量完全相同的分割，把核輪廓與
//  量測環畫出來。因為用的是同一批視野，讀者可以把 S10 直接疊到 Figure 4C 上看。
//
//  ⚠️ 分割參數必須與 PS6_percell.ijm 逐項相同，否則這張圖就不是稽核而是另一次分析：
//     Gaussian σ=2 → Otsu → Fill Holes → Watershed → 30–600 µm²、circ ≥0.35、exclude edge
//     → 每核 Enlarge 5 µm 後與核 XOR = 量測環
//
//  輸出 ~/ps6_suppfig/：每組 <組>_seg_dapi.tif（DAPI 底）與 <組>_seg_ps6.tif（p-S6 底），
//  外加拼版 FigureS10_segmentation.tif（2 列 × 5 欄）與 manifest.csv
// ============================================================================

// ---------------------------- 輸入 ----------------------------
var DIR = getDirectory("Select raw p-S6 LIF input folder");
var PLABEL = newArray("CTL", "PM", "PMLep", "PMLepRapa", "Sac");
var PFILE  = newArray("CTL2.lif", "PM2.lif", "LEP2.lif", "RAPA2.lif", "PS6-Sac.lif");
var PSER   = newArray(1, 5, 12, 13, 0);          // 與 Figure 4C 的 SET C 相同

// ---------------------------- 顯示範圍（與 Figure 4C 一致，不可另訂）----------
var D0MIN = 500;   var D0MAX = 14000;   // DAPI
var D1MIN = 400;   var D1MAX = 4000;    // p-S6

// ---------------------------- 分割參數（與 PS6_percell.ijm 相同）--------------
var CH_DAPI = 1;
var CH_PS6  = 2;
var PX       = 0.325;
var RING_UM  = 5;
var NUC_MIN  = 30;
var NUC_MAX  = 600;
var NUC_CIRC = 0.35;
var GAUSS_SIG = 2;

// ---------------------------- 樣式 ----------------------------
var LW_NUC   = 2;                // 核輪廓線寬
var LW_RING  = 2;                // 環外緣線寬
var LABEL_ON = true;             // 是否標細胞編號
var LABEL_FS = 20;
var BAR_UM   = 20;
var BAR_H    = 8;
var BAR_MARGIN = 18;
var BAR_FONT = 24;
var GUTTER   = 14;
var ASSEMBLE = true;

// ---------------------------- 內部 ----------------------------
var OUT, PATH_MAN, PATH_LOG;
var bufLog = "";

// ============================== 主體 ==============================
print("\\Clear");
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "area mean redirect=None decimal=3");
setBatchMode(false);              // Bio-Formats 在批次模式下會開了圖卻不進影像清單

OUT = getDirectory("Select empty p-S6 overlay output folder");
File.makeDirectory(OUT);
PATH_MAN = OUT + "manifest.csv";
PATH_LOG = OUT + "overlay_log.txt";
File.saveString("panel,file,series,nuclei,ringUM,nucMin,nucMax,circ,D0MIN,D0MAX,D1MIN,D1MAX\n", PATH_MAN);
File.saveString("PS6_suppfig_overlay.ijm\n", PATH_LOG);
logLine("分割參數：Gaussian " + GAUSS_SIG + " / Otsu / Watershed / " + NUC_MIN + "-" + NUC_MAX
        + " um2 / circ >= " + NUC_CIRC + " / ring " + RING_UM + " um");
logLine("顯示範圍：DAPI " + D0MIN + "-" + D0MAX + "  p-S6 " + D1MIN + "-" + D1MAX);

for (i = 0; i < PLABEL.length; i++) {
    path = DIR + PFILE[i];
    if (!File.exists(path)) { logLine("!! 找不到 " + path); continue; }
    run("Close All"); roiManager("reset");
    run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default "
        + "view=Hyperstack stack_order=XYCZT series_" + (PSER[i] + 1));
    if (nImages == 0) { logLine("!! 開檔失敗 " + PFILE[i] + " s" + PSER[i]); continue; }
    idSrc = getImageID();
    getDimensions(W, H, nC, nZ, nT);
    logLine("");
    logLine("[" + PLABEL[i] + "] " + PFILE[i] + " s" + PSER[i] + "  " + W + "x" + H);

    selectImage(idSrc); run("Duplicate...", "title=D16 duplicate channels=" + CH_DAPI);
    idD = getImageID();
    selectImage(idSrc); run("Duplicate...", "title=G16 duplicate channels=" + CH_PS6);
    idG = getImageID();

    // ---- 分割（與定量逐項相同）----
    selectImage(idD);
    run("Duplicate...", "title=MASK");
    idM = getImageID();
    run("Gaussian Blur...", "sigma=" + GAUSS_SIG);
    setAutoThreshold("Otsu dark");
    run("Convert to Mask"); run("Fill Holes"); run("Watershed");
    run("Analyze Particles...", "size=" + NUC_MIN + "-" + NUC_MAX
        + " circularity=" + NUC_CIRC + "-1.00 exclude add");
    nNuc = roiManager("count");
    logLine("    分出核數 " + nNuc);
    closeIfOpen(idM); closeTitle("MASK");

    // ---- 外擴環（索引 nNuc..2nNuc-1 為外緣）----
    for (k = 0; k < nNuc; k++) {
        selectImage(idG);
        roiManager("select", k);
        run("Enlarge...", "enlarge=" + RING_UM);
        roiManager("add");
    }

    drawOverlay(idD, D0MIN, D0MAX, nNuc, OUT + PLABEL[i] + "_seg_dapi.tif");
    drawOverlay(idG, D1MIN, D1MAX, nNuc, OUT + PLABEL[i] + "_seg_ps6.tif");

    File.append(PLABEL[i] + "," + PFILE[i] + "," + PSER[i] + "," + nNuc + "," + RING_UM + ","
        + NUC_MIN + "," + NUC_MAX + "," + NUC_CIRC + "," + D0MIN + "," + D0MAX + ","
        + D1MIN + "," + D1MAX, PATH_MAN);
    run("Close All"); roiManager("reset");
}

if (ASSEMBLE) assemble();
flushLog();
print("完成。輸出：" + OUT);

// ============================== 函式 ==============================

// 把 16-bit 通道依指定顯示範圍轉成灰階 RGB，再畫核輪廓（青）、環外緣（黃）與編號
function drawOverlay(id16, dmin, dmax, nNuc, savePath) {
    selectImage(id16); run("Select None");
    run("Duplicate...", "title=BASE");
    idB = getImageID();
    setMinAndMax(dmin, dmax);
    run("Grays");                    // ⚠️ Bio-Formats 帶進來的 LUT 是綠的，不強制灰階會烘成綠圖
    run("8-bit");
    run("RGB Color");
    idR = getImageID();

    setFont("SansSerif", LABEL_FS, "bold antialiased");
    for (k = 0; k < nNuc; k++) {
        selectImage(idR);
        setForegroundColor(255, 255, 0);           // 環外緣 = 黃
        run("Line Width...", "line=" + LW_RING);
        roiManager("select", nNuc + k);
        run("Draw", "slice");
        setForegroundColor(0, 255, 255);           // 核輪廓 = 青
        run("Line Width...", "line=" + LW_NUC);
        roiManager("select", k);
        run("Draw", "slice");
        if (LABEL_ON) {
            roiManager("select", k);
            Roi.getBounds(bx, by, bw, bh);
            selectImage(idR); run("Select None");
            setColor(255, 255, 255);
            drawString("" + (k + 1), bx + bw / 2 - 6, by - 3);
        }
    }
    selectImage(idR); run("Select None");
    drawBar(idR, round(BAR_UM / PX));
    saveAs("Tiff", savePath);
    idR = getImageID();
    closeIfOpen(idR);
    closeIfOpen(idB);
    closeTitle("BASE");
}

function drawBar(id, barPx) {
    selectImage(id);
    w = getWidth(); h = getHeight();
    bx = w - BAR_MARGIN - barPx; by = h - BAR_MARGIN - BAR_H;
    setForegroundColor(255, 255, 255);
    makeRectangle(bx, by, barPx, BAR_H);
    run("Fill", "slice");
    run("Select None");
    setFont("SansSerif", BAR_FONT, "bold antialiased");
    setColor(255, 255, 255);
    lab = "" + round(barPx * PX) + " " + fromCharCode(181) + "m";
    tw = getStringWidth(lab);
    drawString(lab, bx + barPx / 2 - tw / 2, by - 4);
}

// 2 列（DAPI 底 / p-S6 底）× 5 欄
function assemble() {
    first = OUT + PLABEL[0] + "_seg_dapi.tif";
    if (!File.exists(first)) { logLine("!! 拼版缺檔，跳過"); return; }
    open(first); cw = getWidth(); ch = getHeight(); close();
    n = PLABEL.length;
    tw = n * cw + (n - 1) * GUTTER;
    th = 2 * ch + GUTTER;
    setPasteMode("Copy");
    newImage("ASSEMBLED", "RGB white", tw, th, 1);
    idA = getImageID();
    for (i = 0; i < n; i++) {
        pasteInto(idA, OUT + PLABEL[i] + "_seg_dapi.tif", i * (cw + GUTTER), 0, cw, ch);
        pasteInto(idA, OUT + PLABEL[i] + "_seg_ps6.tif",  i * (cw + GUTTER), ch + GUTTER, cw, ch);
    }
    selectImage(idA); run("Select None");
    saveAs("Tiff", OUT + "FigureS10_segmentation.tif");
    logLine("");
    logLine("拼版 " + tw + " x " + th + " px");
    logLine("置於 170 mm 雙欄寬時 = " + d2s(tw / (170 / 25.4), 0) + " dpi");
}

function pasteInto(idTarget, srcPath, px, py, w, h) {
    if (!File.exists(srcPath)) { logLine("!! 缺檔 " + srcPath); return; }
    open(srcPath);
    idS = getImageID();
    run("Select All"); run("Copy");
    closeIfOpen(idS);
    selectImage(idTarget);
    makeRectangle(px, py, w, h);
    run("Paste");
    run("Select None");
}

function closeIfOpen(id) { if (isOpen(id)) { selectImage(id); close(); } }
function closeTitle(t)   { if (isOpen(t))  { selectWindow(t); close(); } }

function logLine(t) {
    print(t);
    bufLog = bufLog + t + "\n";
}

function flushLog() {
    if (lengthOf(bufLog) > 0) { File.append(trimNL(bufLog), PATH_LOG); bufLog = ""; }
}

function trimNL(t) {
    while (endsWith(t, "\n")) { t = substring(t, 0, lengthOf(t) - 1); }
    return "" + t;
}
