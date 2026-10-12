// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  COLOC_D_supplementary.ijm  —— 兩張補充圖
//
//   S-A  染色對照代表圖   3 列（二抗對照 / CTL / 氯喹）× 3 欄（BODIPY / LC3B / merge）
//   S-B  閾值遮罩示範     1 列 × 6 欄（原圖與二值遮罩並排，含重疊遮罩）
//
//  兩張都沿用 Figure 5F 的顯示範圍、裁切尺寸與比例尺樣式，故可與主圖直接對照。
//  來源：~/coloc_out/blind/F0xx.tif ＋ ~/coloc_out/roi/F0xx.zip
// ============================================================================

// ---- 與 5F 完全相同的樣式（函式要用，必須 var）----
var PX         = 0.1802;
var D1MIN = 0;  var D1MAX = 150;   // DAPI
var D2MIN = 0;  var D2MAX = 200;   // BODIPY
var D3MIN = 0;  var D3MAX = 130;   // LC3B
var BAR_H      = 8;
var BAR_MARGIN = 18;
var BAR_LABEL  = true;
var BAR_FONT   = 24;
var TG = 100;   // BODIPY 閾值
var TR = 90;    // LC3B 閾值

// ---- 加強對比版（只用於另一張獨立的補充圖；三列一律套同一組數值）----
var ENHANCE   = true;
var D2MAX_ENH = 50;    // BODIPY 加強版上限
var D3MAX_ENH = 30;    // LC3B  加強版上限

OUT   = getDirectory("Select colocalization working folder");
BLIND = OUT + "blind/";
ROIDIR= OUT + "roi/";
SUP   = OUT + "supplementary/";

// ---- S-A：對照組（依定量資料挑最具代表性者）----
A_LAB  = newArray("Secondary_only", "CTL",      "Chloroquine");
A_CODE = newArray("F018",           "F016",     "F039");
A_ROI  = newArray(1,                1,          2);
A_SRC  = newArray("Sac s1",         "CTL s3",   "CQs s2");
A_W_UM = 111;  A_H_UM = 60;   // 與 5F 同尺寸同倍率
A_BAR  = 20;

// ---- S-B：閾值遮罩示範（用 PM+Lep 的一個細胞）----
B_CODE = "F017";
B_ROI  = 4;              // 1-based
B_SRC  = "LEP1 s1";
B_UM   = 60;             // 正方形裁切邊長
B_BAR  = 10;

GUTTER = 14;

// ---------------------------------------------------------------------------
setBatchMode(true);
File.makeDirectory(SUP);
print("\\Clear");

aw = round(A_W_UM / PX); ah = round(A_H_UM / PX);
bs = round(B_UM / PX);

// =================== S-A 對照組代表圖 =======================================
for (i = 0; i < lengthOf(A_LAB); i++) {
    code = A_CODE[i]; ri = A_ROI[i] - 1; lab = A_LAB[i];
    rf = roiPath(ROIDIR, code);
    if (!File.exists(BLIND + code + ".tif") || rf == "") {
        print("MISSING: " + code);
    } else {
        close("*");
        open(BLIND + code + ".tif"); rename("SRC");
        run("Properties...", "channels=3 slices=1 frames=1 pixel_width=" + PX + " pixel_height=" + PX + " voxel_depth=1 unit=um");
        roiManager("reset"); roiManager("Open", rf); roiManager("Select", ri);
        Roi.getBounds(bx, by, bw, bh);
        cx = bx + bw/2; cy = by + bh/2;
        run("Select None");

        x0 = clampX(round(cx - aw/2), aw);
        y0 = clampY(round(cy - ah/2), ah);

        chanPanel("SRC", 2, "Green",   D2MIN, D2MAX, x0, y0, aw, ah, SUP + "SA_" + lab + "_BODIPY.tif", 0);
        chanPanel("SRC", 3, "Magenta", D3MIN, D3MAX, x0, y0, aw, ah, SUP + "SA_" + lab + "_LC3B.tif",   0);
        mergePanel("SRC", x0, y0, aw, ah, SUP + "SA_" + lab + "_merge.tif", round(A_BAR/PX), -1);
        if (ENHANCE) {
            chanPanel("SRC", 2, "Green",   D2MIN, D2MAX_ENH, x0, y0, aw, ah, SUP + "SAe_" + lab + "_BODIPY.tif", 0);
            chanPanel("SRC", 3, "Magenta", D3MIN, D3MAX_ENH, x0, y0, aw, ah, SUP + "SAe_" + lab + "_LC3B.tif",   0);
        }
        print("S-A " + lab + " : " + code + " (" + A_SRC[i] + ") roi " + A_ROI[i]);
        roiManager("reset"); close("*");
    }
}

// =================== S-B 閾值遮罩示範 =======================================
rf = roiPath(ROIDIR, B_CODE);
if (!File.exists(BLIND + B_CODE + ".tif") || rf == "") {
    print("MISSING: " + B_CODE);
} else {
    close("*");
    open(BLIND + B_CODE + ".tif"); rename("SRC");
    run("Properties...", "channels=3 slices=1 frames=1 pixel_width=" + PX + " pixel_height=" + PX + " voxel_depth=1 unit=um");
    roiManager("reset"); roiManager("Open", rf); roiManager("Select", B_ROI - 1);
    Roi.getBounds(bx, by, bw, bh);
    cx = bx + bw/2; cy = by + bh/2;
    run("Select None");
    x0 = clampX(round(cx - bs/2), bs);
    y0 = clampY(round(cy - bs/2), bs);

    // 原圖兩通道
    chanPanel("SRC", 2, "Green",   D2MIN, D2MAX, x0, y0, bs, bs, SUP + "SB_1_BODIPY_raw.tif", 0);
    chanPanel("SRC", 3, "Magenta", D3MIN, D3MAX, x0, y0, bs, bs, SUP + "SB_3_LC3B_raw.tif",   0);

    // 二值遮罩（只保留細胞 ROI 之內，與定量口徑一致）
    maskPanel("SRC", 2, TG, "Green",   B_ROI - 1, x0, y0, bs, SUP + "SB_2_BODIPY_mask.tif");
    maskPanel("SRC", 3, TR, "Magenta", B_ROI - 1, x0, y0, bs, SUP + "SB_4_LC3B_mask.tif");

    // 重疊遮罩（兩者皆為陽性的像素 = M2 分子所在）
    overlapPanel("SRC", B_ROI - 1, x0, y0, bs, SUP + "SB_5_overlap.tif");

    // merge ＋ 細胞輪廓 ＋ 比例尺
    mergePanel("SRC", x0, y0, bs, bs, SUP + "SB_6_merge.tif", round(B_BAR/PX), B_ROI - 1);
    print("S-B : " + B_CODE + " (" + B_SRC + ") roi " + B_ROI);
    roiManager("reset"); close("*");
}

// =================== 拼版 ===================================================
assembleA(SUP, A_LAB, aw, ah, GUTTER);
if (ENHANCE) assembleEnh(SUP, A_LAB, aw, ah, GUTTER);
assembleB(SUP, bs, GUTTER);

setBatchMode(false);
print("DONE -> " + SUP);
File.saveString(getInfo("log"), SUP + "supplementary_log.txt");

// ===========================================================================
function roiPath(dir, code) {
    if (File.exists(dir + code + ".zip")) return dir + code + ".zip";
    if (File.exists(dir + code + ".roi")) return dir + code + ".roi";
    return "";
}
function clampX(v, w) { if (v < 0) v = 0; if (v + w > 1024) v = 1024 - w; return v; }
function clampY(v, h) { if (v < 0) v = 0; if (v + h > 1024) v = 1024 - h; return v; }

function drawBar(sw, sh, barPx) {
    if (barPx > 0) {
        setForegroundColor(255, 255, 255);
        bx0 = sw - barPx - BAR_MARGIN;
        by0 = sh - BAR_H - BAR_MARGIN;
        makeRectangle(bx0, by0, barPx, BAR_H);
        run("Fill", "slice");
        run("Select None");
        if (BAR_LABEL) {
            setFont("SansSerif", BAR_FONT, "bold antialiased");
            setColor(255, 255, 255);
            txt = d2s(barPx * PX, 0) + " " + fromCharCode(181) + "m";
            tw = getStringWidth(txt);
            tx = bx0 + barPx/2 - tw/2;
            if (tx < 2) tx = 2;
            drawString(txt, tx, by0 - round(BAR_FONT * 0.25));
        }
    }
}

// 單通道原圖面板
function chanPanel(src, ch, lut, dmin, dmax, x0, y0, w, h, outPath, barPx) {
    selectWindow(src); run("Select None");
    run("Duplicate...", "title=CH duplicate channels=" + ch);
    run(lut); setMinAndMax(dmin, dmax); run("RGB Color");
    makeRectangle(x0, y0, w, h);
    run("Duplicate...", "title=CHC");
    run("Select None");
    drawBar(w, h, barPx);
    saveAs("Tiff", outPath); rename("done");
    close("done"); close("CH");
}

// 二值遮罩面板：閾值 → 只留細胞 ROI 內 → 上色
function maskPanel(src, ch, thr, lut, roiIdx, x0, y0, s, outPath) {
    selectWindow(src); run("Select None");
    run("Duplicate...", "title=MK duplicate channels=" + ch);
    run("32-bit");
    run("Macro...", "code=v=(v>=" + thr + ")*255");
    setMinAndMax(0, 255);
    run("8-bit");
    roiManager("Select", roiIdx);
    setBackgroundColor(0, 0, 0);
    run("Clear Outside", "slice");
    run("Select None");
    run(lut); setMinAndMax(0, 255); run("RGB Color");
    makeRectangle(x0, y0, s, s);
    run("Duplicate...", "title=MKC");
    run("Select None");
    saveAs("Tiff", outPath); rename("done");
    close("done"); close("MK");
}

// 重疊遮罩：兩通道皆陽性的像素，白色
function overlapPanel(src, roiIdx, x0, y0, s, outPath) {
    selectWindow(src); run("Select None");
    run("Duplicate...", "title=OG duplicate channels=2"); run("32-bit");
    run("Macro...", "code=v=(v>=" + TG + ")*1");
    selectWindow(src);
    run("Duplicate...", "title=OR duplicate channels=3"); run("32-bit");
    run("Macro...", "code=v=(v>=" + TR + ")*1");
    imageCalculator("Multiply create 32-bit", "OG", "OR"); rename("OVL");
    run("Multiply...", "value=255");
    setMinAndMax(0, 255); run("8-bit");
    roiManager("Select", roiIdx);
    setBackgroundColor(0, 0, 0);
    run("Clear Outside", "slice");
    run("Select None");
    run("Grays"); setMinAndMax(0, 255); run("RGB Color");
    makeRectangle(x0, y0, s, s);
    run("Duplicate...", "title=OVC");
    run("Select None");
    saveAs("Tiff", outPath); rename("done");
    close("done"); close("OG"); close("OR"); close("OVL");
}

// merge 面板（可選：畫出細胞輪廓）
function mergePanel(src, x0, y0, w, h, outPath, barPx, roiIdx) {
    selectWindow(src); run("Select None");
    run("Duplicate...", "title=mD duplicate channels=1"); run("Blue");    setMinAndMax(D1MIN, D1MAX);
    selectWindow(src); run("Duplicate...", "title=mG duplicate channels=2"); run("Green");   setMinAndMax(D2MIN, D2MAX);
    selectWindow(src); run("Duplicate...", "title=mL duplicate channels=3"); run("Magenta"); setMinAndMax(D3MIN, D3MAX);
    run("Merge Channels...", "c2=[mG] c3=[mD] c6=[mL] create");
    rename("MG"); Stack.setDisplayMode("composite"); run("RGB Color");
    if (roiIdx >= 0) {
        roiManager("Select", roiIdx);
        setForegroundColor(255, 255, 255);
        run("Line Width...", "line=2");
        run("Draw", "slice");
        run("Select None");
    }
    makeRectangle(x0, y0, w, h);
    run("Duplicate...", "title=MGC");
    run("Select None");
    drawBar(w, h, barPx);
    saveAs("Tiff", outPath); rename("done");
    close("done"); close("MG");
}

// 拼版：S-A 3x3
function assembleA(dir, labs, w, h, gut) {
    cols = newArray("BODIPY", "LC3B", "merge");
    totW = 3*w + 2*gut; totH = 3*h + 2*gut;
    close("*");
    newImage("SA", "RGB white", totW, totH, 1);
    for (r = 0; r < lengthOf(labs); r++) {
        for (c = 0; c < lengthOf(cols); c++) {
            fp = dir + "SA_" + labs[r] + "_" + cols[c] + ".tif";
            if (File.exists(fp)) {
                open(fp); rename("P"); run("Select All"); run("Copy"); close("P");
                selectWindow("SA");
                makeRectangle(c*(w+gut), r*(h+gut), w, h);
                run("Paste");
            }
        }
    }
    selectWindow("SA"); run("Select None");
    saveAs("Tiff", dir + "FigureS_controls.tif");
    print("拼版 -> " + dir + "FigureS_controls.tif  (" + totW + " x " + totH + " px)");
    close("*");
}

// 拼版：加強對比版（三列一律相同的加強範圍，僅供輔助檢視）
function assembleEnh(dir, labs, w, h, gut) {
    cols = newArray("BODIPY", "LC3B");
    totW = 2*w + gut; totH = 3*h + 2*gut;
    close("*");
    newImage("SAE", "RGB white", totW, totH, 1);
    for (r = 0; r < lengthOf(labs); r++) {
        for (c = 0; c < 2; c++) {
            fp = dir + "SAe_" + labs[r] + "_" + cols[c] + ".tif";
            if (File.exists(fp)) {
                open(fp); rename("P"); run("Select All"); run("Copy"); close("P");
                selectWindow("SAE");
                makeRectangle(c*(w+gut), r*(h+gut), w, h);
                run("Paste");
            }
        }
    }
    selectWindow("SAE"); run("Select None");
    saveAs("Tiff", dir + "FigureS_controls_enhanced.tif");
    print("拼版 -> " + dir + "FigureS_controls_enhanced.tif  (顯示上限 BODIPY " + D2MAX_ENH + " / LC3B " + D3MAX_ENH + "，三列相同)");
    close("*");
}

// 拼版：S-B 1x6
function assembleB(dir, s, gut) {
    names = newArray("SB_1_BODIPY_raw", "SB_2_BODIPY_mask", "SB_3_LC3B_raw",
                     "SB_4_LC3B_mask",  "SB_5_overlap",     "SB_6_merge");
    totW = 6*s + 5*gut;
    close("*");
    newImage("SB", "RGB white", totW, s, 1);
    for (c = 0; c < lengthOf(names); c++) {
        fp = dir + names[c] + ".tif";
        if (File.exists(fp)) {
            open(fp); rename("P"); run("Select All"); run("Copy"); close("P");
            selectWindow("SB");
            makeRectangle(c*(s+gut), 0, s, s);
            run("Paste");
        }
    }
    selectWindow("SB"); run("Select None");
    saveAs("Tiff", dir + "FigureS_thresholds.tif");
    print("拼版 -> " + dir + "FigureS_thresholds.tif  (" + totW + " x " + s + " px)");
    close("*");
}
