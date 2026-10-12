// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  COLOC_C_representative.ijm  v4  —— Figure 5F 代表圖擷取
//
//  v4 變更：
//   · 裁切為 111x60 µm 橫框（三個面板同尺寸，含 Rapa）。PM 的框對準左下角
//     那兩個細胞，橫框剛好容得下。插圖仍為 20 µm 正方形。
//   · 新增「直接指定裁切中心」CX_PX / CY_PX（-1 = 用該 ROI 的形心）。
//     PM 已填入左下角兩個細胞的中心 (308, 857)。
//   · 位移參數 OFF_*_UM 移除，由 CX/CY 取代（更直接、不會被邊界夾住而失效）。
//
//  來源：~/coloc_out/blind/F0xx.tif ＋ ~/coloc_out/roi/F0xx.zip
// ============================================================================

// ---- 顯示範圍（所有組別一律相同；函式要用到，必須是 var 全域）----
// ---- 比例尺與插圖方框的樣式（全部自己畫，不用 Scale Bar 外掛）----
var PX         = 0.1802;  // µm/px（函式要用，必須是 var）
var BAR_H      = 8;       // 比例尺厚度（px）
var BAR_MARGIN = 18;      // 比例尺距右下角的邊距（px）
var BAR_LABEL  = true;    // 是否在比例尺上方寫出長度
var BAR_FONT   = 24;      // 標籤字級（px）
var BOX_LW     = 3;       // 插圖方框線寬（px）

var D1MIN = 0;  var D1MAX = 150;   // DAPI
var D2MIN = 0;  var D2MAX = 200;   // BODIPY
var D3MIN = 0;  var D3MAX = 130;   // LC3B

OUT    = getDirectory("Select colocalization working folder");
BLIND  = OUT + "blind/";
ROIDIR = OUT + "roi/";
BASE   = OUT + "figure5F/";

// ---- 要輸出哪幾組 ----
SETS = newArray("D4");
//  A  F057 PM1 / F026 LEP2 / F040 Rapa1      分數最佳
//  B  F053 PM2 / F044 LEP2 / F043 Rapa3      次佳
//  C  F053 PM2 / F044 LEP2 / F058 Rapa2      全 rep2
//  D  F057 PM1 / F017 LEP1 / F040 Rapa1      全 rep1
//  D2 F053 PM2 / F017 LEP1 / F003 Rapa2      PM 與 Rapa 換視野
//  D3 F053 PM2 / F017 LEP1 / F058 Rapa2      Rapa 再換視野
//  D4 F053 PM2 / F017 LEP1 / F004 Rapa1      Rapa 再換（現用）
//     F058 roi6 是全 Rapa 組最接近組均值的細胞（z=0.13，次佳 0.41）
//  E  F035 PM3 / F014 LEP3 / F043 Rapa3      全 rep3
//  F  F053 PM2 / F033 LEP3 / F058 Rapa2      每組 6 個細胞
//  M  自訂（見下方 M 區塊）

// ---- 裁切尺寸：橫框，三個面板一律相同（µm）----
//  三個面板必須同尺寸，否則同寬排版時放大倍率不一致、比例尺長度也會不同。
//  故 Rapa 也是 111x60，不能只給 PM 與 LEP 用橫框。
CROP_W_UM = 111;     // 寬 → 616 px
CROP_H_UM = 60;      // 高 → 333 px
INSET_UM  = 30;      // 插圖維持正方形（30 µm → 166 px，拼版時放大 2x 到列高）
BAR_MAIN  = 20;      // 主圖比例尺 µm
BAR_INS   = 10;      // 插圖比例尺 µm

// ---- 裁切中心（視野像素座標；-1 = 用該 ROI 的形心）----
//      順序 = PM / PM+Lep / PM+Lep+Rapa
//  PM 對準左下角那兩個細胞（核在 (204,786) 與 (382,846)）。
//  111x60 的橫框中心 (308,857) → x 0-616 / y 690-1023，兩個細胞都能完整入鏡，
//  且 (353,464) 那個細胞落在框外。
//  LEP 移到視野最右下角：616x333 的框貼齊右下 → 左上角 (408,691)、中心 (716,857)
CX_PX = newArray(308, 716, -1);
CY_PX = newArray(857, 857, -1);

// ---- 在上述中心之上再位移（視野像素；X 負向左、Y 負向上）----
//      LEP ：上移 1.5 個框高 = -500 px；左右累計淨 -226 px → 中心 (490, 357)
//      Rapa：上移一個框高   = -333 px；左移 1/4 框寬 = -154 px
OFF_X_PX = newArray(0, -226, -154);
OFF_Y_PX = newArray(0, -500, -333);

// ---- 模式開關 ----
ASSEMBLE       = true;    // 拼成一張 4 欄 x 3 列的完成圖（不燒文字標籤，標籤留給排版軟體）
GUTTER         = 14;      // 面板間距（px）
MAP            = false;    // 額外輸出整個視野的定位圖（細胞輪廓＋編號＋裁切框）
DUMP_CENTROIDS = false;   // 輸出全部 ROI 的形心座標到 roi_centroids.csv（跑過一次即可）
PREVIEW        = false;   // 只掃綠通道亮度
D2MAX_TEST     = newArray(140, 170, 200, 240);

TG = 100; TR = 90;   // 與定量相同的閾值，僅用於決定插圖中心

// ---------------------------------------------------------------------------
setBatchMode(true);
File.makeDirectory(BASE);
print("\\Clear");

if (DUMP_CENTROIDS) {
    dcsv = OUT + "roi_centroids.csv";
    File.saveString("code,roi,bx,by,bw,bh,cx,cy,area_px\n", dcsv);
    rl = getFileList(ROIDIR);
    for (q = 0; q < lengthOf(rl); q++) {
        fn = rl[q];
        if (endsWith(fn, ".zip") || endsWith(fn, ".roi")) {
            cd = substring(fn, 0, indexOf(fn, "."));
            if (File.exists(BLIND + cd + ".tif")) {
                close("*");
                open(BLIND + cd + ".tif");
                roiManager("reset");
                roiManager("Open", ROIDIR + fn);
                for (k = 0; k < roiManager("count"); k++) {
                    roiManager("Select", k);
                    Roi.getBounds(qbx, qby, qbw, qbh);
                    getRawStatistics(qn);
                    File.append(cd + "," + (k+1) + "," + qbx + "," + qby + "," + qbw + "," + qbh + "," +
                                d2s(qbx + qbw/2, 1) + "," + d2s(qby + qbh/2, 1) + "," + qn, dcsv);
                }
                roiManager("reset"); close("*");
            }
        }
    }
    print("ROI 形心已匯出 -> " + dcsv);
}

cw = round(CROP_W_UM / PX);
ch = round(CROP_H_UM / PX);
insetPx = round(INSET_UM / PX);
print("裁切尺寸 " + cw + " x " + ch + " px  (" + CROP_W_UM + " x " + CROP_H_UM + " µm)");

for (si = 0; si < lengthOf(SETS); si++) {
    SET = SETS[si];
    if (SET == "A")       { codes = newArray("F057","F026","F040"); rois = newArray(1,2,4); srcs = newArray("PM1 s1","LEP2 s4","Rapa1 s4"); }
    else if (SET == "B")  { codes = newArray("F053","F044","F043"); rois = newArray(3,6,2); srcs = newArray("PM2 s3","LEP2 s2","Rapa3 s4"); }
    else if (SET == "C")  { codes = newArray("F053","F044","F058"); rois = newArray(3,6,6); srcs = newArray("PM2 s3","LEP2 s2","Rapa2 s4"); }
    else if (SET == "D")  { codes = newArray("F057","F017","F040"); rois = newArray(1,4,2); srcs = newArray("PM1 s1","LEP1 s1","Rapa1 s4"); }
    else if (SET == "D2") { codes = newArray("F053","F017","F003"); rois = newArray(3,4,2); srcs = newArray("PM2 s3","LEP1 s1","Rapa2 s0"); }
    else if (SET == "D3") { codes = newArray("F053","F017","F058"); rois = newArray(3,4,6); srcs = newArray("PM2 s3","LEP1 s1","Rapa2 s4"); }
    else if (SET == "D4") { codes = newArray("F053","F017","F004"); rois = newArray(3,4,4); srcs = newArray("PM2 s3","LEP1 s1","Rapa1 s0"); }
    else if (SET == "E")  { codes = newArray("F035","F014","F043"); rois = newArray(1,3,2); srcs = newArray("PM3 s0","LEP3 s4","Rapa3 s4"); }
    else if (SET == "F")  { codes = newArray("F053","F033","F058"); rois = newArray(3,6,6); srcs = newArray("PM2 s3","LEP3 s0","Rapa2 s4"); }
    else                  { codes = newArray("F057","F026","F040"); rois = newArray(1,2,4); srcs = newArray("manual","manual","manual"); }
    labels = newArray("PM", "PM_Lep", "PM_Lep_Rapa");

    FIG = BASE + "SET_" + SET + "/";
    File.makeDirectory(FIG);
    PREV = FIG + "preview/";
    if (PREVIEW) File.makeDirectory(PREV);
    else File.saveString("set,panel,code,source,roi,cx_px,cy_px,crop_w,crop_h,inset_cx,inset_cy,D1MAX,D2MAX,D3MAX\n", FIG + "manifest.csv");

    print("===== SET " + SET + " =====");

    for (i = 0; i < lengthOf(labels); i++) {
        lab = labels[i]; code = codes[i]; ri = rois[i] - 1;

        rf = "";
        if (File.exists(ROIDIR + code + ".zip")) rf = ROIDIR + code + ".zip";
        if (rf == "" && File.exists(ROIDIR + code + ".roi")) rf = ROIDIR + code + ".roi";

        if (!File.exists(BLIND + code + ".tif") || rf == "") {
            print("MISSING: " + code);
        } else {
            close("*");
            open(BLIND + code + ".tif");
            rename("SRC");
            run("Properties...", "channels=3 slices=1 frames=1 pixel_width=" + PX + " pixel_height=" + PX + " voxel_depth=1 unit=um");

            roiManager("reset");
            roiManager("Open", rf);
            roiManager("Select", ri);
            Roi.getBounds(bx, by, bw, bh);
            cx = bx + bw/2; cy = by + bh/2;
            if (CX_PX[i] >= 0) cx = CX_PX[i];
            if (CY_PX[i] >= 0) cy = CY_PX[i];
            cx = cx + OFF_X_PX[i];
            cy = cy + OFF_Y_PX[i];

            // ---- 主框位置（與 makeCrop 相同的夾邊規則）----
            x0 = round(cx - cw / 2); y0 = round(cy - ch / 2);
            if (x0 < 0) x0 = 0;
            if (y0 < 0) y0 = 0;
            if (x0 + cw > getWidth())  x0 = getWidth()  - cw;
            if (y0 + ch > getHeight()) y0 = getHeight() - ch;

            // ---- 插圖中心：在「主框之內」找 LC3B∩BODIPY 重疊訊號最強處 ----
            //  必須限制在主框內，否則插圖會放大主圖上根本看不到的區域。
            run("Select None");
            run("Duplicate...", "title=Gm duplicate channels=2"); run("32-bit");
            run("Macro...", "code=v=(v>=" + TG + ")*1");
            selectWindow("SRC");
            run("Duplicate...", "title=Rm duplicate channels=3"); run("32-bit");
            run("Macro...", "code=v=(v>=" + TR + ")*1");
            imageCalculator("Multiply create 32-bit", "Gm", "Rm"); rename("OV");
            run("Gaussian Blur...", "sigma=4");

            half = round(insetPx / 2);
            sx1 = x0 + half; sy1 = y0 + half;
            sw2 = cw - insetPx; sh2 = ch - insetPx;
            if (sw2 < 1) sw2 = 1;
            if (sh2 < 1) sh2 = 1;
            makeRectangle(sx1, sy1, sw2, sh2);
            getRawStatistics(np, mn, mi, mx);
            icx = sx1 + sw2 / 2; icy = sy1 + sh2 / 2;
            if (mx > 0) {
                found = false;
                for (yy = sy1; yy < sy1 + sh2 && !found; yy++) {
                    for (xx = sx1; xx < sx1 + sw2 && !found; xx++) {
                        if (getPixel(xx, yy) >= mx * 0.999) { icx = xx; icy = yy; found = true; }
                    }
                }
            }
            close("Gm"); close("Rm"); close("OV");
            selectWindow("SRC"); run("Select None");

            if (MAP) {
                selectWindow("SRC"); run("Select None");
                run("Duplicate...", "title=mapD duplicate channels=1"); run("Blue");    setMinAndMax(D1MIN, D1MAX);
                selectWindow("SRC"); run("Duplicate...", "title=mapG duplicate channels=2"); run("Green");   setMinAndMax(D2MIN, D2MAX);
                selectWindow("SRC"); run("Duplicate...", "title=mapL duplicate channels=3"); run("Magenta"); setMinAndMax(D3MIN, D3MAX);
                run("Merge Channels...", "c2=[mapG] c3=[mapD] c6=[mapL] create");
                rename("MAPIMG"); Stack.setDisplayMode("composite"); run("RGB Color");
                run("Labels...", "color=white font=24 show use draw");
                roiManager("Show All with labels");
                run("From ROI Manager");
                mx0 = round(cx - cw/2); my0 = round(cy - ch/2);
                if (mx0 < 0) mx0 = 0;
                if (my0 < 0) my0 = 0;
                if (mx0 + cw > getWidth())  mx0 = getWidth()  - cw;
                if (my0 + ch > getHeight()) my0 = getHeight() - ch;
                makeRectangle(mx0, my0, cw, ch);
                run("Add Selection...", "stroke=yellow width=4");
                run("Select None");
                run("Flatten");
                saveAs("Tiff", FIG + lab + "_MAP.tif");
                close("MAPIMG");
                print("  " + lab + " MAP: 框左上角 (" + mx0 + "," + my0 + ")  右下角 (" + (mx0+cw) + "," + (my0+ch) + ")  " + cw + "x" + ch);
            }

            if (PREVIEW) {
                keep = D2MAX;
                for (k = 0; k < lengthOf(D2MAX_TEST); k++) {
                    D2MAX = D2MAX_TEST[k];
                    makeCrop("SRC", lab, cx, cy, cw, ch, round(BAR_MAIN/PX), PREV, "D2max" + D2MAX_TEST[k] + "_", true, 0, 0, 0, 1.0);
                }
                D2MAX = keep;
                print("  " + lab + " : preview x" + lengthOf(D2MAX_TEST));
            } else {
                makeCrop("SRC", lab, cx,  cy,  cw,      ch,      round(BAR_MAIN/PX), FIG, "",       false, icx, icy, insetPx, 1.0);
                makeCrop("SRC", lab, icx, icy, insetPx, insetPx, round(BAR_INS/PX),  FIG, "inset_", false, 0, 0, 0, 0.5);
                File.append(SET + "," + lab + "," + code + "," + srcs[i] + "," + rois[i] + "," +
                            d2s(cx,0) + "," + d2s(cy,0) + "," + cw + "," + ch + "," +
                            d2s(icx,0) + "," + d2s(icy,0) + "," +
                            D1MAX + "," + D2MAX + "," + D3MAX, FIG + "manifest.csv");
                print("  " + lab + " : " + code + " roi " + rois[i] + "  centre(" + d2s(cx,0) + "," + d2s(cy,0) + ")");
            }

            roiManager("reset");
            close("*");
        }
    }
}

// ---- 拼版 ----
if (ASSEMBLE && !PREVIEW) {
    for (si = 0; si < lengthOf(SETS); si++) {
        SET = SETS[si];
        FIG = BASE + "SET_" + SET + "/";
        labels = newArray("PM", "PM_Lep", "PM_Lep_Rapa");
        cols   = newArray("BODIPY", "LC3B", "merge");
        totW = 3 * cw + ch + 3 * GUTTER;
        totH = 3 * ch + 2 * GUTTER;
        close("*");
        newImage("FIG5F", "RGB white", totW, totH, 1);
        okAll = true;
        for (r = 0; r < 3; r++) {
            yy = r * (ch + GUTTER);
            for (cq = 0; cq < 3; cq++) {
                fpath = FIG + labels[r] + "_" + cols[cq] + ".tif";
                if (!File.exists(fpath)) { okAll = false; }
                else {
                    open(fpath); rename("PANEL");
                    run("Select All"); run("Copy"); close("PANEL");
                    selectWindow("FIG5F");
                    makeRectangle(cq * (cw + GUTTER), yy, cw, ch);
                    run("Paste");
                }
            }
            ipath = FIG + labels[r] + "_inset_merge.tif";
            if (!File.exists(ipath)) { okAll = false; }
            else {
                open(ipath); rename("INS");
                run("Size...", "width=" + ch + " height=" + ch + " average interpolation=Bilinear");
                run("Select All"); run("Copy"); close("INS");
                selectWindow("FIG5F");
                makeRectangle(3 * (cw + GUTTER), yy, ch, ch);
                run("Paste");
            }
        }
        selectWindow("FIG5F"); run("Select None");
        saveAs("Tiff", FIG + "Figure5F_assembled.tif");
        close("*");
        if (okAll) print("拼版完成 -> " + FIG + "Figure5F_assembled.tif  (" + totW + " x " + totH + " px)");
        else print("拼版完成但有面板缺檔 -> " + FIG + "Figure5F_assembled.tif");
    }
}

setBatchMode(false);
print("DONE -> " + BASE);
File.saveString(getInfo("log"), BASE + "representative_log.txt");

// ===========================================================================
function makeCrop(src, lab, cx, cy, sw, sh, barPx, dir, prefix, onlyMerge, boxX, boxY, boxS, styleScale) {
    x = round(cx - sw/2); y = round(cy - sh/2);
    selectWindow(src); getDimensions(W, H, nC, nZ, nT);
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x + sw > W) x = W - sw;
    if (y + sh > H) y = H - sh;

    selectWindow(src); run("Select None");
    makeRectangle(x, y, sw, sh);
    run("Duplicate...", "title=CROP duplicate");
    run("Select None");

    if (!onlyMerge) {
        selectWindow("CROP"); run("Duplicate...", "title=B duplicate channels=2");
        run("Green");   setMinAndMax(D2MIN, D2MAX); run("RGB Color");
        saveAs("Tiff", dir + lab + "_" + prefix + "BODIPY.tif"); rename("Brgb");

        selectWindow("CROP"); run("Duplicate...", "title=L duplicate channels=3");
        run("Magenta"); setMinAndMax(D3MIN, D3MAX); run("RGB Color");
        saveAs("Tiff", dir + lab + "_" + prefix + "LC3B.tif"); rename("Lrgb");
    }

    selectWindow("CROP"); run("Duplicate...", "title=mD duplicate channels=1");
    run("Blue");    setMinAndMax(D1MIN, D1MAX);
    selectWindow("CROP"); run("Duplicate...", "title=mG duplicate channels=2");
    run("Green");   setMinAndMax(D2MIN, D2MAX);
    selectWindow("CROP"); run("Duplicate...", "title=mL duplicate channels=3");
    run("Magenta"); setMinAndMax(D3MIN, D3MAX);
    run("Merge Channels...", "c2=[mG] c3=[mD] c6=[mL] create");
    rename("MRG");
    Stack.setDisplayMode("composite");
    run("RGB Color");
    //  styleScale：插圖在拼版時會被放大 2x，故其比例尺與方框先畫成一半，
    //  放大後才會與主圖的粗細、字級一致。
    bh = round(BAR_H * styleScale);      if (bh < 2) bh = 2;
    bm = round(BAR_MARGIN * styleScale); if (bm < 4) bm = 4;
    lw = round(BOX_LW * styleScale);     if (lw < 1) lw = 1;
    fs = round(BAR_FONT * styleScale);   if (fs < 8) fs = 8;

    setForegroundColor(255, 255, 255);
    if (boxS > 0) {
        makeRectangle(round(boxX - boxS / 2 - x), round(boxY - boxS / 2 - y), boxS, boxS);
        run("Line Width...", "line=" + lw);
        run("Draw", "slice");
        run("Select None");
    }
    // 比例尺：白色實心矩形 + 上方標籤，全部用像素繪製，不經 Scale Bar 外掛
    if (barPx > 0) {
        bx0 = sw - barPx - bm;
        by0 = sh - bh - bm;
        makeRectangle(bx0, by0, barPx, bh);
        run("Fill", "slice");
        run("Select None");
        if (BAR_LABEL) {
            setFont("SansSerif", fs, "bold antialiased");
            setColor(255, 255, 255);
            txt = d2s(barPx * PX, 0) + " " + fromCharCode(181) + "m";
            tw = getStringWidth(txt);
            tx = bx0 + barPx / 2 - tw / 2;
            if (tx < 2) tx = 2;
            drawString(txt, tx, by0 - round(fs * 0.25));
        }
    }
    saveAs("Tiff", dir + lab + "_" + prefix + "merge.tif");
    rename("MRGdone");

    close("CROP");
    if (!onlyMerge) { close("Brgb"); close("Lrgb"); }
    close("MRGdone");
}
