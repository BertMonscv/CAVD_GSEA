// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  PS6_percell.ijm  v4  —— p-S6 (Ser235/236) IF 逐細胞定量
//  與 Vimentin / α-SMA 同一口徑（核 Otsu → 5 µm 核周環 → 在 C1 量強度）
//
//  ⚠️ 第一次跑請保持 PROBE = true：只跑第一個 .lif 的前 3 個 series，
//     每一步都印到 Log 視窗。看完 Log 確認沒問題，再改成 false 跑全部。
//
//  v3 相對 v2 的修正（v2 產出空檔）：
//   1. setBatchMode 預設關閉 —— Bio-Formats Importer 在批次模式下常常
//      開了圖卻不進巨集的影像清單，導致每個 series 都被判為「開檔失敗」
//   2. 一律用 **影像 ID**（selectImage）而非視窗標題（selectWindow）
//   3. perfield.csv **每個 series 都會寫一列**，含 status 欄；就算全部跳過
//      也看得出是在哪一關跳的，不會再出現空檔
//
//  v3.1 新增：SAVE_ROI（每個視野存 RoiSet.zip）與 SAVE_OVERLAY（分割稽核疊圖）
//  v4  修正兩件事：
//   (a) percell/perfield 出空表 → **完全不再用字串緩衝**，每產生一列就直接
//       File.append 寫盤；並在開頭做一次「寫得進去嗎」的自我檢測
//   (b) 焦距指標改在 **DAPI(C0) 通道**上量 —— 舊的 focus 用 C1 的 fieldMode
//       當分母，會系統性懲罰訊號強的組（LEP 會被剔掉 57%），不可用於剔除
// ============================================================================

// ---------- 執行模式 ----------
var PROBE        = true;    // true = 只跑一點點並詳細印 Log
var PROBE_FILES  = 1;
var PROBE_SERIES = 3;
var BATCH        = false;   // 確認能跑之後可以改 true 加速

// ---------- 路徑 ----------
var DIR     = getDirectory("Select raw p-S6 LIF input folder");
var OUTNAME = "ps6_percell";
var RECURSE = true;

// ---------- 量測參數 ----------
var CH_DAPI  = 1;     // Bio-Formats 通道編號（1-based）
var CH_FITC  = 2;
var RING_UM  = 5;     // 核周環寬度 µm
var NUC_MIN  = 30;    // 核面積 µm² 下限
var NUC_MAX  = 600;   // 核面積 µm² 上限
var NUC_CIRC = 0.35;
var GAUSS_SIG = 2;
var FLOOR    = 430;   // 相機地板
var SATVAL   = 65535;

// ---------- 分割稽核輸出 ----------
var SAVE_ROI         = true;   // 每個視野存一份 RoiSet.zip 到 ~/ps6_percell/roi/
var SAVE_OVERLAY     = true;   // 輸出「DAPI + 核輪廓 + 環外緣」疊圖到 ~/ps6_percell/overlay/
var OVERLAY_PER_FILE = 2;      // 每個 .lif 輸出幾張疊圖；-1 = 全部視野都出
var OVERLAY_LABEL    = true;   // 在每個細胞旁標上編號（對得回 percell.csv 的 cell 欄）
var DIR_ROI, DIR_OVL;

// 掃描時查出的重複／誤置影像
var EXCLUDE = newArray(
    "PS6-Sac.lif:Image 1",   // ≡ LEP1.lif 的 Image 7
    "LEP3.lif:Image 7"       // ≡ LEP3.lif 的 Image 20
);

// ---------- 內部狀態（IJ 宏一次只能開一個檔案控制代碼 → 緩衝 + File.append）----------
var PATH_PER, PATH_FLD, PATH_LOG, HDR_PER, HDR_FLD;
var bufLog = "";
var nRowPer = 0;
var nRowFld = 0;

// ============================== 主體 ==============================
print("\\Clear");
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "area mean min median redirect=None decimal=3");
if (BATCH) setBatchMode(true);

out = getDirectory("Select empty p-S6 quantification output folder");
File.makeDirectory(out);
if (!File.exists(out)) exit("無法建立輸出資料夾：" + out);

DIR_ROI = out + "roi/";
DIR_OVL = out + "overlay/";
if (SAVE_ROI)     File.makeDirectory(DIR_ROI);
if (SAVE_OVERLAY) File.makeDirectory(DIR_OVL);

PATH_PER = out + "percell.csv";
PATH_FLD = out + "perfield.csv";
PATH_LOG = out + "run_log.txt";
HDR_PER = "file,group,rep,series,seriesName,cell,ringMean,ringMedian,ringMax,ringArea_um2,nucArea_um2,fieldMode,net,pxSize_um";
HDR_FLD = "file,group,rep,series,seriesName,status,nCells,ringMedian_field,fieldMode,fieldMedian,fieldMin,fieldMax,satPct,net_field,dapiNV,dapiEdge,dapiNucEdge,focusC1_legacy,expC0_s,expC1_s,timestamp,posX,posY,pxSize_um";
File.saveString(HDR_PER + "\n", PATH_PER);
File.saveString(HDR_FLD + "\n", PATH_FLD);
File.saveString("PS6_percell.ijm v4   " + gettime() + "\n", PATH_LOG);

// --- 自我檢測：確認這兩個路徑真的寫得進去（上一版出空表就是卡在這裡）---
File.append("__SELFTEST__", PATH_PER);
tLen = lengthOf(File.openAsString(PATH_PER));
logLine("自我檢測：percell.csv 寫入後長度 = " + tLen + " 位元組（應 > " + lengthOf(HDR_PER) + "）");
if (tLen <= lengthOf(HDR_PER) + 2) logLine("!! File.append 對此路徑無效 —— 請把輸出改到別的資料夾，或確認 CSV 沒有被 Excel/Numbers 開著");
File.saveString(HDR_PER + "\n", PATH_PER);   // 清掉測試行

logLine("PROBE=" + PROBE + "  BATCH=" + BATCH + "  RING_UM=" + RING_UM + "  FLOOR=" + FLOOR);
logLine("輸出資料夾 = " + out);
logLine("DIR = " + DIR);
logLine("DIR 是資料夾嗎？ " + File.isDirectory(DIR));

files = listLif(DIR);
logLine("找到 .lif 檔數 = " + files.length);
if (files.length == 0) logLine("!! DIR 裡沒有 .lif —— 請檢查路徑（結尾要有斜線）");

nFile = files.length;
if (PROBE && nFile > PROBE_FILES) nFile = PROBE_FILES;

nTot = 0; nCellTot = 0;
for (f = 0; f < nFile; f++) {
    path = files[f];
    base = File.getName(path);
    grp  = groupOf(base);
    rep  = repOf(base);
    Ext.setId(path);
    Ext.getSeriesCount(nSer);
    logLine("");
    logLine("[檔案] " + base + "   group=" + grp + " rep=" + rep + "   series 總數=" + nSer);

    stem = base;
    if (endsWith(toLowerCase(stem), ".lif")) stem = substring(stem, 0, lengthOf(stem) - 4);
    ovlDone = 0;

    nRun = nSer;
    if (PROBE && nRun > PROBE_SERIES) nRun = PROBE_SERIES;

    for (s = 0; s < nRun; s++) {
        Ext.setSeries(s);
        Ext.getSeriesName(sName);

        if (isProcessed(sName)) { fldSkip(base, grp, rep, s, sName, "svcc");     logLine("  s" + s + " [" + sName + "] → SVCC，跳過"); continue; }
        if (isExcluded(base, sName)) { fldSkip(base, grp, rep, s, sName, "excluded"); logLine("  s" + s + " [" + sName + "] → 在 EXCLUDE 名單，跳過"); continue; }

        run("Close All"); roiManager("reset");
        run("Bio-Formats Importer",
            "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (s+1));
        logLine("  s" + s + " [" + sName + "] 開檔後 nImages=" + nImages);
        if (nImages == 0) { fldSkip(base, grp, rep, s, sName, "open_fail"); continue; }

        idSrc = getImageID();
        getDimensions(iw, ih, nC, nZ, nT);
        bd = bitDepth();
        getPixelSize(unit, pw, ph);
        logLine("       " + iw + "x" + ih + "  nC=" + nC + " nZ=" + nZ + "  " + bd + "-bit  px=" + d4(pw) + " " + unit);

        if (bd != 16) { fldSkip(base, grp, rep, s, sName, "not16bit"); continue; }
        if (nC < 2)   { fldSkip(base, grp, rep, s, sName, "nC<2");     continue; }
        if (pw <= 0) { pw = 1; logLine("       [warn] 無空間標定，面積以像素計"); }

        info = getMetadata("Info");
        e0 = grab(info, "C=0", "Exposure");
        e1 = grab(info, "C=1", "Exposure");
        ts = grab(info, "", "Timestamp");
        px_ = grab(info, "", "PosX");
        py_ = grab(info, "", "PosY");

        // --- 拆通道（用 ID，不靠視窗名）---
        selectImage(idSrc);
        run("Duplicate...", "title=D16 duplicate channels=" + CH_DAPI);
        idD = getImageID();
        selectImage(idSrc);
        run("Duplicate...", "title=G16 duplicate channels=" + CH_FITC);
        idG = getImageID();

        // --- 視野層統計（C1，無選區）---
        selectImage(idG); run("Select None");
        getRawStatistics(nPix, gMean, gMin, gMax, gStd, gHist);
        // 16-bit 的直方圖長度是 max+1，不是 65536 → 索引前必須檢查長度
        fMode = argmax(gHist);
        fMed  = percentileFromHist(gHist, nPix, 0.50);
        satP  = 0;
        if (gHist.length > SATVAL) satP = 100.0 * gHist[SATVAL] / nPix;
        logLine("       C1: min=" + gMin + " mode=" + fMode + " med=" + fMed + " max=" + gMax);

        // --- 焦距指標：在 DAPI(C0) 上量，與 p-S6 讀出通道完全無關 ---
        selectImage(idD); run("Select None");
        getRawStatistics(dPix, dMean, dMin, dMax, dStd);
        dapiNV = 0;                                   // 正規化變異數（經典自動對焦指標）
        if (dMean > 0) dapiNV = (dStd * dStd) / (dMean * dMean);
        selectImage(idD);
        run("Duplicate...", "title=DEDGE");
        idE = getImageID();
        run("32-bit");                                // 先轉 32-bit，避免 Sobel 在 16-bit 被截頂
        run("Find Edges");
        run("Select None");
        getRawStatistics(ePix, eMean, eMin, eMax, eStd);
        dapiEdge = 0;                                 // 全場梯度強度 ÷ 全場亮度
        if (dMean > 0) dapiEdge = eMean / dMean;
        logLine("       DAPI 焦距：NV=" + d4(dapiNV) + "  edge=" + d4(dapiEdge));

        // --- 核分割（DAPI 曝光在檔間不一致 → 逐視野自動閾值）---
        selectImage(idD);
        run("Duplicate...", "title=DMASK");
        idM = getImageID();
        run("Gaussian Blur...", "sigma=" + GAUSS_SIG);
        setAutoThreshold("Otsu dark");
        getThreshold(tLo, tHi);
        run("Convert to Mask");
        run("Fill Holes");
        run("Watershed");
        run("Analyze Particles...", "size=" + NUC_MIN + "-" + NUC_MAX + " circularity=" + NUC_CIRC + "-1.00 exclude add");
        nNuc = roiManager("count");
        logLine("       DAPI Otsu 閾值=" + tLo + "  分出核數=" + nNuc);

        // --- 建環並量測 ---
        for (k = 0; k < nNuc; k++) {
            selectImage(idG);
            roiManager("select", k);
            run("Enlarge...", "enlarge=" + RING_UM);
            roiManager("add");
        }
        ringMeds = newArray(nNuc);
        for (k = 0; k < nNuc; k++) {
            selectImage(idG);
            roiManager("select", k);
            nucA = getValue("Area");
            sel = newArray(2); sel[0] = k; sel[1] = nNuc + k;
            roiManager("select", sel);
            roiManager("XOR");
            rMean = getValue("Mean");
            rMed  = getValue("Median");
            rMax  = getValue("Max");
            rArea = getValue("Area");
            ringMeds[k] = rMed;
            File.append(base + "," + grp + "," + rep + "," + s + "," + sName + "," + (k+1) + ","
                  + d3(rMean) + "," + d3(rMed) + "," + d3(rMax) + "," + d3(rArea) + "," + d3(nucA) + ","
                  + fMode + "," + d3(rMean - fMode) + "," + d4(pw), PATH_PER);
            nRowPer++;
            nCellTot++;
        }

        // --- 逐核銳利度：核內 Sobel 均值 ÷ 核內 DAPI 均值（不隨細胞密度變動）---
        dapiNucEdge = 0;
        if (nNuc > 0) {
            nucEdge = newArray(nNuc);
            for (k = 0; k < nNuc; k++) {
                selectImage(idE); roiManager("select", k); em = getValue("Mean");
                selectImage(idD); roiManager("select", k); dm = getValue("Mean");
                if (dm > 0) nucEdge[k] = em / dm;
            }
            dapiNucEdge = median(nucEdge);
        }

        rmField = 0;
        if (nNuc > 0) rmField = median(ringMeds);
        denom = fMode - FLOOR; if (denom <= 0) denom = 1;
        focusC1 = (rmField - fMode) / denom;   // 舊指標，僅留作對照，**不得用於剔除視野**
        if (nNuc > 0) logLine("       ring 中位=" + d3(rmField) + "  net=" + d3(rmField - fMode) + "  nucEdge=" + d4(dapiNucEdge));

        // ---------- 存 ROI（核 0..n-1、外擴 n..2n-1、環 2n..3n-1）----------
        if (SAVE_ROI && nNuc > 0) {
            for (k = 0; k < nNuc; k++) {
                selectImage(idG);
                sel = newArray(2); sel[0] = k; sel[1] = nNuc + k;
                roiManager("select", sel);
                roiManager("XOR");
                roiManager("add");
            }
            for (k = 0; k < nNuc; k++) {
                roiManager("select", k);              roiManager("rename", "nuc_"  + IJ.pad(k+1, 3));
                roiManager("select", nNuc + k);       roiManager("rename", "out_"  + IJ.pad(k+1, 3));
                roiManager("select", 2 * nNuc + k);   roiManager("rename", "ring_" + IJ.pad(k+1, 3));
            }
            roiManager("deselect");
            roiManager("save", DIR_ROI + stem + "_s" + IJ.pad(s, 3) + "_RoiSet.zip");
        }

        // ---------- 分割稽核疊圖：DAPI 底圖 + 核輪廓(青) + 環外緣(黃) ----------
        wantOvl = false;
        if (SAVE_OVERLAY && nNuc > 0) {
            if (OVERLAY_PER_FILE < 0) wantOvl = true;
            if (ovlDone < OVERLAY_PER_FILE) wantOvl = true;
        }
        if (wantOvl) {
            selectImage(idD);
            run("Duplicate...", "title=OVL");
            idO = getImageID();
            run("Enhance Contrast", "saturated=0.35");   // 僅供目視，非定量
            run("RGB Color");
            run("Line Width...", "line=2");
            setFont("SansSerif", 16, "bold antialiased");
            for (k = 0; k < nNuc; k++) {
                selectImage(idO);
                setForegroundColor(0, 255, 255);            // 核 = 青
                roiManager("select", k);
                run("Draw", "slice");
                setForegroundColor(255, 255, 0);            // 環外緣 = 黃
                roiManager("select", nNuc + k);
                run("Draw", "slice");
                if (OVERLAY_LABEL) {
                    roiManager("select", k);
                    Roi.getBounds(bx, by, bw, bh);
                    selectImage(idO);
                    run("Select None");
                    setForegroundColor(255, 255, 255);
                    drawString("" + (k+1), bx, by);
                }
            }
            selectImage(idO);
            run("Select None");
            saveAs("PNG", DIR_OVL + stem + "_s" + IJ.pad(s, 3) + "_seg.png");
            ovlDone++;
        }

        File.append(base + "," + grp + "," + rep + "," + s + "," + sName + ",ok," + nNuc + ","
              + d3(rmField) + "," + fMode + "," + fMed + "," + gMin + "," + gMax + "," + d4(satP) + ","
              + d3(rmField - fMode) + ","
              + d4(dapiNV) + "," + d4(dapiEdge) + "," + d4(dapiNucEdge) + "," + d3(focusC1) + ","
              + e0 + "," + e1 + "," + ts + "," + px_ + "," + py_ + "," + d4(pw), PATH_FLD);
        nRowFld++;
        nTot++;
        run("Close All"); roiManager("reset");
    }
    Ext.close();
    flush();
}

if (BATCH) setBatchMode(false);
logLine("");
logLine("完成：量測 " + nTot + " 個視野、" + nCellTot + " 個細胞");
logLine("寫出 percell 列數 = " + nRowPer + "（檔案 " + File.length(PATH_PER) + " 位元組）");
logLine("寫出 perfield 列數 = " + nRowFld + "（檔案 " + File.length(PATH_FLD) + " 位元組）");
if (PROBE) logLine("※ 這是 PROBE 模式的結果。確認上面沒問題後，把 PROBE 改成 false 再跑一次全部。");
flush();
print("完成，輸出在：" + out);

// ============================== 函式 ==============================

function logLine(t) {
    print(t);
    bufLog = bufLog + t + "\n";
}

function flush() {
    if (lengthOf(bufLog) > 0) { File.append(trimNL(bufLog), PATH_LOG); bufLog = ""; }
}

function trimNL(t) {
    while (endsWith(t, "\n")) { t = substring(t, 0, lengthOf(t) - 1); }
    return "" + t;
}

// 任何跳過的 series 也寫一列，這樣 perfield.csv 永遠不會是空的
function fldSkip(fn, g, r, si, sn, st) {
    File.append(fn + "," + g + "," + r + "," + si + "," + sn + "," + st
          + ",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA", PATH_FLD);
    nRowFld++;
}

function listLif(dir) {
    res = newArray(0);
    if (!File.isDirectory(dir)) return res;
    l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        p = dir + l[i];
        if (File.isDirectory(p)) {
            if (RECURSE) { sub = listLif(p); res = Array.concat(res, sub); }
        } else {
            if (endsWith(toLowerCase(l[i]), ".lif")) res = Array.concat(res, newArray(p));
        }
    }
    return res;
}

function isProcessed(nm) {
    u = toUpperCase(nm);
    if (indexOf(u, "_LNG_") >= 0)     return true;
    if (indexOf(u, "SVCC") >= 0)      return true;
    if (indexOf(u, "LIGHTNING") >= 0) return true;
    return false;
}

function isExcluded(fn, sn) {
    tag = fn + ":" + sn;
    for (i = 0; i < EXCLUDE.length; i++) { if (EXCLUDE[i] == tag) return true; }
    return false;
}

function groupOf(nm) {
    u = toUpperCase(nm);
    if (indexOf(u, "SAC")  >= 0) return "Sac";
    if (indexOf(u, "RAPA") >= 0) return "RAPA";
    if (indexOf(u, "LEP")  >= 0) return "LEP";
    if (indexOf(u, "PM")   >= 0) return "PM";
    if (indexOf(u, "CTL")  >= 0) return "CTL";
    return "NA";
}

function repOf(nm) {
    for (i = 1; i <= 9; i++) { if (indexOf(nm, "" + i) >= 0) return "" + i; }
    return "NA";
}

function grab(info, tag, key) {
    lines = split(info, "\n");
    for (i = 0; i < lines.length; i++) {
        L = lines[i];
        if (indexOf(L, key) < 0) continue;
        if (tag != "" && indexOf(L, tag) < 0) continue;
        p = indexOf(L, "=");
        if (p < 0) continue;
        v = String.trim(substring(L, p + 1));
        v = replace(v, ",", ";");
        if (lengthOf(v) > 0) return "" + v;
    }
    return "NA";
}

function argmax(h) {
    bi = 0; bv = -1;
    for (i = 0; i < h.length; i++) { if (h[i] > bv) { bv = h[i]; bi = i; } }
    return bi;
}

function percentileFromHist(h, nPix, q) {
    target = nPix * q; c = 0;
    for (i = 0; i < h.length; i++) { c += h[i]; if (c >= target) return i; }
    return h.length - 1;
}

function median(a) {
    if (a.length == 0) return 0;
    b = Array.copy(a); Array.sort(b); n = b.length;
    if (n % 2 == 1) return b[(n - 1) / 2];
    return (b[n / 2 - 1] + b[n / 2]) / 2;
}

function d3(x) { return "" + d2s(x, 3); }
function d4(x) { return "" + d2s(x, 4); }

function gettime() {
    getDateAndTime(y, mo, dw, dd, hh, mm, ss, ms);
    return "" + y + "-" + IJ.pad(mo + 1, 2) + "-" + IJ.pad(dd, 2) + " " + IJ.pad(hh, 2) + ":" + IJ.pad(mm, 2);
}
