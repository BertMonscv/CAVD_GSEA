// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// nucFrac and fracG_inNuc are historical exploratory outputs not used in the manuscript.
// Their DAPI-based interpretation is not validated; do not use them as study endpoints.
// ============================================================================
//  COLOC_B_measure.ijm  v2  —— per-cell 定量（第三步，勾完輪廓後跑）
//
//  v2 修正（v1 的兩個錯誤）：
//   1. ImageJ 宏語言不支援三元運算子 ?: → 遮罩改用 (v>=T)*1
//   2. 不再依賴 Split Channels 產生的視窗名稱 → 改用 Duplicate ... channels=N
//      直接抽通道並自己命名；每張開始前先 close("*") 清乾淨殘留視窗
//
//  讀 blind/F0xx.tif + roi/F0xx.zip，逐細胞輸出：
//    fracG / fracR        BODIPY(閾100) 與 LC3B(閾90) 的閾上面積分數
//    meanD / meanG / meanR 三通道未閾值均值
//    M1                   BODIPY 訊號落在 LC3B 陽性像素的比例
//    M2                   LC3B   訊號落在 BODIPY 陽性像素的比例
//    pearson              未閾值全細胞區的 Pearson 相關
//    enrich               M2 ÷ fracG（隨機共分布時 = 1）
//    nucFrac              細胞 ROI 內 DAPI 核遮罩佔比
//    fracG_inNuc          閾上綠像素中落在核內的比例（405 激發疑慮排查用）
// ============================================================================

OUT    = getDirectory("Select colocalization working folder");
BLIND  = OUT + "blind/";
ROIDIR = OUT + "roi/";
CSV    = OUT + "percell_blind.csv";

TG = 100;              // BODIPY 閾值（8-bit）
TR = 90;               // LC3B 閾值（8-bit）
PXFALLBACK = 0.1802;   // µm/px，僅在 TIFF 無標定時使用

// ---------------------------------------------------------------------------
setBatchMode(true);
print("\\Clear");
close("*");
roiManager("reset");
File.saveString("code,roi,areaPx,areaUm2,meanD,meanG,meanR,fracG,fracR,M1,M2,pearson,enrich,nucFrac,fracG_inNuc\n", CSV);

list = getFileList(BLIND);
nDone = 0; nCells = 0;

for (i = 0; i < lengthOf(list); i++) {
    nm = list[i];
    if (endsWith(nm, ".tif")) {
        code = substring(nm, 0, lengthOf(nm) - 4);

        rf = "";
        if (File.exists(ROIDIR + code + ".zip")) rf = ROIDIR + code + ".zip";
        if (rf == "" && File.exists(ROIDIR + code + ".roi")) rf = ROIDIR + code + ".roi";

        if (rf == "") {
            print("SKIP (no ROI): " + code);
        } else {
            close("*");
            open(BLIND + nm);
            orig = getImageID();
            getPixelSize(unit, pw, ph);
            if (pw < 0.01 || pw > 1) { pw = PXFALLBACK; ph = PXFALLBACK; }
            getDimensions(W, H, nC, nZ, nT);

            if (nC < 3) {
                print("SKIP (channels=" + nC + "): " + code);
                close("*");
            } else {
                // --- 抽通道，自己命名（不用 Split Channels）---
                selectImage(orig); run("Duplicate...", "title=D8 duplicate channels=1");
                selectImage(orig); run("Duplicate...", "title=G8 duplicate channels=2");
                selectImage(orig); run("Duplicate...", "title=R8 duplicate channels=3");
                selectImage(orig); close();

                // --- 32-bit 副本 ---
                selectWindow("G8"); run("Duplicate...", "title=G32"); run("32-bit");
                selectWindow("R8"); run("Duplicate...", "title=R32"); run("32-bit");

                // --- 二值遮罩（1/0）；用 (v>=T)*1，宏語言不支援 ?: ---
                selectWindow("G8"); run("Duplicate...", "title=MG"); run("32-bit");
                run("Macro...", "code=v=(v>=" + TG + ")*1");
                selectWindow("R8"); run("Duplicate...", "title=MR"); run("32-bit");
                run("Macro...", "code=v=(v>=" + TR + ")*1");

                // --- DAPI 核遮罩：逐視野 Otsu ---
                selectWindow("D8"); run("Duplicate...", "title=NUCtmp");
                setAutoThreshold("Otsu dark");
                getThreshold(lo, hi);
                resetThreshold();
                close("NUCtmp");
                if (lo < 1) lo = 1;
                selectWindow("D8"); run("Duplicate...", "title=NUC"); run("32-bit");
                run("Macro...", "code=v=(v>=" + lo + ")*1");

                // --- 交叉乘積 ---
                imageCalculator("Multiply create 32-bit", "G32", "MR"); rename("GinR");
                imageCalculator("Multiply create 32-bit", "R32", "MG"); rename("RinG");
                imageCalculator("Multiply create 32-bit", "G32", "R32"); rename("GR");
                imageCalculator("Multiply create 32-bit", "MG",  "NUC"); rename("MGnuc");

                // --- 逐 ROI 量測 ---
                roiManager("reset");
                roiManager("Open", rf);
                nR = roiManager("count");

                for (k = 0; k < nR; k++) {
                    selectWindow("D8"); roiManager("Select", k);
                    getRawStatistics(nPx, meanD);

                    selectWindow("G8"); roiManager("Select", k);
                    getRawStatistics(np, meanG, mn, mx, sdG);

                    selectWindow("R8"); roiManager("Select", k);
                    getRawStatistics(np, meanR, mn, mx, sdR);

                    selectWindow("MG");    roiManager("Select", k); getRawStatistics(np, fracG);
                    selectWindow("MR");    roiManager("Select", k); getRawStatistics(np, fracR);
                    selectWindow("GinR");  roiManager("Select", k); getRawStatistics(np, mGinR);
                    selectWindow("RinG");  roiManager("Select", k); getRawStatistics(np, mRinG);
                    selectWindow("GR");    roiManager("Select", k); getRawStatistics(np, mGR);
                    selectWindow("NUC");   roiManager("Select", k); getRawStatistics(np, nucFrac);
                    selectWindow("MGnuc"); roiManager("Select", k); getRawStatistics(np, mMGnuc);

                    M1 = -1; if (meanG > 0) M1 = mGinR / meanG;
                    M2 = -1; if (meanR > 0) M2 = mRinG / meanR;
                    pe = -1; if (sdG > 0 && sdR > 0) pe = (mGR - meanG * meanR) / (sdG * sdR);
                    en = -1; if (fracG > 0) en = M2 / fracG;
                    gn = -1; if (fracG > 0) gn = mMGnuc / fracG;

                    File.append(code + "," + (k + 1) + "," + nPx + "," + d2s(nPx * pw * ph, 2) + "," +
                                d2s(meanD, 3) + "," + d2s(meanG, 3) + "," + d2s(meanR, 3) + "," +
                                d2s(fracG, 6) + "," + d2s(fracR, 6) + "," +
                                d2s(M1, 5) + "," + d2s(M2, 5) + "," + d2s(pe, 5) + "," +
                                d2s(en, 4) + "," + d2s(nucFrac, 5) + "," + d2s(gn, 5), CSV);
                    nCells++;
                }

                roiManager("reset");
                close("*");
                nDone++;
                print("measured " + code + "  (" + nR + " cells)");
            }
        }
    }
}

setBatchMode(false);
print("DONE: " + nDone + " fields, " + nCells + " cells -> " + CSV);
File.saveString(getInfo("log"), OUT + "measure_log.txt");
