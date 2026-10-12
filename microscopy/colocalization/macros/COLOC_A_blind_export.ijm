// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  COLOC_A_blind_export.ijm  —— 盲法匯出（第一步）
//
//  把 12 個 .lif 的全部 60 個視野以偽隨機順序匯出成編碼 TIFF（F001…F060），
//  並套用完全相同的顯示範圍。對照組（Sac / CQs / CTL）一併編碼，你在勾輪廓時
//  不會知道任何一張的組別。
//
//  對應 §2.13 的承諾：fixed pseudo-random order under coded identifiers,
//  identical display range applied to every image, key withheld until all
//  outlines are complete.
//
//  ⚠️ 跑完之後不要打開 KEY_DO_NOT_OPEN.csv，直到所有輪廓勾完。
// ============================================================================

DIR   = getDirectory("Select raw LC3B-BODIPY LIF input folder");        // .lif 來源
OUT   = getDirectory("Select colocalization working folder");  // 工作目錄
BLIND = OUT + "blind/";                       // 編碼影像
SEED  = 20260803;                             // 洗牌種子（記在方法裡即可重現）

// 顯示範圍：三個通道各自固定，全部影像一律相同。可以改，但改了就要全部重匯出。
D1_MIN = 0; D1_MAX = 100;   // C1 DAPI
D2_MIN = 0; D2_MAX = 120;   // C2 BODIPY
D3_MIN = 0; D3_MAX =  80;   // C3 LC3B

// ---------------------------------------------------------------------------
setBatchMode(true);
run("Bio-Formats Macro Extensions");
File.makeDirectory(OUT);
File.makeDirectory(BLIND);
File.makeDirectory(OUT + "roi/");
print("\\Clear");

// 1) 收集所有 (檔案, series)
list = getFileList(DIR);
srcF = newArray(0);
srcS = newArray(0);
for (f = 0; f < lengthOf(list); f++) {
    nm = list[f];
    if (endsWith(toLowerCase(nm), ".lif")) {
        nS = 0;
        Ext.setId(DIR + nm);
        Ext.getSeriesCount(nS);
        for (s = 0; s < nS; s++) {
            srcF = Array.concat(srcF, nm);
            srcS = Array.concat(srcS, s);
        }
    }
}
n = lengthOf(srcF);
print("total fields = " + n);

// 2) 偽隨機洗牌
random("seed", SEED);
r = newArray(n);
for (i = 0; i < n; i++) r[i] = random();
ord = Array.rankPositions(r);

// 3) 逐張匯出（Log 只印編碼，不印來源）
key = "code,file,series\n";
for (i = 0; i < n; i++) {
    j    = ord[i];
    code = "F" + IJ.pad(i + 1, 3);

    run("Bio-Formats Importer", "open=[" + DIR + srcF[j] + "] color_mode=Composite view=Hyperstack stack_order=XYCZT series_" + (srcS[j] + 1));

    Stack.setChannel(1); setMinAndMax(D1_MIN, D1_MAX);
    Stack.setChannel(2); setMinAndMax(D2_MIN, D2_MAX);
    Stack.setChannel(3); setMinAndMax(D3_MIN, D3_MAX);
    Stack.setChannel(2);

    saveAs("Tiff", BLIND + code + ".tif");
    close("*");

    key = key + code + "," + srcF[j] + "," + srcS[j] + "\n";
    print("exported " + code);
}

File.saveString(key, OUT + "KEY_DO_NOT_OPEN.csv");
setBatchMode(false);
print("DONE -> " + BLIND);
print("key saved (do not open yet): " + OUT + "KEY_DO_NOT_OPEN.csv");
