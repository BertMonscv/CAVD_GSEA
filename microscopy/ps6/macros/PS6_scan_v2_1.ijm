// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================
//  PS6_scan_v2.ijm  —  p-S6 IF 下機審計宏（THUNDER，不做定量）
//  ① 時間戳與真實採集順序   ② 四組曝光/LED 是否一致
//  ③ 座標分簇 → 哪些視野同一張爬片
//  ④ 逐通道統計 → 焦度、飽和篩選
//  ⑤ raw 與 SVCC 並排掃 → CC 參數核對 + 敏感度分析備料
//  輸出：~/ps6_out/PS6_scan.csv、PS6_scan_log.txt、info_full_*.txt
// ============================================================

var DIR       = getDirectory("Select raw p-S6 LIF input folder");   // p-S6 的 .lif 資料夾
var SCAN_SVCC = true;   // true = raw 與 SVCC 都掃（is_raw 欄區分）；false = 只掃 raw
var RECURSE   = true;   // true = 一併掃子資料夾

var OUT = getDirectory("Select empty p-S6 scan output folder");
var csv = "";
var nLif = 0;
var nRow = 0;

// ---- 路徑正規化 ----
if (!endsWith(DIR, "/")) DIR = DIR + "/";
if (!File.exists(DIR)) exit("找不到資料夾：\n" + DIR + "\n請確認外接碟已掛載，或改 DIR 這一行");
File.makeDirectory(OUT);
if (!File.exists(OUT)) exit("無法建立輸出資料夾：" + OUT);

setBatchMode(true);
run("Bio-Formats Macro Extensions");
print("\\Clear");
print("掃描起點：" + DIR);

csv = "file,series_idx,series_name,is_raw,bitDepth,nC,width,height,nZ,pxW_um,"
    + "timestamp,posX,posY,posZ,expo_all,led_all,cc_params,"
    + "C0_min,C0_med,C0_mode,C0_p99,C0_max,C0_satPct,C0_zeroPct,"
    + "C1_min,C1_med,C1_mode,C1_p99,C1_max,C1_satPct,C1_zeroPct\n";

walk(DIR);

File.saveString(csv, OUT + "PS6_scan.csv");
print("------------------------------------------");
print("完成：" + nLif + " 個 .lif，" + nRow + " 筆 series");
if (nLif == 0) print("!! 沒找到任何 .lif —— 檢查路徑大小寫，或把 RECURSE 設 true");
print("輸出 → " + OUT);
selectWindow("Log");
saveAs("Text", OUT + "PS6_scan_log.txt");
setBatchMode(false);


// ---------- 走訪資料夾（可遞迴）----------
function walk(d) {
    list = getFileList(d);
    for (k = 0; k < list.length; k++) {
        p = d + list[k];
        if (File.isDirectory(p)) {
            if (RECURSE) walk(p);
            continue;
        }
        if (endsWith(toLowerCase(list[k]), ".lif")) scanLif(p, list[k]);
    }
}


// ---------- 掃單一 .lif ----------
function scanLif(path, fname) {
    nLif++;
    print("=== [" + nLif + "] " + fname + " ===");
    Ext.setId(path);
    Ext.getSeriesCount(nSer);
    print("  series count = " + nSer);

    dumped = false;
    for (s = 0; s < nSer; s++) {
        Ext.setSeries(s);
        Ext.getSeriesName(sname);

        // --- 按名稱判 raw / SVCC（不可用 i+=2，交替結構不規則）---
        isRaw = 1;
        if (indexOf(sname, "_Lng_")     >= 0) isRaw = 0;
        if (indexOf(sname, "SVCC")      >= 0) isRaw = 0;
        if (indexOf(sname, "Lightning") >= 0) isRaw = 0;
        if (indexOf(sname, "Decon")     >= 0) isRaw = 0;
        if (isRaw == 0 && !SCAN_SVCC) continue;

        run("Bio-Formats Importer",
            "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (s+1));
        getDimensions(w, h, nC, nZ, nT);
        bd = bitDepth();
        getPixelSize(unit, pw, ph);
        info = getMetadata("Info");

        // 每個 .lif 存一份完整 Info，供我核對鍵名
        if (!dumped) {
            base = replace(fname, ".lif", "");
            base = replace(base, "[^A-Za-z0-9_-]", "_");
            File.saveString(info, OUT + "info_full_" + base + ".txt");
            dumped = true;
        }

        ts = grab(info, "Timestamp", 1);
        if (ts == "NA") ts = grab(info, "AcquiredDate", 1);
        if (ts == "NA") ts = grab(info, "Date", 1);

        px = grab(info, "StagePosX", 1);
        if (px == "NA") px = grab(info, "PositionX", 1);
        py = grab(info, "StagePosY", 1);
        if (py == "NA") py = grab(info, "PositionY", 1);
        pz = grab(info, "PosZ", 1);
        if (pz == "NA") pz = grab(info, "PositionZ", 1);

        ex = grab(info, "ExposureTime", 4);
        if (ex == "NA") ex = grab(info, "Exposure", 4);
        led = grab(info, "IntensityDev", 4);
        if (led == "NA") led = grab(info, "LightSourceIntensity", 4);

        // 注意開頭的 "" ——  IJ macro 若讓字串回傳函式當串接式的第一個運算元，
        // 直譯器會預期數值回傳並報 "Numeric return value expected"
        cc = "" + grab(info, "Clearing", 2) + "~" + grab(info, "FeatureScale", 1)
           + "~" + grab(info, "Strength", 1) + "~" + grab(info, "Decon", 1);

        // --- 守門：非 16-bit 或通道不足直接跳過 ---
        if (bd != 16 || nC < 2) {
            print("  !! 跳過 s" + s + " (" + sname + ")  bitDepth=" + bd + " nC=" + nC);
            close("*");
            continue;
        }

        safeName = replace(sname, ",", ";");
        row = fname + "," + s + "," + safeName + "," + isRaw + "," + bd + ","
            + nC + "," + w + "," + h + "," + nZ + "," + d2s(pw, 5) + ","
            + ts + "," + px + "," + py + "," + pz + "," + ex + "," + led + "," + cc;

        for (c = 1; c <= 2; c++) {
            Stack.setChannel(c);
            row = row + "," + chanStats();
        }
        csv = csv + row + "\n";
        nRow++;

        tag = "raw ";
        if (isRaw == 0) tag = "SVCC";
        print("  s" + s + " [" + tag + "] " + safeName + "  t=" + ts
              + "  X=" + px + " Y=" + py + "  expo=" + ex);
        close("*");
    }
    Ext.close();
}


// ---------- 單通道直方圖統計 ----------
function chanStats() {
    nb = 65536;
    getHistogram(vals, cnt, nb);
    tot = 0;
    for (i = 0; i < nb; i++) tot += cnt[i];
    if (tot == 0) return "0,0,0,0,0,0,0";

    mn = -1; mx = -1; mo = 0; moC = -1;
    for (i = 0; i < nb; i++) {
        if (cnt[i] > 0) {
            if (mn < 0) mn = i;
            mx = i;
            if (cnt[i] > moC) { moC = cnt[i]; mo = i; }
        }
    }
    half = tot * 0.50;  top = tot * 0.99;
    med = -1; p99 = -1; acc = 0;
    for (i = 0; i < nb; i++) {
        acc += cnt[i];
        if (med < 0 && acc >= half) med = i;
        if (p99 < 0 && acc >= top)  p99 = i;
    }
    satP  = 100.0 * cnt[65535] / tot;
    zeroP = 100.0 * cnt[0] / tot;
    return "" + mn + "," + med + "," + mo + "," + p99 + "," + mx + ","
           + d2s(satP, 4) + "," + d2s(zeroP, 4);
}


// ---------- 從 Info 抓含關鍵字的值 ----------
function grab(txt, key, nHit) {
    lines = split(txt, "\n");
    res = ""; n = 0;
    for (i = 0; i < lines.length; i++) {
        if (indexOf(lines[i], key) >= 0) {
            p = indexOf(lines[i], "=");
            if (p < 0) continue;
            v = String.trim(substring(lines[i], p + 1));
            v = replace(v, ",", ";");
            if (v == "") continue;
            if (n > 0) res = res + "|";
            res = res + v;
            n++;
            if (n >= nHit) i = lines.length;
        }
    }
    if (n == 0) return "NA";
    return res;
}
