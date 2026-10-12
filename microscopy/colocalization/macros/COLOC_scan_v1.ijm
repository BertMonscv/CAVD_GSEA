// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// In the manuscript, histograms are post hoc threshold checks; thresholds were fixed in a pilot.
// ============================================================================
//  COLOC_scan_v1.ijm  —— LC3B / BODIPY 正式版 .lif 下机审计（不做定量）
//  只输出元数据与逐通道强度统计，用于：位深确认、饱和检查、采集参数一致性、
//  raw vs LIGHTNING 分辨、视野重复检测、阈值的 pooled 分布。
//
//  设计约束（基于此前在 Leica .lif 上反复失败的经验）：
//   - 不调用 Ext.getPlanePositionX/Y/Z、Ext.getPlaneTiming*、Ext.getPixelsPhysicalSize*
//     （Leica 把座标放在 Tile 节点，OME Plane.PositionX 为 null，必挂）
//   - 座标/激光/增益/针孔一律从 Info 字串抄行，交由后续解析
//   - 每个 series 结束即把 Log 冲进硬盘，挂掉也留得住进度
//   - 不用 arr.length、不用 requires()
//
//  输出（都在 OUT 目录）：
//   scan_coloc.csv     逐 series 逐通道统计
//   hist_coloc.csv     逐文件逐通道 256-bin 汇总直方图（用于定阈值）
//   meta_hits.csv      Info 里的关键行原文
//   info_full__*.txt   每个文件第一个 series 的完整 Info dump
//   scan_log.txt       运行日志
// ============================================================================

// ----------------------------- 只改这两行 -----------------------------------
DIR = getDirectory("Select raw LC3B-BODIPY LIF input folder");     // 放正式版 .lif 的资料夹，结尾必须有 /
OUT = getDirectory("Select colocalization working folder");                                // 留空 = 自动写到 ~/coloc_out/（外接盘常是只读）
// ---------------------------------------------------------------------------

setBatchMode(true);
run("Bio-Formats Macro Extensions");

// ---- 输出目录：确认可写，否则退回 ~/coloc_out/ ----
if (OUT == "") OUT = getDirectory("home") + "coloc_out/";
File.makeDirectory(OUT);
if (!File.exists(OUT)) {
    print("WARN: 无法建立 " + OUT + " -> 改用 home");
    OUT = getDirectory("home") + "coloc_out/";
    File.makeDirectory(OUT);
}
File.saveString("write test", OUT + "_writetest.txt");
if (!File.exists(OUT + "_writetest.txt"))
    exit("输出目录不可写：" + OUT + "\n请手动指定一个可写路径给 OUT。");
File.delete(OUT + "_writetest.txt");
print("OUT = " + OUT);

csv  = OUT + "scan_coloc.csv";
hcsv = OUT + "hist_coloc.csv";
mcsv = OUT + "meta_hits.csv";
logf = OUT + "scan_log.txt";

File.saveString("file,series,seriesName,isProcessed,W,H,nC,nZ,nT,bits,ch,minV,maxV,mean,mode,median,p99,satPct,zeroPct\n", csv);
File.saveString("file,ch,bin,count\n", hcsv);
File.saveString("file,line\n", mcsv);
print("\\Clear");

// Info 里要抄回来的关键字（大小写不敏感）
KEYS = newArray("PosX","PosY","PosZ","Pos ","Objective","NumericalAperture","Zoom",
                "PinholeAiry","Pinhole","LaserName","LaserWavelength","LaserOutputPower",
                "IntensityDev","DetectorName","Gain","Offset","ScanSpeed","Line Average",
                "LineAverage","Frame Average","Accumulation","DimensionDescription",
                "NumberOfElements","Length","Lightning","LIGHTNING","Clearing","Deconv",
                "RefractiveIndex","ImmersionRefraction","Format","BitSize");

list = getFileList(DIR);
nFiles = lengthOf(list);

for (f = 0; f < nFiles; f++) {
    nm = list[f];
    if (endsWith(toLowerCase(nm), ".lif")) {
        path = DIR + nm;
        print("=== FILE: " + nm);
        File.saveString(getInfo("log"), logf);

        nS = 0;
        Ext.setId(path);
        Ext.getSeriesCount(nS);
        print("    series count = " + nS);
        File.saveString(getInfo("log"), logf);

        // 每个文件的 pooled 直方图（最多 4 通道 × 256 bin）
        pool = newArray(4 * 256);
        for (i = 0; i < 4 * 256; i++) pool[i] = 0;
        infoDumped = false;

        for (s = 0; s < nS; s++) {
            sname = "";
            Ext.setSeries(s);
            Ext.getSeriesName(sname);
            print("  S1 open series " + s + " : " + sname);
            File.saveString(getInfo("log"), logf);

            run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (s + 1));

            // ---- Info dump：只做一次，且在任何可能出错的操作之前 ----
            if (!infoDumped) {
                info = getMetadata("Info");
                File.saveString(info, OUT + "info_full__" + nm + ".txt");
                lines = split(info, "\n");
                for (li = 0; li < lengthOf(lines); li++) {
                    L = lines[li];
                    hit = false;
                    for (k = 0; k < lengthOf(KEYS); k++) {
                        if (indexOf(L, KEYS[k]) >= 0) hit = true;
                    }
                    if (hit) File.append(nm + ",\"" + replace(L, "\"", "'") + "\"", mcsv);
                }
                infoDumped = true;
                print("  S2 info dumped");
                File.saveString(getInfo("log"), logf);
            }

            getDimensions(W, H, nC, nZ, nT);
            bits = bitDepth();
            proc = 0;
            up = toUpperCase(sname);
            if (indexOf(up, "LIGHTNING") >= 0) proc = 1;
            if (indexOf(up, "_LNG") >= 0) proc = 1;
            if (indexOf(up, "SVCC") >= 0) proc = 1;
            if (indexOf(up, "DECONV") >= 0) proc = 1;

            nBins = 256;
            if (bits == 16) nBins = 65536;

            for (c = 1; c <= nC; c++) {
                if (nC > 1) Stack.setChannel(c);
                if (nZ > 1) Stack.setSlice(1);   // 单层协议；若真有 Z 只取第 1 层并在日志提示

                values = newArray(1);
                counts = newArray(1);
                getHistogram(values, counts, nBins);

                total = 0; sum = 0.0; modeI = 0; modeC = -1;
                minV = -1; maxV = 0;
                for (i = 0; i < nBins; i++) {
                    ci = counts[i];
                    if (ci > 0) {
                        total = total + ci;
                        sum = sum + i * ci;
                        if (minV < 0) minV = i;
                        maxV = i;
                        if (ci > modeC) { modeC = ci; modeI = i; }
                    }
                }
                if (total == 0) total = 1;
                mean = sum / total;

                cum = 0; med = -1; p99 = -1;
                for (i = 0; i < nBins; i++) {
                    cum = cum + counts[i];
                    if (med < 0 && cum >= 0.50 * total) med = i;
                    if (p99 < 0 && cum >= 0.99 * total) p99 = i;
                }
                zeroPct = 100.0 * counts[0] / total;
                satPct  = 100.0 * counts[nBins - 1] / total;

                File.append(nm + "," + s + ",\"" + sname + "\"," + proc + "," + W + "," + H + "," +
                            nC + "," + nZ + "," + nT + "," + bits + "," + c + "," +
                            minV + "," + maxV + "," + d2s(mean, 2) + "," + modeI + "," +
                            med + "," + p99 + "," + d2s(satPct, 4) + "," + d2s(zeroPct, 4), csv);

                // pooled 直方图：只累加未处理(raw)的 series，压到 256 bin
                if (proc == 0 && c <= 4) {
                    div = nBins / 256;
                    for (i = 0; i < nBins; i++) {
                        b = floor(i / div);
                        if (b > 255) b = 255;
                        pool[(c - 1) * 256 + b] = pool[(c - 1) * 256 + b] + counts[i];
                    }
                }
            }
            close("*");
            print("  S3 done series " + s);
            File.saveString(getInfo("log"), logf);
        }
        Ext.close();

        for (c = 1; c <= 4; c++) {
            any = 0;
            for (b = 0; b < 256; b++) any = any + pool[(c - 1) * 256 + b];
            if (any > 0) {
                for (b = 0; b < 256; b++) {
                    File.append(nm + "," + c + "," + b + "," + pool[(c - 1) * 256 + b], hcsv);
                }
            }
        }
        print("=== FILE DONE: " + nm);
        File.saveString(getInfo("log"), logf);
    }
}

setBatchMode(false);
print("ALL DONE");
File.saveString(getInfo("log"), logf);
