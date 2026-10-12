// Publication adaptation (2026-10-12): portable directory dialogs only.
// Measurement parameters and analysis logic are unchanged; see README and manifest.
// ============================================================================
//  COLOC_A2_outline.ijm  —— 盲法勾輪廓的流程宏（第二步）
//
//  逐張開 blind/F0xx.tif → 你用 freehand 勾細胞、按 T 加入 ROI Manager
//  → 按 OK → 自動存成 roi/F0xx.zip → 自動開下一張。
//
//  特性：
//   - 已經有 roi/F0xx.zip 的自動跳過 → 可以分多次做，隨時接續
//   - ROI Manager 空著就按 OK 會再問一次（防止忘記按 T）；確認「真的沒有合格細胞」
//     才記進 roi/_zero_fields.txt 並跳下一張
//   - 存檔前自動 Deselect（只選一個 ROI 時直接存會只存那一個）
//   - 想中途離開：在對話框按 Cancel，已完成的都已存檔
//   - 勾錯想重來：刪掉那個 F0xx.zip，再跑一次本宏即可
//
//  ⚠️ 全程不要動 Brightness/Contrast。顯示範圍已經烘進 TIFF，
//     調了就破壞 §2.13「identical display range」這個承諾。
// ============================================================================

OUT    = getDirectory("Select colocalization working folder");
BLIND  = OUT + "blind/";
ROIDIR = OUT + "roi/";

// 勾輪廓時用哪個顯示模式：composite = 三通道疊合（建議），grayscale = 單通道
DISPLAY_MODE = "composite";
START_CHANNEL = 2;   // 疊合時把作用中通道停在 C2(BODIPY)，切通道用鍵盤 < >

// ---------------------------------------------------------------------------
setBatchMode(false);
File.makeDirectory(ROIDIR);
setTool("freehand");
run("ROI Manager...");
print("\\Clear");

list = getFileList(BLIND);

// 先數一下進度
nTif = 0; nHave = 0;
for (i = 0; i < lengthOf(list); i++) {
    if (endsWith(list[i], ".tif")) {
        nTif++;
        code = substring(list[i], 0, lengthOf(list[i]) - 4);
        if (File.exists(ROIDIR + code + ".zip")) nHave++;
    }
}
print("共 " + nTif + " 張，已完成 " + nHave + " 張，待處理 " + (nTif - nHave) + " 張");

idx = 0;
for (i = 0; i < lengthOf(list); i++) {
    nm = list[i];
    if (endsWith(nm, ".tif")) {
        idx++;
        code = substring(nm, 0, lengthOf(nm) - 4);

        if (File.exists(ROIDIR + code + ".zip")) {
            print("skip " + code + " (已有 ROI)");
        } else {
            done = false;
            while (!done) {
                open(BLIND + nm);
                Stack.setDisplayMode(DISPLAY_MODE);
                if (DISPLAY_MODE == "composite") Stack.setChannel(START_CHANNEL);
                roiManager("reset");
                roiManager("Show All with labels");

                waitForUser("勾輪廓　" + code + "　（第 " + idx + " / " + nTif + " 張）",
                    "freehand 勾一個細胞 → 按 T 加入 ROI Manager，重複到這張勾完。\n" +
                    "　\n" +
                    "納排標準（§2.13 已寫死，務必照做）：\n" +
                    "　· 觸及影像邊緣的細胞 → 不勾\n" +
                    "　· 邊界不明的重疊細胞 → 不勾\n" +
                    "　· 有絲分裂或形態受損的細胞 → 不勾\n" +
                    "　· 碎片 → 不勾\n" +
                    "　· 其餘每一個合格細胞都要勾（不可挑亮的勾）\n" +
                    "　\n" +
                    "鍵盤：T = 加入 ROI Manager　|　< > = 切換通道\n" +
                    "不要調 Brightness/Contrast。\n" +
                    "　\n" +
                    "勾完按 OK 自動存檔並開下一張；按 Cancel 離開（進度已保留）。");

                nR = roiManager("count");

                if (nR > 0) {
                    roiManager("Deselect");
                    roiManager("Save", ROIDIR + code + ".zip");
                    print(code + " : " + nR + " cells  -> saved");
                    done = true;
                } else {
                    ok = getBoolean("ROI Manager 是空的。\n\n" + code + " 這張真的沒有任何合格細胞嗎？\n\n" +
                                    "是 = 記錄為 0 細胞並跳下一張\n" +
                                    "否 = 回到這張重勾（可能是忘了按 T）");
                    if (ok) {
                        File.append(code, ROIDIR + "_zero_fields.txt");
                        print(code + " : 0 cells (記錄為無合格細胞)");
                        done = true;
                    }
                }
                close("*");
            }
        }
    }
}

print("ALL DONE");
File.saveString(getInfo("log"), OUT + "outline_log.txt");
