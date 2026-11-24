# Tüm deneyleri tek seferde çalıştırma (Windows + conda)

## Kısa yol

`basin_project\run_all.cmd` dosyasına çift tıklayın (veya Miniconda terminalinde çalıştırın).
Betik sırasıyla şunları yapar:

1. conda'yı bulur; `basin-project` ortamı yoksa conda-forge'dan kurar (python 3.12, numpy, scipy, matplotlib, packaging, psutil). Var olan ortama dokunmaz.
2. Ortam kontrolü (`check_environment`).
3. Tüm testler (52 test, ~1 dk). Bir test kırılırsa deneyler **başlamaz**.
4. `python -m basinlab.run_all` → bütün keşifsel deneyler, sırayla, her biri ayrı süreçte.

## Seçenekler (cmd'ye argüman olarak verilir)

| Komut | Ne yapar |
|---|---|
| `run_all.cmd --quick` | Her yeni deneyin küçük (duman) sürümü, ~3 dk. **Önce bunu çalıştırın.** Çıktılar `data/*_quick/` |
| `run_all.cmd` | Tam koşu. Daha önce çıktısı olan eski deneyler atlanır; yeniler çalışır (~3.5–4 saat, laptopa göre değişir) |
| `run_all.cmd --list` | Deney listesi, bağımlılıklar, kaba süre tahmini |
| `run_all.cmd --only regime_map,diagnostic_eval` | Sadece seçilenler |
| `run_all.cmd --skip grid_convergence` | Birini dışarıda bırak |
| `run_all.cmd --stage legacy --force` | Eski 16 deneyi baştan üret (~35 dk) |
| `run_all.cmd --force` | Her şeyi baştan |
| `run_all.cmd --confirmatory` | Ek olarak **kilitli** doğrulayıcı test (aşağıya bakın) |
| `run_all.cmd --dry-run` | Plan, çalıştırmadan |

## Kesilirse

Tekrar çalıştırın. Biten deneyler atlanır (`data/_runs/state.json`); uzun deneyler (`regime_map`, `grid_convergence`, `beam2mode_basins`, `fe_arch_basins`, doğrulayıcı) kaldığı parçadan devam eder (`data/<deney>/parts/`). Bir deneyin kaynak kodu değişirse otomatik yeniden çalışır; bağımlı deneyler (ör. `diagnostic_eval`) da yenilenir.

## Sonuçlar

- `data/_runs/LATEST.txt` → son koşunun klasörü; içinde `REPORT.md` (durum, süre, tepe bellek) ve her deneyin `.log` dosyası.
- Her deneyin sonuçları: `data/<deney>/summary.json` (ayar, sürüm, commit, süre, sonuç).
- Hata olan deney diğerlerini durdurmaz; rapor ve çıkış kodu (1) bildirir.

## Tam koşudaki deneyler

| Deney | İş paketi | Kaba süre | Ne ölçer |
|---|---|---:|---|
| fe_arch_validation | WP5 | 2 dk | FE kirişte Euler burkulma, burkulma sonrası genlik, ağ yakınsaması, eyer araması |
| beam2mode_basins | WP4 | 15 dk | İki modlu burkulmuş kiriş havzaları; q2=0'da SDOF'a indirgeme; bağlaşım etkisi |
| regime_map | WP2 | 70 dk | 5 c × 8 h × 4 aile × ρ ızgarası + BDF2 + discrete-gradient; doygunluk işaretli |
| diagnostic_eval | WP6 | <1 dk | İki-ρ ve iki-h testleri gerçek hatayı yakalıyor mu |
| fe_arch_basins | WP5 | 25 dk | FE burkulmuş kirişte havza uyuşmazlığı (gen-α, 3 Ωh, 6 ρ) |
| grid_convergence | WP1 | 70 dk | 100/200/400 ızgarada mismatch, sınır kayması, sabit kutu entropisi |
| confirmatory_smoke | – | <1 dk | Doğrulayıcı kodun küçük duman testi (doğrulayıcı değildir) |

Süreler bu oturumdaki ölçümlerden kaba tahmindir; laptopunuzda farklı olabilir. İlk `--quick` koşusunun raporundaki süreler fikir verir.

## Kilitli doğrulayıcı test

`confirmatory_family_comparison`, `protocol/PREREG.md` şu koşulları sağlamadan **çalışmaz**:
1. dosyada `Status: APPROVED` satırı var ve `DRAFT` kelimesi yok,
2. dosya git'te commit edilmiş ve değiştirilmemiş.

Önce PREREG'deki açık maddeleri kararlaştırın (kodda "PROPOSED" olarak uygulanan öneriler dosyanın sonuna eklendi), sonra onaylayıp commit edin, sonra `run_all.cmd --confirmatory`. Kural: sonuç ne çıkarsa çıksın eşikler değiştirilmez.
