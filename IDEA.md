# Makale 2 — Fikir ve Yol Haritası

**Çalışma başlığı:** *Numerical damping selects equilibria: how integrator dissipation shifts basin boundaries in nonlinear structural dynamics*

**Kaynak:** ACME review'daki ("When the Solver Becomes the System") açık problem #2, "Numerical damping and basin boundaries" (Sec. 4.x `sec:time:rhoinf`, Fig. 6c, `sections/11_open_problems.tex`).
**Durum:** Ön denemeler tamamlandı (pilot 1–7: SDOF çift kuyu, eğik kuyu, dört integratör ailesi, 2 DOF). Asıl koşular ayrı bir makale klasöründe yapılacak.

---

## 1. Tek cümlelik iddia

İntegratörün yüksek frekansı bastırmak için seçilen sayısal sönümü (ρ∞), çok kararlı bir yapının hangi dengeye yerleşeceğini değiştirir. Etkiyi belirleyen ρ∞ değil, **fiziksel frekanstaki çevrim-başına sayısal enerji kaybı oranı D**'dir (aileler arası ortak koordinat). D = 1 civarında havza hatası minimumdur ve küçük D hatası havza sınırında büyür. Bu nedenle çok kararlı analizlerde ρ∞ yerine D (ya da iki ρ∞ ile havza testi) raporlanmalıdır.

## 2. Motivasyon ve boşluk

- Mühendisler ρ∞'yi yüksek frekanslı gürültüyü bastırmak için seçer, denge seçimi için değil.
- Review'da (Fig. 6c) çift kuyulu osilatörde yanlış kuyuya düşen başlangıç koşulu oranının ρ∞'ya monoton olmayan şekilde bağlı olduğu gösterildi. Bu yalnızca 3 ρ∞ ve 3 h değeriyle, tek sönümle, tek modelle yapıldı. Review'da "böyle bir çalışma bilmiyoruz" diye açık problem olarak bırakıldı.
- Yakın literatür:
  - Çok kararlılık ve havza kararlılığı (Brzeski vd. 2016, impacting oscillator, up/down sweep) integratörü **kusursuz** varsayar.
  - İntegratör kıyas çalışmaları (Kuhl 1999, Lavrenčič 2020, Fernandes 2022) kararlılık, enerji ve doğruluğa bakar, havza sınırına bakmaz.
  - Sönüm–havza ilişkisi fiziksel sönüm için çalışılmış (bistabil Duffing, fraktal sınırlar), **sayısal** sönüm için bulunmadı.
- Hızlı web taraması (Kasım 2026) doğrudan çakışan bir çalışma göstermedi. **Bu tarama kapsamlı değil.** Makale yazımına geçmeden önce Scopus/WoS/Google Scholar'da sistematik tarama şart (bkz. Faz 0).

## 3. Araştırma soruları

- **RQ1 (olgu):** (ρ∞, h, c) uzayında, integratör kaynaklı havza hatası nerede önemli, nerede önemsiz? (rejim haritası)
- **RQ2 (mekanizma):** Havza hatasını minimize eden ρ∞\*, çevrim-başına enerji-kaybı eşleşmesiyle tahmin edilebilir mi?
- **RQ3 (büyütme):** Küçük enerji-kaybı hatası havza sınırında neden kat kat büyüyor? Kazancı öngörebilir miyiz?
- **RQ4 (genellik):** Bulgu SDOF oyuncak modelden gerçek yapısal modellere (snap-through, bistabil kiriş) ve diğer integratör ailelerine taşınıyor mu?
- **RQ5 (pratik):** Mühendis için ucuz bir tanı/kural ne olmalı? (review'daki DT3'ün somutlaştırılması)

## 4. Hipotezler

- **H1:** Optimum ρ∞\* vardır (U şekli). Düşük ρ∞'da sayısal sönüm gerçeği aşar (erken yakalama), ρ∞=1'de eksik kalır (kuyuyu aşma).
- **H2:** ρ∞\*, çevrim-başına enerji kaybı oranı D(h, ρ∞, c) = 1 koşuluyla ±0.1 içinde tahmin edilir. Doğrusal (kuyu-dibi) sönüm eşleşmesi h-bağımlılığını yanlış verir, genlik-ortalamalı enerji kaybı doğru verir.
- **H3 (büyütme):** Havza hatası ≈ (birikmiş sayısal enerji hatası) / (çevrim başına fiziksel enerji kaybı). Yani yakalama çevrim sayısı değişen başlangıç koşullarının oranı. Hafif sönümde payda küçük, kazanç büyük.
- **H4:** Etki ağır sönümde kaybolur (ön deneme: c=0.4'te düz).
- **H5:** Etki dissipatif tüm ailelerde görülür, ama **ρ∞ ortak koordinat değildir**; hata D (çevrim-başına enerji kaybı oranı) ile örgütlenir. Düşük frekanslı sönüm hatası O(h) olan şemalar (dissipatif Newmark) ρ∞ düşükken en kötü sonucu verir. Enerji-koruyan şemalar ρ∞-tipi parametre taşımadığı için referans. *(Ön denemeyle desteklendi: pilot 5, 7; Faz 2'de doygunluk dışı doğrulanacak.)*

## 5. Ön deneme kanıtı (SDOF çift kuyu)

Model: `x'' + c x' − x + x³ = 0`, havzalar ±1, referans RK4 (h=0.01, t=60–80, t=60/120 arasında ayrışma yok). Izgara 100×100, [−2,2]². Kodlar: `pilots/`. İntegratörün yerel kopyası `pilots/genalpha.py` ile özgün `zenodo_package_v1/benchmarks/genalpha.py` aynıdır.

**Düzeltme ve uyarı:** Önceki denklemde kuvvet işaretleri ters yazılmıştı; yukarıdaki denklem pilot koduna ve çift kuyu potansiyeline uygundur. Pilot 4, `round(T/h)` nedeniyle bazı `h` değerlerinde referans ve integratörü farklı zamanlarda karşılaştırır (T=8, h=0.3 için 8 ve 8.1; h=0.6 için 8 ve 7.8). Dolayısıyla ilgili D sonuçları yeniden hesaplanmalıdır. Bu sorun, h=0.4 ve h=0.5'teki sapmayı tek başına açıklamaz.

**(a) Rejim (c=0.15), yanlış kuyu oranı %:**

| h \ ρ∞ | 1.0 | 0.8 | 0.6 | 0.4 | 0.2 | 0.0 |
|---|---|---|---|---|---|---|
| 0.1 | 1.0 | 1.1 | 1.1 | 1.2 | 1.2 | 2.8 |
| 0.2 | 4.5 | 4.6 | 4.5 | 3.8 | 5.7 | 15.1 |
| 0.4 | 18.0 | 17.9 | 14.6 | 11.5 | 19.0 | 40.7 |
| 0.6 | 32.3 | 31.8 | 27.0 | 20.2 | 32.8 | 47.6 |
| 1.0 | 31.8 | 32.3 | 34.7 | 42.7 | 39.7 | 42.5 |

h ≥ 0.8'de harita rastgeleye yakın (≈%30–47), doygun.

**(b) Sönüme bağlılık (argmin ρ∞, h=0.2…0.6):** c=0.05 → ≈0.53–0.8; c=0.15 → ≈0.37–0.40 (h'ye karşı düz); c=0.4 → etki yok (%2–16, düz). c=0.05'te ρ∞ ≤ 0.3 için hata ≈%40 (harita fraktal/rastgele), ρ∞≈0.5–0.6'da %7 (h=0.2).

**(c) Mekanizma testi:**
- *Doğrusal sönüm eşleşmesi* (ω₀=√2'de çarpan analizi): yön doğru, **h-bağımlılığı yanlış** (tahmin 0.18→0.42, gözlenen düz). Hata ≈0.08–0.11. **Yetersiz.**
- *Çevrim-başına enerji-kaybı eşleşmesi* (D = ΔE_num / ΔE_ref, medyan, 2372 IC, T=8): c=0.05 h≤0.4 için hata ≤0.06; c=0.15, h=0.2'de 0.01; ortalama hata ≈0.06–0.07; h-bağımlılığını doğru (düz) veriyor. **Daha iyi, tam değil** (c=0.15, h≥0.3'te tahmin gözlemin ≈0.05–0.09 üstünde).
- *Büyütme:* ρ∞=1'de enerji-kaybı eksikliği yalnızca %2–12, havza hatası %4–37.

**Kanıtın zayıf yanları (pilot 1–4 için):** tek model, tek ailesi, doygun ölçü (yüzde fark), tek ızgara, doygun hücrelerde argmin güvenilmez. Aşağıdaki (d)–(f) bunların bir kısmını kapattı.

**(d) Entegratör aileleri (pilot 5, `pilots/pilot5_families.py`):** her ailede parametre, yüksek-frekans spektral yarıçapı ρ∞'a (sayısal ölçülmüş) eşlenerek karşılaştırıldı. c=0.15, h=0.4, yanlış kuyu oranı %:

| ρ∞ | 1.0 | 0.9 | 0.8 | 0.7 | 0.6 | 0.5 | 0.4 | 0.2 | 0.0 |
|---|---|---|---|---|---|---|---|---|---|
| gen-α | 18.0 | 18.0 | 17.9 | 16.8 | 14.6 | 12.7 | **11.5** | 19.0 | 40.7 |
| HHT (ρ∞≥0.5) | 18.0 | 15.2 | 13.8 | 13.2 | **13.0** | – | – | – | – |
| WBZ | 18.0 | 14.7 | 12.9 | 11.8 | **11.5** | 12.0 | 14.9 | 26.4 | – |
| Newmark (γ>½) | 18.0 | 14.1 | 37.0 | 45.9 | 49.0 | – | – | – | – |

- Gen-α, HHT ve WBZ aynı U şeklini verir, ama **optimum ρ∞ aileye göre kayar** (0.4 / ≥0.6 / 0.6). Yani **ρ∞ evrensel koordinat değil.**
- **Dissipatif Newmark** ρ∞ = 0.8'de %37–49 hata verir (gen-α/WBZ: %13–18). Neden: γ>½ yüzünden düşük frekanslı sönüm O(h), yani fiziksel salınımı da sönümler (gen-α/HHT/WBZ'de O(h³)). Bu, mekanizmanın (düşük frekansta enerji kaybı) doğrudan kanıtı.
- Eğik kuyu (asimetri, e=0.1) sonucu neredeyse değiştirmedi (gen-α 18.0 → 17.8, optimum 0.4 → 0.4). Etki simetriye bağlı değil.

**(e) Çökme (collapse) testi (pilot 7, `pilots/pilot7_collapse.py`):** dört ailede tüm (aile, ρ∞) noktaları için D (çevrim-başına enerji kaybı oranı) ile havza hatası karşılaştırıldı.

| | c=0.15, h=0.4 | c=0.05, h=0.2 |
|---|---|---|
| Spearman (hata, \|D−1\|) | **0.81** | **0.91** |
| Spearman (hata, 1−ρ∞) | 0.09 | 0.46 |

- Hata, **ρ∞'dan çok D'ye** göre örgütleniyor. Evrensel koordinat D. Bu, H5'i güçlendirir ve makalenin ana mesajını değiştirir: *"ρ∞ yerine sayısal enerji kaybı oranı raporlanmalı/kontrol edilmeli."*
- D'ye en yakın ρ∞, havza optimumundan genelde ≤0.1 farklı (gen-α 0.5 vs 0.4 ve 0.6 vs 0.5; HHT aynı veya 0.8 vs 0.6 [aralık sınırlı]; WBZ 0.7 vs 0.6 ve 0.9 vs 0.8). Hâlâ sistematik küçük sapma var (D eşleşmesi optimumun biraz **üstünü** işaret ediyor). Ayrıca D≈0.98 olan ρ∞=1'de bile hata %14–18, yani küçük D hatası büyük havza hatasına çevriliyor (H3 kazancı).
- Dikkat: ortak sıralama doygun (%40 civarı) noktaları da içeriyor, bu korelasyonu şişirebilir. Doygunluk dışı alt kümede ve basin-entropy/sınır-kayması ölçüleriyle yinelenmeli.

**(f) 2 DOF, yüksek frekanslı ikinci modlu (pilot 6, `pilots/pilot6_2dof.py`):** x bistabil mod, y sert mod (ω_y≈3.3, c_y=0.3), g=1 bağlaşım, gen-α, hata %:

| h \ ρ∞ | 1.0 | 0.8 | 0.6 | 0.4 | 0.2 | 0.0 |
|---|---|---|---|---|---|---|
| 0.2 | 5.9 | 6.1 | 6.0 | **4.2** | 5.5 | 15.5 |
| 0.4 | 21.6 | 21.6 | 17.1 | **12.1** | 19.1 | 41.5 |

Aynı U şekli ve aynı optimum (≈0.4); sert modun varlığı etkiyi kaldırmıyor. (Bağlaşım zayıf, sonuç SDOF'a çok yakın. Gerçek yapısal vaka için daha güçlü mod etkileşimi gerekiyor.)

## 6. Önerilen kuramsal çerçeve (geliştirilecek)

1. **Enerji-denge modeli.** Sönümlü çift kuyuda yörünge, enerjisi bariyerin altına düşünce yakalanır. Yakalanan kuyu, enerjinin bariyer altına düştüğü çevrim ve fazla bağlı. Sayısal enerji kaybı hatası δE(E, h, ρ∞) bu eşiği kaydırır.
2. **D(h, ρ∞, c; E):** genlik-ortalamalı sayısal/gerçek enerji kaybı. Yaklaşım yolları: (i) çarpan analizi + yavaş değişen genlikte averaging (harmonik denge) ile analitik, (ii) sayısal ölçüm (pilot 4 gibi).
3. **ρ∞\*(h, c):** D = 1.
4. **Sınır kayması ve kazanç:** Yakalama çevrim sayısı N_c(IC) için, havza hatası ≈ P(ΔN_c ≠ 0). Kazancı G = (havza hatası) / (|δE|/E) olarak tanımla ve c, h ile ölçekle. (H3'ün test edilmesi.)
5. **Genişleme:** Zorlamalı sistem (harmonik), periyodik çekerler. Burada havza çekerlerin arasındadır, kuyuların değil.

## 7. Yöntem

**Modeller (basitten gerçeğe):**
1. SDOF çift kuyu (Duffing, c taraması, asimetrik eğim ε ile).
2. Von Mises kafesi (shallow two-bar truss, gerçek geometri, tek DOF, snap-through).
3. **2 DOF bistabil model:** mod-bağlı indirgenmiş burkulmuş kiriş (iki modlu Galerkin) ya da kütle-yay snap-through (Royal Society 2024 çalışmasındaki gibi). Mod etkileşimi ve homoklinik yapı yüzünden asıl "gerçekçilik" katkısı bu.
4. (Opsiyonel/genişletme) Bistabil eğri kiriş SEM; FE modeli (FEniCSx/OpenSees).

**İntegratörler:** generalized-α (ρ∞), HHT-α, WBZ-α, Newmark (γ/β dissipatif), Bathe (q1/ρ∞-Bathe), BDF2. **Kontroller:** enerji-koruyan (discrete gradient, Simo/Gonzalez), trapezoidal (ρ∞=1).
**Referans:** DOP853 / RK8 sıkı tolerans (≤1e−12), h_ref ve tolerans yakınsama testli. Yakınsama kanıtlanmadan referans kabul edilmez.

**Parametre taraması:** ρ∞ ∈ [0,1] (≥11 değer), h ∈ log-aralık (≥8), c (≥4: hafif→ağır), asimetri ε, ızgara 200²–400² + yakınsama testi.

**Ölçüler (doygunluğa dirençli):**
- Yanlış-kuyu oranı (bilinen ölçü, ön denemeyle sürekli olması için).
- **Basin entropy** Sb ve boundary basin entropy Sbb (Daza vd.) ile fraktallık/ln2 testi.
- Belirsizlik üssü (uncertainty exponent).
- **Sınır kayması:** referans ve integratör sınırı arasındaki ortalama işaretli mesafe (Hausdorff).
- Enerji-kaybı oranı D, kazanç G.
- Her ölçü için ızgara çözünürlüğü yakınsama eğrisi.

**Yinelenebilirlik:** Python + NumPy (+ SciPy), laptopta çalışır, tam zaman/bellek bütçesi raporlanır. Zenodo'ya açık veri/kod (review'daki gibi).

## 8. Yol haritası (kapılı)

| Faz | İçerik | Çıktı | Geçiş kapısı |
|---|---|---|---|
| 0 | Sistematik literatür taraması (Scopus/WoS/arXiv: "numerical damping" × "basin of attraction"; Newmark/HHT × multistability; ayrıca basin stability literatürü) | Boşluk raporu | Doğrudan çakışan çalışma yoksa devam. Varsa konumlandırma/ayrıştırma planı |
| 1 | Altyapı: referans çözücü, havza ölçüleri (Sb, Sbb, sınır kayması), yakınsama testleri | Doğrulanmış kod tabanı | Referans ve ölçü yakınsaması gösterildi |
| 2 | SDOF çift kuyu tam rejim haritası (ρ∞, h, c, ε), tüm integratör aileleri | Rejim haritası (RQ1, RQ4 kısmen) | H1, H4 doğrulandı, H5 en az 2 ailede |
| 3 | Mekanizma: D analitik/yarı-analitik, ρ∞\* tahmini, kazanç G | Tahmin kuralı (RQ2, RQ3) | ρ∞\* hatası ≤0.1 tüm güvenilir hücrelerde; H3 ölçeklenmesi gösterildi |
| 4 | Von Mises kafesi + 2 DOF bistabil model + **gerçek FE vakası (sığ kemer/kabuk snap-through, zorunlu)** | Yapısal doğrulama (RQ4) | Aynı etki ve tahmin kuralı yapısal modelde ve FE'de de geçerli |
| 5 | ~~Zorlamalı (periyodik çeker) genişleme~~ **Kapsam dışı** (devam makalesi; bkz. §11.3-8) | | |
| 6 | Pratik kural, tanı protokolü (DT3 somutlaştırma), vaka çalışması | Öneri bölümü (RQ5) | İki-ρ∞ testi gerçekten sapmaları yakalıyor |
| 7 | Yazım, açık veri/kod, kapak mektubu | Makale | |

**Kırılma noktaları (fikir zayıflarsa):**
- Faz 3 tutmazsa → makale **ölçüm/rejim makalesine** dönüşür (daha düşük etki, hâlâ yayımlanabilir; hedef dergi bir kademe aşağı).
- Faz 4'te yapısal modelde etki yoksa → "toy-model etkisi" riski. Etkinin hangi koşulda ortaya çıktığı için bir rejim sınırı makalesi (negatif sonuç dahil) yazılır.
- Faz 2'de etki yalnızca generalized-α'da çıkıyorsa → başlık ve iddia daraltılır.

## 9. Katkılar (hedeflenen)

1. İlk sistematik ρ∞–h–c havza rejim haritası.
2. Enerji-kaybı eşleşme kuralı ρ∞\*(h, c), doğrusal kuralın neden yetmediğinin açıklamasıyla.
3. Havza sınırında hata büyütme kazancının ölçüsü.
4. Yapısal (snap-through/bistabil) doğrulama.
5. Mühendis için pratik kural: ρ∞ ve h raporlanmalı; çok kararlı sistemde iki ρ∞ ile havza testi.
6. Açık, laptop ölçekli benchmark.

## 10. Riskler ve hakem itirazları

- **"Oyuncak model":** Faz 4 ile karşılanır (von Mises + 2 DOF + mümkünse FE).
- **"Havza ölçüsü ızgaraya bağlı / doygun":** basin entropy, sınır kayması ve çözünürlük yakınsaması.
- **"Referans güvenilir mi":** referans yakınsama testi, iki bağımsız referans çözücü.
- **"Yeni bir şey yok, bilinen sayısal sönüm":** vurgu sayısal sönümün *denge seçimine* etkisi ve nicel tahmin kuralıdır.
- **"Mühendislik önemi":** gerçek yapısal örnek ve pratik kural (ve gerçek problemde ne zaman önemli olduğunu gösteren rejim haritası, ağır sönümde etkinin yok olduğu dahil).
- **Fiziksel çok kararlılık ile karışma:** Zorlamasız sönümlü sistemde havza gerçek bir özellik; integratör etkisi, referans RK ile **aynı** başlangıç koşulunda farktır. Bu ayrım metinde açık olmalı (review'daki "genuine multistability" notu).

## 11. Yayın stratejisi (hedef: SCIE, hibrit dergi)

**Ölçütler:** SCIE/JCR indeksli, hibrit (abonelik + isteğe bağlı açık erişim), Q1–Q2, konu uyumu, makul karar süresi, yazarın APC yükü olmaması.

### 11.1 Dergi sıralaması

Sayılar aggregator sitelerden (LetPub/JournalCompass/Springer sayfaları, Kasım 2026) alındı. **Başvurmadan önce JCR'dan doğrulayın.**

| Sıra | Dergi | Yayıncı | IF (≈) | Model | Neden |
|---|---|---|---|---|---|
| 1 | **Computational Mechanics** | Springer | 4.5 | Hibrit, APC ≈ 3 990 USD | Kapsam: hesaplamalı mekanik. Karar ≈116 gün (medyan). TÜBİTAK ULAKBİM–Springer "Oku & Yayımla" anlaşması, WoS Q1/Q2 Springer hibrit dergilerinde **sorumlu yazar EKUAL üyesiyse APC'siz açık erişim** (yıllık kota var, 2025'te ağustosta bitmesi bekleniyordu → erken başvur/kabul zamanlamasına dikkat). Review'unuz zaten Springer'da (ACME) |
| 2 | **Int. J. Numer. Methods Eng. (IJNME)** | Wiley | 3.3 | Hibrit, APC ≈ 4 250 USD | Zaman integrasyonu literatürünün (Newmark, HHT, gen-α, enerji-momentum) tarihsel evi, hedef okuyucu tam burada. **Wiley–TÜBİTAK anlaşması hibrit dergilerde sorumlu yazara APC'siz açık erişim sağlıyor.** IF daha düşük, hakemler alanın uzmanı (sıkı ama adil) |
| 3 | **Journal of Sound and Vibration (JSV)** | Elsevier | ≈4.9 | Hibrit (abonelik varsayılan, APC isteğe bağlı) | Yapısal dinamik/titreşim okuyucusu. **Dikkat:** editör politikası, katkısı yalnızca ODE sistemlerinin matematiksel işlenmesi olan makaleleri dışlıyor → fiziksel/yapısal içgörü ve gerçek yapı vakası şart. Abonelik yolu APC'siz |
| 4 | **Computers & Structures** | Elsevier | ≈4–5 (doğrulayın) | Hibrit | Yapısal FE odaklı, snap-through/stabilite okuyucusu. FE vakası güçlüyse uygun |
| Hedef üstü | **CMAME** | Elsevier | 7.6 (Q1) | Hibrit, APC ≈ 4 830 USD | Çok prestijli, kabul ≈%20, "yöntem yeniliği" bekler. Yalnızca Faz 3–4 çok güçlü çıkar ve ρ∞ seçimi için yeni bir ilke/algoritma sunulabilirse (örn. havza-güvenli parametre seçimi) denenir |

**Elenen:** *Nonlinear Dynamics* (IF 5.7, Q1, konu uyumu çok iyi) artık **hibrit değil**: Springer, 1 Ocak 2027'den itibaren tamamen açık erişime geçeceğini, 11 Ağustos 2026 sonrası gelen başvurularda APC alınacağını duyurmuş. Hibrit ölçütünü karşılamıyor. İsterseniz APC desteği/waiver koşuluyla yedek olarak değerlendirilebilir.

**Doğrulanacaklar (kütüphane/BAP ile):**
- Hatay MKÜ'nün EKUAL üyeliği ve sorumlu yazar olarak anlaşma kapsamı (Springer ve Wiley).
- Elsevier için güncel TÜBİTAK anlaşması var mı, kapsam nedir? (Araştırmada net sonuç bulamadım.) Yoksa Elsevier dergilerinde abonelik yolu (APC'siz) kullanılabilir.
- Her derginin güncel JCR çeyrekliği (özellikle Computational Mechanics ve IJNME).

### 11.2 Önerilen ana strateji

**Önce Computational Mechanics, ardından IJNME; JSV/Computers & Structures yedek; CMAME yalnızca sonuçlar beklenenin üstündeyse.** Tek bir "ana" okuyucu kitlesi seçiyoruz: **yapısal dinamik/hesaplamalı mekanik analistleri ve FE kodu geliştiricileri**. Çünkü ρ∞ kararını bu kitle veriyor, mesaj onlara ulaşmalı. Dinamik sistem terminolojisi (havza, Sb) araç olarak kullanılacak, makalenin dili değil.

### 11.3 Editör ve hakemi ikna edecek tasarım kararları

1. **Katkı türü "gözlem" değil "ölçüt + doğrulama".** Editör için başlık/özet: *"…and a predictive criterion for damping selection"*. Saf parametre taraması "yeterince yeni değil" riski taşır. Enerji-kaybı eşleşme kuralı (H2) ve kazanç (H3) bu riski azaltan unsurlar.
2. **Oyuncak model itirazını baştan kapat.** Faz 4 (yapısal vaka) **zorunlu**, "opsiyonel" değil. En az: (i) von Mises kafesi, (ii) 2 DOF/iki modlu bistabil kiriş ve (iii) **gerçek bir FE vakası** (sığ kemer/kabuk snap-through; tercihen literatürden yerleşik bir benchmark ve açık kaynak FE kodu ile yeniden üretilebilir). Computational Mechanics/IJNME okuyucusu FE sonucu görmek ister; JSV ise fiziksel içgörü ister.
3. **Okuyucunun tanıdığı integratörler ve benchmarklar.** Newmark, HHT-α, generalized-α, Bathe, enerji-momentum (kontrol). Mümkünse review'da atıf yapılan yerleşik benchmarkların (Kuhl, Lavrenčič, Armero) en az birini yeniden kullan ve onlardan bir adım ileri git (havza ekseni).
4. **"Zaten biliniyor" itirazına cevap:** Sayısal sönümün doğruluğu/kararlılığı etkilediği bilinir; **denge seçimine etkisi ve nicel tahmin kuralı** bilinmiyor. Giriş, bu ayrımı bir tabloyla (bilinen / bilinmeyen) net koysun ve literatür taraması (Faz 0) bunu belgelesin.
5. **Pratik değer:** Çok kararlı problem analisti için üç satırlık kural (iki ρ∞ ile havza testi, h-yakınsaması, rapor edilecekler) + açık, laptop ölçekli kod. "Okuyucu bu makaleden ne alıp gidiyor?" sorusunun cevabı ilk sayfada olmalı.
6. **Hakem itirazı → yerleşik cevap tablosu** (§10'daki liste, makalede kısa bir "Validity" bölümüne dönüşür): referans yakınsaması, ızgara yakınsaması, ölçü seçimi, fiziksel çok kararlılıktan ayrım.
7. **Negatif sonucu da raporla:** Ağır sönümde etki yok, belirli h'de doygunluk var. Dürüst sınırlar hakem güveni kazandırır ve "kapsam abartısı" itirazını önler.
8. **Kapsamı dar tut.** Zorlamalı/periyodik çeker bölümü (Faz 5) ana makaleden **çıkarıldı**; ayrı bir devam makalesi olabilir. Uzunluk hedefi ≈ 25–30 sayfa.
9. **Açık veri/kod (Zenodo)** ve yinelenebilirlik: review'daki gibi. Hibrit dergilerin çoğu bunu olumlu karşılar.
10. **Önerilen hakemler** (editör formu için): review'da atıf yaptığımız zaman integrasyonu ve havza kararlılığı alanlarındaki yazarlar (örn. Kuhl, Erlicher, Lavrenčič, Bathe, Brzeski). İsimlerin güncel ve çıkar çatışması olmadığını doğrulayın.

### 11.4 Dergiye göre vurgu (aynı çalışmanın üç "paketi")

- **Computational Mechanics / IJNME:** zaman integrasyonu + FE vakası + ölçüt. Başlık "…integrators…structural dynamics…". Dinamik sistem terimleri minimum.
- **JSV:** fiziksel içgörü ve yapı vakası önde, ODE ayrıntısı arka planda; "enerji–faz–sönüm dengesi" hikâyesi.
- **CMAME (hedef üstü):** algoritmik katkı: havza-güvenli parametre seçimi / otomatik iki-ρ∞ doğrulama algoritması + FE vakası.

### 11.5 Göndermeden önceki kontrol listesi

- [ ] Faz 0 literatür taraması bitti, çakışma yok ya da konumlandırma net.
- [ ] Referans ve ızgara yakınsaması rapor edildi.
- [ ] En az bir gerçek FE vakası + iki yapısal indirgenmiş model.
- [ ] En az 3 integratör ailesi.
- [ ] Ölçüt H2 hatası ≤0.1 (güvenilir hücrelerde) ve kazanç ölçeklemesi gösterildi.
- [ ] Kod ve veri Zenodo'da.
- [ ] Kapak mektubu: "ne yeni, neden bu dergi, okuyucu ne kazanıyor".
- [ ] EKUAL/APC durumu doğrulandı, JCR çeyrekliği doğrulandı.
- [ ] ACME review ile atıf/özdüzeltme: review'u atıf yap, kendi içinde tekrar yok.

## 12. Bu oturumda olgunlaştırılacak açık noktalar

1. Faz 0 literatür taramasının kapsamı ve anahtar kelimeler (kullanıcı kütüphane erişimiyle Scopus/WoS).
2. D'nin analitik türetimi mümkün mü? (averaging / harmonik denge) Yoksa tamamen sayısal mı?
3. Yapısal vakanın seçimi: von Mises + 2 modlu kiriş Galerkin yeterli mi, FE'ye ihtiyaç var mı?
4. Hangi havza ölçüsü ana ölçü olacak (yüzde fark kolay okunur, basin entropy daha dirençli).
5. İntegratör kapsamı: kaç aile, hangi gruplanma (dissipatif ailelerin ortak bir "ρ∞" ölçeğine eşlenmesi: spektral yarıçap).
6. Zorlamalı/periyodik çeker bölümü kapsamda mı?
7. Başlık ve vurgu: "equilibrium selection" mı, "basin boundary shift" mı, "design rule" mı?

## 13. Dosyalar

- `pilots/pilot1.py`: ilk ρ∞×h ızgarası (c=0.15), tablo 5a.
- `pilots/pilot3.py`: c taraması + doğrusal tahmin karşılaştırması (5b; 5c-doğrusal).
- `pilots/pilot4.py`: enerji-kaybı eşleşmesi (5c-enerji).
- Bağımlılık: `..\zenodo_package_v1\benchmarks\genalpha.py` (pilotlardaki `sys.path` satırı yerel yolu gösterir, makale klasöründe düzeltilmeli).
- Çalıştırma: Python 3.12 + numpy (bu makinede `%LOCALAPPDATA%\Programs\Python\Python312\python.exe`).
