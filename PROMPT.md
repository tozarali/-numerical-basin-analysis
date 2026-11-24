# Yeni oturum başlangıç promptu

Aşağıdaki bloğu yeni oturumun ilk mesajı olarak yapıştırın. Oturumu bu klasörde (`paper2_idea`) veya ondan oluşturacağınız yeni proje klasöründe başlatın.

---

## Kopyalanacak prompt

Ben Ali Tozar (Hatay Mustafa Kemal Üniversitesi, Fizik). Önceki bir oturumda, Computer Methods/ACME'de "under review" olan review makalemden ("When the Solver Becomes the System", açık problem #2: numerical damping and basin boundaries) yeni bir **araştırma makalesi** fikri geliştirdik. Fikir `IDEA.md` içinde, ön deneme kodları `pilots/` içinde. Bu oturumda bu fikri **asıl makale projesine** dönüştüreceğiz.

**Önce şunları oku (sırayla):**
1. `IDEA.md` (tüm fikir, hipotezler H1–H5, ön deneme sonuçları §5a–f, yol haritası Faz 0–7, dergi stratejisi §11, açık noktalar §12)
2. `pilots/pilot1.py … pilot7_collapse.py` (ön denemeler; `pilot5_families.py` ve `pilot7_collapse.py` en güncel olanlar)
3. Gerekirse review kaynakları: `..\submission_v1.5\submission_v1.5\sections\04_time_integration.tex` (Sec. `time:rhoinf`) ve `..\zenodo_package_v1\benchmarks\genalpha.py`, `fig06_b1.py`

**Özet iddia:** Zaman integratörünün sayısal sönümü (ρ∞) çok kararlı yapılarda hangi dengeye yerleşileceğini değiştirir. Etkiyi ρ∞ değil, fiziksel frekanstaki çevrim-başına sayısal enerji kaybı oranı D örgütler (Spearman 0.81–0.91, dört ailede). Küçük D hatası havza sınırında büyür.

**Amaç:** SCIE indeksli, hibrit bir dergide yayın. Öncelik: Computational Mechanics, sonra IJNME; yedek JSV / Computers & Structures; CMAME ancak sonuçlar çok güçlüyse. Editör ilk bakışta ilgilenmeli, hakemler itiraz edemeyecek kadar sağlam doğrulama olmalı. Hedef okuyucu: yapısal dinamik analistleri ve FE geliştiricileri.

**Ortam:** Windows. Python 3.12 `%LOCALAPPDATA%\Programs\Python\Python312\python.exe` konumunda, numpy ve scipy kurulu (matplotlib kurulmamış olabilir, gerekirse kur; pip ağda takılabilir, `--default-timeout 20 --retries 2` kullan). Tüm koşular laptop ölçeğinde olmalı.

**Bilinen zayıf noktalar (öncelikle çöz):**
- Ortak korelasyon doygun (~%40) noktaları içeriyor; doygunluk dışı alt kümede ve basin entropy / sınır kayması ölçüleriyle yinele.
- D eşleşmesi optimumun biraz üstünü işaret ediyor (≤0.1 sapma); nedenini incele.
- 2 DOF modelinde mod bağlaşımı zayıf; güçlü etkileşimli gerçek yapısal vaka (von Mises kafesi, iki modlu bistabil kiriş) ve en az bir FE vakası (sığ kemer/kabuk snap-through, açık kaynak kodla) gerekli.
- Referans çözücü ve ızgara yakınsaması kanıtlanmalı.
- Literatür boşluğu sadece hızlı web taramasıyla kontrol edildi; sistematik tarama gerekli.

**Bu oturumda yapmanı istediğim:**
1. `IDEA.md`'yi oku ve kısa bir durum özeti + kendi eleştirin (hakem gözüyle en zayıf 3 nokta) ver.
2. Yeni proje klasör yapısını öner ve kur (kod, veri, figürler, LaTeX, Zenodo paketi). Pilot kodlarını temiz, test edilebilir bir pakete dönüştür (yerel `sys.path` bağımlılığını kaldır).
3. Faz 0 ve Faz 1'i yürüt: literatür tarama planı/anahtar kelimeler, referans çözücü + havza ölçüleri (yanlış kuyu oranı, basin entropy, sınır kayması) + yakınsama testleri.
4. Her faz sonunda `IDEA.md` §8'deki geçiş kapısını açıkça değerlendir; kapı geçilmezse alternatif yolu (ölçüm/rejim makalesi) öner.
5. Sıradaki büyük algoritmik katkı: D'nin analitik/yarı-analitik türetimi (çarpan analizi + genlik ortalaması), böylece D ücretsiz hesaplanabilsin.

**Çalışma kuralları:** Sonuçları şişirme; negatif sonuçları ve sınırları açıkça raporla. Her iddiayı bir koşuyla destekle. Dosyaları doğrudan oluştur/düzenle, belirsiz kararlarda bana sor. Yanıtlar Türkçe ve kısa olsun, makale metni İngilizce.
