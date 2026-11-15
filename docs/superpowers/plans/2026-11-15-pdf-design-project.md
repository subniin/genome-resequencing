# 설계 프로젝트 조건 맞춤 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `알고리즘응용_설계.pdf` 설계 프로젝트 조건(N 최대 3×10^9, M 2×10^7~2×10^8, L 32~100, D, repeat)에 맞춰 기존 de novo 조립 실험을 확장한다. M=2×10^7까지 이 PC(16GB)에서 실측하고, read 오류와 repeat 실험을 추가하고, README를 발표 슬라이드 구성(62쪽)으로 다시 쓴다.

**Architecture:** 기존 C++17 CLI 도구(`gen_reference`, `gen_reads`, `ssake`, `greedy_overlap`, `evaluate`)와 Python 실험 스크립트 구조를 유지한다. 메모리 병목인 read 저장과 해시 색인을 2-bit 압축 read + 정렬 배열 색인으로 바꾸고, 평가 도구의 문자열 k-mer 맵을 정수 코드 정렬 배열로 바꾼다. SSAKE 알고리즘 자체는 단순화하지 않고 유지한다.

**Tech Stack:** C++17 (g++ MinGW-w64 16.2.0, `-O2 -std=c++17 -static`), Python 3.11 + psutil + matplotlib, `mt19937_64` 고정 seed.

**Spec:** `README.md` 3.1장 + `C:\Users\coolc\Documents\Vault\Areas\학교\족보및교재및교안\3-2\알고리즘\팀플\알고리즘응용_설계.pdf` 60~64쪽. 이전 계획: `docs/superpowers/plans/2026-11-14-denovo-resequencing.md`.

## Global Constraints

- 문제: 길이 L인 short read M개로 그 read가 나온 길이 N의 원본 서열을 복원한다 (PDF 60쪽)
- 조건 범위 (PDF 61쪽): N 최대 3,000,000,000 / M 약 20 million ~ 200 million / L 약 32~100 / D = read당 mismatch 수 / repeat 있을 수 있음 / single-end 또는 paired-end / de novo 또는 reference
- 이 저장소 선택: de novo, single-end, 정방향 read만. paired-end와 역상보는 future work로 README에 적는다
- 실험 L은 32 / 50 / 75 / 100만 쓴다. 기본 L=50, o=⌈L/2⌉=25, coverage 10, SNP 0.1%(PDF 6쪽 "0.1% genetic variants"), D=2
- 실측 최대 규모: M=2×10^7 (N=10^8, L=50, coverage 10). M=2×10^8은 실측 B/read로 외삽해 README에 적는다
- SSAKE 피크 메모리 목표: M=2×10^7, L=50에서 2.0GB 이하
- read id는 `uint32_t`. M ≥ 2^32이면 오류로 중단
- reference 생성 k(K_GEN)=16, 평가 anchor k(K_EVAL)=24
- 난수는 `mt19937_64` + 고정 seed. `rand()` 금지
- 입출력 형식은 기존 유지: read 한 줄에 하나, contig FASTA, 정답 SNP TSV(`pos<TAB>ref<TAB>alt`, 0-based)
- reference는 평가에서만 사용. 조립기는 reads만 입력으로 받는다
- 단일 스레드, Ryzen 5 7500F / DDR5 16GB / Windows 11에서 측정
- 발표 필수 항목(PDF 62쪽): Introduction(문제 정의, 데이터, 입출력), Benchmark, Algorithm, Result(머신 사양, 시간·공간 복잡도 big-O 또는 실측, 벤치마크 비교), Future work(장단점, 개선 방안)
- 실행은 `main`에서 새 브랜치(worktree)로. 커밋은 task마다
- 알고리즘 난이도는 낮추지 않는다: 2-bit 압축, 정렬 배열 색인, 조건부 다수결 연장(Task 8)을 구현한다. de Bruijn graph와 BWT는 future work
- de novo 유지: `ssake`, `greedy_overlap`은 reads만 입력. reference는 `evaluate`에서만 읽는다

## Review Focus

1. reads 파일에 길이가 다른 read, 빈 줄, CRLF가 섞임: 빈 줄과 CR은 무시, 길이가 다르면 줄 번호와 함께 오류로 중단 (Task 1)
2. o > 32 (L=75 → o 38, L=100 → o 50): 색인 key는 앞 32염기만 쓰고 겹침 전체는 mismatch ≤ D로 검증해 L=100 조립이 정상 (Task 4)
3. read 앞쪽 key 구간에 시퀀싱 오류: 그 read는 색인으로 못 찾지만 크래시 없이 contig가 끊길 뿐 (Task 4)
4. repeat 길이 ≥ L: 조립기는 멈추거나 무한 반복하지 않고, 평가는 repeat k-mer를 anchor에서 빼고 misassembly를 센다 (Task 2, 5)
5. read 0개 / 1개: 빈 FASTA / read 그대로 contig 1개 (Task 4)

---

### Task 1: 2-bit 압축 read 저장 (`common.hpp`)

**Files:**
- Modify: `src/common.hpp` (`base_index`, `kmer_code`를 `src/gen_reference.cpp`에서 옮기고 추가 구조체 정의)
- Modify: `src/gen_reference.cpp` (옮긴 두 함수 삭제)
- Test: `tests/test_common.cpp`

**Interfaces:**
- Produces:
  - `inline int base_index(char c)` (A0 C1 G2 T3, 기존과 동일)
  - `inline uint64_t kmer_code(const std::string& seq, size_t end, int k)` — `end`에서 끝나는 k-mer, k ∈ [1, 32]
  - `struct PackedReads { int L = 0; std::vector<uint64_t> data; size_t count = 0; size_t size() const; void push(const std::string& r); char base(size_t id, int i) const; uint64_t code(size_t id, int start, int k) const; std::string str(size_t id) const; }` — read 하나에 `ceil(2L/64)`개 word. 첫 `push`가 L을 정하고, 다른 길이는 `std::invalid_argument`, ACGT 외 문자는 `std::runtime_error`
  - `inline PackedReads load_reads_packed(const std::string& path)` — 한 줄씩 스트리밍. 빈 줄 무시, 끝 CR 제거, 대문자화. 길이 불일치 시 메시지에 줄 번호 포함
  - `inline void dedupe_packed(PackedReads& r)` — 첫 사본만 남기고 입력 순서 유지 (`dedupe_reads`와 같은 의미)

- [ ] **Step 1: Write the failing test** (`tests/test_common.cpp`에 추가)

```cpp
// PackedReads round trip and k-mer codes match the string versions.
PackedReads p;
std::string a = "ACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAAC"; // 68 bases, crosses a word
std::string b = std::string(68, 'G');
p.push(a); p.push(b); p.push(a);
assert(p.size() == 3 && p.L == 68);
assert(p.str(0) == a && p.str(1) == b);
assert(p.base(0, 67) == a[67]);
assert(p.code(0, 30, 32) == kmer_code(a, 61, 32));
bool threw = false;
try { p.push("ACGT"); } catch (const std::invalid_argument&) { threw = true; }
assert(threw);
dedupe_packed(p);
assert(p.size() == 2 && p.str(0) == a && p.str(1) == b);

// load_reads_packed: skips blank lines, strips CR, rejects mixed lengths.
{ std::ofstream f("build/tmp_reads.txt", std::ios::binary); f << "acgt\r\n\r\nTTGG\n"; }
PackedReads q = load_reads_packed("build/tmp_reads.txt");
assert(q.size() == 2 && q.str(0) == "ACGT" && q.str(1) == "TTGG");
{ std::ofstream f("build/tmp_reads.txt", std::ios::binary); f << "ACGT\nACG\n"; }
threw = false;
try { load_reads_packed("build/tmp_reads.txt"); } catch (const std::invalid_argument&) { threw = true; }
assert(threw);
std::remove("build/tmp_reads.txt");
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_common.cpp -o build/test_common.exe`
Expected: 컴파일 오류 `'PackedReads' was not declared`

- [ ] **Step 3: Implement** `PackedReads`, `load_reads_packed`, `dedupe_packed`, 옮긴 `base_index`/`kmer_code` in `src/common.hpp`

`dedupe_packed`는 id 배열을 (read word 비교) 기준 정렬해 중복을 찾고, 각 그룹의 최소 id만 남긴 뒤 원래 순서로 압축한다. 임시 메모리는 id 배열(4B/read)만 쓴다.

- [ ] **Step 4: Run all tests to verify they pass**

Run: `for t in common gen_reference gen_reads greedy_overlap ssake evaluate; do g++ -O2 -std=c++17 -static tests/test_$t.cpp -o build/test_$t.exe && ./build/test_$t.exe || echo FAIL $t; done`
Expected: `FAIL` 출력 없음

- [ ] **Step 5: Commit**

```bash
git add src/common.hpp src/gen_reference.cpp tests/test_common.cpp
git commit -m "feat: add 2-bit packed read storage"
```

---

### Task 2: reference 생성기 k=16 지원과 repeat 삽입

**Files:**
- Modify: `src/gen_reference.cpp`
- Test: `tests/test_gen_reference.cpp`

**Interfaces:**
- Consumes: Task 1 `kmer_code`
- Produces:
  - `gen_unique_kmer_seq(size_t N, int k, std::mt19937_64& rng)` — 허용 k를 [1, 16]으로 넓힘 (4^16 bit = 512MB)
  - `inline void insert_repeats(std::string& seq, int count, int len, std::mt19937_64& rng)` — `count`번, 무작위 원본 구간 `[s, s+len)`을 겹치지 않는 무작위 위치 `[t, t+len)`에 덮어쓴다(|s−t| ≥ len). `len > seq.size()/2`이면 `std::invalid_argument`
  - CLI 플래그 추가: `--repeats C --repeat-len R` (기본 0, 0)

- [ ] **Step 1: Write the failing test**

```cpp
auto rng = make_rng(9);
std::string s = gen_unique_kmer_seq(20000, 16, rng);
std::string before = s;
insert_repeats(s, 0, 100, rng);
assert(s == before);
insert_repeats(s, 3, 100, rng);
assert(s.size() == 20000);
for (char c : s) assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
// Each copy makes at least one 100-mer appear twice.
std::unordered_map<std::string, int> cnt;
for (size_t i = 0; i + 100 <= s.size(); ++i) cnt[s.substr(i, 100)]++;
int dup = 0;
for (auto& kv : cnt) dup += kv.second > 1;
assert(dup >= 3);
bool threw = false;
try { insert_repeats(s, 1, 15000, rng); } catch (const std::invalid_argument&) { threw = true; }
assert(threw);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_gen_reference.cpp -o build/test_gen_reference.exe`
Expected: 컴파일 오류 `'insert_repeats' was not declared` (k=16은 런타임 `k must be in [1, 15]`)

- [ ] **Step 3: Implement** `insert_repeats`, k 상한 16, CLI 플래그 in `src/gen_reference.cpp`

- [ ] **Step 4: Run test to verify it passes**

Run: `./build/test_gen_reference.exe; echo $?`
Expected: `0`

- [ ] **Step 5: Commit**

```bash
git add src/gen_reference.cpp tests/test_gen_reference.cpp
git commit -m "feat: allow k=16 and insert repeats in reference"
```

---

### Task 3: read 생성기 스트리밍, 희소 SNP 샘플링, 시퀀싱 오류

**Files:**
- Modify: `src/gen_reads.cpp`
- Test: `tests/test_gen_reads.cpp`

**Interfaces:**
- Produces:
  - `insert_snps(std::string& seq, int count, std::mt19937_64& rng)` — 시그니처 유지. 길이 N 위치 배열 대신 `unordered_set<size_t>` 거절 샘플링 (N=10^8에서 800MB 배열 제거)
  - `template <class Emit> void sample_reads(const std::string& sample, int L, double coverage, double err_rate, std::mt19937_64& rng, Emit emit)` — M = llround(N·coverage/L) 개, 시작 위치 균등 [0, N−L]. 각 염기를 확률 `err_rate`로 다른 염기로 치환. `emit(size_t start, const std::string& read)` 호출
  - `make_dataset(ref, L, coverage, snp_count, seed, double err_rate = 0.0)` — 내부에서 `insert_snps` 후 `sample_reads`로 `reads` 채움
  - CLI: `--err P` 추가 (기본 0). main은 `Dataset.reads`에 모으지 않고 `sample_reads`로 파일에 바로 쓴다. `make_dataset`과 같은 rng 호출 순서
  - `err_rate`가 [0, 1) 밖이면 `std::invalid_argument`

- [ ] **Step 1: Write the failing test** (기존 테스트 유지, 아래 추가)

```cpp
// err_rate 0: every read equals the sample at its start.
{
    auto r = make_rng(4);
    std::string sample = ref; // 10000 bases from the existing test
    size_t n = 0;
    sample_reads(sample, 50, 10, 0.0, r, [&](size_t st, const std::string& rd) {
        assert(rd == sample.substr(st, 50)); ++n; });
    assert(n == 2000);
}
// err_rate 0.05: mismatch fraction within [0.04, 0.06].
{
    auto r = make_rng(4);
    size_t mm = 0, total = 0;
    sample_reads(ref, 50, 10, 0.05, r, [&](size_t st, const std::string& rd) {
        for (int i = 0; i < 50; ++i) mm += rd[i] != ref[st + i];
        total += 50; });
    double f = (double)mm / total;
    assert(f > 0.04 && f < 0.06);
}
bool thr = false;
try { make_dataset(ref, 50, 10, 0, 1, 1.5); } catch (const std::invalid_argument&) { thr = true; }
assert(thr);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_gen_reads.cpp -o build/test_gen_reads.exe`
Expected: 컴파일 오류 `'sample_reads' was not declared`

- [ ] **Step 3: Implement** `sample_reads`, 희소 `insert_snps`, `err_rate` 인자와 `--err` in `src/gen_reads.cpp`

- [ ] **Step 4: Run test to verify it passes**

Run: `./build/test_gen_reads.exe; echo $?`
Expected: `0`

- [ ] **Step 5: Commit**

```bash
git add src/gen_reads.cpp tests/test_gen_reads.cpp
git commit -m "feat: stream reads, sparse SNP sampling, read errors"
```

---

### Task 4: SSAKE를 압축 read + 정렬 배열 색인으로 교체

**Files:**
- Modify: `src/ssake.cpp`
- Test: `tests/test_ssake.cpp`

**Interfaces:**
- Consumes: Task 1 `PackedReads`, `load_reads_packed`, `dedupe_packed`, `kmer_code`
- Produces:
  - `struct KeyIndex { std::vector<uint64_t> keys; std::vector<uint32_t> ids; std::pair<size_t, size_t> range(uint64_t key) const; }` — keys 오름차순, `range`는 `equal_range` 결과 인덱스 구간
  - `inline KeyIndex build_prefix_index(const PackedReads& r, int kk)` — key = read 앞 kk염기 코드
  - `inline KeyIndex build_suffix_index(const PackedReads& r, int kk)` — key = read 뒤 kk염기 코드
  - `inline std::vector<std::string> assemble_ssake(PackedReads reads, int o, int D)` — 내부에서 `dedupe_packed`. kk = min(o, 32)
  - `inline std::vector<std::string> assemble_ssake(const std::vector<std::string>& input, int o, int D)` — 기존 시그니처 유지, PackedReads로 바꿔 위 함수 호출. 빈 입력은 빈 결과
  - CLI: `load_reads_packed` → `assemble_ssake(std::move(reads), o, D)` → `write_fasta`

동작은 기존 `assemble_ssake`와 같다: seed는 입력 순서, 3' 연장 후 5' 연장, 겹침 k를 `min(contig 길이, L−1)`부터 o까지 내려가며 key로 후보를 찾고 겹침 전체 mismatch ≤ D인 첫 미사용 read를 붙인다. 3' key는 contig의 마지막 k염기 중 앞 kk염기(`kmer_code(contig, contig.size()-k+kk-1, kk)`), 5' key는 contig 앞 k염기 중 뒤 kk염기(`kmer_code(contig, k-1, kk)`). read 수가 `UINT32_MAX` 이상이면 `std::length_error`.

- [ ] **Step 1: Write the failing test** (기존 테스트 유지, 아래 추가)

```cpp
// L=100, o=50 (> 32): key is capped, full overlap still checked.
{
    auto r = make_rng(21);
    std::string g;
    for (int i = 0; i < 1000; ++i) g.push_back("ACGT"[r() % 4]);
    std::vector<std::string> rd;
    for (int st = 0; st + 100 <= 1000; st += 30) rd.push_back(g.substr(st, 100));
    rd.push_back(g.substr(900, 100));
    std::vector<std::string> c = assemble_ssake(rd, 50, 2);
    assert(c.size() == 1 && c[0] == g);
}
// Error inside a read's key region: no crash, contigs stay substrings of the genome
// except at most D mismatches each.
{
    auto r = make_rng(22);
    std::string g;
    for (int i = 0; i < 400; ++i) g.push_back("ACGT"[r() % 4]);
    std::vector<std::string> rd;
    for (int st = 0; st + 50 <= 400; st += 10) rd.push_back(g.substr(st, 50));
    rd[5][3] = rd[5][3] == 'A' ? 'C' : 'A';
    std::vector<std::string> c = assemble_ssake(rd, 25, 2);
    assert(!c.empty());
}
// Empty and single-read inputs.
assert(assemble_ssake(std::vector<std::string>{}, 15, 2).empty());
{
    std::vector<std::string> one = {std::string(30, 'A')};
    std::vector<std::string> c = assemble_ssake(one, 15, 2);
    assert(c.size() == 1 && c[0] == one[0]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_ssake.cpp -o build/test_ssake.exe && ./build/test_ssake.exe`
Expected: L=100 케이스에서 `std::invalid_argument` 없이 통과할 수도 있으므로, 컴파일 오류(`'build_prefix_index' was not declared`)가 나도록 아래 한 줄도 추가한다:
```cpp
{ PackedReads p; p.push("ACGTACGT"); p.push("CGTACGTT"); KeyIndex ix = build_prefix_index(p, 4); auto rg = ix.range(kmer_code("CGTA", 3, 4)); assert(rg.second - rg.first == 1 && ix.ids[rg.first] == 1); }
```

- [ ] **Step 3: Implement** `KeyIndex`, 두 build 함수, 두 `assemble_ssake` 오버로드, CLI in `src/ssake.cpp`

- [ ] **Step 4: Run all tests**

Run: Task 1 Step 4 명령
Expected: `FAIL` 출력 없음

- [ ] **Step 5: 메모리 확인 (M=2×10^6, L=50)**

Run: `./build/gen_reference.exe --N 10000000 --k 16 --seed 1 --out build/ref10m.txt && ./build/gen_reads.exe --ref build/ref10m.txt --L 50 --coverage 10 --snps 10000 --seed 1 --out build/r10m && python -c "import psutil,subprocess,time;p=subprocess.Popen(['build/ssake.exe','--reads','build/r10m/reads.txt','--o','25','--D','2','--out','build/r10m/c.fa']);q=psutil.Process(p.pid);m=0
while p.poll() is None: m=max(m,q.memory_info().rss);time.sleep(0.02)
print(m/2e6,'B/read')"`
Expected: 100 B/read 이하 (이전 구현 357 B/read)

- [ ] **Step 6: Commit**

```bash
git add src/ssake.cpp tests/test_ssake.cpp
git commit -m "perf: packed reads and sorted key index in SSAKE"
```

---

### Task 5: 평가 도구 정수 k-mer anchor (k=24)와 10M 이상치 확인

**Files:**
- Modify: `src/evaluate.cpp`
- Test: `tests/test_evaluate.cpp`

**Interfaces:**
- Consumes: Task 1 `kmer_code`
- Produces: `evaluate_contigs(contigs, ref, sample, truth, int k)` — 시그니처 유지, k ∈ [1, 32] 밖이면 `std::invalid_argument`. reference k-mer는 `(code, pos)` 정렬 배열, 두 번 이상 나온 code는 제외. contig k-mer는 rolling code + 이진 탐색. CLI `--k` 기본값 24

가설: N=10M에서 misassembly 8757과 recall 0.905는 평가 artifact다. k=14에서는 SNP가 바뀐 14-mer가 reference의 다른 고유 14-mer와 우연히 같을 확률이 약 N/4^14 ≈ 3.7%라 연속 spurious anchor가 생긴다.

- [ ] **Step 1: Write the failing test** (기존 테스트 유지, 아래 추가)

```cpp
// SNP turns "X C G" into "X A G" where "X A G" exists elsewhere in ref (X = 13 bases).
// k=14 sees two consecutive spurious anchors; k=24 does not.
{
    auto r = make_rng(31);
    std::string X, left, mid, right;
    for (int i = 0; i < 13; ++i) X.push_back("ACGT"[r() % 4]);
    for (int i = 0; i < 300; ++i) left.push_back("ACGT"[r() % 4]);
    for (int i = 0; i < 300; ++i) mid.push_back("ACGT"[r() % 4]);
    for (int i = 0; i < 300; ++i) right.push_back("ACGT"[r() % 4]);
    std::string ref = left + X + "AG" + mid + X + "CG" + right;
    size_t p = left.size() + 13 + 2 + mid.size() + 13; // the C
    std::string sample = ref;
    sample[p] = 'A';
    std::vector<TruthSnp> truth = {{p, 'A'}};
    Metrics m14 = evaluate_contigs({sample}, ref, sample, truth, 14);
    Metrics m24 = evaluate_contigs({sample}, ref, sample, truth, 24);
    assert(m14.misassemblies >= 1);
    assert(m24.misassemblies == 0);
    assert(std::fabs(m24.snp_recall - 1.0) < 1e-9);
    bool thr = false;
    try { evaluate_contigs({sample}, ref, sample, truth, 33); } catch (const std::invalid_argument&) { thr = true; }
    assert(thr);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_evaluate.cpp -o build/test_evaluate.exe && ./build/test_evaluate.exe`
Expected: `k=33` 케이스 assertion 실패 (현재는 예외 없음)

- [ ] **Step 3: Implement** 정수 anchor in `src/evaluate.cpp`

- [ ] **Step 4: Run all tests**

Run: Task 1 Step 4 명령
Expected: `FAIL` 출력 없음

- [ ] **Step 5: 10M 이상치 재확인**

이전 조건(N=10M, K_GEN=14, L=30, coverage 10, SNP 0.5%, o=15)으로 contig를 다시 만들고 `--k 14`와 `--k 24`로 각각 평가한다.
Run: `./build/gen_reference.exe --N 10000000 --k 14 --seed 1 --out build/ref10m14.txt && ./build/gen_reads.exe --ref build/ref10m14.txt --L 30 --coverage 10 --snps 50000 --seed 1 --out build/old10m && ./build/ssake.exe --reads build/old10m/reads.txt --o 15 --D 2 --out build/old10m/c.fa && for k in 14 24; do ./build/evaluate.exe --contigs build/old10m/c.fa --ref build/ref10m14.txt --sample build/old10m/sample.txt --truth build/old10m/truth_snps.tsv --k $k; done`
Expected: k=24에서 misassemblies가 k=14보다 크게 줄고 recall이 올라감. 결과 두 줄을 커밋 메시지 본문과 Task 7 README "평가 수정 이력"에 그대로 적는다. 줄지 않으면 가설 기각으로 적고 원인 조사는 future work로 남긴다

- [ ] **Step 6: Commit**

```bash
git add src/evaluate.cpp tests/test_evaluate.cpp
git commit -m "fix: integer k-mer anchors in evaluator, default k=24"
```

---

### Task 6: 실험 스크립트와 그래프를 설계 조건으로 갱신

**Files:**
- Modify: `scripts/run_experiments.py`
- Modify: `scripts/plot.py`

**Interfaces:**
- Consumes: Task 2 `--repeats/--repeat-len`, Task 3 `--err`, Task 5 `--k`
- Produces:
  - 상수 `SEED=1`, `K_GEN=16`, `K_EVAL=24`, 기본 `L=50`, `COV=10`, `D=2`
  - `run_case(step, variable, value, tool, N, L, coverage, snp_count, o, D=2, err=0.0, repeats=0, repeat_len=0)` — reference 파일명에 `N`, `repeats`, `repeat_len`, run 태그에 `err`, `D` 포함
  - `COLUMNS`에 `"D", "err", "repeats", "repeat_len"` 추가
  - `plot.py`: `plot_step1`(x축 M, log-log 시간·메모리), `plot_step2`(LW 곡선, L=50 기준), `plot_step3`, `plot_step4`(L별 contig 수·N50), `plot_step5`(err별 accuracy·precision), `plot_step6`(repeat_len별 misassembly)

단계 정의 (SNP 수 = round(N × 비율)):

| 단계 | 고정 | 변경 | 도구 |
|---|---|---|---|
| 0 | N 100,000, L 50, cov 10, SNP 0.1%, o 25, D 2 | - | 둘 다 |
| 1 | L 50, cov 10, SNP 0.1%, o 25 | M = 2×10^4, 2×10^5, 2×10^6, 2×10^7 (N = 5M) | ssake 전부, greedy는 M=2×10^4만 |
| 2 | 단계 0 | coverage 5, 10, 15 | 둘 다 |
| 3 | 단계 0 | SNP 0.1%, 0.5%, 1.0% | 둘 다 |
| 4 | 단계 0 | L 32, 50, 75, 100 (o=⌈L/2⌉); L 50에서 o 20, 25, 35 | 둘 다 |
| 5 | 단계 0 | err 0, 0.005, 0.01, 0.02 (D 2); err 0.01에서 D 0, 1, 2, 3 | 둘 다 |
| 6 | 단계 0 | repeats 20개, repeat_len 25, 50, 100, 200 | 둘 다 |

- [ ] **Step 1: 스크립트 수정 후 단계 0만 실행**

Run: `python scripts/run_experiments.py --steps 0`
Expected: `results/step0.csv` 2행, 두 도구 모두 `snp_recall` ≥ 0.99, `misassemblies` 0

- [ ] **Step 2: 그래프 함수 확인**

Run: `python scripts/plot.py`
Expected: 오류 없이 종료. 없는 step CSV는 건너뛴다는 메시지만 출력

- [ ] **Step 3: Commit**

```bash
git add scripts/run_experiments.py scripts/plot.py
git commit -m "feat: experiment steps 0-6 for design-project conditions"
```

---

### Task 7: 전체 실험 실행과 README 재작성 (발표 구성)

**Files:**
- Modify: `README.md`
- Create: `results/step0.csv` ~ `results/step6.csv`, `results/*.png` (`git add -f`)

- [ ] **Step 1: 전체 실험 실행**

Run: `python -u scripts/run_experiments.py --steps 0 1 2 3 4 5 6 && python scripts/plot.py`
Expected: `wrote results/step6.csv`까지 출력. 단계 1 M=2×10^7 행 `peak_mem_mb` ≤ 2000

- [ ] **Step 2: README를 PDF 62쪽 구성으로 다시 쓰기**

장 순서와 각 장에 반드시 들어갈 것:
1. Introduction — 문제 정의(PDF 60쪽 문장), 조건 표(현재 3.1장 유지), 데이터 생성 방법(고유 16-mer reference, repeat 삽입, SNP, read 오류), 입력·출력 형식
2. Benchmark — Greedy Overlap 동작과 선택 이유
3. Algorithm — SSAKE 계열 seed-and-extend, 2-bit 압축 read, 정렬 배열 색인(key 최대 32염기), 양방향 연장, (Task 8 실행 시) 다수결 연장
4. Result — 머신 사양, 복잡도 표(아래), 단계 0~6 표와 그래프, 벤치마크 비교, M=2×10^8 메모리 외삽(실측 B/read × 2×10^8), 평가 수정 이력(Task 5 Step 5 결과)
5. Future work — 장단점, 개선 방안(paired-end, 역상보, de Bruijn graph 비교, 디스크 기반 색인)

복잡도 표 (M = read 수, L = read 길이, o = 최소 겹침):

| 도구 | 시간 | 공간 |
|---|---|---|
| SSAKE (이 구현) | O(M log M + M·(L−o)·(log M + L)) | O(M·L/4 + 24M) bytes |
| Greedy Overlap | O(M²·L²) | O(M²) (후보 쌍) |

결과 해석은 CSV 수치만 근거로 쓰고, 확인하지 않은 원인은 "가설"로 표시한다.

- [ ] **Step 3: Commit**

```bash
git add README.md
git add -f results/
git commit -m "docs: rewrite README for design-project presentation; add results"
```


---

### Task 8 (조건부): 다수결 연장

**실행 조건:** Task 7 단계 5에서 err 0.01, D 2의 SSAKE `accuracy` < 0.999 또는 `snp_precision` < 0.9. 조건을 만족하지 않으면 건너뛰고 README future work에 "필요 없었음"과 수치를 적는다.

**Files:**
- Modify: `src/ssake.cpp`
- Test: `tests/test_ssake.cpp`

**Interfaces:**
- Produces: `assemble_ssake(PackedReads reads, int o, int D, bool vote)` — `vote=false`는 Task 4 동작. `vote=true`면 3' 연장 한 번에 겹침 ≥ o, mismatch ≤ D인 미사용 후보를 모두 모으고, 겹침 바로 다음 염기의 다수결(후보 ≥ 2, 최다 염기 비율 ≥ 2/3)로 한 염기씩 늘린다. 후보가 그 위치를 지나면 사용 처리. 조건을 못 채우면 연장 중단. 5'도 대칭. CLI `--vote 0|1`

- [ ] **Step 1: Write the failing test**

```cpp
// Coverage 6 reads of a 400-base genome; one read carries an error at its 3' end.
// vote=true must reproduce the genome exactly.
{
    auto r = make_rng(41);
    std::string g;
    for (int i = 0; i < 400; ++i) g.push_back("ACGT"[r() % 4]);
    PackedReads p;
    for (int st = 0; st + 50 <= 400; st += 8) {
        std::string rd = g.substr(st, 50);
        if (st == 80) rd[49] = rd[49] == 'A' ? 'C' : 'A';
        p.push(rd);
    }
    p.push(g.substr(350, 50));
    std::vector<std::string> c = assemble_ssake(p, 25, 2, true);
    assert(c.size() == 1 && c[0] == g);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `g++ -O2 -std=c++17 -static tests/test_ssake.cpp -o build/test_ssake.exe`
Expected: 컴파일 오류 (4인자 오버로드 없음)

- [ ] **Step 3: Implement** 다수결 연장 in `src/ssake.cpp`, `run_experiments.py` 단계 5에 `vote=1` 행 추가

- [ ] **Step 4: Run all tests and step 5**

Run: Task 1 Step 4 명령 후 `python -u scripts/run_experiments.py --steps 5 && python scripts/plot.py`
Expected: `FAIL` 없음. step5.csv에 vote 행 추가

- [ ] **Step 5: README 단계 5 표 갱신 후 Commit**

```bash
git add src/ssake.cpp tests/test_ssake.cpp scripts/run_experiments.py README.md
git add -f results/step5.csv results/step5_*.png
git commit -m "feat: majority-vote extension for reads with errors"
```
