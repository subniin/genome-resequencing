# De novo Re-sequencing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 실습 자료(`1029실습_프로젝트준비.pdf`) 조건의 랜덤 ACGT 시퀀스에서 short read만으로 SNP가 포함된 sample genome을 contig로 복원하는 SSAKE 계열 조립기를 만들고, Greedy Overlap 벤치마크와 비교 실험한다.

**Architecture:** C++17 CLI 도구 4개(`gen_reads`, `greedy_overlap`, `ssake`, `evaluate`)가 텍스트 파일(reads, FASTA, TSV)로 연결된다. Python은 실험 반복 실행과 그래프만 맡는다.

**Tech Stack:** C++17 (g++ MinGW-w64 UCRT 16.2.0, `-O2 -std=c++17 -static`), Python 3.11 + matplotlib, 난수 `mt19937_64` 고정 seed.

**Spec:** `README.md` (섹션 2~7)

## Global Constraints

- 데이터 조건: `README.md` 3장에 정의됨 (A/T/G/C 랜덤 시퀀스, 길이 10K/100K/1M/10M, read 길이 15/30/45 중 하나, coverage 5 이상, 반복 줄이는 방법을 보고서에 기록)
- 기본 설정: N=10,000, L=30, coverage 10, M=3,333, SNP 50개(0.5%, 치환만), D=2, 최소 겹침 o=⌈L/2⌉=15
- reference: 외부 데이터 없이 ACGT 랜덤 생성. 길이 14인 부분 문자열(14-mer)이 모두 고유하도록 생성 (반복 제거)
- read는 정방향만, 시작 위치 균등 무작위. paired-end·역상보 없음. 기본 설정에서는 시퀀싱 오류 없음
- 난수는 `mt19937_64` + 고정 seed. `rand()` 금지
- 입출력: read는 한 줄에 하나(`reads.txt`), contig는 FASTA, 정답 SNP는 TSV(`pos<TAB>ref<TAB>alt`, 0-based)
- reference는 평가 단계에서만 사용. `ssake`, `greedy_overlap`은 reads만 입력으로 받는다
- 단일 스레드, 같은 PC에서 시간·메모리 측정 (Ryzen 5 7500F / DDR5 16GB)
- 알고리즘 설명보다 실측 결과 기록을 우선

## Review Focus

1. 직접 만든 reference 파일에 소문자·CRLF·ACGT 외 문자가 있음: 대문자로 정규화하고 ACGT 외 문자는 오류로 중단 (Task 1)
2. 중복 read(같은 위치에서 두 번 뽑힘)와 완전히 동일한 read: 조립이 멈추거나 무한 루프에 빠지지 않는다 (Task 3)
3. 최소 겹침 o가 L 이상이거나 read 수가 0/1: 오류 메시지를 내거나 read 자체를 contig로 낸다 (Task 2, 3)
4. contig가 reference에서 anchor를 못 찾음(너무 짧거나 SNP 밀집): 정렬 실패로 세고 크래시하지 않는다 (Task 4)
5. SNP 위치가 서로 겹치거나 시퀀스 끝(첫/마지막 염기)에 있음: 중복 없이 요청한 개수를 뽑는다 (Task 1)
6. N이 너무 커서 14-mer를 고유하게 만들 수 없음(N > 4^14): 오류로 중단 (Task 0)

## File Structure

| 파일 | 책임 |
|---|---|
| `src/common.hpp` | 텍스트/FASTA 읽기·쓰기, 시드 있는 RNG, mismatch 세기 (header-only) |
| `src/gen_reads.cpp` | reference 구간 선택 → SNP 삽입 → read 샘플링 |
| `src/greedy_overlap.cpp` | 벤치마크: all-pairs 겹침 탐색 |
| `src/ssake.cpp` | 본 알고리즘: prefix/suffix 해시 + 양방향 연장 |
| `src/evaluate.cpp` | contig→reference 정렬, SNP 호출, 지표 계산 |
| `src/gen_reference.cpp` | 14-mer가 모두 고유한 ACGT 랜덤 시퀀스 생성 |
| `scripts/run_experiments.py` | 단계 0~4 반복 실행, 시간·메모리 측정, `results/*.csv` 저장 |
| `scripts/plot.py` | `results/*.csv` → 그래프 |
| `tests/test_*.cpp` | 도구별 `assert` 기반 자체 검사 (프레임워크 없음) |

`evaluate`를 C++로 둔 이유: N=10M 단계에서 Python 정렬은 느리다. 문제가 되지 않으면 Python으로 바꿔도 된다.

---

### Task 0: 환경 스캐폴딩과 reference 생성

**Files:**
- Create: `src/common.hpp`, `src/gen_reference.cpp`, `.gitignore`, `tests/test_common.cpp`, `tests/test_gen_reference.cpp`

**Interfaces:**
- Produces (`common.hpp`):
  - `std::string read_fasta_seq(const std::string& path)`: 헤더 제외, 대문자, 개행·CR 제거
  - `std::vector<std::string> read_lines(const std::string& path)`: 빈 줄 제외
  - `void write_fasta(const std::string& path, const std::vector<std::string>& contigs)`: 이름 `contig_<i>`
  - `int mismatches_bounded(const std::string& a, size_t ai, const std::string& b, size_t bi, size_t len, int limit)`: mismatch 수, limit 초과 시 `limit+1`에서 중단
  - `std::mt19937_64 make_rng(uint64_t seed)`
- Produces CLI: `gen_reference --N 10000 --k 14 --seed 1 --out data/ref_10000.txt`
  - `std::string gen_unique_kmer_seq(size_t N, int k, std::mt19937_64& rng)`: 염기를 하나씩 붙이되 새 k-mer가 이미 나왔으면 다른 염기를 뽑고, 4개 모두 막히면 직전 염기부터 되돌려 다시 뽑는다. 사용한 k-mer는 4^k 비트 bitset으로 관리

- [ ] **Step 1: 테스트 작성** `tests/test_common.cpp`: `mismatches_bounded("ACGT",0,"ACCT",0,4,2)==1`, `mismatches_bounded("AAAA",0,"TTTT",0,4,2)==3`, 임시 FASTA(소문자+CRLF)를 `read_fasta_seq`로 읽으면 대문자·개행 없는 문자열
- [ ] **Step 2: 실패 확인** `g++ -O2 -std=c++17 tests/test_common.cpp -o build/test_common.exe` → 컴파일 실패 (`common.hpp` 없음)
- [ ] **Step 3: `common.hpp` 구현**
- [ ] **Step 4: 통과 확인** `build\test_common.exe` → 종료 코드 0
- [ ] **Step 5: `gen_reference` 테스트 작성** (`tests/test_gen_reference.cpp`): N=10000, k=14 → 길이 10000, ACGT만, 모든 14-mer 고유, 같은 seed 두 번 결과 동일. `gen_unique_kmer_seq(5, 1, rng)`처럼 N > 4^k면 예외 (Review Focus 6)
- [ ] **Step 6: `gen_reference.cpp` 구현.** `.gitignore`에 `build/`, `data/`, `results/`. 실행해 N = 10K, 100K, 1M, 10M 파일 4개 생성 확인
- [ ] **Step 7: 반복 비교 기록.** 같은 N=10K에서 `mt19937_64`로 그냥 뽑은 시퀀스의 중복 14-mer 수와 고유 생성 시퀀스(0개)를 출력하는 한 줄 스크립트 결과를 README 3장 "반복을 줄이는 방법" 아래에 적는다 (실습 자료의 보고서 요구)
- [ ] **Step 8: Commit** `chore: add common helpers and unique-kmer reference generator`

### Task 1: read 생성기

**Files:**
- Create: `src/gen_reads.cpp`, `tests/test_gen_reads.cpp`

**Interfaces:**
- Consumes: `read_fasta_seq`, `make_rng`
- Produces CLI: `gen_reads --ref data/ref_10000.txt --L 30 --coverage 10 --snps 50 --seed 1 --out data/run0/` (M = round(N·coverage/L), N은 reference 길이)
  - `sample.txt` (SNP 삽입 후), `truth_snps.tsv`, `reads.txt`. reference는 입력 파일이 곧 정답 원본
  - 내부 함수 `std::vector<Snp> insert_snps(std::string& seq, int count, std::mt19937_64& rng)` (`struct Snp {size_t pos; char ref, alt;}`)

- [ ] **Step 1: 테스트 작성** (`insert_snps`와 구간 선택 함수 단위): 
  - `insert_snps(seq, 50, rng)` 후 `sample`과 원본의 서로 다른 위치 수 == 50 (기본 설정 SNP 개수), 모든 `alt != ref`, 위치 중복 없음 (Review Focus 5)
  - 소문자·CRLF reference는 정규화되어 읽히고, ACGT 외 문자는 예외 (Review Focus 1)
  - 같은 seed 두 번 실행 결과 동일
- [ ] **Step 2: 실패 확인**
- [ ] **Step 3: 구현.** read 시작 위치는 `[0, N-L]` 균등.
- [ ] **Step 4: 통과 확인** + 기본 설정 실행, `reads.txt` 줄 수 3333, 각 줄 길이 30
- [ ] **Step 5: Commit** `feat: add read generator with SNP insertion`

### Task 2: 벤치마크 Greedy Overlap

**Files:**
- Create: `src/greedy_overlap.cpp`, `tests/test_greedy_overlap.cpp`

**Interfaces:**
- Consumes: `read_lines`, `write_fasta`, `mismatches_bounded`
- Produces CLI: `greedy_overlap --reads reads.txt --o 15 --D 2 --out contigs.fa`
  - `std::vector<std::string> assemble_greedy(const std::vector<std::string>& reads, int o, int D)`

알고리즘: 모든 read 쌍에 대해 겹침 길이 k ∈ [o, L-1]의 최대 겹침을 구해 가장 긴 겹침부터 병합(탐욕). mismatch ≤ D.

- [ ] **Step 1: 테스트 작성**: 길이 300 문자열에서 L=100, 간격 50으로 만든 read 5개 → `assemble_greedy` 결과 contig 1개, 원본과 동일. read 1개 → contig 1개(그 read). `o >= L` → 예외 (Review Focus 3)
- [ ] **Step 2: 실패 확인**
- [ ] **Step 3: 구현** O(M²·L). N=100K 기본 설정만 돌린다 (5M은 대상 아님)
- [ ] **Step 4: 통과 확인**
- [ ] **Step 5: Commit** `feat: add greedy overlap baseline`

### Task 3: SSAKE 계열 조립기

**Files:**
- Create: `src/ssake.cpp`, `tests/test_ssake.cpp`

**Interfaces:**
- Consumes: `read_lines`, `write_fasta`, `mismatches_bounded`
- Produces CLI: `ssake --reads reads.txt --o 15 --D 2 --out contigs.fa`
  - `std::vector<std::string> assemble_ssake(const std::vector<std::string>& reads, int o, int D)`

알고리즘 결정 사항:
- 색인: read의 앞 o개 염기(prefix) → read id 목록, 뒤 o개 염기(suffix) → read id 목록
- seed 순서: 입력 순서 고정 (재현성)
- 3' 연장: contig 끝에서 겹침 길이 k를 L-1부터 o까지 줄이며, `contig.substr(len-k, o)`로 prefix 색인을 조회하고 후보를 `mismatches_bounded(..., k, D)`로 검증. 가장 긴 k의 후보를 붙이고 사용 처리. 없으면 중단
- 5' 연장: suffix 색인으로 대칭 구현
- 사용된 read는 색인 조회에서 제외. 후보 여러 개면 입력 순서상 첫 번째 (다수결은 단계 5, 지금은 구현 안 함)

- [ ] **Step 1: 테스트 작성**: Task 2와 같은 5-read 입력 → contig 1개. read를 섞어도 같은 contig 집합. 완전히 동일한 read 3개 포함 입력 → 무한 루프 없이 종료하고 contig 하나 (Review Focus 2). 겹침 없는 두 read → contig 2개
- [ ] **Step 2: 실패 확인**
- [ ] **Step 3: 구현**
- [ ] **Step 4: 통과 확인**, 기본 설정 데이터에서 `greedy_overlap`과 contig 수 비교
- [ ] **Step 5: Commit** `feat: add SSAKE-style assembler`

### Task 4: 평가 도구

**Files:**
- Create: `src/evaluate.cpp`, `tests/test_evaluate.cpp`

**Interfaces:**
- Consumes: `read_fasta_seq`, `read_lines`, `mismatches_bounded`
- Produces CLI: `evaluate --contigs contigs.fa --ref ref_10000.txt --sample sample.txt --truth truth_snps.tsv --k 14`
  - 출력: `key=value` 줄. 키 `accuracy`, `snp_recall`, `snp_precision`, `n_contigs`, `n50`, `covered_frac`, `misassemblies`, `unaligned_contigs`

방법: reference에서 고유한 k-mer(k=12~20)를 anchor로 contig 위치를 잡고 염기 단위로 비교해 SNP를 호출. 호출한 SNP를 truth와 대조해 recall/precision. 오조립은 한 contig 안에서 anchor 위치가 비단조이거나 점프하는 경우로 센다.

- [ ] **Step 1: 테스트 작성**: 합성 reference 1000염기 + SNP 3개 sample, contig 1개 = sample 그대로 → `snp_recall=1`, `snp_precision=1`, `accuracy=1`. contig를 앞뒤 뒤바꿔 이은 것 → `misassemblies>=1`. reference와 무관한 무작위 contig → `unaligned_contigs=1`, 크래시 없음 (Review Focus 4)
- [ ] **Step 2: 실패 확인**
- [ ] **Step 3: 구현**
- [ ] **Step 4: 통과 확인**, Task 3 결과 contigs.fa에 실행해 지표 확인
- [ ] **Step 5: Commit** `feat: add evaluator for contig accuracy and SNP recovery`

### Task 5: 실험 스크립트와 결과 정리

**Files:**
- Create: `scripts/run_experiments.py`, `scripts/plot.py`
- Modify: `README.md` (섹션 8 진행 상황, 결과 표·그래프, 측정 PC 사양)

**Interfaces:**
- Consumes: Task 1~4의 CLI 4개
- Produces: `results/step{1..4}.csv` (열: 변수값, 도구, seconds, peak_mem_mb, 평가 지표 전부)

실험 단계 (한 번에 한 변수): 

| 단계 | 변경 | 범위 |
|---|---|---|
| 0 | 기본 설정 | `ssake`, `greedy_overlap` 모두 |
| 1 | N | 10K, 100K, 1M, 10M (`greedy_overlap`은 100K까지만) |
| 2 | coverage | 5, 10, 15 → README의 Lander–Waterman 예상 contig 수(≈137, 22, 3)와 비교 |
| 3 | SNP 비율 | 0.1%, 0.5%, 1.0% |
| 4 | L, o | L = 15 / 30 / 45 (o = ⌈L/2⌉), L=30에서 o 변경 |

시간·메모리는 `psutil`로 자식 프로세스의 wall time과 peak working set을 측정한다. 5~6단계(오류·다수결, 역상보·repeat)는 선택이며 이 계획에서는 다루지 않는다.

- [ ] **Step 1: `run_experiments.py` 작성** 각 실행의 seed 고정, 도구 컴파일은 사전에 `g++ -O2 -std=c++17 -static` 으로 `build/`에 생성
- [ ] **Step 2: 단계 0 실행**, CSV 행 2개 확인
- [ ] **Step 3: 단계 1~4 실행**
- [ ] **Step 4: `plot.py`** 시간/메모리 vs N(로그축), contig 수 vs coverage(예상 곡선 겹침), SNP recall vs SNP 비율
- [ ] **Step 5: README 갱신**: 진행 상황 체크박스, 결과 표, 그래프 경로, 원논문과 달라진 구현 지점, 미정 사항 4개 확정값
- [ ] **Step 6: Commit** `docs: add experiment results`

## 미정 사항 (구현 전에 기본값으로 둔 것)

| 항목 | 이 계획의 기본값 | 바꾸려면 |
|---|---|---|
| D의 의미 | 겹침 검증에서 허용하는 mismatch 수 | Task 2, 3의 검증 호출만 수정 |
| 역상보 가닥 | 기본 설정에서 제외 (선택 단계 6) | Task 1에 `--rc` 옵션과 Task 3에 역상보 색인 추가 |
| 입출력 형식 | 위 Global Constraints | 모든 CLI 영향 |
