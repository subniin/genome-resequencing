# genome-resequencing

짧은 read(short read)만으로 유전체를 복원하는 **De novo 조립**을 다룬 프로젝트입니다. 설계 프로젝트 조건(`알고리즘응용_설계.pdf` 60~64쪽)에 맞춰 실험했고, 발표 구성(Introduction, Benchmark, Algorithm, Result, Future work)에 따라 정리합니다.

## 1. Introduction

### 1.1 문제 정의

길이 L인 short read M개가 주어졌을 때, 그 read가 나온 길이 N의 원본 서열을 복원합니다 (PDF 60쪽). 이 저장소는 다음 선택을 했습니다.

- **de novo**: 조립기는 read만 입력으로 받는다. reference는 평가 단계에서만 읽는다.
- **single-end**, 정방향 read만 사용한다. paired-end와 역상보 가닥은 future work다.
- SNP가 있는 sample genome을 복원하고, 정답 SNP를 얼마나 찾았는지 평가한다.

### 1.2 조건과 실험 범위

| 항목 | PDF 조건 (61쪽) | 이 저장소 실험 |
|---|---|---|
| N (유전체 길이) | 최대 3×10^9 | 10^5 ~ 10^8 (실측 최대 N=10^8) |
| M (read 개수) | 약 2×10^7 ~ 2×10^8 | 2×10^4 ~ 2×10^7 (실측). 2×10^8은 외삽 |
| L (read 길이) | 약 32 ~ 100 | 32 / 50 / 75 / 100 |
| D (read당 mismatch 수) | 허용치로 정함 | D = 2 (기본), 0 ~ 3 비교 |
| repeat | 있을 수 있음 | 단계 6에서 20개 삽입 |
| read 형태 | single-end 또는 paired-end | single-end |
| 방식 | de novo 또는 reference | de novo |

기본 설정은 N = 10^5, L = 50, coverage 10, SNP 0.1%, o = ⌈L/2⌉ = 25, D = 2입니다. 한 번에 한 변수만 바꿉니다.

### 1.3 데이터 생성

- **reference**: 모든 16-mer가 서로 다른 ACGT 랜덤 서열을 만든다. 한 염기씩 붙이면서 새로 생긴 16-mer가 이미 나왔으면 다른 염기를 고르고, 막히면 한 칸 되돌아간다. 4^16 ≈ 4.3×10^9 공간에서 N=10^8까지 가능하다. 난수는 `mt19937_64`에 고정 seed를 쓴다.
- **repeat**: 단계 6에서만 쓴다. 길이 R 구간을 고른 뒤 |s−t| ≥ R인 다른 위치에 덮어써 20개 복사본을 만든다.
- **sample**: reference에 SNP를 무작위 위치에 치환으로 넣는다 (SNP 0.1%가 기본).
- **read**: sample에서 길이 L인 구간을 시작 위치 균등 분포로 뽑는다. M = round(N·coverage/L). 단계 5에서는 염기마다 확률 `err`로 다른 염기로 바꾼다.

### 1.4 입출력 형식

- 입력 read 파일: 한 줄에 read 하나. 빈 줄과 CR은 무시하고, 길이가 다르면 줄 번호와 함께 중단한다.
- 출력: contig FASTA (`>contig_0`, ...).
- 정답 SNP: `pos<TAB>ref<TAB>alt` TSV, 0-based.
- 평가 출력: accuracy, snp_recall, snp_precision, n_contigs, n50, covered_frac, misassemblies, unaligned_contigs.

## 2. Benchmark

비교 대상은 **Greedy Overlap**입니다. 모든 read 쌍의 최대 겹침을 한 번씩 계산하고, 겹침이 긴 순서대로 연결한다. 양쪽 끝이 비어 있고 사이클을 만들지 않을 때만 잇는다(union-find). 구현이 단순해서 정확성 기준선으로 쓰기 좋지만, 쌍 비교가 M²에 비례해 큰 M에서는 돌릴 수 없다. 그래서 실측은 M ≤ 3.1×10^4 범위에서만 했다.

## 3. Algorithm

### 3.1 SSAKE 계열 seed-and-extend

read 하나를 seed로 잡고, 그 3' 끝과 겹치는 read를 찾아 반복해서 붙입니다. 더 붙일 read가 없으면 contig를 닫고 남은 read로 새 seed를 시작합니다. 3' 연장 후에 5' 연장을 같은 방식으로 한다. 한쪽만 연장하면 seed 왼쪽 read가 별도 contig로 남아 조각난다.

### 3.2 이 구현의 자료구조

- **2-bit 압축 read** (`PackedReads`): 염기 하나를 2비트로 저장한다. read 하나는 2L 비트(L=50이면 13바이트 정도). 문자열(염기당 1바이트) 저장 대비 메모리가 약 1/4이다.
- **정렬 배열 색인**: read마다 key를 만든다. key는 앞 kk = min(o, 32)염기(prefix) 또는 뒤 kk염기(suffix)다. (key, read id) 쌍을 key 순으로 정렬해 두 배열(keys, ids)로 둔다. 조회는 `equal_range`로 한다.
- **연장 조건**: contig 끝 k염기(k는 o에서 L−1까지 큰 값부터)에 대응하는 key로 후보를 찾고, 겹침 전체 길이 k에 대해 mismatch ≤ D를 확인한다. kk가 32로 잘리므로 o > 32 (L=75, 100)에서도 key는 32염기만 쓰고, 검증은 겹침 전체로 한다.
- **중복 제거**: 같은 read는 첫 번째 사본만 남긴다. 정렬 후 같은 내용끼리 묶어 id가 가장 작은 것만 남긴다.
- **다수결 연장**: 오류가 있는 read 때문에 연장 후보가 충돌하는 경우를 다룬다. 단계 5 결과에 따라 Task 8에서 추가하며, 아래 3.4에 기록한다.

### 3.3 복잡도

M = read 수, L = read 길이, o = 최소 겹침.

| 도구 | 시간 | 공간 |
|---|---|---|
| SSAKE (이 구현) | O(M log M + M·(L−o)·(log M + L)) | O(M·L/4 + 24M) bytes |
| Greedy Overlap | O(M²·L²) (쌍마다 겹침 계산) | O(M²) (후보 쌍) |

SSAKE 공간의 24M은 prefix와 suffix 색인(key 8바이트 + id 4바이트, 두 색인)에서 온다. 실측 메모리는 이 공식보다 크다 (단계 1 결과 참고). 차이는 벡터 용량 여유와 색인 정렬 중 임시 버퍼 때문일 것으로 보이나(가설), 확인하지 않았다.

### 3.4 다수결 연장 (단계 5 결과에 따른 추가)

단계 5(4.7절)에서 SSAKE의 오류율 0.01, D=2 결과가 실행 조건(accuracy < 0.999 또는 precision < 0.9)을 만족했다. 다수결 연장 구현과 결과는 Task 8에서 이 절에 추가한다.

### 3.5 알려진 한계

- 탐욕 방식이라 repeat에서 연장이 멈추거나 잘못 이어질 수 있다 (단계 6에서 확인).
- read 오류가 있으면 contig에 오류 염기가 남고, 그 염기가 SNP로 잘못 호출된다 (단계 5에서 precision이 크게 떨어짐).
- 최소 겹침 o가 너무 크면 연결 가능한 쌍이 줄어 contig가 늘어난다 (단계 4 o=35).

## 4. Result

### 4.1 측정 환경

- CPU: AMD Ryzen 5 7500F, RAM: DDR5 16GB, OS: Windows 11
- 빌드: MinGW-w64 g++ 16.2.0, `g++ -O2 -std=c++17 -static` (동적 링크는 segfault가 나서 정적 링크 사용)
- 단일 스레드, 실험 스크립트는 Python 3 (psutil로 메모리 측정)
- 시간은 wall time. 메모리는 자식 프로세스 포함 RSS를 20ms 간격으로 샘플링한 최댓값이다.

원본: `results/step0.csv` ~ `results/step6.csv`. 그래프: `results/step1_time_memory.png`, `results/step2_contigs_vs_coverage.png`, `results/step3_recall_vs_snp.png`, `results/step4_L_contigs_n50.png`, `results/step5_err_accuracy.png`, `results/step6_repeat_misassembly.png`.

### 4.2 단계 0: 기본 설정

| 도구 | 시간 (s) | 메모리 (MB) | contig | N50 | recall | misassembly |
|---|---|---|---|---|---|---|
| SSAKE | 0.138 | 4.1 | 135 | 1266 | 1.000 | 0 |
| Greedy | 43.148 | 11.8 | 135 | 1266 | 1.000 | 0 |

### 4.3 단계 1: M에 따른 시간·메모리 (L 50, coverage 10, SNP 0.1%)

N = M·L/coverage로 두었다. M을 2×10^4 ~ 2×10^7로 올리면 N은 10^5 ~ 10^8이다.

| M | N | SSAKE 시간 (s) | SSAKE 메모리 (MB) | contig | N50 | recall | misassembly |
|---|---|---|---|---|---|---|---|
| 2×10^4 | 10^5 | 0.031 | 4.2 | 135 | 1266 | 1.000 | 0 |
| 2×10^5 | 10^6 | 0.291 | 12.8 | 1213 | 1350 | 1.000 | 0 |
| 2×10^6 | 10^7 | 3.91 | 88.3 | 12189 | 1349 | 1.000 | 0 |
| **2×10^7** | **10^8** | **72.47** | **1001.4** | 122100 | 1335 | 0.99993 | 2 |

- 실측 최대(M=2×10^7, N=10^8)에서 SSAKE 피크 메모리는 1001 MB로 목표 2.0GB 이하를 만족했다.
- read 1개당 메모리는 M=2×10^6에서 약 46 B, M=2×10^7에서 약 52 B다.
- Greedy Overlap은 M=2×10^4에서만 돌렸다: 42.62초, 11.8 MB. 같은 입력에서 SSAKE는 0.031초라 약 1,400배 차이다. Greedy는 쌍 비교가 M²에 비례하므로 M=2×10^7은 실행하지 않았다.

**M=2×10^8 외삽 (실측 아님).** 선형 가정으로 메모리는 52 B × 2×10^8 ≈ 10.5 GB다. 16GB 램에서 이론상 들어가지만 OS와 read 파일 버퍼를 빼면 여유가 작다. 시간은 72초에서 선형으로 늘린다고 가정하면 약 12~15분이다. 두 값 모두 외삽이다.

### 4.4 단계 2: coverage에 따른 contig 수 (N 10^5, L 50)

| coverage | M | Lander-Waterman 예상 | SSAKE contig | Greedy contig | SNP recall |
|---|---|---|---|---|---|
| 5 | 10^4 | 821 | 770 | 770 | 1.000 (둘 다) |
| 10 | 2×10^4 | 135 | 135 | 135 | 1.000 (둘 다) |
| 15 | 3×10^4 | 16.6 | 15 | 15 | 1.000 (둘 다) |

예상값은 `N·(c/L)·exp(−(c/L)(L−o))`로 계산한 것이다. 두 조립기 결과가 같다. Greedy는 coverage 15에서 87.6초였다.

### 4.5 단계 3: SNP 비율 (N 10^5, coverage 10)

| SNP 비율 | SNP 수 | SSAKE contig | SSAKE recall | Greedy recall | misassembly |
|---|---|---|---|---|---|
| 0.1% | 100 | 135 | 1.000 | 1.000 | 0 |
| 0.5% | 500 | 133 | 1.000 | 1.000 | 0 |
| 1.0% | 1000 | 132 | 1.000 | 1.000 | 0 |

SNP 비율이 1%까지 올라도 recall은 1.000으로 같다. 이 설계에서 read는 sample에서 뽑으므로 SNP가 read에 그대로 들어 있다. 그래서 SNP는 조립을 방해하지 않고 평가에서만 기록된다.

### 4.6 단계 4: read 길이 L과 최소 겹침 o

| L | o | M | SSAKE contig | SSAKE N50 | SSAKE 시간 (s) | Greedy 시간 (s) | recall (둘 다) |
|---|---|---|---|---|---|---|---|
| 32 | 16 | 31250 | 194 | 806 | 0.058 | 57.40 | 1.000 |
| 50 | 25 | 20000 | 135 | 1266 | 0.030 | 42.57 | 1.000 |
| 75 | 38 | 13333 | 85 | 2013 | 0.032 | 30.44 | 1.000 |
| 100 | 50 | 10000 | 76 | 1976 | 0.032 | 24.89 | 1.000 |

L=100(o=50)은 o가 32를 넘는 경로다. key는 32염기만 쓰고 겹침 전체를 검증했는데도 recall 1.000, misassembly 0으로 정상 조립되었다.

L=50에서 o만 바꾼 결과:

| o | SSAKE contig | SSAKE N50 | Greedy contig |
|---|---|---|---|
| 20 | 42 | 3732 | 42 |
| 25 | 135 | 1266 | 135 |
| 35 | 922 | 164 | 922 |

o를 크게 잡을수록 contig 수가 늘었다. o=35에서 contig 922개는 o가 L에 가까워져 연결 가능한 쌍이 줄기 때문으로 보인다(가설). 확인하지 않았다. 모든 설정에서 misassembly는 0이었다.

### 4.7 단계 5: read 오류와 mismatch 허용치 D (N 10^5, L 50, coverage 10)

오류 실험은 D = 2로 고정하고 오류율을 바꾼 것과, 오류율 0.01에서 D를 바꾼 것 두 가지다.

| 오류율 | SSAKE accuracy | SSAKE precision | SSAKE recall | SSAKE contig | Greedy accuracy | Greedy precision | Greedy recall | Greedy contig |
|---|---|---|---|---|---|---|---|---|
| 0 | 1.000000 | 1.000 | 1.000 | 135 | 1.000000 | 1.000 | 1.000 | 135 |
| 0.005 | 0.956229 | 0.032 | 0.990 | 2859 | 0.981652 | 0.078 | 0.990 | 870 |
| 0.01 | 0.917205 | 0.016 | 1.000 | 4844 | 0.945312 | 0.032 | 1.000 | 1625 |
| 0.02 | 0.813749 | 0.009 | 0.990 | 8502 | 0.838960 | 0.013 | 0.970 | 3656 |

오류율 0.01에서 D를 바꾼 결과:

| D | SSAKE accuracy | SSAKE recall | SSAKE contig | Greedy accuracy | Greedy recall | Greedy contig |
|---|---|---|---|---|---|---|
| 0 | 0.929540 | 0.990 | 7443 | 0.929246 | 1.000 | 7487 |
| 1 | 0.920993 | 1.000 | 5474 | 0.915677 | 1.000 | 3370 |
| 2 | 0.917205 | 1.000 | 4844 | 0.945312 | 1.000 | 1625 |
| 3 | 0.916982 | 1.000 | 4753 | 0.970636 | 1.000 | 1135 |

판정: 오류율 0.01, D=2에서 SSAKE accuracy는 0.917로 0.999보다 낮고, SNP precision은 0.016으로 0.9보다 낮다. 이 조건은 다수결 연장(Task 8) 실행 조건을 만족한다. precision이 낮은 이유는 contig에 남은 시퀀싱 오류 염기가 reference와 다른 염기로 집계되어 정답 없는 SNP로 호출되기 때문으로 보인다(가설). 오류 염기 수를 직접 세지 않았다.

### 4.8 단계 6: repeat (N 10^5, L 50, 반복 20개)

| repeat 길이 | SSAKE misassembly | SSAKE recall | SSAKE precision | SSAKE accuracy | Greedy misassembly | Greedy recall | Greedy precision |
|---|---|---|---|---|---|---|---|
| 25 | 0 | 1.000 | 1.000 | 1.000000 | 0 | 1.000 | 1.000 |
| 50 | 7 | 0.980 | 0.560 | 0.998954 | 21 | 1.000 | 0.952 |
| 100 | 11 | 0.980 | 0.551 | 0.992458 | 15 | 0.970 | 0.858 |
| 200 | 10 | 0.990 | 0.723 | 0.983675 | 13 | 1.000 | 0.980 |

repeat 길이가 read 길이(50)보다 짧은 25에서는 오조립이 없었다. 길이가 50 이상이 되면 misassembly가 생겼다. 200에서 SSAKE의 covered_frac은 0.963, unaligned는 1이었다. 반복이 있으면 precision이 크게 떨어진다. 즉 repeat에서 잘못 이어진 contig가 SNP 호출을 오염시킨다.

### 4.9 평가 수정 이력

초기 평가 도구는 문자열 k-mer 맵을 썼고, anchor k=14를 기본으로 했다. N=10^7(L=30, SNP 0.5%, o=15)에서 두 값을 비교한 결과다.

| k | accuracy | SNP recall | SNP precision | misassembly | unaligned |
|---|---|---|---|---|---|
| 14 | 0.999761 | 0.904060 | 0.985996 | 8786 | 0 |
| 24 | 0.999846 | 0.999340 | 0.996252 | 11 | 40 |

k=14에서 보이던 높은 오조립 수와 낮은 recall은 평가 방식의 문제로 보인다. k=24로 올리면 misassembly가 8786에서 11로 줄고 recall이 0.904에서 0.999로 오른다. 14-mer가 reference에서 우연히 여러 번 나오면서 가짜 anchor가 생긴다는 설명은 가설이며, 가짜 anchor를 직접 세지는 않았다. 평가를 k=24로 고쳐 단계 0~6을 다시 실행했다. 이전 README의 수치(N=10^4)는 이 실험 설계와 다르므로 사용하지 않는다.

### 4.10 구현에서 정한 것 (원 논문과 다른 점)

- 중복 read는 두 조립기 모두 입력 단계에서 첫 번째 사본만 남긴다.
- SSAKE 연장 후보가 여러 개면 입력 순서상 첫 번째를 쓴다. 다수결은 Task 8에서 추가한다.
- accuracy는 contig 전체 염기 대비 sample과 일치하는 염기 비율이다. 정렬되지 않은 contig는 전부 불일치로 센다.
- 평가는 reference 고유 k-mer(k=24) anchor로 대각선 구간을 잡고 꼬리까지 채점한다.

## 5. Future work

### 장점과 단점

- **장점**: 메모리가 read당 약 50바이트로 작고, M=2×10^7을 단일 PC에서 1GB 이하로 돌렸다. 정확한 read에서는 recall 1.000을 보였다.
- **단점**: 오류 read에 약하다 (단계 5: 오류율 0.01에서 precision 0.016). repeat 길이가 read 길이 이상이면 오조립이 생긴다 (단계 6). 역상보와 paired-end를 쓰지 않아 정보가 부족하다.

### 개선 방안

- paired-end와 역상보 가닥을 입력에 추가해 repeat 연결을 돕는다.
- de Bruijn graph 기반 조립과 비교한다. 오류 k-mer를 빈도로 걸러내는 방식이 오류 read에 효과적일 것으로 기대하나(가설) 아직 구현하지 않았다.
- BWT 기반 색인으로 read 저장과 조회 메모리를 더 줄인다.
- 색인을 디스크에 두어 M=2×10^8 실측을 가능하게 한다.

## 6. 진행 상황

- [x] 문제 정의와 데이터 조건 (PDF 60~64쪽)
- [x] reference 생성기 (k=16, repeat 삽입)
- [x] read 생성기 (희소 SNP, 시퀀싱 오류)
- [x] 2-bit 압축 read와 정렬 배열 색인
- [x] 평가 도구 (정수 k-mer anchor, k=24)
- [x] 벤치마크 Greedy Overlap
- [x] 실험 단계 0~6 실행과 결과 정리
- [ ] 다수결 연장 (단계 5 조건 충족, Task 8)

## 7. 참고

- Warren et al., SSAKE, Bioinformatics, 2007 (원문 확인 필요)
- Lander & Waterman, 1988 (끊김 수 근사)
- Eulerian path / de Bruijn 그래프 (비교 대상 알고리즘의 배경)
