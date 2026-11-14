# genome-resequencing

short read로 유전체를 복원하는 **De novo 방식**을 다루는 프로젝트입니다. 문제 정의, 데이터 설정, 선택한 알고리즘, 평가 계획, 진행 상황을 기록합니다.


## 1. 프로젝트 개요

길이 L의 short read M개가 주어졌을 때 길이 N의 유전체를 복원하는 re-sequencing 문제를 직접 정의하고, 알고리즘을 설계·구현해 실험합니다. 목표는 SNP가 포함된 sample genome을 복원하는 것입니다.

- 문제 조건(N, M, L, D, repeat, single/paired, reference 유무)은 직접 정해 쓴다.
- 정리 항목: 문제 정의, 데이터, 입출력, 벤치마크, 알고리즘, 결과(실행 환경 사양, 복잡도, 벤치마크 비교), future work.
- 알고리즘 설명보다 실제 실험 결과를 남기는 것을 우선한다.

## 2. 문제 정의

두 가지 방식이 있습니다.

| 방식 | 입력 | 출력 |
|---|---|---|
| Reference 방식 | reference + read | read의 위치와 SNP 목록 |
| **De novo 방식 (이 저장소)** | read만 | 겹침으로 이어 붙인 contig 집합 |

이 저장소는 De novo 방식을 다룹니다.

- 입력: sample genome에서 뽑은 short read 집합. reference는 주지 않는다.
- 출력: sample genome을 복원한 contig 집합(FASTA).
- 평가 단계에서만 reference를 사용한다. contig를 reference에 정렬해 SNP를 찾고, 정답과 비교한다.

## 3. 데이터 설정

데이터 조건은 다음과 같이 정의합니다. A, T, G, C 4개 염기로 시퀀스를 랜덤 생성합니다(길이 10,000 / 100,000 / 1,000,000 / 10,000,000). 일반 난수로 만들면 생기는 반복은 줄이는 방법을 쓰고 보고서에 기록합니다. read 길이는 15 / 30 / 45 중 하나로 하고, coverage는 5 이상으로 생성합니다.

ACGT 랜덤 시퀀스를 reference로 생성합니다. reference에 SNP를 넣어 sample genome을 만들고, sample에서 길이 L인 read를 무작위 위치로 뽑습니다.

기본 설정:

| 항목 | 값 |
|---|---|
| reference 길이 N | 10,000 (실습 자료의 시퀀스 길이 중 가장 작은 것) |
| read 길이 L | 30 (자료의 15 / 30 / 45 중 가운데) |
| coverage | 10 (= M×L/N, 자료의 요구는 5 이상) |
| read 개수 M | 3,333 |
| SNP | 50개 (0.5%), 치환만 |
| 허용 mismatch D | 2 |
| 최소 겹침 o | 15 (= ⌈L/2⌉) |
| repeat | 없음 (생성 시 제거) |

반복을 줄이는 방법: 시퀀스를 한 염기씩 늘리면서, 새로 생기는 길이 14짜리 부분 문자열(14-mer)이 이미 나온 적 있으면 그 염기를 다시 뽑는다. 결과적으로 시퀀스 안의 모든 14-mer가 서로 다르다. 4^14 ≈ 2.7억이라 N=10,000,000까지 가능하다. 일반 난수만 쓴 시퀀스와 반복 14-mer 수를 비교해 보고서에 적는다.

| N | 일반 난수(`mt19937_64`) 중복 14-mer 수 | 고유 생성 시퀀스 중복 14-mer 수 |
|---|---|---|
| 10,000 | 0 | 0 |
| 1,000,000 | 1,780 | 0 |
| 10,000,000 | 184,685 | 0 |

N=10K에서는 일반 난수도 중복이 거의 없어 차이가 보이지 않는다. 기대 충돌 수가 약 N²/(2·4^14)이라 N이 커야 차이가 드러난다.

가정:

1. read는 sample의 정방향에서만, 시작 위치를 균등하게 뽑는다. 역상보 가닥과 paired-end는 쓰지 않는다.
2. read에는 SNP만 있고 시퀀싱 오류는 기본 설정에서 넣지 않는다. 오류 추가는 선택 확장이다.
3. reference 안에는 repeat가 없다. repeat 포함 시퀀스는 선택 실험(단계 6)에서만 쓴다.

## 4. 선택한 알고리즘: SSAKE 계열 Seed-and-Extend

### 아이디어

read 하나를 seed로 잡고, 그 3' 끝과 겹치는 read를 찾아 반복해서 붙입니다. 겹치는 read가 없으면 contig를 끝냅니다. 처리하지 않은 read가 남아 있으면 새 seed로 다시 시작합니다.

### 이 문제에 맞춘 설계

- 색인: 모든 read의 앞 o개 염기(prefix)를 해시에 넣는다. o는 최소 겹침 길이다. 5' 쪽 연장을 위해 뒤 o개 염기(suffix) 색인도 둔다.
- 연장: contig 끝 o개와 같은 prefix를 가진 read를 찾고, 전체 겹침을 검사(mismatch ≤ D)한 뒤 가장 긴 겹침의 read를 붙인다. 붙인 read는 사용 처리한다.
- 양방향 연장: 정방향 read만 쓰므로 3' 연장 후 5' 연장을 대칭으로 구현한다. 한쪽만 하면 seed 왼쪽 read들이 별도 contig로 남아 contig가 겹쳐서 조각난다.
- 후보가 여러 개일 때(오류가 있는 확장 설정): 다음 염기에 대한 다수결로 연장 여부를 정한다. 기본 설정에서는 사용하지 않는다.

### 알려진 한계

- 탐욕 방식이라 repeat에서 연장이 멈추거나 오조립이 생길 수 있다.
- seed 순서에 따라 contig 경계가 달라질 수 있다. 고정 길이 read에서는 별도 정렬 기준이 없으므로 순서를 고정해 재현성을 확보한다.
- 최소 겹침 o가 너무 작으면 우연한 일치로 잘못 병합된다. 우연한 prefix 일치 기대 횟수는 약 M/4^o 이므로 N이 클수록 o를 키워야 한다.
- coverage나 L이 작으면 끊김이 늘어난다.

원 논문의 세부 규칙(자료구조, 동률 처리, 파라미터)은 아직 원문으로 확인하지 않았습니다. 구현이 원 알고리즘과 다른 부분은 결과를 정리할 때 함께 적습니다.

### 계산 기준값 (예상, 실측 전)

Lander–Waterman 근사로 read 시작 밀도 λ = M/N = cov/L, 연결 조건 길이 L−o일 때 끊김 수는 약 N·λ·e^(−λ(L−o)) 입니다. N=10,000, L=30, o=15 기준:

| coverage | 예상 contig 수 |
|---|---|
| 5 | 약 137 |
| 10 | 약 22 |
| 15 | 약 3 |

이 값을 구현 결과와 비교하는 기준선으로 씁니다.

## 5. 평가 지표

| 지표 | 정의 |
|---|---|
| Reconstruction Accuracy | 복원한 서열 중 sample과 일치하는 비율 |
| SNP Recovery Rate | contig를 reference에 정렬해 얻은 SNP 중 정답 SNP를 찾은 비율(recall), 오검출은 precision으로 따로 본다 |
| 실행 시간, 메모리 | 같은 PC에서 측정, 사양 기재 |
| contig 수, N50, 복원율 | 조립의 연속성 |
| 오조립 수 | 서로 다른 구간을 잘못 잇는 경우 |

벤치마크는 순진한 all-pairs 겹침 탐색(Greedy Overlap)으로 두고, 속도와 정확도를 비교합니다.

contig→reference 정렬은 reference에서 고유한 짧은 k-mer(12~20)를 anchor로 위치를 잡은 뒤 염기 단위로 비교하는 간단한 방식을 씁니다. SNP 비율이 1%이면 100염기마다 SNP가 하나라 anchor를 짧게 잡아야 합니다.

## 6. 실험 계획

기본 설정으로 정상 동작을 확인한 뒤 한 번에 한 변수씩 바꿉니다.

| 단계 | 변경 | 관찰 |
|---|---|---|
| 0 | 기본 설정 | 알고리즘 동작 확인 |
| 1 | N = 10K, 100K, 1M, 10M | 시간·메모리 증가 곡선 (Greedy Overlap은 100K까지만) |
| 2 | coverage = 5, 10, 15 | 끊김 수, 복원율 |
| 3 | SNP = 0.1%, 0.5%, 1.0% | SNP 복원율 (de novo는 SNP가 조립을 방해하지 않으므로 거의 평평할 것으로 예상) |
| 4 | L = 15 / 30 / 45 (o = ⌈L/2⌉), L=30에서 o 변경 | 파라미터 민감도 |
| 5 (선택) | read 오류 최대 D개, 다수결 추가 | 오류 내성, 다수결 비용 |
| 6 (선택) | 역상보 가닥, repeat 포함 구간 | 알고리즘 한계 |

reference는 14-mer가 모두 고유하도록 만들기 때문에, 길이가 커져도 repeat 때문에 끊기지 않습니다. repeat에 의한 한계는 단계 6에서 repeat를 일부러 넣어 따로 관찰합니다.

## 7. 개발 환경과 구현 계획

### 환경

| 항목 | 내용 |
|---|---|
| 핵심 구현 언어 | C++17 (g++ / MinGW-w64) |
| 보조 스크립트 | Python 3 (데이터 생성 보조, 결과 그래프는 matplotlib) |
| 빌드 | `g++ -O2 -std=c++17 -static` (Windows MinGW 실행 파일은 `-static` 필요) |
| OS | Windows 11 |
| 실행 PC | 데스크탑, GPU RTX 4060 Ti (GPU는 사용하지 않음, CPU 단일 스레드) |
| CPU / RAM | AMD Ryzen 5 7500F (6코어 12스레드, 최대 3.7GHz) / DDR5 16GB (8GB×2, 5200MT/s 동작) |
| 난수 | `mt19937_64`, 고정 seed (재현성). Windows `rand()`는 주기가 짧아 반복이 생기므로 쓰지 않는다 |

시간과 메모리는 같은 PC에서 측정하고, 결과를 적을 때 CPU/RAM 사양을 함께 기록합니다. 언어와 빌드 옵션은 구현을 시작하면서 바뀔 수 있습니다.

### 구현 계획

0. reference 생성기: 14-mer가 모두 고유한 ACGT 랜덤 시퀀스 (N = 10K, 100K, 1M, 10M).
1. read 생성기: reference에 SNP 삽입 → read 샘플링, 정답 SNP 목록 저장.
2. 벤치마크: 순진한 all-pairs Greedy Overlap.
3. 본 알고리즘: prefix/suffix 해시 색인, 양방향 연장.
4. 평가 도구: contig→reference 정렬, SNP 호출, 지표 계산.
5. 실험 스크립트와 결과 그래프.

## 8. 진행 상황

- [x] 문제 정의와 데이터 설정
- [x] 알고리즘 후보 비교와 선택
- [x] reference 생성기
- [x] read 생성기
- [x] 벤치마크(Greedy Overlap)
- [x] SSAKE 계열 구현 (단계 0)
- [x] 평가 도구
- [x] 실험 단계 0~4 (단계 5~6 제외)
- [x] 결과 정리

## 9. 결과

측정 환경: AMD Ryzen 5 7500F, DDR5 16GB, Windows 11, `g++ -O2 -std=c++17 -static`, 단일 스레드. 시간은 wall time, 메모리는 자식 포함 RSS를 20ms 간격으로 샘플링한 최댓값이다. 0.1초보다 짧은 실행은 샘플링이 안 되어 메모리가 0.0으로 찍히며, 이 값은 측정 불가를 뜻한다.

결과 원본: `results/step0.csv` ~ `results/step4.csv`. 그래프: `results/step1_time_memory.png`, `results/step2_contigs_vs_coverage.png`, `results/step3_recall_vs_snp.png`.

평가 수정 이력: 초기 평가 도구는 contig 끝 꼬리와 spurious anchor 때문에 recall과 accuracy를 낮게 찍었다. 수정 후(커밋 25d3ad5 이후) 값만 아래에 적는다.

### 단계 1: N에 따른 시간·메모리 (coverage 10, L 30, SNP 0.5%)

| N | SSAKE 시간(s) | SSAKE 메모리(MB) | Greedy 시간(s) | Greedy 메모리(MB) | contig 수 | misassembly | SNP recall | accuracy |
|---|---|---|---|---|---|---|---|---|
| 10K | 0.11 | 측정 불가 | 0.76 | 5.5 | 20 | 0 | 1.000 | 1.000000 |
| 100K | 0.07 | 16.9 | 59.3 | 12.2 | 215 | 4 | 1.000 | 1.000000 |
| 1M | 0.76 | 124.3 | - | - | 1920 | 83 | 0.991 | 0.999983 |
| 10M | 11.08 | 1189.1 | - | - | 19328 | 8757 | 0.905 | 0.999771 |

Greedy는 100K까지만 돌렸다(all-pairs 비용). N이 커지면 SSAKE 시간과 메모리는 거의 선형으로 늘었다. 10M에서 SNP recall이 0.905로 떨어진 원인은 아직 확인하지 않았다. 후보는 contig가 짧아져(N50 836) 정렬 anchor가 줄어드는 것이다. misassembly 8757은 1M 이상에서 늘었으며, 이것도 원인 확인이 필요하다.

### 단계 2: coverage에 따른 contig 수 (N 10K)

| coverage | 예상(LW) | SSAKE contig | Greedy contig | SNP recall (둘 다) |
|---|---|---|---|---|
| 5 | 137 | 130 | 130 | 1.000 |
| 10 | 22 | 20 | 20 | 1.000 |
| 15 | 3 | 3 | 3 | 1.000 |

contig 수는 LW 예상과 가깝다. 두 조립기 결과가 같다.

### 단계 3: SNP 비율 (N 10K, coverage 10)

| SNP 비율 | SSAKE recall | Greedy recall | misassembly |
|---|---|---|---|
| 0.1% | 1.000 | 1.000 | 0 |
| 0.5% | 1.000 | 1.000 | 0 |
| 1.0% | 1.000 | 1.000 | 0 |

SNP 비율을 1%까지 올려도 recall이 변하지 않았다. 예전 0.98~0.99 값은 평가 버그 때문이었다. 이 결과는 de novo 조립이 SNP에 영향을 받지 않는다는 것을 보인다. 단, 이 실험의 read는 sample에서 뽑으므로 SNP가 read에 그대로 들어 있다. 실험은 SNP 수에 따른 평가 결과만 보여준다.

### 단계 4: 파라미터 민감도 (N 10K, coverage 10)

| 변경 | 값 | SSAKE contig | SSAKE misassembly | SSAKE recall | Greedy contig | Greedy misassembly | Greedy recall |
|---|---|---|---|---|---|---|---|
| L (o = ⌈L/2⌉) | 15 (o 8) | 100 | 12 | 0.98 | 55 | 210 | 0.52 |
| L | 30 (o 15) | 20 | 0 | 1.00 | 20 | 0 | 1.00 |
| L | 45 (o 23) | 15 | 0 | 1.00 | 15 | 0 | 1.00 |
| o (L 30) | 10 | 4 | 0 | 1.00 | 4 | 0 | 1.00 |
| o (L 30) | 15 | 20 | 0 | 1.00 | 20 | 0 | 1.00 |
| o (L 30) | 20 | 101 | 0 | 1.00 | 101 | 0 | 1.00 |

L=15에서 Greedy는 misassembly 210, recall 0.52로 SSAKE(12, 0.98)보다 크게 나빴다. 가설: Greedy는 o 길이 겹침에서도 D=2 mismatch를 허용해서 우연한 겹침을 잘못 잇는다. SSAKE는 o 염기 seed를 정확히 맞춰야 연장한다. 이 가설은 아직 확인하지 않았다.

### 구현에서 정한 것 (원 논문과 다른 점)

- 중복 read는 두 조립기 모두 입력 단계에서 첫 번째 사본만 남긴다. 중복을 남기면 같은 위치 read가 독립 contig를 만들어 N=10K에서 contig 수가 301까지 부풀었다(제거 후 20).
- SSAKE 연장 후보가 여러 개면 입력 순서상 첫 번째를 쓴다. 다수결은 구현하지 않았다.
- 평가 accuracy는 contig 전체 염기 대비 sample과 일치하는 염기 비율이다. 정렬되지 않은 contig는 전부 불일치로 센다.
- 평가는 contig 끝까지 대각선 구간을 확장해 채점한다(꼬리 포함). 대각선 anchor가 하나뿐인 고립 anchor는 spurious로 보고 무시한다. 두 개 이상 연속된 spurious anchor는 아직 걸러내지 못한다.

## 10. 미정 사항 (확정)


- D의 의미: read와 reference 사이 mismatch 허용치로 쓴다. 겹침 검증에만 사용.
- 역상보 가닥: 기본 설정에서 제외. 선택 단계 6에서만 다룬다.
- 입출력 형식: read는 한 줄에 하나, contig는 FASTA, 정답 SNP는 `pos\tref\talt` TSV(0-based).
- 측정 PC 사양: 위 결과 섹션의 환경과 같다.

## 11. 참고

- Warren et al., SSAKE, Bioinformatics, 2007 (원문 확인 필요)
- Lander & Waterman, 1988 (끊김 수 근사)
- Eulerian path / De Bruijn 그래프 (비교 대상 알고리즘의 배경)
