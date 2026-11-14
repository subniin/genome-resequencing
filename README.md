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

E. coli K-12 MG1655 유전체의 일부를 reference로 사용합니다. reference에 SNP를 넣어 sample genome을 만들고, sample에서 read를 무작위 위치로 뽑습니다.

기본 설정:

| 항목 | 값 |
|---|---|
| reference 길이 N | 100,000 |
| read 개수 M | 10,000 |
| read 길이 L | 100 |
| coverage | 10 (= M×L/N) |
| SNP | 100개 (0.1%), 치환만 |
| 허용 mismatch D | 2 |
| repeat | 고려하지 않음 |

가정:

1. read는 sample의 정방향에서만, 시작 위치를 균등하게 뽑는다. 역상보 가닥과 paired-end는 쓰지 않는다.
2. read에는 SNP만 있고 시퀀싱 오류는 기본 설정에서 넣지 않는다. 오류 추가는 선택 확장이다.
3. 선택한 구간에 repeat가 있으면 구간을 다시 뽑는다.

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

Lander–Waterman 근사로 read 시작 밀도 λ = M/N = cov/L, 연결 조건 길이 L−o일 때 끊김 수는 약 N·λ·e^(−λ(L−o)) 입니다. N=100,000, L=100, o=30 기준:

| coverage | 예상 contig 수 |
|---|---|
| 1 | 거의 조립되지 않음 (염기의 약 63%만 덮임) |
| 5 | 약 150 |
| 10 | 약 9 |
| 15 | 약 1 |

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
| 1 | N = 10K, 100K, 1M, 약 5M | 시간·메모리 증가 곡선 |
| 2 | coverage = 1, 5, 10, 15 | 끊김 수, 복원율 |
| 3 | SNP = 0.1%, 0.5%, 1.0% | SNP 복원율 (de novo는 SNP가 조립을 방해하지 않으므로 거의 평평할 것으로 예상) |
| 4 | o(최소 겹침) 변경, L = 50 / 100 | 파라미터 민감도 |
| 5 (선택) | read 오류 최대 D개, 다수결 추가 | 오류 내성, 다수결 비용 |
| 6 (선택) | 역상보 가닥, repeat 포함 구간 | 알고리즘 한계 |

N=1M 이상에서는 E. coli에 repeat(rRNA operon 등)가 있어 greedy 방식이 그 지점에서 끊길 수 있습니다. 이는 알고리즘의 한계로 보고합니다.

## 7. 구현 계획

1. read 생성기: reference 구간 추출 → SNP 삽입 → read 샘플링, 정답 SNP 목록 저장.
2. 벤치마크: 순진한 all-pairs Greedy Overlap.
3. 본 알고리즘: prefix/suffix 해시 색인, 양방향 연장.
4. 평가 도구: contig→reference 정렬, SNP 호출, 지표 계산.
5. 실험 스크립트와 결과 그래프.

## 8. 진행 상황

- [x] 문제 정의와 데이터 설정
- [x] 알고리즘 후보 비교와 선택
- [ ] read 생성기
- [ ] 벤치마크(Greedy Overlap)
- [ ] SSAKE 계열 구현 (단계 0)
- [ ] 평가 도구
- [ ] 실험 단계 1~4
- [ ] 결과 정리

## 9. 미정 사항

- De novo에서 D의 의미: 현재는 read와 reference의 mismatch 허용치로 보고 겹침 검증에 쓴다.
- 역상보 가닥을 기본 설정에 넣을지.
- 입출력 형식: read는 한 줄에 하나, contig는 FASTA, 정답 SNP는 별도 목록.
- 측정 PC 사양.

## 10. 참고

- Warren et al., SSAKE, Bioinformatics, 2007 (원문 확인 필요)
- Lander & Waterman, 1988 (끊김 수 근사)
- Eulerian path / De Bruijn 그래프 (비교 대상 알고리즘의 배경)
