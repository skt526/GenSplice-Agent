# GenSplice-Agent 🧬

> **통합 전사체 발현량(DEG) 및 대립적 스플라이싱(Alternative Splicing) 시각화·AI 해석 대시보드**

`GenSplice-Agent`는 유전자의 발현량 변화(Quantity, DEG)와 스플라이싱 형태의 질적 변화(Quality, Alternative Splicing)를 단일 사분면 좌표계로 통합하고, Google Gemini AI를 통해 분자생물학적 기전을 해설하는 전사체학 전문 통합 플랫폼입니다.

---

## 🚦 명령 워크플로우 (Command Workflow)

### 🧪 1. 시스템 동작 검증용 (Testing)
설치 직후 예제 데이터셋(GSE52778 Dexamethasone 모델)으로 파이프라인 전 과정이 이상 없이 작동하는지 검증합니다:
```bash
bash install.sh ➔ ./test
```
*(또는 `./install` ➔ `./test`)*

### 🧬 2. 실제 전사체 샘플 분석용 (Real Analysis)
원하는 레퍼런스 게놈(Human, Mouse, Rice 등)을 세팅하고 실제 파이프라인 분석을 구동합니다:
```bash
bash install.sh ➔ ./ref human ➔ ./GenSplice
```
*(또는 `./install` ➔ `./ref human` ➔ `./GenSplice`)*

---

## 📦 필수 의존성 환경 (Conda/Bioconda Ecosystem)

| 구분 | 도구/패키지명 | 추천 버전 | 용도 및 역할 |
| :--- | :--- | :--- | :--- |
| **CLI 분석 도구** | `fastp` | `>= 0.23.4` | FASTQ 시퀀싱 데이터 품질 관리(QC) 및 어댑터 트리밍 |
| | `STAR` | `>= 2.7.11a` | RNA-seq 리드 게놈 정렬 및 2-pass 정렬 엔진 |
| | `rmats` | `>= 4.3.0` | Alternative Splicing (SE, RI, A5SS, A3SS, MXE) 분석 |
| | `subread` | `>= 2.0.6` | `featureCounts` 전사체 정량 산출 |
| | `R` (r-base) | `>= 4.2` | rMATS 및 DESeq2 통계 검증 백엔드 |
| **파이썬 GUI** | `streamlit` | `>= 1.35.0` | 웹 대시보드 인터페이스 |
| | `polars` | `>= 0.20.0` | 고속 데이터 결합 및 4사분면 파서 |
| | `plotly` | `>= 5.20.0` | 4사분면 cross-plot & Dual Volcano 시각화 |
| | `google-genai` | `>= 0.1.0` | Gemini API 기반 생물학적 기전 해석 에이전트 |

---

## 🧪 검증용 데이터셋 명세 (GSE52778 Airway Smooth Muscle)

`./test` 구동 시 사용되는 예제 데이터셋은 학계 검증 표준 세트입니다:
- **GEO Accession**: GSE52778 (SRA: SRP033346) (Himes et al., 2014)
- **연구 모델**: Human Airway Smooth Muscle Cell (Control vs Dexamethasone 처리군)
- **검증 타깃 유전자**: 면역/스플라이싱 타깃 유전자 (*DUSP1*, *CRISPLD2*)의 DEG 및 5대 스플라이싱 변이 자동 검출 확인.

---

## 📁 디렉토리 구조 (Directory Architecture)

```text
GenSplice-Agent/
├── install.sh              # 1-Click 환경 설치 스크립트 (Conda/Bioconda)
├── test.sh                 # 시스템 작동 검증 테스트 스크립트 (GSE52778)
├── ref.sh                  # 레퍼런스 게놈 (FASTA, GTF) 원클릭 다운로더
├── GenSplice               # 메인 파이프라인 실행 래퍼 명령 스크립트
├── run_pipeline.py         # 전체 파이프라인 오케스트레이터 (Python 메인 실행기)
├── config.yaml             # 참조 유전체 경로 및 파이프라인 스레드 설정
├── environment.yml         # Conda/Bioconda 환경 패키지 명세서
├── inputs/                 # [사용자 FASTQ 데이터 투입 폴더]
│   ├── control/            # 대조군 FASTQ 파일들 (.fastq / .fq.gz)
│   └── treatment/          # 실험군/노화 FASTQ 파일들 (.fastq / .fq.gz)
├── {organism}-ref/         # 참조 유전체(FASTA, GTF) 및 STAR 인덱스 저장 폴더
├── outputs/                # 중간 결과 및 최종 결과 자동 저장 폴더
│   ├── 01_clean_fq/        # fastp QC/트리밍 결과 FASTQ
│   ├── 02_aligned_bam/     # STAR 정렬 결과 BAM 파일
│   ├── 03_deg/             # DESeq2 / featureCounts 정량 결과
│   └── 04_rmats/           # rMATS 5대 이벤트 분석 결과
└── app.py                  # GenSplice-Agent Streamlit 대시보드 앱
```

---

## 📥 입력 데이터 포맷 명세 (Input FASTQ Specification)

실제 RNA-seq 데이터를 분석할 때, `inputs/` 디렉터리에 대조군(Control)과 실험군(Treatment) 샘플 파일들을 배치합니다:

```text
inputs/
├── control/
│   ├── sampleA_1.fq.gz
│   └── sampleA_2.fq.gz
└── treatment/
    ├── sampleB_1.fq.gz
    └── sampleB_2.fq.gz
```

### 1. 지원 파일 확장자 (File Extensions)
- **압축 파일 (권장)**: `.fq.gz`, `.fastq.gz`
- **비압축 파일**: `.fq`, `.fastq`

### 2. 페어드 엔드 (Paired-End) 파일 명명 규칙
자동 샘플 쌍 매핑을 위해 아래 접미사 패턴 중 하나를 사용해야 합니다:
- **패턴 1**: `*_1.fq.gz` / `*_2.fq.gz` 또는 `*_1.fastq.gz` / `*_2.fastq.gz`
- **패턴 2**: `*_R1.fastq.gz` / `*_R2.fastq.gz` 또는 `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`

### 3. 싱글 엔드 (Single-End) 지원
- `sampleA.fq.gz` 와 같이 단일 파일로 투입된 경우 자동으로 Single-End 데이터로 감지되어 QC 및 Alignment가 수행됩니다.

---

## 📄 문서
- [BLUEPRINT.md](BLUEPRINT.md): GenSplice-Agent 시스템 설계 청사진 및 알고리즘 명세
- [environment.yml](environment.yml): Conda/Bioconda 환경 명세서
