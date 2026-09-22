# GenSplice-Agent 🧬

> **통합 전사체 발현량(DEG) 및 대립적 스플라이싱(Alternative Splicing) 시각화·AI 해석 대시보드**

`GenSplice-Agent`는 유전자의 발현량 변화(Quantity, DEG)와 스플라이싱 형태의 질적 변화(Quality, Alternative Splicing)를 단일 사분면 좌표계로 통합하고, Google Gemini AI를 통해 분자생물학적 기전을 해설하는 전사체학 전문 통합 대시보드입니다.

---

## 📦 예상 의존성 목록 (Conda/Bioconda Ecosystem)

서버 환경 간 라이브러리 충돌을 방지하기 위해 **Conda / Mamba** 기반 완전 환경 격리를 제공합니다.

### 1. 생물정보학 CLI 분석 도구 (Bioconda)
| 도구명 | 추천 버전 | 용도 및 역할 |
| :--- | :--- | :--- |
| **`fastp`** | `>= 0.23.4` | FASTQ 시퀀싱 데이터 품질 관리(QC) 및 어댑터 트리밍 |
| **`STAR`** | `>= 2.7.11a` | RNA-seq 리드 게놈 정렬 및 스플라이스 정렬 엔진 |
| **`rmats`** | `>= 4.3.0` | Alternative Splicing (SE, RI, A5SS, A3SS, MXE) 분석 도구 |
| **`subread`** | `>= 2.0.6` | `featureCounts` 전사체 정량 산출 CLI 도구 |
| **`R` (r-base)** | `>= 4.2` | rMATS 및 DESeq2 통계 검증 R 백엔드 연동 |

### 2. 파이썬 웹 대시보드 라이브러리 (`environment.yml`)
| 패키지명 | 최소 버전 | 용도 |
| :--- | :--- | :--- |
| `streamlit` | `>= 1.35.0` | 대시보드 웹 인터페이스 GUI |
| `polars` | `>= 0.20.0` | 수만 건 데이터 고속 결합(Join) 및 필터링 엔진 |
| `pyarrow` | `>= 15.0.0` | Polars 고속 인메모리 파케/테이블 처리 |
| `plotly` | `>= 5.20.0` | 4사분면 cross-plot 및 듀얼 볼케이노 인터랙티브 시각화 |
| `google-genai` | `>= 0.1.0` | Gemini API 연동 생물학적 기전 자동 해석 에이전트 |
| `pandas` | `>= 2.0.0` | 시각화 및 호환용 데이터 프레임 연동 |
| `openpyxl` | `>= 3.1.0` | 엑셀 데이터 파일 입출력 지원 |

---

## 🚀 1. 원클릭 환경 구축 (`install.sh`)

`bash install.sh` (또는 `./install`) 단 한 줄만 실행하면 **Conda/Mamba 감지, `environment.yml` 기반 가상환경 생성, 생물정보학 CLI 도구 설치 및 무결성 검증까지 원클릭으로 완료**됩니다.

```bash
git clone https://github.com/skt526/GenSplice-Agent.git
cd GenSplice-Agent
bash install.sh
```

---

## 🧬 2. 표준 레퍼런스 게놈 다운로더 (`ref`)

명령어 `./ref` (또는 `bash ref.sh`)를 실행하면 인터랙티브 메뉴를 통해 **최신 버전의 Reference Genome (FASTA) 및 유전자 주석(GTF) 데이터**를 즉각 다운로드하고 자동으로 압축 해제하여 `./{선택한 종}-ref/` 폴더에 세팅합니다.

### 인터랙티브 메뉴 실행
```bash
./ref
```
```text
Please select a Reference Genome to download:
  1) Human (Homo sapiens - GRCh38)
  2) Mouse (Mus musculus - GRCm39)
  3) Rice (Oryza sativa - IRGSP-1.0)
  4) Arabidopsis thaliana (TAIR10)
  5) Maize / Corn (Zea mays - B73)
  6) Exit
```

### 파라미터 직접 지정 실행
```bash
./ref human        # ./human-ref/ 생성 (GRCh38 FASTA + GTF)
./ref mouse        # ./mouse-ref/ 생성 (GRCm39 FASTA + GTF)
./ref rice         # ./rice-ref/ 생성 (IRGSP-1.0 FASTA + GTF)
./ref arabidopsis  # ./arabidopsis-ref/ 생성 (TAIR10 FASTA + GTF)
./ref maize        # ./maize-ref/ 생성 (Zm-B73-NAM-5.0 FASTA + GTF)
```

## 📁 디렉토리 구조 (Directory Architecture)

```text
GenSplice-Agent/
├── install.sh              # 1-Click 환경 설치 스크립트 (Conda/Bioconda)
├── ref.sh                  # 레퍼런스 게놈 (FASTA, GTF) 원클릭 다운로더
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

## 🖥 3. 대시보드 실행

환경 구축 후 아래 명령어로 Conda 환경을 활성화하고 대시보드를 구동합니다:

```bash
conda activate gensplice-agent
streamlit run app.py
```

---

## 📄 문서
- [BLUEPRINT.md](BLUEPRINT.md): GenSplice-Agent 시스템 설계 청사진 및 알고리즘 명세
- [environment.yml](environment.yml): Conda/Bioconda 환경 명세서
