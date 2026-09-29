# GenSplice-Agent AntiGravity 통합 개발 청사진 (Blueprint)

유전자의 발현 총량(Quantity, DEG)뿐만 아니라 스플라이싱 형태의 질적 변화(Quality, Alternative Splicing)가 생체 조절 기전에서 갖는 중요성을 직관적으로 입증하기 위한 `GenSplice-Agent`의 AntiGravity 전용 통합 개발 청사진(Blueprint)입니다.

---

### 1. 핵심 생물학적 가설 및 시각화 컨셉

기존 유전체 연구의 대다수는 총 발현량(DEG)에만 매몰되어, "전체 유전자 발현량에는 유의미한 변화가 없으나, 기능적 엑손이 탈락되거나 인트론이 삽입되어 단백질 기능이 완전히 바뀌는 핵심 조절자(Alternative Splicing)"를 놓치는 중대한 한계를 지닙니다.

실제로 가뭄 스트레스 전사체 연구에서 총 발현량(DEG)에는 유의차가 없었으나 스플라이싱 변이체(Isoform)만 통계적으로 유의미하게 변화한 유전자(예: *GRMZM2G140355*)가 발견된 사례처럼, 발현량과 발현 양상을 동시에 추적해야 질환 및 노화의 진짜 기전을 포착할 수 있습니다.

이 도구는 DEG 결과와 rMATS 스플라이싱 결과를 단일 좌표계로 결합하여 전체 전사체를 다음 4개 사분면(Quadrant)으로 즉각 분류합니다:

```text
                  ▲ ΔPSI (Percent Spliced In 변화)
                  │
     [Q2] Splicing-Driven Only   │   [Q1] Dual Responders
   (총 발현량 불변, 스플라이싱만 극변)  │  (발현량도 변하고, 형태도 바뀜)
   ★ 숨겨진 핵심 스플라이싱 조절자 ★   │
──────────────────┼──────────────────▶ Log2 Fold Change (DEG 발현량)
     [Q3] Invariant Background   │   [Q4] Expression-Driven Only
        (변화 없는 유전자군)        │   (형태는 유지, 발현 총량만 증감)
                  │
```

---

### 2. Zero-Server-Cost 로컬 리눅스 CLI 및 단독 실행형 HTML 리포트 아키텍처

대용량 NGS 전사체 데이터(FASTQ 파일 당 수십 GB~수백 GB)를 웹 서버에 업로드하고 클라우드 인프라에서 분석/시각화하는 방식은 심각한 서버 유지비(AWS/GCP 용량 및 컴퓨팅 비용)와 데이터 전송 병목 현상을 유발합니다.

`GenSplice-Agent`는 **100% 로컬 리눅스/macOS CLI 기반**으로 구동되며, 데이터 정렬부터 DEG/AS 분석 및 최종 시각화 리포트까지 사용자의 로컬 컴퓨터에서 처리합니다.

#### 선행 연구 및 펍메드(PubMed) 논문 사례
학계 선행 연구들 역시 대용량 유전체 데이터의 웹 서버 비용 문제를 극복하기 위해 **로컬 CLI + 단독 실행형(Standalone) HTML 리포트 내보내기 방식**을 주요 패러다임으로 채택하여 논문을 게재해 왔습니다:

1. **MultiQC** (Ewels et al., *Bioinformatics* 2016, DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)):
   - 다양한 NGS 툴 결과를 단일 독립형 HTML 대시보드로 요약. 서버 설치 없이 로컬에서 즉시 열람 가능.
2. **fastp** (Chen et al., *Bioinformatics* 2018, DOI: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)):
   - Ultra-fast FASTQ preprocessor로, 자바스크립트/그래프가 임베디드된 로컬 독립형 HTML 리포트 출력.
3. **Degust** (Powell et al., Monash Univ.):
   - 로컬 RNA-Seq DEG 인터랙티브 시각화 툴로, 단일 HTML/JS 파일 추출 지원.
4. **maser** (Kinser et al., *Bioinformatics* 2022, DOI: [10.1093/bioinformatics/btac655](https://doi.org/10.1093/bioinformatics/btac655)):
   - rMATS 스플라이싱 이벤트 시각화 R 패키지로, 리포트 생성 및 비주얼 추출 지원.
5. **RNA-Seq-Pop** (Incorvaia et al., *GigaScience* 2021):
   - 집단 전사체 파이프라인으로 로컬 실행 후 독립형 HTML 리포트 내보내기.

`GenSplice-Agent`는 이 검증된 학술적 패러다임을 계승하여 **서버 비용 0원(Zero Cloud Server Cost)**으로 동작하는 단독 실행형 HTML 분석 엔진을 완성합니다.

---

### 3. 검증된 오픈소스 자산 활용 명세

바퀴를 다시 발명하지 않고, 검증된 학계 표준의 데이터 규격과 오픈소스 UI/통계 로직을 결합합니다.

1. **입력 데이터 규격**:
   - **DEG 자산**: DESeq2 / edgeR의 표준 출력 포맷 (`gene_id`, `log2FoldChange`, `pvalue`, `padj`).
   - **Alternative Splicing 자산**: rMATS v4.x 표준 출력 포맷 (`SE.MATS.JC.txt`, `RI.MATS.JC.txt` 등 5대 이벤트 파일).

2. **분석 및 UI 인터페이스 차용**:
   - **maser (Bioconductor)**: rMATS의 5개 텍스트 파일을 묶어 읽어 들이는 파싱 로직 및 FDR/$\Delta\text{PSI}$ 필터링 기준 차용.
   - **Degust (Monash Univ.)**: 수만 개 행을 버벅임 없이 슬라이더로 실시간 필터링하는 반응형 데이터 연동 방식 차용.
   - **BioChatter (Helmholtz Institute)**: 사용자가 자신의 API Key를 직접 넣는 BYOK(Bring Your Own Key) 보안 아키텍처 및 생체 지식 프롬프트 템플릿 구조 차용.

3. **고속 데이터 및 벤치마크 엔진**:
   - Pandas/R 대비 **Polars**를 사용하여 수만 개의 유전자 카운트 및 스플라이싱 좌표를 0.05초 내에 `gene_id` 기준으로 Inner/Outer Join.
   - 메모리(RAM Peak RSS) 및 실행 속도 벤치마크 평가 모듈 탑재 (`core/benchmark.py`).

---

### 4. 확장 신규 모듈 명세

`GenSplice-Agent`의 학술적 가치와 wet-lab 검증 연계성을 극대화하기 위해 다음 3대 신규 엔진을 추가 구성합니다.

#### 모듈 1: 초고속 벤치마크 평가 엔진 (`core/benchmark.py`)
- **목적**: Polars vs Pandas vs R 간의 데이터 처리 속도(Execution Time) 및 메모리 점유율(Peak Memory Usage, MB)을 정량적 벤치마크 측정.
- **기능**:
  - 수천~수십만 행의 DEG 및 rMATS 대용량 TSV 데이터를 대상으로 조인/필터링 성능 측정.
  - 리포트 대시보드에 벤치마크 결과 그래프 및 요약 표 제공.

#### 모듈 2: 엑손-인트론 구조 & Sashimi Visualizer (`visualizer/exon_structure.py`)
- **목적**: Q1/Q2 핵심 타깃 유전자의 스플라이싱 이벤트(SE, RI, A5SS, A3SS, MXE) 엑손 구조와 접합부(Junction) 카운트를 직관적인 베터 그래픽으로 시각화.
- **기능**:
  - GTF 좌표 기반 엑손(Exon) - 인트론(Intron) 영역 자동 매핑.
  - 대조군(Control) vs 처리군(Treatment)의 Junction Read Count 및 Inclusion/Exclusion 엑손 스킵 다이어그램 생성 (Plotly/HTML 호환).

#### 모듈 3: Isoform-Specific RT-qPCR Primer Design Engine (`core/primer_designer.py`)
- **목적**: Q1 및 Q2 사분면에 위치한 타깃 유전자의 스플라이싱 변이(Isoform)를 Wet-lab에서 즉시 검증할 수 있는 맞춤형 RT-qPCR 프라이머 자동 설계.
- **기능**:
  - Isoform-specific Exon-Exon Junction을 교차하는 프라이머 쌍(Forward/Reverse) 자동 생성.
  - $T_m$ (녹는 온도: $58\sim 62^\circ\text{C}$), GC content ($40\sim 60\%$), Self-dimerization, Hairpin risk 및 PCR 생성물 크기(Amplicon size: $80\sim 200\text{ bp}$) 자동 검증.
  - HTML 리포트에 Primer Table 및 주문용 Sequence 복사 기능 제공.

---

### 5. AntiGravity 기반 시스템 아키텍처 및 디렉토리 구조

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
├── outputs/                # 중간 결과 및 최종 결과 자동 저장 폴더
│   ├── 01_clean_fq/        # fastp QC/트리밍 결과 FASTQ
│   ├── 02_aligned_bam/     # STAR 정렬 결과 BAM 파일
│   ├── 03_deg/             # DESeq2 / featureCounts 정량 결과
│   └── 04_rmats/           # rMATS 5대 이벤트 분석 결과
├── app.py                  # GenSplice-Agent Streamlit 대시보드 앱
├── config.py               # 기본 임계값(FDR, ΔPSI, Log2FC) 및 색상 테마
├── core/
│   ├── __init__.py
│   ├── deg_loader.py       # DESeq2/edgeR 결과 CSV/TSV 파서 (Polars)
│   ├── rmats_loader.py     # rMATS 5대 이벤트 파일 병합 로더 (Polars)
│   ├── merger.py           # DEG + AS 테이블 간의 고속 Join 및 사분면 라벨러
│   ├── benchmark.py        # [신규] Polars 대용량 데이터 벤치마크 평가 엔진
│   └── primer_designer.py  # [신규] Isoform-Specific RT-qPCR 프라이머 설계 엔진
├── visualizer/
│   ├── __init__.py
│   ├── quadrant_plot.py    # Log2FC vs ΔPSI 4사분면 인터랙티브 산점도 (Plotly)
│   ├── dual_volcano.py     # DEG Volcano와 AS Volcano를 나란히 배치한 듀얼 뷰
│   ├── exon_structure.py  # [신규] 엑손-인트론 구조 & Sashimi-style 비주얼 엔진
│   └── report_exporter.py  # 단독 실행형 HTML 독립 리포트 생성기
├── ai/
│   ├── __init__.py
│   └── gemini_evaluator.py # google-genai SDK 기반 양적/질적 변화 종합 해석기
├── requirements.txt        # 파이썬 의존성 패키지 명세서
└── test_data/              # 검증용 Mock 데이터셋 (DEG 1개, rMATS 5개)
```

---

### 6. 데이터 통합 로직 및 사분면 매핑 알고리즘 (`core/merger.py`)

DEG와 Alternative Splicing 데이터는 유전자 심볼(`geneSymbol` 또는 `GeneID`)을 키로 결합됩니다.

#### (1) 조인 및 결측치 처리 규칙
- 하나의 유전자에 여러 개의 스플라이싱 이벤트가 존재할 경우: $\vert{}\Delta\text{PSI}\vert{}$ 값이 가장 크거나 FDR이 가장 낮은 대표 이벤트를 유전자 단위의 대표 스플라이싱 이벤트로 선정.
- DEG에는 존재하나 AS 이벤트가 검출되지 않은 경우: $\Delta\text{PSI} = 0, \text{Event} = \text{'None'}$으로 처리.
- AS에는 존재하나 DEG에 없는 경우: $\text{Log}_2\text{FC} = 0, \text{FDR}_{\text{DEG}} = 1.0$으로 처리.

#### (2) 4사분면 자동 분류 기준 수식
- **유의미 임계치**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge \theta_{\text{DEG}}$ (기본값 $0.5$ 또는 $1.0$), $\vert{}\Delta\text{PSI}\vert{} \ge \theta_{\text{AS}}$ (기본값 $0.1$), $\text{FDR} \le 0.05$
- **분류 체계**:
  - **Q1 (Dual Impact)**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge \theta_{\text{DEG}} \land \vert{}\Delta\text{PSI}\vert{} \ge \theta_{\text{AS}}$
  - **Q2 (Splicing-Driven / Masked)**: $\vert{}\text{Log}_2\text{FC}\vert{} < \theta_{\text{DEG}} \land \vert{}\Delta\text{PSI}\vert{} \ge \theta_{\text{AS}}$ $\rightarrow$ **이 연구의 핵심 어필 유전자군**
  - **Q4 (Abundance-Driven)**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge \theta_{\text{DEG}} \land \vert{}\Delta\text{PSI}\vert{} < \theta_{\text{AS}}$
  - **Q3 (Background / Static)**: 나머지 전체 (기본 시각화 Off 처리로 리소스 최적화)

---

### 7. BMC Bioinformatics 논문 출판 블루프린트 (Publication Roadmap)

#### 목표 저널
- **BMC Bioinformatics** (IF: ~3.0+, Software Article Section) 또는 **Bioinformatics (Oxford)**.

#### 논문 제목 제안
> **GenSplice-Agent: A local-first transcriptomic pipeline and interactive standalone exporter for dual-impact gene discovery and isoform-specific RT-qPCR validation**

#### 핵심 학술적 소구 포인트 (Novelty & Value Proposition)
1. **Zero-Server-Cost Standalone Exporter**:
   - expensive cloud web server 없이 100% 로컬 리눅스에서 동작하며, 단일 `.html` 리포트 내보내기로 대용량 전사체 데이터 공유 및 브라우저 탐색 가능.
2. **Dual-Impact 4-Quadrant Transcriptomic Cross-Plot**:
   - DEG(양적)와 Alternative Splicing(질적) 분석을 최초로 직관적 4사분면 매핑으로 통합하여 hidden splicing regulators (Q2) 발견.
3. **Integrated Wet-lab Validation Engine**:
   - Q1/Q2 타깃 유전자에 대한 Isoform-specific RT-qPCR 프라이머 자동 설계 및 Sashimi-like 엑손 구조 시각화 연동.
4. **High-Performance Polars Benchmarking**:
   - Pandas 대비 최대 10x 이상 빠른 Polars 기반 데이터 변환으로 수만 개 유전자 실시간 4사분면 재분류 지원.

---

### 8. 의존성 및 패키지 명세 (`requirements.txt`)

```text
streamlit>=1.35.0
polars>=0.20.0
pyarrow>=15.0.0
plotly>=5.20.0
google-genai>=0.1.0
pandas>=2.0.0
openpyxl>=3.1.0
```
