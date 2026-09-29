# GenSplice-Agent AntiGravity 통합 개발 청사진 (Blueprint)

유전자의 발현 총량(Quantity, DEG)뿐만 아니라 스플라이싱 형태의 질적 변화(Quality, Alternative Splicing)가 생체 조절 기전에서 갖는 중요성을 직관적으로 입증하기 위한 `GenSplice-Agent`의 SCIE 저널 게재용 통합 개발 청사진(Blueprint)입니다.

---

### 1. 핵심 생물학적 가설 및 시각화 컨셉

기존 유전체 연구의 대다수는 총 발현량(DEG)에만 매몰되어, "전체 유전자 발현량에는 유의미한 변화가 없으나, 기능적 엑손이 탈락되거나 인트론이 삽입되어 단백질 기능이 완전히 바뀌는 핵심 조절자(Alternative Splicing)"를 놓치는 중대한 한계를 지닙니다.

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

### 2. SCIE 저널 게재용 핵심 차별화 모듈 (Technical Novelty)

`GenSplice-Agent`는 단순 파이프라인 래퍼를 넘어, 학술적 독창성(Novelty)을 확보하기 위해 다음 4대 핵심 기술 모듈을 확장 탑재합니다.

#### 🧬 모듈 A: NMD & Frameshift 고위험군(High-Risk) 자동 판정 엔진 (`core/risk_evaluator.py`)
스플라이싱 이벤트가 단백질 기능에 미치는 치명적 영향을 게놈 좌표 레벨에서 정밀 분석하여 위험도를 분류합니다.
- 🚨 **HIGH RISK**: 
  - **Frame-shift**: 포함/탈락 엑손 길이($L$)가 3의 배수가 아님 ($\Delta L \pmod 3 \neq 0$).
  - **NMD (Nonsense-Mediated Decay)**: 프레임시프트로 인해 조기 종단 코돈(PTC, Premature Termination Codon)이 생성되어 mRNA가 자동 파괴되는 경우 (50-nt rule 적용).
- ⚠️ **MEDIUM RISK**:
  - **In-frame Domain Loss**: $3$의 배수 변이로 프레임은 유지되나, 주요 단백질 도메인(Pfam/UniProt CDS 영역)이 소실된 경우.
- ℹ️ **LOW RISK**:
  - 5'-UTR / 3'-UTR 스플라이싱 또는 비암호화 영역 변형.

#### 🎨 모듈 B: 정상 vs 변이 Isoform 구조 비교 시각화 (`visualizer/isoform_switch.py`)
- Control 대표 Isoform과 Treatment에서 전환된 Isoform의 엑손-인트론 구조를 Plotly 기반 듀얼 트랙으로 시각화.
- 변이 엑손(Exon Skipping 등)을 빨간색/주황색으로 강조(Highlight)하고 소실된 도메인 및 예상되는 단백질 기능 변화를 직관적으로 제시.

#### 🧪 모듈 C: rMAPS 기반 상위 RBP (RNA-Binding Protein) 모티프 역추적 엔진 (`core/rmaps_wrapper.py`)
- **목적**: Q1/Q2 사분면 유전자들의 스플라이싱 변화를 유발한 상위 조절 인자(Master Splicing Factor) 역추적.
- **기능**:
  - rMATS 결과(SE, RI, A5SS, A3SS, MXE)를 입력으로 하여 **rMAPS (rMAPS2)** 백엔드 모듈 연동.
  - 변이 엑손 및 양측 250bp 인트론 영역에서 RBP 결합 모티프(CIS-BP-RNA / RBPDB 기준)의 통계적 농축(Enrichment) 분석 수행.
  - SRSF1, PTBP1, hnRNPA1 등 어떤 스플라이싱 인자가 해당 전사체 변동의 주범인지 예측 결과 제공.

#### 🤖 모듈 D: Multi-LLM (Gemini/ChatGPT) & No-API Dual 에이전트 섹션 (`ai/multi_llm_evaluator.py`)
대시보드 상에 **AI 해석 섹션**을 구성하고 듀얼 엑세스 모드를 제공합니다.
1. **API 입력 모드 (Gemini API Key 또는 ChatGPT API Key)**:
   - 사용자가 Google Gemini Key 또는 OpenAI ChatGPT Key를 입력하면, DEG/AS 수치 및 NMD/Frameshift 위험도 결과를 종합 분석.
   - 해당 약물/실험 조건이 특정 신호 전달 경로를 **어떤 기전(Mechanism of Action, MoA)으로 작용 또는 망가뜨리고 있는지** 생물학적 보고서 자동 생성.
2. **No-API 버튼 (Fallback Mode)**:
   - API 키가 없거나 사용을 원치 않을 경우 **"No-API (PubMed 논문 검색)"** 버튼 제공.
   - 버튼 클릭 시 타깃 유전자별 **PubMed REST API 기반 실시간 논문 검색 및 문헌 링크 섹션**으로 즉시 이동.

---

### 3. 검증된 오픈소스 자산 및 기술 스택

1. **입력 데이터 규격**:
   - **DEG 자산**: DESeq2 / edgeR의 표준 출력 포맷 (`gene_id`, `log2FoldChange`, `pvalue`, `padj`).
   - **Alternative Splicing 자산**: rMATS v4.x 표준 출력 포맷 (`SE.MATS.JC.txt`, `RI.MATS.JC.txt` 등 5대 이벤트 파일).
2. **백엔드 분석 & RBP 엔진**:
   - **rMATS / rMAPS**: Differential Alternative Splicing 및 RBP motif enrichment 계산 백엔드.
   - **Polars**: 수만 개의 유전자 카운트 및 스플라이싱 좌표를 0.05초 내 고속 Join/Filter.
3. **인터페이스 & 에이전트**:
   - **Streamlit**: 파이썬 기반 인터랙티브 대시보드 UI.
   - **google-genai / openai SDK**: Gemini 2.5/3.0 및 GPT-4o API 호출 기반 분자 기전 해석.

---

### 4. AntiGravity 기반 시스템 아키텍처 및 디렉토리 구조

```text
GenSplice-Agent/
├── install.sh              # 1-Click 환경 설치 스크립트 (Conda/Bioconda: rMATS, rMAPS 포함)
├── ref.sh                  # 레퍼런스 게놈 (FASTA, GTF) 원클릭 다운로더
├── run_pipeline.py         # 전체 파이프라인 오케스트레이터 (Python 메인 실행기)
├── config.yaml             # 참조 유전체 경로 및 파이프라인 스레드 설정
├── environment.yml         # Conda/Bioconda 환경 패키지 명세서 (rMAPS2 포함)
├── inputs/                 # [사용자 FASTQ 데이터 투입 폴더]
│   ├── control/            # 대조군 FASTQ 파일들 (.fastq / .fq.gz)
│   └── treatment/          # 실험군 FASTQ 파일들 (.fastq / .fq.gz)
├── outputs/                # 결과 자동 저장 폴더
│   ├── 01_clean_fq/        # fastp QC/트리밍 결과 FASTQ
│   ├── 02_aligned_bam/     # STAR 정렬 결과 BAM 파일
│   ├── 03_deg/             # DESeq2 / featureCounts 정량 결과
│   ├── 04_rmats/           # rMATS 5대 이벤트 분석 결과
│   └── 05_rmaps/           # [신규] rMAPS RBP 모티프 농축 분석 결과
├── app.py                  # GenSplice-Agent Streamlit 대시보드 앱 (AI & PubMed 섹션 통합)
├── config.py               # 기본 임계값(FDR, ΔPSI, Log2FC) 및 위험도 라벨 설정
├── core/
│   ├── __init__.py
│   ├── deg_loader.py       # DESeq2/edgeR 결과 TSV 파서 (Polars)
│   ├── rmats_loader.py     # rMATS 5대 이벤트 파일 병합 로더 (Polars)
│   ├── merger.py           # DEG + AS 테이블 간 고속 Join 및 4사분면 라벨러
│   ├── risk_evaluator.py   # [신규] NMD & Frameshift 고위험군 자동 판정 모듈
│   └── rmaps_wrapper.py    # [신규] rMAPS RBP 모티프 역추적 파이프라인 래퍼
├── visualizer/
│   ├── __init__.py
│   ├── quadrant_plot.py    # Log2FC vs ΔPSI 4사분면 인터랙티브 산점도 (Plotly)
│   ├── dual_volcano.py     # DEG Volcano & AS Volcano 듀얼 뷰
│   ├── isoform_switch.py   # [신규] 정상 vs 변이 Isoform 구조 비교 시각화 엔진
│   └── rbp_motif_plot.py   # [신규] RBP 모티프 농축 P-value 워터폴 그래프
├── ai/
│   ├── __init__.py
│   ├── multi_llm_evaluator.py # [신규] Gemini / ChatGPT / No-API Dual 해석 에이전트
│   └── pubmed_fetcher.py   # [신규] No-API 모드용 PubMed REST API 논문 검색기
├── requirements.txt        # 파이썬 의존성 패키지 명세서
└── test_data/              # 검증용 Mock 데이터셋 (DEG 1개, rMATS 5개)
```

---

### 5. AI & No-API 섹션 인터랙션 로직 명세 (`ai/multi_llm_evaluator.py`)

#### (1) UI 렌더링 구성
```text
┌─────────────────────────────────────────────────────────────────────────┐
│ 🤖 GenSplice AI Mechanism Interpretation & Literature Explorer          │
├─────────────────────────────────────────────────────────────────────────┤
│ Select Provider: [ Gemini (google-genai) | ChatGPT (openai) ]           │
│ API Key: [ Enter your API Key here...                              ]     │
│                                                                         │
│ 🔘 [ Run AI MoA Analysis ]      OR      🔘 [ No-API: PubMed Search ]    │
└─────────────────────────────────────────────────────────────────────────┘
```

#### (2) 분기 처리 로직
- **[Run AI MoA Analysis] 클릭 시**:
  - 선택된 유전자(Q1/Q2/Q4)의 $\text{Log}_2\text{FC}$, $\Delta\text{PSI}$, NMD/Frameshift 판정 결과, rMAPS 상위 RBP 모티프 결과를 프롬프트로 구성.
  - Gemini / ChatGPT API를 호출하여 약물의 **기전(MoA) 및 세포 신호전달망 교란 영향** 보고서 출력.
- **[No-API: PubMed Search] 클릭 시**:
  - LLM API 호출 없이 아래 **PubMed Paper Explorer** 섹션으로 즉시 스크롤 이동.
  - 선택 유전자 심볼(e.g., `CRISPLD2 alternative splicing`) 기반 NCBI PubMed 검색 결과 및 바로가기 URL 카드 생성.

---

### 6. SCIE 저널 출판 로드맵 (Publication Strategy)

#### 목표 저널
- **BMC Bioinformatics** (SCIE, IF ~3.0)
- **Frontiers in Genetics / Frontiers in Bioinformatics** (SCIE, IF ~2.8~3.2)
- **Genes** (MDPI, SCIE, IF ~2.8)
- **PLOS ONE** (SCIE, IF ~2.9)

#### 논문 제목 제안
> **GenSplice-Agent: An integrative pipeline for dual-impact transcriptomic cross-plotting, NMD risk assessment, and RBP motif tracking with multi-LLM interpretation**

#### 학술적 논문 소구점 (Novelty Summary)
1. **Dual-Impact 4-Quadrant Mapping**: DEG(양적)와 AS(질적) 통합을 통한 hidden splicing regulator (Q2) 발견.
2. **Automated NMD & Frameshift Risk Profiling**: 게놈 좌표 기반 엑손 스킵의 치명적 도메인/NMD 위험도 자동 분류.
3. **rMAPS-Integrated Upstream RBP Discovery**: 변이 엑손 상위 조절 RBP 모티프의 통계적 농축 파이프라인 탑재.
4. **Flexible Multi-LLM / No-API Hybrid Architecture**: Gemini/ChatGPT 키 입력 기반 분자 기전 자동 요약 및 No-API PubMed 문헌 검색 듀얼 지원.

---

### 7. 파이썬 의존성 명세 (`requirements.txt`)

```text
streamlit>=1.35.0
polars>=0.20.0
pyarrow>=15.0.0
plotly>=5.20.0
google-genai>=0.1.0
openai>=1.12.0
requests>=2.31.0
pandas>=2.0.0
openpyxl>=3.1.0
biopython>=1.81
```
