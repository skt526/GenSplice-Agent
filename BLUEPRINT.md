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

### 2. 검증된 오픈소스 자산 활용("긴빠이") 명세

바퀴를 다시 발명하지 않고, 검증된 학계 표준의 데이터 규격과 오픈소스 UI/통계 로직을 결합합니다.

1. **입력 데이터 규격**:
* **DEG 자산**: DESeq2 / edgeR의 표준 출력 포맷 (`gene_id`, `log2FoldChange`, `pvalue`, `padj`).
* **Alternative Splicing 자산**: rMATS v4.x 표준 출력 포맷 (`SE.MATS.JC.txt`, `RI.MATS.JC.txt` 등 5대 이벤트 파일).


2. **분석 및 UI 인터페이스 차용**:
* **maser (Bioconductor)**: rMATS의 5개 텍스트 파일을 묶어 읽어 들이는 파싱 로직 및 FDR/$\Delta\text{PSI}$ 필터링 기준 차용.
* **Degust (Monash Univ.)**: 웹 브라우저 단에서 수만 개 행을 버벅임 없이 슬라이더로 실시간 필터링하는 반응형 데이터 연동 방식 차용.
* **BioChatter (Helmholtz Institute)**: 사용자가 자신의 API Key를 직접 넣는 BYOK(Bring Your Own Key) 보안 아키텍처 및 생체 지식 프롬프트 템플릿 구조 차용.


3. **고속 데이터 엔진**:
* Pandas 대신 **Polars**를 사용하여 수만 개의 유전자 카운트 및 스플라이싱 좌표를 0.05초 내에 `gene_id` 기준으로 Inner/Outer Join.



---

### 3. AntiGravity 기반 시스템 아키텍처 및 디렉토리 구조

로컬 머신의 AntiGravity 환경에서 원클릭으로 구동할 수 있는 경량 모듈형 구조입니다.

```text
GenSplice-Agent/
├── app.py                      # Streamlit 메인 대시보드 UI 컨트롤러
├── config.py                   # 기본 임계값(FDR, ΔPSI, Log2FC) 및 색상 테마
├── core/
│   ├── __init__.py
│   ├── deg_loader.py           # DESeq2/edgeR 결과 CSV/TSV 파서 (Polars)
│   ├── rmats_loader.py         # rMATS 5대 이벤트 파일 병합 로더 (Polars)
│   └── merger.py               # DEG + AS 테이블 간의 고속 Join 및 사분면 라벨러
├── visualizer/
│   ├── __init__.py
│   ├── quadrant_plot.py        # Log2FC vs ΔPSI 4사분면 인터랙티브 산점도 (Plotly)
│   └── dual_volcano.py         # DEG Volcano와 AS Volcano를 나란히 배치한 듀얼 뷰
├── ai/
│   ├── __init__.py
│   └── gemini_evaluator.py     # google-genai SDK 기반 양적/질적 변화 종합 해석기
├── requirements.txt            # streamlit, polars, pyarrow, plotly, google-genai
└── test_data/                  # 검증용 Mock 데이터셋 (DEG 1개, rMATS 5개)
```

---

### 4. 데이터 통합 로직 및 사분면 매핑 알고리즘 (`core/merger.py`)

DEG와 Alternative Splicing 데이터는 유전자 심볼(`geneSymbol` 또는 `GeneID`)을 키로 결합됩니다.

#### (1) 조인 및 결측치 처리 규칙

* 하나의 유전자에 여러 개의 스플라이싱 이벤트(예: 엑손 2번 스킵, 5번 스킵)가 존재할 경우:
* $\vert{}\Delta\text{PSI}\vert{}$ 값이 가장 크거나 FDR이 가장 낮은 대표 이벤트를 유전자 단위의 대표 스플라이싱 이벤트로 선정.


* DEG에는 존재하나 AS 이벤트가 검출되지 않은 경우: $\Delta\text{PSI} = 0, \text{Event} = \text{'None'}$으로 처리.
* AS에는 존재하나 DEG에 없는 경우: $\text{Log}_2\text{FC} = 0, \text{FDR}_{\text{DEG}} = 1.0$으로 처리.

#### (2) 4사분면 자동 분류 기준 수식

* **유의미 임계치**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge \theta_{\text{DEG}}$ (기본값 $1.0$), $\vert{}\Delta\text{PSI}\vert{} \ge \theta_{\text{AS}}$ (기본값 $0.1$), $\text{FDR} \le 0.05$
* **분류 체계**:
* **Q1 (Dual Impact)**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge 1.0 \land \vert{}\Delta\text{PSI}\vert{} \ge 0.1$
* **Q2 (Splicing-Driven / Masked)**: $\vert{}\text{Log}_2\text{FC}\vert{} < 1.0 \land \vert{}\Delta\text{PSI}\vert{} \ge 0.1$ $\rightarrow$ **이 연구의 핵심 어필 유전자군**

* **Q4 (Abundance-Driven)**: $\vert{}\text{Log}_2\text{FC}\vert{} \ge 1.0 \land \vert{}\Delta\text{PSI}\vert{} < 0.1$
* **Q3 (Background / Static)**: 나머지 전체



---

### 5. UI/UX 화면 구성 설계 (`app.py`)

1. **사이드바 (Control Panel)**:
* **데이터 디렉토리 설정**:
* DEG 결과 파일 경로 (`deg_result.csv`)
* rMATS 결과 디렉토리 경로 (`rmats_output/`)


* **필터 슬라이더**:
* $\text{DEG FDR Cutoff}$ (0.001 ~ 0.1)
* $\text{AS FDR Cutoff}$ (0.001 ~ 0.1)
* $\vert{}\text{Log}_2\text{FC}\vert{} \text{ Cutoff}$ (0.5 ~ 3.0)
* $\vert{}\Delta\text{PSI}\vert{} \text{ Cutoff}$ (0.05 ~ 0.5)


* **API 입력창**: `Google Gemini API Key` (Password 필드 처리)


2. **메인 패널 - 상단 지표 (Summary KPI)**:
* 전체 분석 유전자 수 | Q1(동시 변화) 수 | **Q2(스플라이싱 전용 변화) 수** | Q4(발현량 전용 변화) 수


3. **메인 패널 - 메인 시각화 (Interactive Visualizer)**:
* **탭 1: Quadrant Cross-Plot (중앙 배치)**:
* X축: $\text{Log}_2\text{FC}$ (발현량 변화)
* Y축: $\Delta\text{PSI}$ (스플라이싱 비율 변화)
* 점 색상: 5대 이벤트 유형(SE, RI, A5SS, A3SS, MXE) 또는 사분면 그룹
* 호버 정보: 유전자명, 유전자 설명, $\text{Log}_2\text{FC}$, $\Delta\text{PSI}$, 각 FDR 값


* **탭 2: Dual Volcano Parallel View**:
* 좌측: DEG Volcano ($\text{Log}_2\text{FC}$ vs $-\log_{10}(\text{FDR}_{\text{DEG}})$)
* 우측: AS Volcano ($\Delta\text{PSI}$ vs $-\log_{10}(\text{FDR}_{\text{AS}})$)




4. **메인 패널 - 하단 세부 데이터 테이블 & AI 해석**:
* Q2 영역(발현량 변화 없이 스플라이싱만 극적으로 바뀐 유전자)을 기본 정렬하여 표시.
* 유전자 선택 체크박스 제공 $\rightarrow$ `[Gemini 심층 기전 분석]` 버튼 클릭 시 우측 AI 패널 활성화.



---

### 6. Gemini 기반 AI 에이전트 프롬프트 및 구조 (`ai/gemini_evaluator.py`)

선택된 유전자가 "왜 발현량 변화만으로는 설명될 수 없으며, 스플라이싱 패턴의 변화가 단백질 기능에 어떤 치명적 영향을 미치는지"를 분자생물학적 관점에서 해설하도록 설계합니다.

```python
SYSTEM_INSTRUCTION = """
당신은 전사체학(Transcriptomics) 및 분자유전학 전문 생물정보학 연구원입니다.
제공되는 데이터는 특정 조건(노화, 조직 손상, 스트레스)에서 산출된 DEG(발현량) 및 Alternative Splicing(형태 변화) 통합 분석 결과입니다.
단백질 발현 총량(Log2FC)과 엑손 스킵/인트론 보존(ΔPSI)의 변화를 결합하여, 
단순 발현량 분석만으로는 놓칠 수 있었던 '스플라이싱 주도형 기능 조절 메커니즘'을 전문적이고 설득력 있게 해석하십시오.
"""

USER_PROMPT_TEMPLATE = """
다음은 상위 분석 대상 유전자 정보입니다:
- Gene Symbol: {gene_symbol}
- Gene Description: {description}
- Log2 Fold Change (발현량 변화): {log2fc:.3f} (FDR: {deg_fdr:.2e})
- ΔPSI (스플라이싱 변위): {delta_psi:.3f} (FDR: {as_fdr:.2e})
- Splicing Event Type: {event_type} (좌표: {coordinates})

[질의 사항]
1. 이 유전자의 총 발현량 변화(Log2FC)와 스플라이싱 변화(ΔPSI)의 관계를 분석해 주십시오. (예: 발현량은 미미하나 기능적 도메인을 포함하는 엑손이 소실되었을 가능성 등)
2. 해당 스플라이싱 변이가 단백질 도메인, NMD(Nonsense-mediated decay), 또는 세포 운명/재생/노화 경로에 미칠 수 있는 분자생물학적 가설을 제시하십시오.
3. 이 결과를 실험적으로 증명하기 위해 Isoform-specific qRT-PCR 프라이머를 설계할 때 타깃해야 하는 구체적인 엑손 접합부 검증 전략을 제안하십시오.
"""
```

---

### 7. 즉시 개발 실행을 위한 `requirements.txt` 및 의존성

AntiGravity 터미널에서 아래 패키지를 설치하여 실행을 준비합니다.

```text
streamlit>=1.35.0
polars>=0.20.0
pyarrow>=15.0.0
plotly>=5.20.0
google-genai>=0.1.0
pandas>=2.0.0
openpyxl>=3.1.0
```

이 청사진은 rMATS의 5개 TSV 파일과 DESeq2의 단일 CSV 파일만 있으면 **복잡한 재계산 없이 1초 만에 4사분면 통합 그래프를 로드**하며, "발현 총량과 발현 양상의 동시 분석"이라는 연구 철학을 완벽히 시각화합니다.
