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

## 🚀 원클릭 환경 구축 (`install.sh`)

`bash install.sh` (또는 `./install`) 단 한 줄만 실행하면 **Conda/Mamba 감지, `environment.yml` 기반 가상환경 생성, 생물정보학 CLI 도구 설치 및 무결성 검증까지 원클릭으로 완료**됩니다.

### 1. 저장소 클론
```bash
git clone https://github.com/skt526/GenSplice-Agent.git
cd GenSplice-Agent
```

### 2. 원클릭 의존성 설치
```bash
bash install.sh
```
*(또는 `./install`)*

#### 💡 설치 프로세스 진행 단계 (터미널 실시간 표시)
1. **[1/5] 패키지 관리자 감지**: `mamba` 또는 `conda` 자동 체크 (미설치 시 python venv 대체 폴백)
2. **[2/5] Conda 환경 확인**: `gensplice-agent` 환경 존재 여부 확인
3. **[3/5] Conda 패키지 설치**: `environment.yml` 기반 Bioconda + Python 통합 패키지 설치
4. **[4/5] 생물정보학 CLI 검증**: `fastp`, `STAR`, `rmats.py`, `featureCounts`, `R` 실행 권한 및 검증
5. **[5/5] Python 파이프라인 검증**: `streamlit`, `polars`, `google.genai` 등 정상 파싱 로드 테스트

---

## 🖥 대시보드 실행

설치가 완료되면 아래 명령어로 Conda 환경을 활성화하고 대시보드를 구동합니다:

```bash
conda activate gensplice-agent
streamlit run app.py
```

---

## 📄 문서
- [BLUEPRINT.md](BLUEPRINT.md): GenSplice-Agent 시스템 설계 청사진 및 알고리즘 명세
- [environment.yml](environment.yml): Conda/Bioconda 명세서
