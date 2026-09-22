# GenSplice-Agent 🧬

> **통합 전사체 발현량(DEG) 및 대립적 스플라이싱(Alternative Splicing) 시각화·AI 해석 대시보드**

`GenSplice-Agent`는 유전자의 발현량 변화(Quantity, DEG)와 스플라이싱 형태의 질적 변화(Quality, Alternative Splicing)를 단일 사분면 좌표계로 통합하고, Google Gemini AI를 통해 분자생물학적 기전을 해설하는 전사체학 전문 통합 대시보드입니다.

---

## 📦 예상 의존성 목록 (Required Dependencies)

### 1. 시스템 요구사항
- **Python**: 3.9 이상 (Python 3.10+ 권장)
- **Git**: 2.0 이상

### 2. 핵심 파이썬 패키지 (`requirements.txt`)
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

## 🚀 원클릭 설치 및 실행 방법 (One-Click Installation)

저장소를 클론한 후 `./install` (또는 `bash install.sh`) 명령어를 실행하면, **가상 환경 생성부터 의존성 설치 및 무결성 검증까지 전 과정을 시각적 진행 상태(Progress)와 함께 원클릭으로 완료**할 수 있습니다.

### 1. 저장소 클론
```bash
git clone https://github.com/skt526/GenSplice-Agent.git
cd GenSplice-Agent
```

### 2. 원클릭 의존성 설치
```bash
./install
```
*(또는 `bash install.sh`)*

#### 💡 설치 프로세스 진행 단계 (터미널 실시간 표시)
1. **[1/5] Python 3 환경 확인**: Python 3.9+ 및 경로 체크
2. **[2/5] 가상 환경 생성**: `.venv` 파이썬 가상 환경 자동 생성
3. **[3/5] 패키지 관리자 업그레이드**: `pip`, `setuptools`, `wheel` 최신화
4. **[4/5] 의존성 설치**: `requirements.txt` 패키지 자동 설치 (프로그레스 바 출력)
5. **[5/5] 무결성 검증**: 핵심 모듈(`polars`, `streamlit`, `google.genai` 등) 정상 로드 테스트

---

## 🖥 대시보드 실행

설치가 완료된 후 아래 명령어로 대시보드를 구동합니다:

```bash
source .venv/bin/activate
streamlit run app.py
```

---

## 📄 문서
- [BLUEPRINT.md](BLUEPRINT.md): GenSplice-Agent 시스템 설계 청사진 및 알골리즘 명세
