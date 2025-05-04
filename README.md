# ECBSA: Elliptic Curve Based Symmetric Algorithm

ECBSA는 타원 곡선 기반의 키 공유와 동적 S-box를 결합한 블록 암호 알고리즘입니다. 공개키 암호 방식을 통해 공유된 비밀로부터 S-box와 라운드 키를 동적으로 생성하여 대칭 암호화를 수행합니다.

## 🔐 개요

- **키 교환:** 송신자와 수신자가 각각 일회용 ECC 키 쌍을 생성하고 공개키를 교환하여 공유 비밀 \( S = k_A * k_B * G \)을 생성합니다.
- **S-box 생성:** 공유 비밀 \( x_S \)의 해시로부터 Fisher-Yates 셔플을 이용한 동적 S-box \( σ \)를 생성합니다.
- **키 스케줄:** 각 라운드 키는 \( d_r = SHA-256(x_S) + r mod p \)에서 생성된 점 \( d_r * G \)의 x좌표를 직렬화해 구성됩니다.
- **암호화/복호화:** S-box 변환과 XOR 연산을 반복하여 입력 블록을 암호화합니다.

## ⚙️ 필요 라이브러리

- [GMP](https://gmplib.org/) (GNU 다중 정밀도 산술 라이브러리)
- [OpenSSL](https://www.openssl.org/) (SHA-256 해시 및 난수 생성용)

## 🛠️ 빌드 방법

다음 명령어로 실행 파일을 컴파일할 수 있습니다:

```bash
gcc *.c -o ECBSA -lgmp -lcrypto
```

## ▶️ 실행 예시
```bash
./ECBSA
```