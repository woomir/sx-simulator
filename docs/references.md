### 📖 참고 문헌 (Literature References)

* **참조 프레임워크 및 문헌 앵커 (직접 구현된 full MSE solver를 뜻하지 않음)**:
  * Wang et al. (2002). "A thermodynamic framework for predicting the phase behavior of aqueous and mixed-solvent extraction systems."
  * ALTA 2024 Metallurgical Conference Materials (Mixer-Settler SX Simulation)
* **종분화 (Speciation) 및 수화 이온 반응 상수 ($K_{MOH}$, $K_{MSO_4}$)**:
  * Baes, C. F., & Mesmer, R. E. (1976). *The Hydrolysis of Cations*. Wiley.
  * Smith, R. M., & Martell, A. E. (1976). *Critical Stability Constants*. Plenum Press.
* **추출제 파라미터 ($\text{pH}_{50}$, $k$, $\alpha$, $E_{\max}$ 등)**:
  * **Cyanex 272 / D2EHPA**:
    * Mohapatra, M. et al. (2007). "Solvent extraction of heavy metals from aqueous solutions." *Hydrometallurgy*.
    * Pereira, D. D. et al. (2014). "Separation of nickel and cobalt from sulfate leach liquor." *Minerals Engineering*.
    * 내부 피팅 데이터 및 국내 배터리 재활용 센터 운전 기준 보정.
* **고로딩(High-Loading) 다핵 착물 형성 관련 보정 아이디어 (v2.0)**:
  * D2EHPA 및 Cyanex 272의 고농도 금속 추출 환경에서 보고된 다핵 올리고머(Multi-nuclear oligomeric complexes) 경향을 참고해, 현재 앱은 유효 화학양론 결합수($n_{\text{eff}}$)를 줄이는 경량 보정으로 반영합니다. 이는 exact oligomeric equilibrium solver를 의미하지 않습니다.
  * Rydberg, J. et al. (2004). *Solvent Extraction Principles and Practice*. Marcel Dekker.
  * Neuman, R. D., et al. (1990). "Formation of polynuclear complexes in acidic organophosphorus extraction systems".
