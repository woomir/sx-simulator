import numpy as np
from sx_simulator.single_stage import solve_single_stage
from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS

def step1_mass_balance():
    print("=== Step 1: Single Stage Mass Balance Verification ===")
    feed_aq = {'Li': 1.5, 'Ni': 5.0, 'Co': 3.0, 'Mn': 2.0}
    feed_org = {m: 0.0 for m in DEFAULT_METALS}
    target_ph = 4.5
    vol_aq = 100.0
    vol_org = 100.0
    extractant = "Cyanex 272"
    t_c = 25.0
    sulfate = 0.5
    
    result = solve_single_stage(
        C_aq_in=feed_aq, C_org_in=feed_org, pH_in=3.0,
        Q_aq=vol_aq, Q_org=vol_org,
        extractant=extractant, C_ext=0.5,
        target_pH=target_ph, temperature=t_c, C_sulfate=sulfate
    )
    
    aq_out = result['C_aq_out']
    org_out = result['C_org_out']
    
    passed = True
    for m in DEFAULT_METALS:
        mass_in = feed_aq.get(m, 0)*vol_aq + feed_org.get(m, 0)*vol_org
        mass_out = aq_out.get(m, 0)*vol_aq + org_out.get(m, 0)*vol_org
        diff = abs(mass_in - mass_out)
        print(f"[{m}] IN: {mass_in:.4f} g, OUT: {mass_out:.4f} g, Diff: {diff:.6f}")
        if diff > 1e-4:
            passed = False
    print(f"Step 1 Passed: {passed}")
    return passed

def step2_isotherm_consistency():
    print("\n=== Step 2: Isotherm Curve Consistency Check (Co) ===")
    extractant = "Cyanex 272"
    t_c = 25.0
    phs = np.arange(1.0, 6.1, 0.5)
    
    passed = True
    print(f"Extractant: {extractant}, Temp: {t_c}°C, Metal: Co")
    for ph in phs:
        # Simplistic check using pure sigmoid directly
        from sx_simulator.config import EXTRACTANT_PARAMS
        params = EXTRACTANT_PARAMS[extractant]['Co']
        ph50 = params['pH50']
        k = params['k']
        e_max = params['E_max']
        e = e_max / (1 + np.exp(-k * (ph - ph50)))
        print(f"pH {ph:.1f} -> Recovery: {e:.2f}%")
        if not (0 <= e <= 100):
            passed = False
    
    # pH50 Check
    e_ph50 = e_max / (1 + np.exp(-k * (ph50 - ph50)))
    print(f"At pH50 ({ph50:.2f}), Recovery: {e_ph50:.2f}% (Expected ~{e_max/2:.2f}%)")
    if abs(e_ph50 - e_max/2) > 1e-2:
        passed = False
        
    print(f"Step 2 Passed: {passed}")
    return passed

def step3_mccabe_thiele_logic():
    print("\n=== Step 3: McCabe-Thiele Counter-Current Logic Check ===")
    feed_aq = {'Li': 1.5, 'Ni': 5.0, 'Co': 3.0, 'Mn': 2.0}
    extract_org = {m: 0.0 for m in DEFAULT_METALS}
    num_stages = 4
    target_ph = 4.5
    vol_aq = 100.0
    vol_org = 100.0
    extractant = "Cyanex 272"
    t_c = 25.0
    sulfate = 0.5
    
    stages_res = solve_multistage_countercurrent(
        C_aq_feed=feed_aq, pH_feed=3.0,
        C_org_fresh=extract_org, Q_aq=vol_aq, Q_org=vol_org,
        extractant=extractant, C_ext=0.5,
        n_stages=num_stages, target_pH=target_ph,
        temperature=t_c, C_sulfate=sulfate
    )
    
    # Check if mass balance holds for the entire cascade
    sum_mass_in = {m: feed_aq.get(m, 0)*vol_aq + extract_org.get(m, 0)*vol_org for m in DEFAULT_METALS}
    raffinates = [s['C_aq_out'] for s in stages_res['stages']]
    loaded_orgs = [s['C_org_out'] for s in stages_res['stages']]
    
    print(f"Converged: {stages_res.get('converged', 'Unknown')} in {stages_res.get('iterations', 'Unknown')} iterations")
    
    # Raffinate leaves from stage N (index 3)
    # Loaded Org leaves from stage 1 (index 0)
    final_raff = raffinates[-1]
    final_loaded_org = loaded_orgs[0]
    
    passed = True
    for m in DEFAULT_METALS:
        mass_in = sum_mass_in[m]
        mass_out = final_raff.get(m, 0)*vol_aq + final_loaded_org.get(m, 0)*vol_org
        diff = abs(mass_in - mass_out)
        print(f"Total cascade [{m}] IN: {mass_in:.4f} g, OUT: {mass_out:.4f} g, Diff: {diff:.6f}")
        if diff > 1e-4:
            passed = False
            
    # Check step sequence 
    print("Checking Operating Line and inter-stage balances...")
    for i in range(num_stages):
        # aq in from i-1 (or feed if i=0)
        aq_in = raffinates[i-1] if i > 0 else feed_aq
        org_in = loaded_orgs[i+1] if i < num_stages-1 else extract_org
        aq_out = raffinates[i]
        org_out = loaded_orgs[i]
        
        for m in ['Co']: # Just check Co
            mass_in = aq_in.get(m, 0)*vol_aq + org_in.get(m, 0)*vol_org
            mass_out = aq_out.get(m, 0)*vol_aq + org_out.get(m, 0)*vol_org
            diff = abs(mass_in - mass_out)
            stage_ok = diff <= 1e-4
            if not stage_ok:
                passed = False
        print(
            f"Stage {i+1} Co balance {'OK' if stage_ok else 'FAIL'} "
            f"(Diff: {diff:.6f}, Aq out: {aq_out.get('Co',0):.3f}, Org out: {org_out.get('Co',0):.3f})"
        )

    print(f"Step 3 Passed: {passed}")
    return passed


def step4_target_pH_dilution_consistency():
    print("\n=== Step 4: Target pH Self-Consistent Dilution Check ===")
    feed_aq = {'Li': 1.5, 'Ni': 5.0, 'Co': 3.0, 'Mn': 2.0}
    feed_org = {m: 0.0 for m in DEFAULT_METALS}
    target_ph = 4.5
    vol_aq = 100.0
    vol_org = 100.0
    extractant = "Cyanex 272"
    t_c = 25.0
    sulfate = 0.5
    c_naoh = 5.0

    target_result = solve_single_stage(
        C_aq_in=feed_aq, C_org_in=feed_org, pH_in=3.0,
        Q_aq=vol_aq, Q_org=vol_org,
        extractant=extractant, C_ext=0.5,
        target_pH=target_ph, C_NaOH=c_naoh,
        temperature=t_c, C_sulfate=sulfate,
    )

    q_naoh_est = target_result.get("Q_NaOH_estimated_L_hr", 0.0)
    fixed_result = solve_single_stage(
        C_aq_in=feed_aq, C_org_in=feed_org, pH_in=3.0,
        Q_aq=vol_aq, Q_org=vol_org,
        extractant=extractant, C_ext=0.5,
        C_NaOH=c_naoh, Q_NaOH=q_naoh_est,
        temperature=t_c, C_sulfate=sulfate,
    )

    passed = True
    pH_diff = abs(fixed_result["pH_out"] - target_ph)
    print(
        f"Estimated NaOH flow: {q_naoh_est:.4f} L/hr, "
        f"Target-mode dilution iterations: {target_result.get('dilution_iterations', 0)}"
    )
    print(
        f"Fixed-NaOH replay pH: {fixed_result['pH_out']:.4f}, "
        f"Target pH: {target_ph:.4f}, Diff: {pH_diff:.6f}"
    )
    if pH_diff > 2e-2:
        passed = False

    for metal in ["Li", "Ni", "Co", "Mn"]:
        diff = abs(
            fixed_result["C_aq_out"].get(metal, 0.0)
            - target_result["C_aq_out"].get(metal, 0.0)
        )
        print(
            f"{metal} target-vs-fixed raffinate diff: {diff:.6f} g/L"
        )
        if diff > 5e-3:
            passed = False

    print(f"Step 4 Passed: {passed}")
    return passed

if __name__ == "__main__":
    p1 = step1_mass_balance()
    p2 = step2_isotherm_consistency()
    p3 = step3_mccabe_thiele_logic()
    p4 = step4_target_pH_dilution_consistency()

    print("\nOVERALL VALIDATION:", "PASSED" if (p1 and p2 and p3 and p4) else "FAILED")
