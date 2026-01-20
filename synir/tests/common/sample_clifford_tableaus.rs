use bitvec::bitvec;
use bitvec::prelude::Lsb0;
use synir::{
    data_structures::{CliffordTableau, PauliString, PropagateClifford}
};

pub fn setup_sample_ct() -> CliffordTableau {
    // Stab: ZZZ, -YIY, XIX
    // Destab: -IXI, XXI, IYY
    // qubit 1x: ZYI
    // qubit 1z: IZZ
    let pauli_1 = PauliString::from_text("ZYIIZZ");

    // qubit 2x: ZIX
    // qubit 2z: XII
    let pauli_2 = PauliString::from_text("ZIXXII");

    // qubit 3x: ZYY
    // qubit 3z: IIZ
    let pauli_3 = PauliString::from_text("ZYYIIZ");

    let signs = bitvec![0, 1, 0, 1, 0, 0];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3], signs)
}

pub fn setup_sample_inverse_ct() -> CliffordTableau {
    // Stab: -ZIYZ, -ZZYZ, -XZXI, IZXX
    // Destab: -YYIZ, -YYXZ, ZIXX, -XZXZ
    // qubit 1x: ZZXI
    // qubit 1z: YYZX
    let pauli_1 = PauliString::from_text("ZZXIYYZX");

    // qubit 2x: IZZZ
    // qubit 2z: YYIZ
    let pauli_2 = PauliString::from_text("IZZZYYIZ");

    // qubit 3x: YYXX
    // qubit 3z: IXXX
    let pauli_3 = PauliString::from_text("YYXXIXXX");

    // qubit 3x: ZZIX
    // qubit 3z: ZZXZ
    let pauli_4 = PauliString::from_text("ZZIXZZXZ");

    let signs = bitvec![1, 1, 1, 0, 1, 1, 0, 1];
    CliffordTableau::from_parts(vec![pauli_1, pauli_2, pauli_3, pauli_4], signs)
}

pub fn identity_2qb_ct() -> CliffordTableau {
    return CliffordTableau::new(2);
}

pub fn sample_swap_ct() -> CliffordTableau {
    let mut clifford_tableau = CliffordTableau::new(2);

    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 0);
    clifford_tableau.cx(0, 1);
    return clifford_tableau;
}

pub fn half_swap_0_1() -> CliffordTableau {
    let mut clifford_tableau = CliffordTableau::new(2);
    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 0);
    return clifford_tableau;
}

pub fn half_swap_1_0() -> CliffordTableau {
    let mut clifford_tableau = CliffordTableau::new(2);
    clifford_tableau.cx(1, 0);
    clifford_tableau.cx(0, 1);
    return clifford_tableau;
}

pub fn sample_2cnot_ladder() -> CliffordTableau {
    let mut clifford_tableau = CliffordTableau::new(3);
    clifford_tableau.cx(0, 1);
    clifford_tableau.cx(1, 2);
    return clifford_tableau;
}

pub fn sample_s_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(1);
    ct.s(0);
    return ct;
}

pub fn sample_s_dgr_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(1);
    ct.s_dgr(0);
    return ct;
}

pub fn sample_v_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(1);
    ct.v(0);
    return ct;
}

pub fn sample_v_dgr_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(1);
    ct.v_dgr(0);
    return ct;
}

pub fn sample_cnot_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    ct.cx(0, 1);
    return ct;
}

pub fn sample_cnot_reverse_gate() -> CliffordTableau {
    let mut ct = CliffordTableau::new(2);
    ct.cx(1, 0);
    return ct;
}
