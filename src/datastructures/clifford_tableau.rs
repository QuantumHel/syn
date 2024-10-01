use bitvec::{prelude::BitVec, slice::BitSlice};
use std::iter;

use super::{
    pauli_string::{cx, PauliString},
    IndexType, PropagateClifford,
};

#[derive(PartialEq, Eq, Debug)]
pub struct CliffordTableau {
    // We keep track of the pauli letters per qubit not per stabilizer
    pauli_columns: Vec<PauliString>,
    signs: BitVec,
    size: usize, // https://quantumcomputing.stackexchange.com/questions/28740/tracking-the-signs-of-the-inverse-tableau
}

impl CliffordTableau {
    /// Constructs a Clifford Tableau of `n` qubits initialized to the identity operation
    pub fn new(n: usize) -> Self {
        CliffordTableau {
            pauli_columns: { (0..n).map(|i| PauliString::from_basis_int(i, n)).collect() },
            signs: BitVec::from_iter(iter::repeat(false).take(2 * n)),
            size: n,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub(crate) fn z_signs(&self) -> &BitSlice {
        let n = self.size();
        &self.signs[0..n]
    }

    pub(crate) fn x_signs(&self) -> &BitSlice {
        let n = self.size();
        &self.signs[n..]
    }
}

impl PropagateClifford for CliffordTableau {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let n = self.size();
        let [control, target] = self.pauli_columns.get_many_mut([control, target]).unwrap();
        self.signs ^= (control.x.clone() & target.z.clone())
            & (control.z.clone() ^ target.x.clone() ^ BitVec::repeat(true, 2 * n));
        cx(control, target);
        self
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.pauli_columns.get_mut(target).unwrap();
        // Verified: SXS^dag = Y
        //           SYS^dag = -X
        //           SZS^dag = Z
        self.signs ^= chains_target.y_bitmask();
        // Defined for Phase gate in https://arxiv.org/pdf/quant-ph/0406196
        chains_target.s();
        self
    }

    fn v(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.pauli_columns.get_mut(target).unwrap();
        // Verified: VXV^dag = X
        //           VYV^dag = Z
        //           VZV^dag = -Y
        chains_target.v();
        self.signs ^= chains_target.y_bitmask();
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::bitvec;
    use bitvec::prelude::Lsb0;

    #[test]
    fn test_clifford_tableau_constructor() {
        let ct_size = 3;
        let ct = CliffordTableau::new(ct_size);
        let x_1 = bitvec![1, 0, 0, 0, 0, 0];
        let z_1 = bitvec![0, 0, 0, 1, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        let x_2 = bitvec![0, 1, 0, 0, 0, 0];
        let z_2 = bitvec![0, 0, 0, 0, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };
        let x_3 = bitvec![0, 0, 1, 0, 0, 0];
        let z_3 = bitvec![0, 0, 0, 0, 0, 1];
        let pauli_3 = PauliString { x: x_3, z: z_3 };
        let signs = bitvec![0, 0, 0, 0, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3],
            signs,
            size: ct_size,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    fn setup_sample_ct() -> CliffordTableau {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        // qubit 1x: ZYX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3 = PauliString { x: x_3, z: z_3 };

        let signs = bitvec![0, 0, 0, 0, 0, 0];
        CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3],
            signs,
            size: ct_size,
        }
    }

    #[test]
    fn test_clifford_tableau_s() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S on qubit 0
        ct.s(0);

        // Stab: ZZZ, (-X)IY, YIX
        // Destab: IXI, YXI, IYY

        // qubit 1x: ZXY
        // qubit 1z: IYI
        let z_1 = bitvec![1, 0, 1, 0, 1, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![0, 1, 0, 0, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_v() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.v(0);

        // Stab: (-Y)ZZ, ZIY, XIX
        // Destab: IXI, XXI, IYY

        // qubit 1x: YZX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![1, 0, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![1, 0, 0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_sdag() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S_dgr on qubit 0, 1
        ct.s_dgr(0);
        ct.s_dgr(1);

        // Stab: ZZZ, XIY, (-Y)IX
        // Destab: I(-Y)I, (-Y)(-Y)I, IXY

        // qubit 1x: ZXY
        // qubit 1z: IYI
        let z_1 = bitvec![1, 0, 1, 0, 1, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: YYX
        let z_2 = bitvec![1, 0, 0, 1, 1, 0];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![0, 0, 1, 1, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_vdag() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Vdgr on qubit 0, 1
        ct.v_dgr(0);
        ct.v_dgr(1);

        // Stab: YYZ, (-Z)IY, XIX
        // Destab: IXI, XXI, I(-Z)Y

        // qubit 1x: YZX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![1, 0, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: YII
        // qubit 2z: XXZ
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![1, 0, 0, 1, 1, 0];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![0, 1, 0, 0, 0, 1];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_h() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.h(0);
        ct.h(1);

        // Stab: XXZ, (-Y)IY, ZIX
        // Destab: IZI, ZZI, I(-Y)Y

        // qubit 1x: XYZ
        // qubit 1z: IZI
        let z_1 = bitvec![0, 1, 1, 0, 1, 0];
        let x_1 = bitvec![1, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: XII
        // qubit 2z: ZZY
        let z_2 = bitvec![0, 0, 0, 1, 1, 1];
        let x_2 = bitvec![1, 0, 0, 0, 0, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![0, 1, 0, 0, 0, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_x() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.x(0);

        // Stab: (-Z)ZZ, (-Y)IY, XIX
        // Destab: IXI, XXI, IYY

        // qubit 1x: ZYX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![1, 1, 0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_y() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Y on qubit 0, 1
        ct.y(0);

        // Stab: (-Z)ZZ, YIY, (-X)IX
        // Destab: IXI, (-X)XI, IYY

        // qubit 1x: ZYX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![1, 0, 1, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_z() {
        // Stab: ZZZ, YIY, XIX
        // Destab: IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Z on qubit 0
        ct.z(0);

        // Stab: ZZZ, (-Y)IY, (-X)IX
        // Destab: IXI, (-X)XI, IYY

        // qubit 1x: ZYX
        // qubit 1z: IXI
        let z_1 = bitvec![1, 1, 0, 0, 0, 0];
        let x_1 = bitvec![0, 1, 1, 0, 1, 0];
        let pauli_1_ref = PauliString { x: x_1, z: z_1 };

        // qubit 2x: ZII
        // qubit 2z: XXY
        let z_2 = bitvec![1, 0, 0, 0, 0, 1];
        let x_2 = bitvec![0, 0, 0, 1, 1, 1];
        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        // qubit 3x: ZYX
        // qubit 3z: IIY
        let z_3 = bitvec![1, 1, 0, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 1];
        let pauli_3_ref = PauliString { x: x_3, z: z_3 };

        let signs_ref = bitvec![0, 1, 1, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    /// This does not generate a valid Clifford Tableau. Only used to check commutation relations
    fn setup_sample_two_qubit_ct(pauli_control: &str) -> CliffordTableau {
        let pauli_1_ref = match pauli_control {
            "i" => {
                // qubit 1x: II
                // qubit 1z: II
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![0, 0, 0, 0];

                PauliString { x: x_1, z: z_1 }
            }
            "x" => {
                // qubit 1x: XX
                // qubit 1z: XX
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString { x: x_1, z: z_1 }
            }
            "y" => {
                // qubit 1x: YY
                // qubit 1z: YY
                let z_1 = bitvec![1, 1, 1, 1];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString { x: x_1, z: z_1 }
            }
            "z" => {
                // qubit 1x: ZZ
                // qubit 1z: ZZ
                let z_1 = bitvec![1, 1, 1, 1];
                let x_1 = bitvec![0, 0, 0, 0];

                PauliString { x: x_1, z: z_1 }
            }
            _ => panic!("Pauli letter not recognized"),
        };

        // qubit 1x: IX
        // qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];

        let pauli_2_ref = PauliString { x: x_2, z: z_2 };

        let signs_ref = bitvec![0, 0, 0, 0];
        CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref],
            signs: signs_ref,
            size: 2,
        }
    }

    #[test]
    fn test_clifford_tableau_cx_i() {
        // Stab: II, IX
        // Destab: IY, IZ
        let mut ct = setup_sample_two_qubit_ct("i");

        // Apply CX to 0 -> 1.
        ct.cx(0, 1);

        // Stab: II, IX
        // Destab: ZY, ZZ

        //qubit 1x: II
        //qubit 1z: ZZ
        let z_1 = bitvec![0, 0, 1, 1];
        let x_1 = bitvec![0, 0, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 1x: IX
        //qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cx_x() {
        // Stab: XI, XX
        // Destab: XY, XZ
        let mut ct = setup_sample_two_qubit_ct("x");
        // Apply CX to 0 -> 1.
        ct.cx(0, 1);

        // Stab: XX, XI
        // Destab: YZ, -YY

        //qubit 1x: XX
        //qubit 1z: YY
        let z_1 = bitvec![0, 0, 1, 1];
        let x_1 = bitvec![1, 1, 1, 1];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 2x: XI
        //qubit 2z: ZY
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![1, 0, 0, 1];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 0, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cx_y() {
        // Stab: YI, YX
        // Destab: YY, YZ
        let mut ct = setup_sample_two_qubit_ct("y");
        // Apply CX to 0 -> 1.
        ct.cx(0, 1);

        // Stab: YX, YI
        // Destab: -XZ, XY

        //qubit 1x: YY
        //qubit 1z: XX
        let z_1 = bitvec![1, 1, 0, 0];
        let x_1 = bitvec![1, 1, 1, 1];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 2x: XI
        //qubit 2z: ZY
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![1, 0, 0, 1];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cx_z() {
        // Stab: ZI, ZX
        // Destab: ZY, ZZ
        let mut ct = setup_sample_two_qubit_ct("z");
        // Apply CX to 0 -> 1.
        ct.cx(0, 1);

        // Stab: ZI, ZX
        // Destab: IY, IZ

        //qubit 1x: ZZ
        //qubit 1z: II
        let z_1 = bitvec![1, 1, 0, 0];
        let x_1 = bitvec![0, 0, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 2x: IX
        //qubit 2z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cz_i() {
        // Stab: II, ZX
        // Destab: ZY, IZ
        let mut ct = setup_sample_two_qubit_ct("i");

        // Apply CX to 0 -> 1.
        ct.cz(0, 1);

        // Stab: II, ZX
        // Destab: ZY, IZ

        //qubit 1x: IZ
        //qubit 1z: ZI
        let z_1 = bitvec![0, 1, 1, 0];
        let x_1 = bitvec![0, 0, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 1x: IX
        //qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cz_x() {
        // Stab: XI, XX
        // Destab: XY, XZ
        let mut ct = setup_sample_two_qubit_ct("x");
        // Apply CX to 0 -> 1.
        ct.cz(0, 1);

        // Stab: XZ, YY
        // Destab: -YX, XI

        //qubit 1x: XY
        //qubit 1z: YX
        let z_1 = bitvec![0, 1, 1, 0];
        let x_1 = bitvec![1, 1, 1, 1];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 1x: ZY
        //qubit 1z: XI
        let z_2 = bitvec![1, 1, 0, 0];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cz_y() {
        // Stab: YI, YX
        // Destab: YY, YZ
        let mut ct = setup_sample_two_qubit_ct("y");
        // Apply CX to 0 -> 1.
        ct.cz(0, 1);

        // Stab: YZ, -XY
        // Destab: XX, YI

        //qubit 1x: YX
        //qubit 1z: XY
        let z_1 = bitvec![1, 0, 0, 1];
        let x_1 = bitvec![1, 1, 1, 1];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 2x: ZY
        //qubit 2z: XI
        let z_2 = bitvec![1, 1, 0, 0];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_cz_z() {
        // Stab: ZI, ZX
        // Destab: ZY, ZZ
        let mut ct = setup_sample_two_qubit_ct("z");
        // Apply CX to 0 -> 1.
        ct.cz(0, 1);

        // Stab: ZI, IX
        // Destab: IY, ZZ

        //qubit 1x: ZI
        //qubit 1z: IZ
        let z_1 = bitvec![1, 0, 0, 1];
        let x_1 = bitvec![0, 0, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };
        //qubit 2x: IX
        //qubit 2z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(clifford_tableau_ref, ct);
    }
}
