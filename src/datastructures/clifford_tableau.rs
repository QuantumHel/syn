use bitvec::{prelude::BitVec, slice::BitSlice};
use itertools::{izip, Itertools};
use std::fmt;
use std::iter::{self, zip};
use std::ops::Mul;

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

    pub(crate) fn compose(&self, rhs: &Self) -> Self {
        rhs.prepend(self)
    }

    /// Implements algorithms from https://doi.org/10.22331/q-2022-06-13-734 and Qiskit Clifford implementation
    pub(crate) fn prepend(&self, lhs: &Self) -> Self {
        let size = self.size();
        let mut pauli_columns = vec![PauliString::from_text(&"I".repeat(2 * size)); size];

        // Matrix-multiplication for M(rhs o self) = M(self) * M(rhs) as this is a row-permutation.
        // Loop re-order to be (k, i, j) as j ordering is contiguous.
        for (k, rhs_pauli_column) in self.pauli_columns.iter().enumerate() {
            for i in 0..size {
                pauli_columns[k].x ^=
                    BitVec::repeat(rhs_pauli_column.x[i], 2 * size) & &lhs.pauli_columns[i].x;
                pauli_columns[k].x ^= BitVec::repeat(rhs_pauli_column.x[i + size], 2 * size)
                    & &lhs.pauli_columns[i].z;
                pauli_columns[k].z ^=
                    BitVec::repeat(rhs_pauli_column.z[i], 2 * size) & &lhs.pauli_columns[i].x;
                pauli_columns[k].z ^= BitVec::repeat(rhs_pauli_column.z[i + size], 2 * size)
                    & &lhs.pauli_columns[i].z;
            }
        }

        let mut i_factors = vec![0_usize; 2 * size];
        // Keep track of the inherent i factors of left-hand tableau (where there are Y's in tableau rows)
        for lhs_pauli_column in lhs.pauli_columns.iter() {
            let local_sign = lhs_pauli_column.x.clone() & &lhs_pauli_column.z;
            for (fact, sign) in zip(i_factors.iter_mut(), local_sign) {
                *fact += sign as usize;
            }
        }

        // Accumulate the i factors when lhs basis is aggregated per rows in rhs tableau.
        // Indices reflect a (i,j) × (j, k) matrix multiplication.
        // Loop re-order to be (i, k, j).
        for (i, i_factor) in i_factors.iter_mut().enumerate() {
            for rhs_pauli_column in self.pauli_columns.iter() {
                let mut x1_select = Vec::new();
                let mut z1_select = Vec::new();
                for (j, lhs_pauli_column) in lhs.pauli_columns.iter().enumerate() {
                    if lhs_pauli_column.x[i] {
                        x1_select.push(rhs_pauli_column.x[j]);
                        z1_select.push(rhs_pauli_column.z[j])
                    }
                    if lhs_pauli_column.z[i] {
                        x1_select.push(rhs_pauli_column.x[j + size]);
                        z1_select.push(rhs_pauli_column.z[j + size])
                    }
                }
                let x1_accumulator = x1_select
                    .iter()
                    .scan(false, |state, x| {
                        *state ^= x;
                        Some(*state)
                    })
                    .collect_vec();

                let z1_accumulator = z1_select
                    .iter()
                    .scan(false, |state, x| {
                        *state ^= x;
                        Some(*state)
                    })
                    .collect_vec();

                let indexer = izip!(
                    x1_select.iter().skip(1),
                    z1_select.iter().skip(1),
                    x1_accumulator.iter(),
                    z1_accumulator.iter()
                )
                .map(lookup)
                .sum::<usize>();
                *i_factor += indexer;
            }
        }

        let mut new_signs = BitVec::repeat(false, 2 * size);

        // Contribution of combination of signs in rhs basis.
        // Calculate matrix vector M(lhs) * sign(rhs)
        for (j, lhs_pauli_column) in lhs.pauli_columns.iter().enumerate() {
            new_signs ^= BitVec::repeat(self.signs[j], 2 * size) & &lhs_pauli_column.x;
            new_signs ^= BitVec::repeat(self.signs[j + size], 2 * size) & &lhs_pauli_column.z;
        }

        // Get rid of `i` factors and convert to sign flips
        let p = i_factors
            .iter()
            .map(|sign| ((sign % 4) / 2) > 0)
            .collect::<BitVec>();

        new_signs ^= p;
        new_signs ^= lhs.signs.clone();

        CliffordTableau {
            pauli_columns,
            signs: new_signs,
            size,
        }
    }

    pub fn adjoint(&self) -> Self {
        // Algorithm taken from https://algassert.com/post/2002
        let size = self.size();
        // Create new CliffordTableau entries

        let mut new_columns = vec![PauliString::from_text(&"I".repeat(2 * size)); size];
        (0..size).for_each(|i| {
            for (j, pauli_column) in self.pauli_columns.iter().enumerate() {
                let (x1, z1, x2, z2) = match (
                    *pauli_column.x.get(i).unwrap(),
                    *pauli_column.z.get(i).unwrap(),
                    *pauli_column.x.get(i + size).unwrap(),
                    *pauli_column.z.get(i + size).unwrap(),
                ) {
                    // II -> II
                    (false, false, false, false) => (false, false, false, false),
                    // IX -> IX
                    (false, false, true, false) => (false, false, true, false),
                    // IY -> XX
                    (false, false, true, true) => (true, false, true, false),
                    // IZ -> XI
                    (false, false, false, true) => (true, false, false, false),
                    // XI -> IZ
                    (true, false, false, false) => (false, false, false, true),
                    // XX -> IY
                    (true, false, true, false) => (false, false, true, true),
                    // XY -> XY
                    (true, false, true, true) => (true, false, true, true),
                    // XZ -> XZ
                    (true, false, false, true) => (true, false, false, true),
                    // YI -> ZZ
                    (true, true, false, false) => (false, true, false, true),
                    // YX -> ZY
                    (true, true, true, false) => (false, true, true, true),
                    // YY -> YY
                    (true, true, true, true) => (true, true, true, true),
                    // YZ -> YZ
                    (true, true, false, true) => (true, true, false, true),
                    // ZI -> ZI
                    (false, true, false, false) => (false, true, true, false),
                    // ZX -> ZX
                    (false, true, true, false) => (false, true, true, false),
                    // ZY -> YX
                    (false, true, true, true) => (true, true, true, false),
                    // ZZ -> YI
                    (false, true, false, true) => (true, true, false, false),
                };

                new_columns[i].x.replace(j, x1);
                new_columns[i].z.replace(j, z1);
                new_columns[i].x.replace(j + size, x2);
                new_columns[i].z.replace(j + size, z2);
            }
        });
        let mut adjoint_table = CliffordTableau {
            pauli_columns: new_columns,
            signs: BitVec::repeat(false, 2 * size),
            size,
        };

        adjoint_table.signs ^= (adjoint_table.compose(self)).signs;
        adjoint_table
    }
}

impl PropagateClifford for CliffordTableau {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let n = self.size();
        let [control, target] = self.pauli_columns.get_many_mut([control, target]).unwrap();
        let mut scratch = BitVec::repeat(true, 2 * n);
        scratch ^= &target.x;
        scratch ^= &control.z;
        scratch &= &control.x;
        scratch &= &target.z;
        self.signs ^= scratch;

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

/// Lookup table that determines additional `i` factors when Pauli matrices are multiplied
fn lookup(accum: (&bool, &bool, &bool, &bool)) -> usize {
    match accum {
        (true, false, true, true) | (true, true, false, true) | (false, true, true, false) => 3,
        (true, true, true, false) | (false, true, true, true) | (true, false, false, true) => 1,
        _ => 0,
    }
}

impl Mul for CliffordTableau {
    type Output = Self;

    fn mul(self, lhs: Self) -> Self {
        self.prepend(&lhs)
    }
}

impl fmt::Display for CliffordTableau {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "CliffordTableau({})", self.size())?;
        for pauli_column in self.pauli_columns.iter() {
            writeln!(f, "{}", pauli_column)?;
        }
        let mut sign_str = String::new();
        for bit in self.signs.iter() {
            match *bit {
                true => sign_str.push('-'),
                false => sign_str.push('+'),
            }
            sign_str.push(' ')
        }
        sign_str.pop();
        write!(f, "{}", sign_str)
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
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
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

        let signs = bitvec![0, 1, 0, 1, 0, 0];
        CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3],
            signs,
            size: ct_size,
        }
    }

    #[test]
    fn test_clifford_tableau_s() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S on qubit 0
        ct.s(0);

        // Stab: ZZZ, -(-X)IY, YIX
        // Destab: -IXI, YXI, IYY

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

        let signs_ref = bitvec![0, 0, 0, 1, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_v() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.v(0);

        // Stab: (-Y)ZZ, -ZIY, XIX
        // Destab: -IXI, XXI, IYY

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

        let signs_ref = bitvec![1, 1, 0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_sdag() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S_dgr on qubit 0, 1
        ct.s_dgr(0);
        ct.s_dgr(1);

        // Stab: ZZZ, -XIY, (-Y)IX
        // Destab: -I(-Y)I, (-Y)(-Y)I, IXY

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

        let signs_ref = bitvec![0, 1, 1, 0, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_vdag() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Vdgr on qubit 0, 1
        ct.v_dgr(0);
        ct.v_dgr(1);

        // Stab: YYZ, -(-Z)IY, XIX
        // Destab: -IXI, XXI, I(-Z)Y

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

        let signs_ref = bitvec![0, 0, 0, 1, 0, 1];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_h() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.h(0);
        ct.h(1);

        // Stab: XXZ, -(-Y)IY, ZIX
        // Destab: -IZI, ZZI, I(-Y)Y

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

        let signs_ref = bitvec![0, 0, 0, 1, 0, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_x() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.x(0);

        // Stab: (-Z)ZZ, -(-Y)IY, XIX
        // Destab: -IXI, XXI, IYY

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

        let signs_ref = bitvec![1, 0, 0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_y() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Y on qubit 0, 1
        ct.y(0);

        // Stab: (-Z)ZZ, -YIY, (-X)IX
        // Destab: -IXI, (-X)XI, IYY

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

        let signs_ref = bitvec![1, 1, 1, 1, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    #[test]
    fn test_clifford_tableau_z() {
        // Stab: ZZZ, -YIY, XIX
        // Destab: -IXI, XXI, IYY
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Z on qubit 0
        ct.z(0);

        // Stab: ZZZ, -(-Y)IY, (-X)IX
        // Destab: -IXI, (-X)XI, IYY

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

        let signs_ref = bitvec![0, 0, 1, 1, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(clifford_tableau_ref, ct);
    }

    /// This does not generate a valid Clifford Tableau. Only used to check commutation relations
    fn setup_sample_two_qubit_ct(pauli_control: char) -> CliffordTableau {
        let pauli_1_ref = match pauli_control {
            'i' => {
                // qubit 1x: II
                // qubit 1z: II
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![0, 0, 0, 0];

                PauliString { x: x_1, z: z_1 }
            }
            'x' => {
                // qubit 1x: XX
                // qubit 1z: XX
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString { x: x_1, z: z_1 }
            }
            'y' => {
                // qubit 1x: YY
                // qubit 1z: YY
                let z_1 = bitvec![1, 1, 1, 1];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString { x: x_1, z: z_1 }
            }
            'z' => {
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
        let mut ct = setup_sample_two_qubit_ct('i');

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
        let mut ct = setup_sample_two_qubit_ct('x');
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
        let mut ct = setup_sample_two_qubit_ct('y');
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
        let mut ct = setup_sample_two_qubit_ct('z');
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
        // Stab: II, IX
        // Destab: IY, IZ
        let mut ct = setup_sample_two_qubit_ct('i');

        // Apply CZ to 0 -> 1.
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
        let mut ct = setup_sample_two_qubit_ct('x');
        // Apply CZ to 0 -> 1.
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
        let mut ct = setup_sample_two_qubit_ct('y');
        // Apply CZ to 0 -> 1.
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
        let mut ct = setup_sample_two_qubit_ct('z');
        // Apply CZ to 0 -> 1.
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

    #[test]
    fn test_clifford_tableau_compose() {
        let mut first_ct = setup_sample_ct();
        first_ct.x(0);
        first_ct.h(0);
        first_ct.cx(0, 1);

        let mut second_ct = CliffordTableau::new(3);
        second_ct.s(1);
        second_ct.h(1);
        second_ct.cx(1, 0);

        let third = first_ct.compose(&second_ct);

        let mut ref_ct = setup_sample_ct();
        ref_ct.x(0);
        ref_ct.h(0);
        ref_ct.cx(0, 1);

        ref_ct.s(1);
        ref_ct.h(1);
        ref_ct.cx(1, 0);

        assert_eq!(ref_ct, third);
    }

    #[test]
    fn test_clifford_tableau_prepend() {
        let mut first_ct = setup_sample_ct();
        first_ct.x(0);
        first_ct.h(0);
        first_ct.cx(0, 1);

        let mut second_ct = CliffordTableau::new(3);
        second_ct.s(1);
        second_ct.h(1);
        second_ct.cx(1, 0);

        let third = second_ct.prepend(&first_ct);

        let mut ref_ct = setup_sample_ct();
        ref_ct.x(0);
        ref_ct.h(0);
        ref_ct.cx(0, 1);

        ref_ct.s(1);
        ref_ct.h(1);
        ref_ct.cx(1, 0);

        assert_eq!(ref_ct, third);
        assert_eq!(ref_ct, second_ct * first_ct);
    }

    #[test]
    fn test_clifford_tableau_inverse() {
        let mut ct = CliffordTableau::new(2);
        ct.x(0);
        ct.h(0);
        ct.cx(0, 1);

        let adjoint_ct = ct.adjoint();

        let ct_size = 2;

        let x_1 = bitvec![0, 0, 1, 1];
        let z_1 = bitvec![1, 0, 0, 0];
        let pauli_1 = PauliString { x: x_1, z: z_1 };

        let x_2 = bitvec![1, 1, 0, 0];
        let z_2 = bitvec![0, 0, 0, 1];
        let pauli_2 = PauliString { x: x_2, z: z_2 };

        let ct_signs = bitvec![1, 0, 0, 0];
        let ct_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs: ct_signs,
            size: ct_size,
        };

        assert_eq!(ct_ref, adjoint_ct);
        let identity = CliffordTableau::new(2);
        assert_eq!(identity, ct * adjoint_ct);
    }

    #[test]
    fn test_clifford_tableau_display() {
        let ct = setup_sample_ct();
        assert_eq!(
            ct.to_string(),
            "CliffordTableau(3)\nZ Y X I X I\nZ I I X X Y\nZ Y X I I Y\n+ - + - + +"
        );
    }
}
