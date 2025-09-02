use bitvec::prelude::BitVec;
use itertools::{izip, Itertools};
use std::fmt;
use std::iter::zip;
use std::ops::Mul;

use crate::data_structures::PauliLetter;

use super::HasAdjoint;
use super::{
    pauli_string::{cx, PauliString},
    IndexType, PropagateClifford,
};

#[derive(PartialEq, Eq, Debug, Clone, Default)]
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
            signs: BitVec::repeat(false, 2 * n),
            size: n,
        }
    }

    pub fn from_parts(pauli_columns: Vec<PauliString>, signs: BitVec) -> Self {
        let size = pauli_columns[0].len() / 2;
        CliffordTableau {
            pauli_columns,
            signs,
            size,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn signs(&self) -> &BitVec {
        &self.signs
    }

    pub(crate) fn x_signs(&self) -> BitVec {
        let n = self.size();
        self.signs[0..n].to_bitvec()
    }

    pub(crate) fn z_signs(&self) -> BitVec {
        let n = self.size();
        self.signs[n..].to_bitvec()
    }

    pub(crate) fn column(&self, i: usize) -> &PauliString {
        &self.pauli_columns[i]
    }

    pub fn compose(&self, rhs: &Self) -> Self {
        rhs.prepend(self)
    }

    /// Implements algorithms from https://doi.org/10.22331/q-2022-06-13-734 and Qiskit Clifford implementation
    pub(crate) fn prepend(&self, lhs: &Self) -> Self {
        let size = self.size();
        let pauli_columns = vec![PauliString::from_text(&"I".repeat(2 * size)); size];
        // Matrix-multiplication for M(rhs o self) = M(self) * M(rhs) as this is a row-permutation.
        // Loop re-order to be (k, i, j) as j ordering is contiguous.
        for (k, rhs_pauli_column) in self.pauli_columns.iter().enumerate() {
            for i in 0..size {
                let mut x = pauli_columns[k].x.write().unwrap();
                let mut z = pauli_columns[k].z.write().unwrap();
                *x ^= BitVec::repeat(rhs_pauli_column.x(i), 2 * size)
                    & lhs.pauli_columns[i].x.read().unwrap().as_bitslice();
                *x ^= BitVec::repeat(rhs_pauli_column.x(i + size), 2 * size)
                    & lhs.pauli_columns[i].z.read().unwrap().as_bitslice();
                *z ^= BitVec::repeat(rhs_pauli_column.z(i), 2 * size)
                    & lhs.pauli_columns[i].x.read().unwrap().as_bitslice();
                *z ^= BitVec::repeat(rhs_pauli_column.z(i + size), 2 * size)
                    & lhs.pauli_columns[i].z.read().unwrap().as_bitslice();
            }
        }

        let mut i_factors = vec![0_usize; 2 * size];
        // Keep track of the inherent i factors of left-hand tableau (where there are Y's in tableau rows)
        for lhs_pauli_column in lhs.pauli_columns.iter() {
            let local_sign = lhs_pauli_column.y_bitmask();
            for (fact, sign) in zip(i_factors.iter_mut(), local_sign) {
                *fact += sign as usize;
            }
        }

        // Accumulate the i factors when lhs basis is aggregated per rows in rhs tableau.
        // Indices reflect a (i, j) × (j, k) matrix multiplication.
        // Loop re-order to be (i, k, j).
        for (i, i_factor) in i_factors.iter_mut().enumerate() {
            for rhs_pauli_column in self.pauli_columns.iter() {
                let mut x1_select = Vec::new();
                let mut z1_select = Vec::new();
                for (j, lhs_pauli_column) in lhs.pauli_columns.iter().enumerate() {
                    if lhs_pauli_column.x(i) {
                        x1_select.push(rhs_pauli_column.x(j));
                        z1_select.push(rhs_pauli_column.z(j))
                    }
                    if lhs_pauli_column.z(i) {
                        x1_select.push(rhs_pauli_column.x(j + size));
                        z1_select.push(rhs_pauli_column.z(j + size));
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
                    .scan(false, |state, z| {
                        *state ^= z;
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
            new_signs ^= BitVec::repeat(self.signs[j], 2 * size)
                & lhs_pauli_column.x.read().unwrap().as_bitslice();
            new_signs ^= BitVec::repeat(self.signs[j + size], 2 * size)
                & lhs_pauli_column.z.read().unwrap().as_bitslice();
        }

        // Get rid of `i` factors and convert to sign flips
        let p = i_factors
            .iter()
            .map(|sign| ((sign % 4) / 2) > 0)
            .collect::<BitVec>();

        new_signs ^= p;
        new_signs ^= lhs.signs.as_bitslice();

        CliffordTableau {
            pauli_columns,
            signs: new_signs,
            size,
        }
    }

    pub fn permute(&mut self, permutation_vector: &[usize]) {
        assert_eq!(
            permutation_vector
                .iter()
                .copied()
                .sorted_unstable()
                .collect::<Vec<_>>(),
            (0..self.size()).collect::<Vec<_>>()
        );
        let pauli_columns = std::mem::take(&mut self.pauli_columns);
        let sorted_pauli_columns = zip(pauli_columns, permutation_vector)
            .sorted_unstable_by_key(|a| a.1)
            .map(|a| a.0)
            .collect::<Vec<_>>();
        self.pauli_columns = sorted_pauli_columns;
    }
}

impl HasAdjoint for CliffordTableau {
    fn adjoint(&self) -> Self {
        // Algorithm taken from https://algassert.com/post/2002
        let size = self.size();
        // Create new CliffordTableau entries

        let new_columns = vec![PauliString::from_text(&"I".repeat(2 * size)); size];
        (0..size).for_each(|i| {
            for (j, pauli_column) in self.pauli_columns.iter().enumerate() {
                let ((x1, z1), (x2, z2)) = reverse_flow(
                    pauli_column.x(i),
                    pauli_column.z(i),
                    pauli_column.x(i + size),
                    pauli_column.z(i + size),
                );

                let mut x = new_columns[i].x.write().unwrap();
                let mut z = new_columns[i].z.write().unwrap();
                x.replace(j, x1);
                z.replace(j, z1);
                x.replace(j + size, x2);
                z.replace(j + size, z2);
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

const I: (bool, bool) = (false, false);
const X: (bool, bool) = (true, false);
const Y: (bool, bool) = (true, true);
const Z: (bool, bool) = (false, true);

fn reverse_flow(x1: bool, z1: bool, x2: bool, z2: bool) -> ((bool, bool), (bool, bool)) {
    match ((x1, z1), (x2, z2)) {
        (I, I) => (I, I),
        (I, X) => (I, X),
        (I, Y) => (X, X),
        (I, Z) => (X, I),
        (X, I) => (I, Z),
        (X, X) => (I, Y),
        (X, Y) => (X, Y),
        (X, Z) => (X, Z),
        (Y, I) => (Z, Z),
        (Y, X) => (Z, Y),
        (Y, Y) => (Y, Y),
        (Y, Z) => (Y, Z),
        (Z, I) => (Z, I),
        (Z, X) => (Z, X),
        (Z, Y) => (Y, X),
        (Z, Z) => (Y, I),
    }
}

impl PropagateClifford for CliffordTableau {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let n = self.size();

        let [control, target] = self
            .pauli_columns
            .get_disjoint_mut([control, target])
            .unwrap();

        let mut scratch = BitVec::repeat(true, 2 * n);
        scratch ^= target.x.read().unwrap().as_bitslice();
        scratch ^= control.z.read().unwrap().as_bitslice();
        scratch &= control.x.read().unwrap().as_bitslice();
        scratch &= target.z.read().unwrap().as_bitslice();
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
    // fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    //     write!(f, "    ||")?;
    //     for i in 0..self.size() {
    //         write!(f, " X{} Z{}|", i + 1, i + 1)?;
    //     }
    //     writeln!(f)?;
    //     write!(f, "+/- ||")?;
    //     for (i, sign) in self.signs().iter().enumerate() {
    //         if *sign {
    //             write!(f, " - ")?;
    //         } else {
    //             write!(f, " + ")?;
    //         }
    //         if i % 2 != 0 {
    //             write!(f, "|")?;
    //         }
    //     }
    //     writeln!(f)?;

    //     for (i, column) in self.pauli_columns.iter().enumerate() {
    //         write!(f, "QB{} ||", i)?;
    //         let mut out = String::new();
    //         let mut letter_count = 0;
    //         for j in 0..column.len() {
    //             // let letter = column.pauli(j);
    //             match column.pauli(j) {
    //                 PauliLetter::I => {
    //                     out.push('I');
    //                     letter_count += 1;
    //                 }
    //                 PauliLetter::X => {
    //                     out.push('X');
    //                     letter_count += 1;
    //                 }
    //                 PauliLetter::Z => {
    //                     out.push('Z');
    //                     letter_count += 1;
    //                 }
    //                 PauliLetter::Y => {
    //                     out.push('Y');
    //                     letter_count += 1;
    //                 }
    //             }
    //             if letter_count % 2 == 1 {
    //                 out.push(' ');
    //             } else {
    //                 out.push_str(" |");
    //             }
    //             out.push(' ');
    //         }
    //         out.push('\n');
    //         write!(f, " {}", out)?;
    //     }
    //     writeln!(f)
    // }

    // build a function that push the correct string character for each Pauli letter
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "    || Stabilizers | Destabilizers |\n")?;
        let column0 = self.pauli_columns[0].len();
        for i in 0..column0 / 2 {
            write!(f, "QB{} || ", i)?;
            let sign = self.signs[i];
            if sign {
                write!(f, "- ")?;
            } else {
                write!(f, "+ ")?;
            }
            for column in self.pauli_columns.iter() {
                let mut out = String::new();
                let ch = get_pauli_char(&column.pauli(i));
                out.push(ch);
                write!(f, "{} ", out)?;
            }
            let space_left = 10 - 2 * self.pauli_columns.len();
            for _ in 0..space_left {
                write!(f, " ")?;
            }
            write!(f, "| ")?;
            let sign = self.signs[i + column0 / 2];
            if sign {
                write!(f, "- ")?;
            } else {
                write!(f, "+ ")?;
            }
            for column in self.pauli_columns.iter() {
                let mut out = String::new();
                let ch = get_pauli_char(&column.pauli(i + column0 / 2));
                out.push(ch);
                write!(f, "{} ", out)?;
            }
            let space_left = 12 - 2 * self.pauli_columns.len();
            for _ in 0..space_left {
                write!(f, " ")?;
            }
            writeln!(f, "|")?;
        }
        writeln!(f)
    }
}
pub fn get_pauli_char(letter: &PauliLetter) -> char {
    match letter {
        PauliLetter::I => 'I',
        PauliLetter::X => 'X',
        PauliLetter::Y => 'Y',
        PauliLetter::Z => 'Z',
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
        let pauli_1 = PauliString::new(x_1, z_1);
        let x_2 = bitvec![0, 1, 0, 0, 0, 0];
        let z_2 = bitvec![0, 0, 0, 0, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);
        let x_3 = bitvec![0, 0, 1, 0, 0, 0];
        let z_3 = bitvec![0, 0, 0, 0, 0, 1];
        let pauli_3 = PauliString::new(x_3, z_3);
        let signs = bitvec![0, 0, 0, 0, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3],
            signs,
            size: ct_size,
        };
        assert_eq!(ct, clifford_tableau_ref);
    }

    fn setup_sample_ct() -> CliffordTableau {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        // qubit 1x: ZYI
        // qubit 1z: IZZ
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];

        let pauli_1 = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII

        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];

        let pauli_2 = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ

        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let pauli_3 = PauliString::new(x_3, z_3);

        let signs = bitvec![0, 1, 0, 1, 0, 0];
        CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3],
            signs,
            size: ct_size,
        }
    }

    fn setup_sample_inverse_ct() -> CliffordTableau {
        // Stab: -ZIYZ, -ZZYZ, -XZXI, IZXX
        // Destab: -YYIZ, -YYXZ, ZIXX, -XZXZ
        let ct_size = 4;
        // qubit 1x: ZZXI
        // qubit 1z: YYZX
        let x_1 = bitvec![0, 0, 1, 0, 1, 1, 0, 1];
        let z_1 = bitvec![1, 1, 0, 0, 1, 1, 1, 0];
        let pauli_1 = PauliString::new(x_1, z_1);

        // qubit 2x: IZZZ
        // qubit 2z: YYIZ
        let x_2 = bitvec![0, 0, 0, 0, 1, 1, 0, 0];
        let z_2 = bitvec![0, 1, 1, 1, 1, 1, 0, 1];
        let pauli_2 = PauliString::new(x_2, z_2);

        // qubit 3x: YYXX
        // qubit 3z: IXXX
        let x_3 = bitvec![1, 1, 1, 1, 0, 1, 1, 1];
        let z_3 = bitvec![1, 1, 0, 0, 0, 0, 0, 0];
        let pauli_3 = PauliString::new(x_3, z_3);

        // qubit 3x: ZZIX
        // qubit 3z: ZZXZ
        let x_4 = bitvec![0, 0, 0, 1, 0, 0, 1, 0];
        let z_4 = bitvec![1, 1, 0, 0, 1, 1, 0, 1];

        let pauli_4 = PauliString::new(x_4, z_4);

        let signs = bitvec![1, 1, 1, 0, 1, 1, 0, 1];
        CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2, pauli_3, pauli_4],
            signs,
            size: ct_size,
        }
    }

    #[test]
    fn test_clifford_tableau_s() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S on qubit 0
        ct.s(0);

        // Stab: ZZZ, -(-X)IY, IXY
        // Destab: -IXI, ZII, ZIZ

        // qubit 1x: ZXI
        // qubit 1z: IZZ
        let z_1 = bitvec![1, 0, 0, 0, 1, 1];
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![0, 0, 0, 1, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_v() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.v(0);

        // Stab: (-Y)ZZ, -ZIY, IXY
        // Destab: -IXI, (-Y)II, (-Y)IZ

        // qubit 1x: YZI
        // qubit 1z: IYY
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];
        let x_1 = bitvec![1, 0, 0, 0, 1, 1];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![1, 1, 0, 1, 1, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_sdag() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply S_dgr on qubit 0, 1
        ct.s_dgr(0);
        ct.s_dgr(1);

        // Stab: ZZZ, -XIY, I(-Y)Y
        // Destab: -I(-Y)I, ZII, ZIZ

        // qubit 1x: ZXI
        // qubit 1z: IZZ
        let z_1 = bitvec![1, 0, 0, 0, 1, 1];
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIY
        // qubit 2z: YII
        let z_2 = bitvec![1, 0, 1, 1, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![0, 1, 1, 0, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_vdag() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Vdgr on qubit 0, 1
        ct.v_dgr(0);
        ct.v_dgr(1);

        // Stab: YYZ, -(-Z)IY, IXY
        // Destab: -IXI, YII, YIZ

        // qubit 1x: YZI
        // qubit 1z: IYY
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];
        let x_1 = bitvec![1, 0, 0, 0, 1, 1];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: YIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![1, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![0, 0, 0, 1, 0, 0];

        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_h() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.h(0);
        ct.h(1);

        // Stab: XXZ, -(-Y)IY, IZY
        // Destab: -IZI, XII, XIZ

        // qubit 1x: XYI
        // qubit 1z: IXX
        let z_1 = bitvec![0, 1, 0, 0, 0, 0];
        let x_1 = bitvec![1, 1, 0, 0, 1, 1];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: XIZ
        // qubit 2z: ZII
        let z_2 = bitvec![0, 0, 1, 1, 0, 0];
        let x_2 = bitvec![1, 0, 0, 0, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![0, 0, 0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_x() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply V on qubit 0
        ct.x(0);

        // Stab: (-Z)ZZ, -(-Y)IY, IXY
        // Destab: -IXI, (-Z)II, (-Z)IZ

        // qubit 1x: ZYI
        // qubit 1z: IZZ
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![1, 0, 0, 1, 1, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_y() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Y on qubit 0, 1
        ct.y(0);

        // Stab: (-Z)ZZ, -YIY, IXY
        // Destab: -IXI, (-Z)II, (-Z)IZ

        // qubit 1x: ZYI
        // qubit 1z: IZZ
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![1, 1, 0, 1, 1, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    #[test]
    fn test_clifford_tableau_z() {
        // Stab: ZZZ, -YIY, IXY
        // Destab: -IXI, ZII, ZIZ
        let ct_size = 3;
        let mut ct = setup_sample_ct();

        // Apply Y on qubit 0, 1
        ct.z(0);

        // Stab: ZZZ, -(-Y)IY, IXY
        // Destab: -IXI, ZII, ZIZ

        // qubit 1x: ZYI
        // qubit 1z: IZZ
        let z_1 = bitvec![1, 1, 0, 0, 1, 1];
        let x_1 = bitvec![0, 1, 0, 0, 0, 0];
        let pauli_1_ref = PauliString::new(x_1, z_1);

        // qubit 2x: ZIX
        // qubit 2z: XII
        let z_2 = bitvec![1, 0, 0, 0, 0, 0];
        let x_2 = bitvec![0, 0, 1, 1, 0, 0];
        let pauli_2_ref = PauliString::new(x_2, z_2);

        // qubit 3x: ZYY
        // qubit 3z: IIZ
        let z_3 = bitvec![1, 1, 1, 0, 0, 1];
        let x_3 = bitvec![0, 1, 1, 0, 0, 0];
        let pauli_3_ref = PauliString::new(x_3, z_3);

        let signs_ref = bitvec![0, 0, 0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1_ref, pauli_2_ref, pauli_3_ref],
            signs: signs_ref,
            size: ct_size,
        };

        assert_eq!(ct, clifford_tableau_ref);
    }

    /// This does not generate a valid Clifford Tableau. Only used to check commutation relations
    fn setup_sample_two_qubit_ct(pauli_control: char) -> CliffordTableau {
        let pauli_1_ref = match pauli_control {
            'i' => {
                // qubit 1x: II
                // qubit 1z: II
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![0, 0, 0, 0];

                PauliString::new(x_1, z_1)
            }
            'x' => {
                // qubit 1x: XX
                // qubit 1z: XX
                let z_1 = bitvec![0, 0, 0, 0];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString::new(x_1, z_1)
            }
            'y' => {
                // qubit 1x: YY
                // qubit 1z: YY
                let z_1 = bitvec![1, 1, 1, 1];
                let x_1 = bitvec![1, 1, 1, 1];

                PauliString::new(x_1, z_1)
            }
            'z' => {
                // qubit 1x: ZZ
                // qubit 1z: ZZ
                let z_1 = bitvec![1, 1, 1, 1];
                let x_1 = bitvec![0, 0, 0, 0];

                PauliString::new(x_1, z_1)
            }
            _ => panic!("Pauli letter not recognized"),
        };

        // qubit 1x: IX
        // qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];

        let pauli_2_ref = PauliString::new(x_2, z_2);

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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 1x: IX
        //qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 2x: XI
        //qubit 2z: ZY
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![1, 0, 0, 1];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 0, 1];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 2x: XI
        //qubit 2z: ZY
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![1, 0, 0, 1];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 2x: IX
        //qubit 2z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 1x: IX
        //qubit 1z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 1x: ZY
        //qubit 1z: XI
        let z_2 = bitvec![1, 1, 0, 0];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 1, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 2x: ZY
        //qubit 2z: XI
        let z_2 = bitvec![1, 1, 0, 0];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 1, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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
        let pauli_1 = PauliString::new(x_1, z_1);
        //qubit 2x: IX
        //qubit 2z: YZ
        let z_2 = bitvec![0, 0, 1, 1];
        let x_2 = bitvec![0, 1, 1, 0];
        let pauli_2 = PauliString::new(x_2, z_2);

        let signs = bitvec![0, 0, 0, 0];
        let clifford_tableau_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs,
            size: 2,
        };
        assert_eq!(ct, clifford_tableau_ref);
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

        assert_eq!(third, ref_ct);
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

        assert_eq!(third, ref_ct);
        assert_eq!(second_ct * first_ct, ref_ct);
    }

    #[test]
    fn test_reverse_flow() {
        let mut output = Vec::new();
        let ordered_ref = (0..16)
            .map(|i| {
                (
                    (i >> 3 & 1 == 1, i >> 2 & 1 == 1),
                    (i >> 1 & 1 == 1, i & 1 == 1),
                )
            })
            .collect::<Vec<_>>();

        for ((xx, xz), (zx, zz)) in ordered_ref.clone() {
            output.push(reverse_flow(xx, xz, zx, zz));
        }
        let mut sorted_output = output.clone();
        sorted_output.sort();

        for (i, j) in zip(&sorted_output, &ordered_ref) {
            assert_eq!(i, j);
        }

        let mut ordered_output = Vec::new();
        for ((xx, xz), (zx, zz)) in output {
            ordered_output.push(reverse_flow(xx, xz, zx, zz));
        }
        for (i, j) in zip(&ordered_output, &ordered_ref) {
            assert_eq!(i, j);
        }
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
        let pauli_1 = PauliString::new(x_1, z_1);

        let x_2 = bitvec![1, 1, 0, 0];
        let z_2 = bitvec![0, 0, 0, 1];
        let pauli_2 = PauliString::new(x_2, z_2);

        let ct_signs = bitvec![1, 0, 0, 0];
        let ct_ref = CliffordTableau {
            pauli_columns: vec![pauli_1, pauli_2],
            signs: ct_signs,
            size: ct_size,
        };

        assert_eq!(adjoint_ct, ct_ref);
        let identity = CliffordTableau::new(2);
        assert_eq!(ct * adjoint_ct, identity);
    }

    #[test]
    fn test_clifford_tableau_simple_inverse() {
        let mut ct = CliffordTableau::new(2);
        ct.x(0);
        ct.x(0);

        let identity = CliffordTableau::new(2);
        assert_eq!(ct, identity);
    }

    #[test]
    fn test_clifford_tableau_inverse_complex() {
        let ct = setup_sample_inverse_ct();
        let adjoint_ct = ct.adjoint();
        let identity = CliffordTableau::new(4);

        assert_eq!(ct * adjoint_ct, identity);
    }

    #[test]
    fn test_clifford_tableau_display() {
        let ct = setup_sample_ct();
        assert_eq!(
            ct.to_string(),
            "    || X1 Z1| X2 Z2| X3 Z3|\n+/- || +  - | +  - | +  + |\nQB0 || Z  Y | I  I | Z  Z | \nQB1 || Z  I | X  X | I  I | \nQB2 || Z  Y | Y  I | I  Z | \n\n"
        );
    }
}
