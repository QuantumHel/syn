use std::iter::zip;

use bitvec::vec::BitVec;

use super::{
    pauli_string::{cx, PauliString},
    IndexType, PropagateClifford,
};

// todo: Make this into a union / type Angle
type Angle = f64;

#[derive(Debug, Clone)]
pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: Vec<Angle>,
    size: usize,
}

impl PauliPolynomial {
    pub fn from_hamiltonian(hamiltonian_representation: Vec<(&str, Angle)>) -> Self {
        assert!(!hamiltonian_representation.is_empty());
        let num_qubits = hamiltonian_representation[0].0.len();
        let mut angles = Vec::<Angle>::new();
        let mut chain_strings = vec![String::new(); num_qubits];
        //let chains = vec![PauliString::new(); num_qubits];
        for (pauli_string, angle) in hamiltonian_representation {
            assert!(pauli_string.len() == chain_strings.len());
            zip(chain_strings.iter_mut(), pauli_string.chars()).for_each(
                |(chain, pauli_letter)| {
                    chain.push(pauli_letter);
                },
            );
            angles.push(angle);
        }
        let chains = chain_strings
            .iter()
            .map(|gadget| PauliString::from_text(gadget))
            .collect::<Vec<PauliString>>();

        PauliPolynomial {
            chains,
            angles,
            size: num_qubits,
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }
}

impl PropagateClifford for PauliPolynomial {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let mut bit_mask = BitVec::repeat(true, self.angles.len());
        let [control, target] = self.chains.get_many_mut([control, target]).unwrap();

        bit_mask ^= &control.z;
        bit_mask ^= &target.x;
        bit_mask &= &control.x;
        bit_mask &= &target.z;

        cx(control, target);
        for (angle, flip) in zip(self.angles.iter_mut(), bit_mask.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }

        self
    }

    fn s(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        // Update angles
        let y_vec = chains_target.y_bitmask();
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }
        chains_target.s();
        self
    }

    fn v(&mut self, target: IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        chains_target.v();
        // Update angles
        let y_vec = chains_target.to_owned().y_bitmask();
        for (angle, flip) in zip(self.angles.iter_mut(), y_vec.iter()) {
            if *flip {
                *angle *= -1.0;
            }
        }

        self
    }
}
