use std::iter::zip;

use super::{
    pauli_string::{cx, PauliString},
    IndexType, PropagateClifford,
};

// todo: Make this into a union / type Angle
type Angle = f64;

pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: Vec<Angle>,
}

impl PauliPolynomial {
    pub fn from_hamiltonian(hamiltonian_repr: Vec<(String, Angle)>) -> Self {
        assert!(!hamiltonian_repr.is_empty());
        let num_qubits = hamiltonian_repr[0].0.len();
        let mut angles = Vec::<Angle>::new();
        let mut chain_strings = vec![Vec::new(); num_qubits];
        //let chains = vec![PauliString::new(); num_qubits];
        for (mut pauli_string, angle) in hamiltonian_repr {
            (0..num_qubits).for_each(|i| {
                chain_strings[i].push(pauli_string.pop().unwrap());
            });
            angles.push(angle);
        }
        let mut chains = Vec::new();
        for chain_string in chain_strings.into_iter() {
            chains.push(PauliString::from_text_string(String::from_iter(
                chain_string.iter(),
            )));
        }
        PauliPolynomial { chains, angles }
    }
}

impl PropagateClifford for PauliPolynomial {
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        let [control, target] = self.chains.get_many_mut([control, target]).unwrap();
        cx(control, target);
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
