use std::{borrow::{Borrow, BorrowMut}, cell::RefCell, rc::Rc};

use super::{pauli_string::{cx, PauliString}, PropagateClifford};

// todo: Make this into a union / type Angle
type Angle = f64;

pub struct PauliPolynomial{
    chains: Vec<PauliString>,
    angles: Vec<Angle>
}

impl PauliPolynomial{
    pub fn from_hamiltonian(hamiltonian_repr: Vec<(String, Angle)>) -> Self{
        assert!(hamiltonian_repr.len() != 0);
        let num_qubits = hamiltonian_repr[0].0.len();
        let mut angles = Vec::<Angle>::new();
        let mut chain_strings = vec![Vec::new(); num_qubits];
        //let chains = vec![PauliString::new(); num_qubits];
        for (mut pauli_string, angle) in hamiltonian_repr{
            for i in 0..num_qubits{
                chain_strings[i].push(pauli_string.pop().unwrap());
            }
            angles.push(angle);
        }
        let mut chains =  Vec::new();
        for chain_string in chain_strings.into_iter(){
            chains.push(PauliString::from_text_string(String::from_iter(chain_string.iter())));
        }
        PauliPolynomial{
            chains,
            angles
        }
    }
}

impl PropagateClifford for PauliPolynomial{
    fn cx(&mut self, control: super::IndexType, target: super::IndexType) -> &mut Self {
        match control < target {
            true => {
                let split = self.chains.split_at_mut(target); 
                cx(split.1.get_mut(0).unwrap(), split.0.get_mut(control).unwrap())
            },
            false => {
                let split = self.chains.split_at_mut(control);
                cx(split.0.get_mut(target).unwrap(), split.1.get_mut(0).unwrap())
            },
        };
        self
    }

    fn s(&mut self, target: super::IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();
        chains_target.s();

        // Update angles
        let y_vec =  chains_target.to_owned().y_bitmask();
        let _ = self.angles.iter_mut().enumerate().map(|(i, a)| match y_vec[i]{
            true => *a *= -1.0,
            false => (),
        });
        self
    }

    fn v(&mut self, target: super::IndexType) -> &mut Self {
        let chains_target = self.chains.get_mut(target).unwrap();

        // Update angles
        let y_vec =  chains_target.to_owned().y_bitmask();
        let _ = self.angles.iter_mut().enumerate().map(|(i, a)| match y_vec[i]{
            true => *a *= -1.0,
            false => (),
        });
        chains_target.s();
        self
    }

}