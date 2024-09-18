use std::{borrow::BorrowMut, cell::RefCell, rc::Rc};

use super::{pauli_string::{cx, PauliString}, PropagateClifford};

type Angle = f64;

pub struct PauliPolynomial{
    chains: Vec<Rc<RefCell<PauliString>>>,
    // todo: Make this into a union / type Angle
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
            chains.push(Rc::new(RefCell::new(PauliString::from_text_string(String::from_iter(chain_string.iter())))));
        }
        PauliPolynomial{
            chains,
            angles
        }
    }
}

impl PropagateClifford for PauliPolynomial{
    fn cx(&mut self, control: super::IndexType, target: super::IndexType) -> &mut Self {
        cx(self.chains[control].borrow().clone(), self.chains[target].borrow().clone());
        self
    }

    fn s(&mut self, target: super::IndexType) -> &mut Self {
        self.chains[target].borrow().clone().s();
        // Update angles
        let y_vec = self.chains[target].borrow().clone().y_bitmask();
        let _ = self.angles.iter_mut().enumerate().map(|(i, a)| match y_vec[i]{
            true => -1.0 * *a,
            false => *a,
        });
        self
    }

    fn v(&mut self, target: super::IndexType) -> &mut Self {
        // Update angles
        let y_vec = self.chains[target].borrow().clone().y_bitmask();
        let _ = self.angles.iter_mut().enumerate().map(|(i, a)| match y_vec[i]{
            true => -1.0 * *a,
            false => *a,
        });
        self.chains[target].borrow().clone().v();
        self
    }

}