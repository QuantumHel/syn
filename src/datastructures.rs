mod clifford_tableau;
mod pauli_polynomial;
mod pauli_string;

type IndexType = usize;

pub trait PropagateClifford where Self:Sized{
    
    fn cx(&mut self, control: IndexType, target: IndexType) -> &mut Self;
    fn s(&mut self, target: IndexType) -> &mut Self;
    fn v(&mut self, target: IndexType) -> &mut Self;
    fn y(&mut self, target: IndexType) -> &mut Self{
        self.s_dgr(target).x(target).s(target)
    }

    fn x(&mut self, target: IndexType) -> &mut Self{
        self.v(target).v(target)
    }

    fn z(&mut self, target: IndexType) -> &mut Self{
        self.s(target).s(target)
    }

    fn s_dgr(&mut self, target: IndexType) -> &mut Self{
        self.z(target).s(target)
    }

    fn v_dgr(&mut self, target: IndexType) -> &mut Self{
        self.x(target).v(target)
    }

    fn h(&mut self, target: IndexType) -> &mut Self {
        self.s(target).v_dgr(target).s(target)
    }

    fn cz(&mut self, control: IndexType, target: IndexType) -> &mut Self {
        self.h(target);
        self.cx(control, target);
        self.h(target)
    }
}