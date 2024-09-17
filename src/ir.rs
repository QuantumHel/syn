mod pauli_exponential;

trait Synthesizer<Program> {
    fn synthesize(program: Program);
}

struct MySynthesizer {
    // e.g. heuristic function
}

struct MyProgram {
    // e.g. pauli polynomial and clifford tableu
}

impl Synthesizer<MyProgram> for MySynthesizer {
    fn synthesize(program: MyProgram) {
        todo!()
    }
}
