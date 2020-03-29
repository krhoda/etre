use dna::{Cat, Nucl, Nucleotide};

#[derive(Debug)]
struct Helix {
    val: Vec<Nucl>,
    src: Cat,
}

impl Helix {
    pub fn new(c: Cat) -> Helix {
        Helix {
            val: Vec::<Nucl>::new(),
            src: c,
        }
    }

    pub fn push(&mut self, n: Nucleotide) {
        let next = self.wrap(n);
        self.val.push(next);
    }

    pub fn concat(&mut self, other: &mut Helix) {
        self.val.append(&mut other.val)
    }

    fn wrap(&self, n: Nucleotide) -> Nucl {
        match n {
            Nucleotide::A => self.src.A(),
            Nucleotide::T => self.src.T(),
            Nucleotide::G => self.src.G(),
            Nucleotide::C => self.src.C(),
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
