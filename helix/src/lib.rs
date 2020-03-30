use dna::{Cat, Nucl, Nucleotide};

type Strand = Vec<Nucl>;

#[derive(Debug)]
pub struct Helix {
    val: Strand,
    src: Cat,
}

impl Helix {
    pub fn new(c: Cat) -> Helix {
        Helix {
            val: Vec::<Nucl>::new(),
            src: c,
        }
    }

    pub fn strands(h: Helix) -> (Helix, Helix) {
        let x = h.inverse();
        (h, x)
    }

    pub fn push(&mut self, n: Nucleotide) {
        let next = self.wrap(n);
        self.val.push(next);
    }

    pub fn concat(&mut self, other: &mut Helix) {
        self.val.append(&mut other.val)
    }

    // TODO: Optimize?
    // TODO: Bubble err, not None
    pub fn from_string(s: String, cat: Cat) -> Option<Helix> {
        let mut h = Helix::new(cat);
        let mut ok = true;
        for c in s.chars() {
            match dna::from_char(c) {
                Some(x) => h.push(x),
                None => {
                    ok = false;
                    break;
                }
            }
        }

        match ok {
            true => Some(h),
            _ => None,
        }
    }

    pub fn inverse(&self) -> Helix {
        let mut next = Helix::new(self.src.clone());

        let mut next_val = Vec::<Nucl>::new();

        for x in &self.val {
            next_val.push(self.src.compliment(x.clone()));
        }

        next.val = next_val.into_iter().rev().collect();
        next
    }

    pub fn pairs(&self) -> Vec<(Nucl, Nucl)> {
        let mut next = Vec::<(Nucl, Nucl)>::new();
        for x in &self.val {
            next.push((x.clone(), self.src.compliment(x.clone())));
        }
        next
    }

    // TODO: better result number then tuple fraction
    pub fn gc_content(&self) -> (u64, u64) {
        let mut n: u64 = 0;
        let mut d: u64 = 0;
        for x in &self.val {
            d = d + 1;
            if Helix::is_g_or_c(x.read().as_ref().unwrap()) {
                n = n + 1;
            }
        }

        (n, d)
    }

    fn is_g_or_c(n: &Nucleotide) -> bool {
        match n {
            Nucleotide::G => true,
            Nucleotide::C => true,
            _ => false,
        }
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
    use super::*;

    #[test]
    fn helix_new() {
        let cat = Cat::new();
        let mut h = Helix::new(cat.clone());
        let mut tvec = Vec::<Nucl>::new();
        assert_eq!(&tvec, &h.val);
        h.push(Nucleotide::A);
        tvec.push(cat.A());

        assert_eq!(&tvec, &h.val);

        h.push(Nucleotide::T);
        tvec.push(cat.T());

        assert_eq!(&tvec, &h.val);

        h.push(Nucleotide::C);
        tvec.push(cat.C());

        assert_eq!(&tvec, &h.val);

        h.push(Nucleotide::G);
        tvec.push(cat.G());

        assert_eq!(&tvec, &h.val);
    }

    fn helix_from_string() {
        let cat = Cat::new();

        let hxstr = String::from("gattaca");
        let badstr = String::from("a123g");

        let mut controlH = Helix::new(cat.clone());

        controlH.push(Nucleotide::G);
        controlH.push(Nucleotide::A);
        controlH.push(Nucleotide::T);
        controlH.push(Nucleotide::T);
        controlH.push(Nucleotide::A);
        controlH.push(Nucleotide::C);
        controlH.push(Nucleotide::A);

        let maybeNone = Helix::from_string(badstr, cat.clone());
        match maybeNone {
            None => println!(""),
            Some(x) => panic!("Should've recieved nothing, got: {:?}", x)
        }

        let maybeH = Helix::from_string(hxstr, cat.clone());
        match maybeH {
            None => panic!("Failed in Helix from_string with good string"),
            Some(h) => assert_eq!(h.val, controlH.val),
        }
    }
}
