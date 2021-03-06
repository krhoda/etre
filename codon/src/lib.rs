use rna::{RNA, RNACell};
use amino::Amino;
use polymer::Strand;

type RNACodon = (RNACell, RNACell, RNACell);

#[derive(Debug)]
pub struct Codon {
    pub amino: Amino,
    pub rna: RNACodon,
}

impl PartialEq for Codon {
    fn eq(&self, other: &Self) -> bool {
        self.amino == other.amino
    }
}

// NOTE: could be a more general -- but this is to test quartz.
impl Codon {
    // deep_equal compares the underlying RNA strand as well.
    pub fn deep_equal(&self, other: &Self) -> bool {
        if self == other {
            self.rna == other.rna 
        } else {
            false
        }
    }

    pub fn from_strand(s: &Strand<RNACell>) -> Option<Codon> {
        let mut x = None;
        let y = Self::tuple_from_strand(&s);
        if s.contents.len() == 3 {
            x = Some(Codon {
                amino: Self::to_amino_encoding(&y),
                rna: y,
            });
        }
        x
    }

    fn tuple_from_strand(s: &Strand<RNACell>) -> (RNACell, RNACell, RNACell) {
        (
            s.contents[0].clone(),
            s.contents[1].clone(),
            s.contents[2].clone(),
        )
    }

    // NOTE: Ignoring the AUG => START command. Will use seperate function for this.
    pub fn to_amino_encoding(r: &RNACodon) -> Amino {
        //
        let (fst, snd, thd) = r;

        let x = fst.read();
        let y = snd.read();
        let z = thd.read();

        let plain_codon = (
            x.as_ref().unwrap(),
            y.as_ref().unwrap(),
            z.as_ref().unwrap(),
        );

        match plain_codon {
            // Reference: https://cnx.org/contents/GFy_h8cu@9.87:QEibhJMi@8/The-Genetic-Code
            // Square 1:
            (&RNA::U, &RNA::U, &RNA::A) | (&RNA::U, &RNA::U, &RNA::G) => Amino::Leu,
            (&RNA::U, &RNA::U, &RNA::U) | (&RNA::U, &RNA::U, &RNA::C) => Amino::Phe,

            // Square 2:
            (&RNA::U, &RNA::C, &RNA::U)
            | (&RNA::U, &RNA::C, &RNA::C)
            | (&RNA::U, &RNA::C, &RNA::A)
            | (&RNA::U, &RNA::C, &RNA::G) => Amino::Ser,

            // Square 3:
            (&RNA::U, &RNA::A, &RNA::U) | (&RNA::U, &RNA::A, &RNA::C) => Amino::Tyr,
            (&RNA::U, &RNA::A, &RNA::A) | (&RNA::U, &RNA::A, &RNA::G) => Amino::STOP,

            // Square 4:
            (&RNA::U, &RNA::G, &RNA::U) | (&RNA::U, &RNA::G, &RNA::C) => Amino::Cys,
            (&RNA::U, &RNA::G, &RNA::A) => Amino::STOP,
            (&RNA::U, &RNA::G, &RNA::G) => Amino::Trp,

            // Square 5:
            (&RNA::C, &RNA::U, &RNA::U)
            | (&RNA::C, &RNA::U, &RNA::C)
            | (&RNA::C, &RNA::U, &RNA::A)
            | (&RNA::C, &RNA::U, &RNA::G) => Amino::Leu,

            // Square 6:
            (&RNA::C, &RNA::C, &RNA::U)
            | (&RNA::C, &RNA::C, &RNA::C)
            | (&RNA::C, &RNA::C, &RNA::A)
            | (&RNA::C, &RNA::C, &RNA::G) => Amino::Pro,

            // Square 7:
            (&RNA::C, &RNA::A, &RNA::U) | (&RNA::C, &RNA::A, &RNA::C) => Amino::His,
            (&RNA::C, &RNA::A, &RNA::A) | (&RNA::C, &RNA::A, &RNA::G) => Amino::Gln,

            // Square 8:
            (&RNA::C, &RNA::G, &RNA::U)
            | (&RNA::C, &RNA::G, &RNA::C)
            | (&RNA::C, &RNA::G, &RNA::A)
            | (&RNA::C, &RNA::G, &RNA::G) => Amino::Arg,

            // Square 9:
            (&RNA::A, &RNA::U, &RNA::U)
            | (&RNA::A, &RNA::U, &RNA::C)
            | (&RNA::A, &RNA::U, &RNA::A) => Amino::Ile,
            // NOTE: Not start, assumed start has been dealt with
            (&RNA::A, &RNA::U, &RNA::G) => Amino::Met,

            // Square 10:
            (&RNA::A, &RNA::C, &RNA::U)
            | (&RNA::A, &RNA::C, &RNA::C)
            | (&RNA::A, &RNA::C, &RNA::A)
            | (&RNA::A, &RNA::C, &RNA::G) => Amino::Thr,

            // Square 11:
            (&RNA::A, &RNA::A, &RNA::U) | (&RNA::A, &RNA::A, &RNA::C) => Amino::Asn,
            (&RNA::A, &RNA::A, &RNA::A) | (&RNA::A, &RNA::A, &RNA::G) => Amino::Lys,

            // Square 12:
            (&RNA::A, &RNA::G, &RNA::U) | (&RNA::A, &RNA::G, &RNA::C) => Amino::Ser,
            (&RNA::A, &RNA::G, &RNA::A) | (&RNA::A, &RNA::G, &RNA::G) => Amino::Arg,

            // Square 13:
            (&RNA::G, &RNA::U, &RNA::U)
            | (&RNA::G, &RNA::U, &RNA::C)
            | (&RNA::G, &RNA::U, &RNA::A)
            | (&RNA::G, &RNA::U, &RNA::G) => Amino::Val,

            // Square 14:
            (&RNA::G, &RNA::C, &RNA::U)
            | (&RNA::G, &RNA::C, &RNA::C)
            | (&RNA::G, &RNA::C, &RNA::A)
            | (&RNA::G, &RNA::C, &RNA::G) => Amino::Ala,

            // Square 15:
            (&RNA::G, &RNA::A, &RNA::U) | (&RNA::G, &RNA::A, &RNA::C) => Amino::Asp,
            (&RNA::G, &RNA::A, &RNA::A) | (&RNA::G, &RNA::A, &RNA::G) => Amino::Glu,

            // Square 16:
            (&RNA::G, &RNA::G, &RNA::U)
            | (&RNA::G, &RNA::G, &RNA::C)
            | (&RNA::G, &RNA::G, &RNA::A)
            | (&RNA::G, &RNA::G, &RNA::G) => Amino::Gly,
        }
    }
}

mod tests {
    use super::*;
    use category::Cat;
    use rna::RNACat;

    #[test]
    fn codon_from_strand() {
        let cat = RNACat::new();
        let strand = cat.polymer_from_string(String::from("aug")).unwrap();
        let codon = Codon::from_strand(&strand).unwrap();
        let test = (
            cat.monomer_from_string(String::from("a")).unwrap(),
            cat.monomer_from_string(String::from("u")).unwrap(),
            cat.monomer_from_string(String::from("g")).unwrap(),
        );
        assert_eq!(codon.rna, test);
    }

    #[test]
    fn amino_from_codon() {
        // TODO: MORE RIGOROUS ACC. TESTS.
        let cat = RNACat::new();
        let strand = cat.polymer_from_string(String::from("aug")).unwrap();
        let codon = Codon::from_strand(&strand).unwrap();

        let amino = Amino::Met;

        assert_eq!(amino, codon.amino);
        assert_eq!((strand.contents[0].clone(), strand.contents[1].clone(), strand.contents[2].clone()), codon.rna);

    }
}
