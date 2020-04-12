#[macro_use]
use rna::{RNA, RNACat, RNACell};
use amino::Amino;
use category::{Cat, ICat, NCat};
use monomer::{IMono, Mono, NucleicAcid};
use once_mono::IMonomer;
use polymer::{Polymer, Strand};

pub struct Codon {
    pub rna: (RNACell, RNACell, RNACell),
}

// NOTE: could be a more general -- but this is to test quartz.
impl Codon {
    pub fn from_strand(s: Strand<RNACell>) -> Option<Codon> {
        let mut x = None;
        if s.contents.len() == 3 {
            x = Some(Codon {
                rna: Self::tuple_from_strand(s),
            });
        }
        x
    }

    fn tuple_from_strand(s: Strand<RNACell>) -> (RNACell, RNACell, RNACell) {
        (
            s.contents[0].clone(),
            s.contents[1].clone(),
            s.contents[2].clone(),
        )
    }

    // NOTE: Ignoring the AUG => START command. Will use seperate function for this.
    pub fn to_amino_encoding(&self) -> Amino {
        let (fst, snd, thd) = &self.rna;

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

// BRUTE FORCE:
// pub enum Codon {
//     AAA,
//     AAU,
//     AAC,
//     AAG,
//     AUA,
//     AUU,
//     AUC,
//     AUG,
//     ACA,
//     ACU,
//     ACC,
//     ACG,
//     AGA,
//     AGU,
//     AGC,
//     AGG,
//     UAA,
//     UAU,
//     UAC,
//     UAG,
//     UUA,
//     UUU,
//     UUC,
//     UUG,
//     UCA,
//     UCU,
//     UCC,
//     UCG,
//     UGA,
//     UGU,
//     UGC,
//     UGG,
//     CAA,
//     CAU,
//     CAC,
//     CAG,
//     CUA,
//     CUU,
//     CUC,
//     CUG,
//     CCA,
//     CCU,
//     CCC,
//     CCG,
//     CGA,
//     CGU,
//     CGC,
//     CGG,
//     GAA,
//     GAU,
//     GAC,
//     GAG,
//     GUA,
//     GUU,
//     GUC,
//     GUG,
//     GCA,
//     GCU,
//     GCC,
//     GCG,
//     GGA,
//     GGU,
//     GGC,
//     GGG,
// }
