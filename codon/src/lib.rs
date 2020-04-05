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
    pub fn to_amino_encoding(&self, first: bool) -> Amino {
        let (fst, snd, thd) = &self.rna;
        let x = fst.read();
        let y = snd.read();
        let z = thd.read();
        let plain_codon = (x.as_ref().unwrap(), y.as_ref().unwrap(), z.as_ref().unwrap());
        match plain_codon {
            // Reference: https://cnx.org/contents/GFy_h8cu@9.87:QEibhJMi@8/The-Genetic-Code
            // Square 1:
            (&RNA::U, &RNA::U, &RNA::A) | (&RNA::U, &RNA::U, &RNA::G) => Amino::Leu,
            (&RNA::U, &RNA::U, &RNA::U) | (&RNA::U, &RNA::U, &RNA::C) => Amino::Phe,

            // Square 2:
            // TODO: COMPLETE:
            (&RNA::U, &RNA::C, &RNA::U) | (&RNA::U, &RNA::C, &RNA::C) => Amino::Ser,

            // Used to compile.
            // TODO: REMOVE:
            _ => Amino::STOP
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
