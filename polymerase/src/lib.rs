#[macro_use]
use rna;
use dna;

pub fn dna_to_mrna(d: dna::DNA) -> rna::RNA {
    match d {
        dna::DNA::A => rna::RNA::A,
        dna::DNA::T => rna::RNA::U,
        dna::DNA::G => rna::RNA::G,
        dna::DNA::C => rna::RNA::C,
    }
}

pub fn mrna_to_dna(r: rna::RNA) -> dna::DNA {
    match r {
        rna::RNA::A => dna::DNA::A,
        rna::RNA::U => dna::DNA::T,
        rna::RNA::G => dna::DNA::G,
        rna::RNA::C => dna::DNA::C,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn t_mrna_to_dna() {
        assert_eq!(mrna_to_dna(rna::RNA::A), dna::DNA::A);
        assert_eq!(mrna_to_dna(rna::RNA::U), dna::DNA::T);
        assert_eq!(mrna_to_dna(rna::RNA::C), dna::DNA::C);
        assert_eq!(mrna_to_dna(rna::RNA::G), dna::DNA::G);
    }

    #[test]
    fn t_dna_to_mrna() {
        assert_eq!(dna_to_mrna(dna::DNA::A), rna::RNA::A);
        assert_eq!(dna_to_mrna(dna::DNA::T), rna::RNA::U);
        assert_eq!(dna_to_mrna(dna::DNA::C), rna::RNA::C);
        assert_eq!(dna_to_mrna(dna::DNA::G), rna::RNA::G);
    }
}