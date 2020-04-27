use monomer::{Mono};
use once_mono::{Monomer};
use category::{Cat};
use polymer::{Polymer, Strand};

// TODO MAKE MONOMER BY CHANGING FROM CHAR TO FROM STR
#[derive(PartialEq, Debug, Clone)]
pub enum Amino {
    START,
    STOP,
    Ala,
    Arg,
    Asn,
    Asp,
    Cys,
    Gln,
    Glu,
    Gly,
    His,
    Ile,
    Leu,
    Lys,
    Met,
    Phe,
    Pro,
    Ser,
    Thr,
    Trp,
    Tyr,
    Val,
}

// TODO: MAKE COMPLIANT WITH GAF -> Gene Association Files.
// For now, a poor folk's implementation imagining Protiens as CSV.
impl Mono for Amino {
    fn from_string(s: String) -> Option<Amino> {
        match s.to_uppercase().as_ref() {
            "START" => Some(Amino::START),
            "STOP" => Some(Amino::STOP),
            "ALA" => Some(Amino::Ala),
            "ARG" => Some(Amino::Arg),
            "ASN" => Some(Amino::Asn),
            "ASP" => Some(Amino::Asp),
            "CYS" => Some(Amino::Cys),
            "GLN" => Some(Amino::Gln),
            "GLU" => Some(Amino::Glu),
            "GLY" => Some(Amino::Gly),
            "HIS" => Some(Amino::His),
            "ILE" => Some(Amino::Ile),
            "LEU" => Some(Amino::Leu),
            "LYS" => Some(Amino::Lys),
            "MET" => Some(Amino::Met),
            "PHE" => Some(Amino::Phe),
            "PRO" => Some(Amino::Pro),
            "SER" => Some(Amino::Ser),
            "THR" => Some(Amino::Thr),
            "TRP" => Some(Amino::Trp),
            "TYR" => Some(Amino::Tyr),
            "VAL" => Some(Amino::Val),
            _ => None,
        }
    }
}

pub type AminoCell = Monomer<Amino>;

pub struct AminoMorphisms {
    pub start: AminoCell,
    pub stop: AminoCell,
    pub ala: AminoCell,
    pub arg: AminoCell,
    pub asn: AminoCell,
    pub asp: AminoCell,
    pub cys: AminoCell,
    pub gln: AminoCell,
    pub glu: AminoCell,
    pub gly: AminoCell,
    pub his: AminoCell,
    pub ile: AminoCell,
    pub leu: AminoCell,
    pub lys: AminoCell,
    pub met: AminoCell,
    pub phe: AminoCell,
    pub pro: AminoCell,
    pub ser: AminoCell,
    pub thr: AminoCell,
    pub trp: AminoCell,
    pub tyr: AminoCell,
    pub val: AminoCell,
}

impl AminoMorphisms {
    pub fn new() -> Self {
        AminoMorphisms {
            start: Monomer::from_string(String::from("start")).unwrap(),
            stop: Monomer::from_string(String::from("stop")).unwrap(),
            ala: Monomer::from_string(String::from("ala")).unwrap(),
            arg: Monomer::from_string(String::from("arg")).unwrap(),
            asn: Monomer::from_string(String::from("asn")).unwrap(),
            asp: Monomer::from_string(String::from("asp")).unwrap(),
            cys: Monomer::from_string(String::from("cys")).unwrap(),
            gln: Monomer::from_string(String::from("gln")).unwrap(),
            glu: Monomer::from_string(String::from("glu")).unwrap(),
            gly: Monomer::from_string(String::from("gly")).unwrap(),
            his: Monomer::from_string(String::from("his")).unwrap(),
            ile: Monomer::from_string(String::from("ile")).unwrap(),
            leu: Monomer::from_string(String::from("leu")).unwrap(),
            lys: Monomer::from_string(String::from("lys")).unwrap(),
            met: Monomer::from_string(String::from("met")).unwrap(),
            phe: Monomer::from_string(String::from("phe")).unwrap(),
            pro: Monomer::from_string(String::from("pro")).unwrap(),
            ser: Monomer::from_string(String::from("ser")).unwrap(),
            thr: Monomer::from_string(String::from("thr")).unwrap(),
            trp: Monomer::from_string(String::from("trp")).unwrap(),
            tyr: Monomer::from_string(String::from("tyr")).unwrap(),
            val: Monomer::from_string(String::from("val")).unwrap(),
        }
    }
}

pub struct AminoCat {
    pub morphisms: AminoMorphisms,
}

impl Cat<AminoCell, Strand<AminoCell>> for AminoCat {
    fn new() -> Self {
        AminoCat {
            morphisms: AminoMorphisms::new(),
        }
    }

    fn monomer_from_string(&self, s: String) -> Option<AminoCell> {
        match Amino::from_string(s) {
            None => None,
            Some(x) => match x {
                Amino::START => Some(self.morphisms.start.clone()),
                Amino::STOP => Some(self.morphisms.stop.clone()),
                Amino::Ala => Some(self.morphisms.ala.clone()),
                Amino::Arg => Some(self.morphisms.arg.clone()),
                Amino::Asn => Some(self.morphisms.asn.clone()),
                Amino::Asp => Some(self.morphisms.asp.clone()),
                Amino::Cys => Some(self.morphisms.cys.clone()),
                Amino::Gln => Some(self.morphisms.gln.clone()),
                Amino::Glu => Some(self.morphisms.glu.clone()),
                Amino::Gly => Some(self.morphisms.gly.clone()),
                Amino::His => Some(self.morphisms.his.clone()),
                Amino::Ile => Some(self.morphisms.ile.clone()),
                Amino::Leu => Some(self.morphisms.leu.clone()),
                Amino::Lys => Some(self.morphisms.lys.clone()),
                Amino::Met => Some(self.morphisms.met.clone()),
                Amino::Phe => Some(self.morphisms.phe.clone()),
                Amino::Pro => Some(self.morphisms.pro.clone()),
                Amino::Ser => Some(self.morphisms.ser.clone()),
                Amino::Thr => Some(self.morphisms.thr.clone()),
                Amino::Trp => Some(self.morphisms.trp.clone()),
                Amino::Tyr => Some(self.morphisms.tyr.clone()),
                Amino::Val => Some(self.morphisms.val.clone()),
            }
        }
    }

    // NOTE: Not in final form just to satisfy the interface.
    // Assumes no START or STOP codons, as they will be infered via RIBOSOME.
    fn polymer_from_string(&self, s: String) -> Option<Strand<AminoCell>> {
        let mut x = Strand::<AminoCell>::new();
        match s.len() % 3 == 0 {
            false => None,
            _ => {
                let mut ok = true;
                let mut z = s.chars().peekable();
                while z.peek().is_some() {
                    let chunk: String = z.by_ref().take(3).collect();
                    match self.monomer_from_string(chunk) {
                        Some(y) => x.push(y),
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }

                match ok {
                    true => Some(x),
                    _ => None,
                }
            }
        }
    }
}