pub trait Mono: std::marker::Sized + PartialEq + Clone {
    fn from_string(c: String) -> Option<Self>;
}

pub trait IMono: Mono {
    fn inverse(c: &Self) -> Self;
}

pub trait NucleicAcid: IMono {
    fn is_g_or_c(&self) -> bool;
}
