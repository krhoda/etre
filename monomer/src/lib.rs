pub trait Monomer {
    fn from_char(c: char) -> Option<Self> where Self: std::marker::Sized;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
