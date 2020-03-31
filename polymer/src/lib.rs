use monomer::{Monomer};
use once_cell::{OnceVal};

// TODO: See how to implement this thing as an iterator.
pub trait Polymer<T> where T: Monomer + PartialEq {
    fn new() -> Self;
    fn push(&mut self, t: OnceVal<T>);
    fn concat(&mut self, other: &mut Self);
}