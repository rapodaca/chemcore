use super::Atom;
use gamma::graph::{Error, Graph};

pub trait Molecule: Graph {
    /// Returns the atomic attributes associated with id,
    /// or Error if id not found.
    fn atom(&self, id: usize) -> Result<&Atom, Error>;

    /// Returns the charge computation associated with id,
    /// or Error if id not found.
    fn charge(&self, id: usize) -> Result<f32, Error>;

    /// Returns the bond order computation associated with the source
    /// and target ids, or Error if either sid or tid not found.
    fn bond_order(&self, sid: usize, tid: usize) -> Result<f32, Error>;
}
