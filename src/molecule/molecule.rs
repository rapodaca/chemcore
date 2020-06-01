use gamma::graph::Graph;
use gamma::graph::Error;

use super::element::Element;
use super::parity::Parity;

/// Exposes services required by a minimal molecule API.
pub trait Molecule<'a, N: 'a>: Graph<'a, usize> {
    fn element(&self, id: &usize) -> Result<Element, Error>;

    fn isotope(&self, id: &usize) -> Result<Option<u16>, Error>;

    fn electrons(&self, id: &usize) -> Result<u8, Error>;

    fn hydrogens(&self, id: &usize) -> Result<u8, Error>;

    fn charge(&self, id: &usize) -> Result<i8, Error>;

    fn atom_parity(&self, id: &usize) -> Result<Option<Parity>, Error>;

    fn bond_order(&self, sid: &usize, tid: &usize) -> Result<u8, Error>;
    
    fn bond_parity(
        &self, sid: &usize, tid: &usize
    ) -> Result<Option<Parity>, Error>;
}