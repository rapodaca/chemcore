use gamma::graph::{ Error, Graph };
use super::Element;
use super::Parity;

pub trait Molecule: Graph {
    fn element(&self, id: usize) -> Result<Element, Error>;

    fn isotope(&self, id: usize) -> Result<Option<u16>, Error>;

    fn electrons(&self, id: usize) -> Result<u8, Error>;

    fn hydrogens(&self, id: usize) -> Result<u8, Error>;

    fn charge(&self, id: usize) -> Result<i8, Error>;

    fn atom_parity(&self, id: usize) -> Result<Option<Parity>, Error>;

    fn bond_order(&self, sid: usize, tid: usize) -> Result<f32, Error>;
    
    fn bond_parity(
        &self, sid: usize, tid: usize
    ) -> Result<Option<Parity>, Error>;
}

