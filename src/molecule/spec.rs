use super::element::Element;
use super::parity::Parity;
use super::bond_order::BondOrder;

#[derive(Default)]
pub struct Atom {
    pub element: Element,
    pub hydrogens: u8,
    pub ion: i8,
    pub isotope: Option<u16>,
    pub parity: Option<Parity>
}

#[derive(Default)]
pub struct Bond {
    pub sid: usize,
    pub tid: usize,
    pub order: BondOrder,
    pub parity: Option<Parity>
}

pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>
}