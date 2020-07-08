use super::{ Element, Parity };

#[derive(Debug, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub hydrogens: u8,
    pub ion: i8,
    pub isotope: Option<u16>,
    pub parity: Option<Parity>
}

impl Atom {
    pub fn new(
        element: Element, hydrogens: u8, ion: i8, isotope: Option<u16>,
        parity: Option<Parity>
    ) -> Self {
        Atom { element, hydrogens, ion, isotope, parity }
    }
}