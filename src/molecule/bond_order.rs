#[derive(Clone, Copy, Eq, Hash, PartialEq, Debug)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple
}

impl Default for BondOrder {
    fn default() -> Self {
        Self::Single
    }
}

impl BondOrder {
    pub fn multiplicity(&self) -> u8 {
        match self {
            BondOrder::Zero =>      0,
            BondOrder::Single =>    1,
            BondOrder::Double =>    2,
            BondOrder::Triple =>    3,
            BondOrder::Quadruple => 4
        }
    }
}