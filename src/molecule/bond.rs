use super::Parity;

#[derive(Debug,PartialEq)]
pub struct Bond {
    pub electrons: u8,
    pub parity: Option<Parity>,
    pub tid: usize
}

impl Bond {
    pub fn new(electrons: u8, parity: Option<Parity>, tid: usize) -> Self {
        Self {
            electrons,
            parity,
            tid
        }
    }

    pub fn order(&self) -> f32 {
        self.electrons as f32 / 2 as f32
    }
}