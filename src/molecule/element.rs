use purr::parts;

#[rustfmt::skip]
#[derive(PartialEq, Debug, Clone)]
pub enum Element {
//  0   1   2   3   4   5   6   7   8   9
        H,  He, Li, Be, B,  C,  N,  O,  F,  //  0
    Ne, Na, Mg, Al, Si, P,  S,  Cl, Ar, K,  //  1
    Ca, Sc, Ti, V, Cr,  Mn, Fe, Co, Ni, Cu, //  2
    Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y,  //  3
    Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, //  4
    Sn, Sb, Te, I,  Xe, Cs, Ba, La, Ce, Pr, //  5
    Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, //  6
    Yb, Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, //  7
    Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, //  8
    Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, //  9
    Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, // 10
    Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og      // 11
}

impl Into<Element> for &parts::Aliphatic {
    fn into(self) -> Element {
        match self {
            parts::Aliphatic::B => Element::B,
            parts::Aliphatic::C => Element::C,
            parts::Aliphatic::N => Element::N,
            parts::Aliphatic::O => Element::O,
            parts::Aliphatic::P => Element::P,
            parts::Aliphatic::S => Element::S,
            parts::Aliphatic::F => Element::F,
            parts::Aliphatic::Cl => Element::Cl,
            parts::Aliphatic::Br => Element::Br,
            parts::Aliphatic::I => Element::I,
            parts::Aliphatic::At => Element::At,
            parts::Aliphatic::Ts => Element::Ts,
        }
    }
}

impl Into<Element> for &parts::Aromatic {
    fn into(self) -> Element {
        match self {
            parts::Aromatic::B => Element::B,
            parts::Aromatic::C => Element::C,
            parts::Aromatic::N => Element::N,
            parts::Aromatic::O => Element::O,
            parts::Aromatic::P => Element::P,
            parts::Aromatic::S => Element::S,
        }
    }
}

impl Into<Element> for &parts::BracketAromatic {
    fn into(self) -> Element {
        match self {
            parts::BracketAromatic::As => Element::As,
            parts::BracketAromatic::B => Element::B,
            parts::BracketAromatic::C => Element::C,
            parts::BracketAromatic::N => Element::N,
            parts::BracketAromatic::O => Element::O,
            parts::BracketAromatic::P => Element::P,
            parts::BracketAromatic::S => Element::S,
            parts::BracketAromatic::Se => Element::Se,
        }
    }
}

impl Into<Element> for &parts::Element {
    fn into(self) -> Element {
        match self {
            parts::Element::H => Element::H,
            parts::Element::He => Element::He,
            parts::Element::Li => Element::Li,
            parts::Element::Be => Element::Be,
            parts::Element::B => Element::B,
            parts::Element::C => Element::C,
            parts::Element::N => Element::N,
            parts::Element::O => Element::O,
            parts::Element::F => Element::F,
            parts::Element::Ne => Element::Ne,
            parts::Element::Na => Element::Na,
            parts::Element::Mg => Element::Mg,
            parts::Element::Al => Element::Al,
            parts::Element::Si => Element::Si,
            parts::Element::P => Element::P,
            parts::Element::S => Element::S,
            parts::Element::Cl => Element::Cl,
            parts::Element::Ar => Element::Ar,
            parts::Element::K => Element::K,
            parts::Element::Ca => Element::Ca,
            parts::Element::Sc => Element::Sc,
            parts::Element::Ti => Element::Ti,
            parts::Element::V => Element::V,
            parts::Element::Cr => Element::Cr,
            parts::Element::Mn => Element::Mn,
            parts::Element::Fe => Element::Fe,
            parts::Element::Co => Element::Co,
            parts::Element::Ni => Element::Ni,
            parts::Element::Cu => Element::Cu,
            parts::Element::Zn => Element::Zn,
            parts::Element::Ga => Element::Ga,
            parts::Element::Ge => Element::Ge,
            parts::Element::As => Element::As,
            parts::Element::Se => Element::Se,
            parts::Element::Br => Element::Br,
            parts::Element::Kr => Element::Kr,
            parts::Element::Rb => Element::Rb,
            parts::Element::Sr => Element::Sr,
            parts::Element::Y => Element::Y,
            parts::Element::Zr => Element::Zr,
            parts::Element::Nb => Element::Nb,
            parts::Element::Mo => Element::Mo,
            parts::Element::Tc => Element::Tc,
            parts::Element::Ru => Element::Ru,
            parts::Element::Rh => Element::Rh,
            parts::Element::Pd => Element::Pd,
            parts::Element::Ag => Element::Ag,
            parts::Element::Cd => Element::Cd,
            parts::Element::In => Element::In,
            parts::Element::Sn => Element::Sn,
            parts::Element::Sb => Element::Sb,
            parts::Element::Te => Element::Te,
            parts::Element::I => Element::I,
            parts::Element::Xe => Element::Xe,
            parts::Element::Cs => Element::Cs,
            parts::Element::Ba => Element::Ba,
            parts::Element::La => Element::La,
            parts::Element::Ce => Element::Ce,
            parts::Element::Pr => Element::Pr,
            parts::Element::Nd => Element::Nd,
            parts::Element::Pm => Element::Pm,
            parts::Element::Sm => Element::Sm,
            parts::Element::Eu => Element::Eu,
            parts::Element::Gd => Element::Gd,
            parts::Element::Tb => Element::Tb,
            parts::Element::Dy => Element::Dy,
            parts::Element::Ho => Element::Ho,
            parts::Element::Er => Element::Er,
            parts::Element::Tm => Element::Tm,
            parts::Element::Yb => Element::Yb,
            parts::Element::Lu => Element::Lu,
            parts::Element::Hf => Element::Hf,
            parts::Element::Ta => Element::Ta,
            parts::Element::W => Element::W,
            parts::Element::Re => Element::Re,
            parts::Element::Os => Element::Os,
            parts::Element::Ir => Element::Ir,
            parts::Element::Pt => Element::Pt,
            parts::Element::Au => Element::Au,
            parts::Element::Hg => Element::Hg,
            parts::Element::Tl => Element::Tl,
            parts::Element::Pb => Element::Pb,
            parts::Element::Bi => Element::Bi,
            parts::Element::Po => Element::Po,
            parts::Element::At => Element::At,
            parts::Element::Rn => Element::Rn,
            parts::Element::Fr => Element::Fr,
            parts::Element::Ra => Element::Ra,
            parts::Element::Ac => Element::Ac,
            parts::Element::Th => Element::Th,
            parts::Element::Pa => Element::Pa,
            parts::Element::U => Element::U,
            parts::Element::Np => Element::Np,
            parts::Element::Pu => Element::Pu,
            parts::Element::Am => Element::Am,
            parts::Element::Cm => Element::Cm,
            parts::Element::Bk => Element::Bk,
            parts::Element::Cf => Element::Cf,
            parts::Element::Es => Element::Es,
            parts::Element::Fm => Element::Fm,
            parts::Element::Md => Element::Md,
            parts::Element::No => Element::No,
            parts::Element::Lr => Element::Lr,
            parts::Element::Rf => Element::Rf,
            parts::Element::Db => Element::Db,
            parts::Element::Sg => Element::Sg,
            parts::Element::Bh => Element::Bh,
            parts::Element::Hs => Element::Hs,
            parts::Element::Mt => Element::Mt,
            parts::Element::Ds => Element::Ds,
            parts::Element::Rg => Element::Rg,
            parts::Element::Cn => Element::Cn,
            parts::Element::Nh => Element::Nh,
            parts::Element::Fl => Element::Fl,
            parts::Element::Mc => Element::Mc,
            parts::Element::Lv => Element::Lv,
            parts::Element::Ts => Element::Ts,
            parts::Element::Og => Element::Og,
        }
    }
}

impl Element {
    pub fn valence_electrons(&self) -> u8 {
        let mut result = self.atomic_number();

        if let Some(core) = self.core() {
            result -= core.atomic_number();
        }

        result.into()
    }

    pub fn atomic_number(&self) -> u8 {
        match self {
            Element::H => 1,
            Element::He => 2,
            Element::Li => 3,
            Element::Be => 4,
            Element::B => 5,
            Element::C => 6,
            Element::N => 7,
            Element::O => 8,
            Element::F => 9,
            Element::Ne => 10,
            Element::Na => 11,
            Element::Mg => 12,
            Element::Al => 13,
            Element::Si => 14,
            Element::P => 15,
            Element::S => 16,
            Element::Cl => 17,
            Element::Ar => 18,
            Element::K => 19,
            Element::Ca => 20,
            Element::Sc => 21,
            Element::Ti => 22,
            Element::V => 23,
            Element::Cr => 24,
            Element::Mn => 25,
            Element::Fe => 26,
            Element::Co => 27,
            Element::Ni => 28,
            Element::Cu => 29,
            Element::Zn => 30,
            Element::Ga => 31,
            Element::Ge => 32,
            Element::As => 33,
            Element::Se => 34,
            Element::Br => 35,
            Element::Kr => 36,
            Element::Rb => 37,
            Element::Sr => 38,
            Element::Y => 39,
            Element::Zr => 40,
            Element::Nb => 41,
            Element::Mo => 42,
            Element::Tc => 43,
            Element::Ru => 44,
            Element::Rh => 45,
            Element::Pd => 46,
            Element::Ag => 47,
            Element::Cd => 48,
            Element::In => 49,
            Element::Sn => 50,
            Element::Sb => 51,
            Element::Te => 52,
            Element::I => 53,
            Element::Xe => 54,
            Element::Cs => 55,
            Element::Ba => 56,
            Element::La => 57,
            Element::Ce => 58,
            Element::Pr => 59,
            Element::Nd => 60,
            Element::Pm => 61,
            Element::Sm => 62,
            Element::Eu => 63,
            Element::Gd => 64,
            Element::Tb => 65,
            Element::Dy => 66,
            Element::Ho => 67,
            Element::Er => 68,
            Element::Tm => 69,
            Element::Yb => 70,
            Element::Lu => 71,
            Element::Hf => 72,
            Element::Ta => 73,
            Element::W => 74,
            Element::Re => 75,
            Element::Os => 76,
            Element::Ir => 77,
            Element::Pt => 78,
            Element::Au => 79,
            Element::Hg => 80,
            Element::Tl => 81,
            Element::Pb => 82,
            Element::Bi => 83,
            Element::Po => 84,
            Element::At => 85,
            Element::Rn => 86,
            Element::Fr => 87,
            Element::Ra => 88,
            Element::Ac => 89,
            Element::Th => 90,
            Element::Pa => 91,
            Element::U => 92,
            Element::Np => 93,
            Element::Pu => 94,
            Element::Am => 95,
            Element::Cm => 96,
            Element::Bk => 97,
            Element::Cf => 98,
            Element::Es => 99,
            Element::Fm => 100,
            Element::Md => 101,
            Element::No => 102,
            Element::Lr => 103,
            Element::Rf => 104,
            Element::Db => 105,
            Element::Sg => 106,
            Element::Bh => 107,
            Element::Hs => 108,
            Element::Mt => 109,
            Element::Ds => 110,
            Element::Rg => 111,
            Element::Cn => 112,
            Element::Nh => 113,
            Element::Fl => 114,
            Element::Mc => 115,
            Element::Lv => 116,
            Element::Ts => 117,
            Element::Og => 118,
        }
    }

    fn core(&self) -> Option<Self> {
        if self.atomic_number() < 3 {
            None
        } else if self.atomic_number() < 11 {
            Some(Element::He)
        } else if self.atomic_number() < 19 {
            Some(Element::Ne)
        } else if self.atomic_number() < 37 {
            Some(Element::Ar)
        } else if self.atomic_number() < 55 {
            Some(Element::Kr)
        } else if self.atomic_number() < 87 {
            Some(Element::Xe)
        } else {
            Some(Element::Rn)
        }
    }
}

#[cfg(test)]
mod valence_electrons {
    use super::*;

    #[test]
    fn bromine() {
        let bromine = Element::Br;

        assert_eq!(bromine.valence_electrons(), 17)
    }

    #[test]
    fn iodine() {
        let iodine = Element::I;

        assert_eq!(iodine.valence_electrons(), 17)
    }

    #[test]
    fn astatine() {
        let astatine = Element::At;

        assert_eq!(astatine.valence_electrons(), 17 + 14)
    }

    #[test]
    fn tennesine() {
        let tennesine = Element::Ts;

        assert_eq!(tennesine.valence_electrons(), 17 + 14)
    }
}
