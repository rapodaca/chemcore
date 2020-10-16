use purr::mol;

#[derive(Clone, Copy, Eq, Hash, PartialEq, Debug)]
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

impl Default for Element {
    fn default() -> Self {
        Self::C
    }
}

impl From<mol::Element> for Element {
    fn from(element: mol::Element) -> Self {
        match element {
            mol::Element::H  => Element::H,
            mol::Element::He => Element::He,
            mol::Element::Li => Element::Li,
            mol::Element::Be => Element::Be,
            mol::Element::B =>  Element::B,
            mol::Element::C =>  Element::C,
            mol::Element::N =>  Element::N,
            mol::Element::O =>  Element::O,
            mol::Element::F =>  Element::F,
            mol::Element::Ne => Element::Ne,
            mol::Element::Na => Element::Na,
            mol::Element::Mg => Element::Mg,
            mol::Element::Al => Element::Al,
            mol::Element::Si => Element::Si,
            mol::Element::P =>  Element::P,
            mol::Element::S =>  Element::S,
            mol::Element::Cl => Element::Cl,
            mol::Element::Ar => Element::Ar,
            mol::Element::K =>  Element::K,
            mol::Element::Ca => Element::Ca,
            mol::Element::Sc => Element::Sc,
            mol::Element::Ti => Element::Ti,
            mol::Element::V =>  Element::V,
            mol::Element::Cr => Element::Cr,
            mol::Element::Mn => Element::Mn,
            mol::Element::Fe => Element::Fe,
            mol::Element::Co => Element::Co,
            mol::Element::Ni => Element::Ni,
            mol::Element::Cu => Element::Cu,
            mol::Element::Zn => Element::Zn,
            mol::Element::Ga => Element::Ga,
            mol::Element::Ge => Element::Ge,
            mol::Element::As => Element::As,
            mol::Element::Se => Element::Se,
            mol::Element::Br => Element::Br,
            mol::Element::Kr => Element::Kr,
            mol::Element::Rb => Element::Rb,
            mol::Element::Sr => Element::Sr,
            mol::Element::Y =>  Element::Y,
            mol::Element::Zr => Element::Zr,
            mol::Element::Nb => Element::Nb,
            mol::Element::Mo => Element::Mo,
            mol::Element::Tc => Element::Tc,
            mol::Element::Ru => Element::Ru,
            mol::Element::Rh => Element::Rh,
            mol::Element::Pd => Element::Pd,
            mol::Element::Ag => Element::Ag,
            mol::Element::Cd => Element::Cd,
            mol::Element::In => Element::In,
            mol::Element::Sn => Element::Sn,
            mol::Element::Sb => Element::Sb,
            mol::Element::Te => Element::Te,
            mol::Element::I =>  Element::I,
            mol::Element::Xe => Element::Xe,
            mol::Element::Cs => Element::Cs,
            mol::Element::Ba => Element::Ba,
            mol::Element::La => Element::La,
            mol::Element::Ce => Element::Ce,
            mol::Element::Pr => Element::Pr,
            mol::Element::Nd => Element::Nd,
            mol::Element::Pm => Element::Pm,
            mol::Element::Sm => Element::Sm,
            mol::Element::Eu => Element::Eu,
            mol::Element::Gd => Element::Gd,
            mol::Element::Tb => Element::Tb,
            mol::Element::Dy => Element::Dy,
            mol::Element::Ho => Element::Ho,
            mol::Element::Er => Element::Er,
            mol::Element::Tm => Element::Tm,
            mol::Element::Yb => Element::Yb,
            mol::Element::Lu => Element::Lu,
            mol::Element::Hf => Element::Hf,
            mol::Element::Ta => Element::Ta,
            mol::Element::W =>  Element::W,
            mol::Element::Re => Element::Re,
            mol::Element::Os => Element::Os,
            mol::Element::Ir => Element::Ir,
            mol::Element::Pt => Element::Pt,
            mol::Element::Au => Element::Au,
            mol::Element::Hg => Element::Hg,
            mol::Element::Tl => Element::Tl,
            mol::Element::Pb => Element::Pb,
            mol::Element::Bi => Element::Bi,
            mol::Element::Po => Element::Po,
            mol::Element::At => Element::At,
            mol::Element::Rn => Element::Rn,
            mol::Element::Fr => Element::Fr,
            mol::Element::Ra => Element::Ra,
            mol::Element::Ac => Element::Ac,
            mol::Element::Th => Element::Th,
            mol::Element::Pa => Element::Pa,
            mol::Element::U =>  Element::U,
            mol::Element::Np => Element::Np,
            mol::Element::Pu => Element::Pu,
            mol::Element::Am => Element::Am,
            mol::Element::Cm => Element::Cm,
            mol::Element::Bk => Element::Bk,
            mol::Element::Cf => Element::Cf,
            mol::Element::Es => Element::Es,
            mol::Element::Fm => Element::Fm,
            mol::Element::Md => Element::Md,
            mol::Element::No => Element::No,
            mol::Element::Lr => Element::Lr,
            mol::Element::Rf => Element::Rf,
            mol::Element::Db => Element::Db,
            mol::Element::Sg => Element::Sg,
            mol::Element::Bh => Element::Bh,
            mol::Element::Hs => Element::Hs,
            mol::Element::Mt => Element::Mt,
            mol::Element::Ds => Element::Ds,
            mol::Element::Rg => Element::Rg,
            mol::Element::Cn => Element::Cn,
            mol::Element::Nh => Element::Nh,
            mol::Element::Fl => Element::Fl,
            mol::Element::Mc => Element::Mc,
            mol::Element::Lv => Element::Lv,
            mol::Element::Ts => Element::Ts,
            mol::Element::Og => Element::Og
        }
    }
}

impl Element {
    pub fn valence_electrons(&self) -> u16 {
        let mut result = self.atomic_number();

        if let Some(core) = self.core() {
            result -= core.atomic_number();
        }

        result
    }

    pub fn atomic_number(&self) -> u16 {
        match self {
            Element::H  => 1,
            Element::He => 2,
            Element::Li => 3,
            Element::Be => 4,
            Element::B =>  5,
            Element::C =>  6,
            Element::N =>  7,
            Element::O =>  8,
            Element::F =>  9,
            Element::Ne => 10,
            Element::Na => 11,
            Element::Mg => 12,
            Element::Al => 13,
            Element::Si => 14,
            Element::P =>  15,
            Element::S =>  16,
            Element::Cl => 17,
            Element::Ar => 18,
            Element::K =>  19,
            Element::Ca => 20,
            Element::Sc => 21,
            Element::Ti => 22,
            Element::V =>  23,
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
            Element::Y =>  39,
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
            Element::I =>  53,
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
            Element::W =>  74,
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
            Element::U =>  92,
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
            Element::Og => 118
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

        assert_eq!(tennesine.valence_electrons(),  17 + 14)
    }
}