use itertools::Itertools;

const HYDROPHOBIC: &'static str = "QNHSTYC";
// const POLAR: &'static str = "AILMFVPG";

#[derive(PartialEq)]
pub enum Amino {
    Hydrophobic,
    Polar
}

pub struct Residue {
    pub kind: Amino,
    pub loc: [i16; 3],
}

impl Residue {
    pub fn is_adjacent(&self, other: &Residue) -> bool {
        ((self.loc[0]-other.loc[0]).abs() == 1) ||
        ((self.loc[1]-other.loc[1]).abs() == 1) ||
        ((self.loc[2]-other.loc[2]).abs() == 1)
    }
}

pub struct Protein {
    pub residues: Vec<Residue>,
}

impl Protein {
    pub fn new(seq: &str) -> Self {
        Self {
            residues: seq.chars().enumerate().map(|(i, c)| Residue {
                kind: if HYDROPHOBIC.contains(c) { Amino::Hydrophobic } else { Amino::Polar },
                loc: [0, i as i16, 0]
            }).collect()
        }
    }

    pub fn energy(&self) -> f64 {
        let candidates = self.residues.iter().filter(|r| r.kind == Amino::Hydrophobic);
        candidates.combinations(2).map(|c| if c[0].is_adjacent(c[1]) {1.} else {0.}).sum()
    }
}

pub trait Strategy {
    fn fold(prot: &mut Protein) -> f64;
}

pub struct MonteCarlo;

impl MonteCarlo {
    fn attempt_move(prot: &mut Protein) {
        todo!()
    }
}

impl Strategy for MonteCarlo {
    fn fold(prot: &mut Protein) -> f64 {
        
        prot.energy()
    }
}

