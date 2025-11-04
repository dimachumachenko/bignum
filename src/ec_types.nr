use crate::field::FieldElt;

pub struct ECParams<let N: u32, let P: u32> {
    pub a: FieldElt<N, P>,
    pub b: FieldElt<N, P>,
}

pub struct AffinePoint<let N: u32, let P: u32> {
    pub x: FieldElt<N, P>,
    pub y: FieldElt<N, P>,
    pub is_inf: u128, 
}

pub struct JacobianPoint<let N: u32, let P: u32> {
    pub X: FieldElt<N, P>,
    pub Y: FieldElt<N, P>,
    pub Z: FieldElt<N, P>,
}

