use crate::montgomery::MontCtx;
use crate::limbs::{MASK_120, LIMB_BITS};
use crate::bignum::BigNum;
use crate::field::FieldElt;
use crate::ec_types::{ECParams, AffinePoint, JacobianPoint};
use crate::ec_jacobian::{affine_to_jacobian, point_double, point_add_mixed, jacobian_ct_select};

fn fe_zero<let N: u32, let P: u32>() -> FieldElt<N, P> {
    FieldElt::<N, P> { mont: BigNum::<N>::from_limbs([0u128; N]) }
}
fn fe_one<let N: u32, let P: u32>() -> FieldElt<N, P> {
    let mut limbs: [u128; N] = [0u128; N];
    limbs[0] = 1u128;
    FieldElt::<N, P> { mont: BigNum::<N>::from_limbs(limbs) }
}

pub fn scalar_mul<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, params: ECParams<N, P>,
    q: AffinePoint<N, P>, k: BigNum<N>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    // _qj is not used directly, kept for symmetry
    let (mut ctx1, _qj) = affine_to_jacobian::<N, P>(ctx, q);

    let mut R = JacobianPoint::<N, P> { X: fe_zero::<N, P>(), Y: fe_one::<N, P>(), Z: fe_zero::<N, P>() };

    for i in 0..N {
        let limb = k.limbs[i] & MASK_120;
        for bitpos in 0..LIMB_BITS {
            let (c1, R2) = point_double::<N, P>(ctx1, params, R);
            let bit: u128 = (limb >> (bitpos as u128)) & 1u128;
            let (c2, R_plus) = point_add_mixed::<N, P>(c1, params, R2, q);
            R = jacobian_ct_select::<N, P>(R2, R_plus, bit);
            ctx1 = c2;
        }
    }
    (ctx1, R)
}

