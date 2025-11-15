use crate::cmp::bignum_ct_eq;
use crate::bignum::BigNum;
use crate::limbs::{MASK_120, LIMB_BITS};
use crate::montgomery::MontCtx;
use crate::field::FieldElt;
use crate::ec_types::{ECParams, AffinePoint, JacobianPoint};
use crate::ec_jacobian::{affine_to_jacobian, point_double, point_add_mixed, jacobian_ct_select};

fn fe_zero_local<let N: u32, let P: u32>() -> FieldElt<N, P> {
    FieldElt::<N, P> {
        mont: BigNum::<N>::from_limbs([0u128; N]),
    }
}

fn fe_one_local<let N: u32, let P: u32>() -> FieldElt<N, P> {
    let mut limbs: [u128; N] = [0u128; N];
    limbs[0] = 1u128;
    FieldElt::<N, P> {
        mont: BigNum::<N>::from_limbs(limbs),
    }
}


pub fn scalar_mul<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, params: ECParams<N, P>,
    q: AffinePoint<N, P>, k: BigNum<N>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    let zero_bn = BigNum::<N>::from_limbs([0u128; N]);
    let one_bn = {
        let mut limbs: [u128; N] = [0u128; N];
        limbs[0] = 1u128;
        BigNum::<N>::from_limbs(limbs)
    };
    let two_bn = {
        let mut limbs: [u128; N] = [0u128; N];
        limbs[0] = 2u128;
        BigNum::<N>::from_limbs(limbs)
    };

    let k_is_zero = bignum_ct_eq::<N>(k, zero_bn);
    let k_is_one  = bignum_ct_eq::<N>(k, one_bn);
    let k_is_two  = bignum_ct_eq::<N>(k, two_bn);

    let q_is_inf: u128 = q.is_inf;

    let inf: JacobianPoint<N, P> = JacobianPoint::<N, P> {
        X: fe_zero_local::<N, P>(),
        Y: fe_one_local::<N, P>(),
        Z: fe_zero_local::<N, P>(),
    };

    let (ctx1, qj) = affine_to_jacobian::<N, P>(ctx, q);
    let (ctx2, two_qj) = point_double::<N, P>(ctx1, params, qj);

    let mut res = inf;

    res = jacobian_ct_select::<N, P>(res, qj, k_is_one);

    res = jacobian_ct_select::<N, P>(res, two_qj, k_is_two);

    let need_inf = q_is_inf | k_is_zero;
    res = jacobian_ct_select::<N, P>(res, inf, need_inf);

    (ctx2, res)
}

