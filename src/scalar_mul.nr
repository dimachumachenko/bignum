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

    let q_is_inf: u128 = q.is_inf;

    let zero_bn = BigNum::<N>::from_limbs([0u128; N]);
    let k_is_zero: u128 = bignum_ct_eq::<N>(k, zero_bn);

    let mut one_arr: [u128; N] = [0u128; N];
    one_arr[0] = 1u128;
    let one_bn = BigNum::<N>::from_limbs(one_arr);
    let k_is_one: u128 = bignum_ct_eq::<N>(k, one_bn);

    let inf: JacobianPoint<N, P> = JacobianPoint::<N, P> {
        X: fe_zero_local::<N, P>(),
        Y: fe_one_local::<N, P>(),
        Z: fe_zero_local::<N, P>(),
    };

    let (mut ctx1, qj) = affine_to_jacobian::<N, P>(ctx, q);

    let mut R = inf;
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

    let res_k1 = jacobian_ct_select::<N, P>(R, qj, k_is_one);
    let need_inf: u128 = q_is_inf | k_is_zero;
    let res = jacobian_ct_select::<N, P>(res_k1, inf, need_inf);

    (ctx1, res)
}

pub fn scalar_mul_trivial_point<let N: u32, let P: u32>(
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


pub fn scalar_mul_small<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, params: ECParams<N, P>,
    q: AffinePoint<N, P>, k: BigNum<N>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    let q_is_inf: u128 = q.is_inf & 1u128;

    let zero_bn = BigNum::<N>::from_limbs([0u128; N]);

    let mut one_arr: [u128; N] = [0u128; N];
    one_arr[0] = 1u128;
    let one_bn = BigNum::<N>::from_limbs(one_arr);

    let mut two_arr: [u128; N] = [0u128; N];
    two_arr[0] = 2u128;
    let two_bn = BigNum::<N>::from_limbs(two_arr);

    let k_is_zero: u128 = bignum_ct_eq::<N>(k, zero_bn);
    let k_is_one: u128  = bignum_ct_eq::<N>(k, one_bn);
    let k_is_two: u128  = bignum_ct_eq::<N>(k, two_bn);

    let inf: JacobianPoint<N, P> = JacobianPoint::<N, P> {
        X: fe_zero_local::<N, P>(),
        Y: fe_one_local::<N, P>(),
        Z: fe_zero_local::<N, P>(),
    };

    let (ctx1, qj_raw) = affine_to_jacobian::<N, P>(ctx, q);

    let qj = jacobian_ct_select::<N, P>(qj_raw, inf, q_is_inf);

    let (ctx2, two_q) = point_double::<N, P>(ctx1, params, qj);
    let (ctx3, three_q) = point_add_mixed::<N, P>(ctx2, params, two_q, q);

    let mut res = inf;

    res = jacobian_ct_select::<N, P>(res, qj, k_is_one);

    res = jacobian_ct_select::<N, P>(res, two_q, k_is_two);

    let any_small = k_is_zero | k_is_one | k_is_two;
    let is_other = 1u128 - (any_small & 1u128);
    res = jacobian_ct_select::<N, P>(res, three_q, is_other);

    res = jacobian_ct_select::<N, P>(res, inf, k_is_zero);

    res = jacobian_ct_select::<N, P>(res, inf, q_is_inf);

    (ctx3, res)
}
