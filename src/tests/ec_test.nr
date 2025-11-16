use crate::montgomery::MontCtx;
use crate::bignum::BigNum;
use crate::field::{FieldElt, field_from_big};
use crate::ec_types::{ECParams, AffinePoint, JacobianPoint};
use crate::ec_jacobian::{point_add_mixed, point_double, jacobian_equal_proj, affine_to_jacobian};
use crate::scalar_mul::{scalar_mul, scalar_mul_trivial_point};
use crate::cmp::bignum_ct_eq;
use crate::scalar_mul::scalar_mul_small;


use crate::params::{test_p_N2, test_r2_N2, test_ninv};
fn make_fe<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>,
    lo: u128,
    hi: u128,
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let mut limbs: [u128; N] = [0u128; N];
    limbs[0] = lo;
    if N > 1u32 {
        limbs[1] = hi;
    }
    let bn = BigNum::<N>::from_limbs(limbs);
    field_from_big::<N, P>(ctx, bn)
}

fn fe_zero<let N: u32, let P: u32>() -> FieldElt<N, P> {
    FieldElt::<N, P> {
        mont: BigNum::<N>::from_limbs([0u128; N]),
    }
}

fn fe_one<let N: u32, let P: u32>() -> FieldElt<N, P> {
    let mut limbs: [u128; N] = [0u128; N];
    limbs[0] = 1u128;
    FieldElt::<N, P> {
        mont: BigNum::<N>::from_limbs(limbs),
    }
}

#[test]
fn ec_add_identity_mixed() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };

    let (ctx1, a_fe) = make_fe::<2, 4>(ctx0, 0u128, 0u128);
    let (ctx2, b_fe) = make_fe::<2, 4>(ctx1, 7u128, 0u128);
    let params = ECParams::<2, 4> { a: a_fe, b: b_fe };

    let inf = JacobianPoint::<2, 4> { X: fe_zero::<2, 4>(), Y: fe_one::<2, 4>(), Z: fe_zero::<2, 4>() };

    let q_inf = AffinePoint::<2, 4> { x: fe_zero::<2, 4>(), y: fe_one::<2, 4>(), is_inf: 1u128 };

    // inf + Q_inf = inf
    let (ctx3, sum1) = point_add_mixed::<2, 4>(ctx2, params, inf, q_inf);
    let (_ctx4, eq1) = jacobian_equal_proj::<2, 4>(ctx3, sum1, inf);
    assert(eq1 == 1u128);
}

#[test]
fn ec_double_of_infinity_is_infinity() {
    let ctx0 = MontCtx::<2, 4> {
        p: test_p_N2(),
        r2: test_r2_N2(),
        ninv: test_ninv(),
    };

    let (ctx1, a_fe) = make_fe::<2, 4>(ctx0, 0u128, 0u128);
    let (ctx2, b_fe) = make_fe::<2, 4>(ctx1, 7u128, 0u128);
    let params = ECParams::<2, 4> { a: a_fe, b: b_fe };

    let inf = JacobianPoint::<2, 4> {
        X: fe_zero::<2, 4>(),
        Y: fe_one::<2, 4>(),
        Z: fe_zero::<2, 4>(),
    };

    let (_ctx3, twoO) = point_double::<2, 4>(ctx2, params, inf);
    let (_ctx4, eq) = jacobian_equal_proj::<2, 4>(_ctx3, twoO, inf);
    assert(eq == 1u128);
}

#[test]
fn ec_scalar_mul_trivial_ks() {
    println("Start");
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };

    let (ctx1, a_fe) = make_fe::<2, 4>(ctx0, 0u128, 0u128);
    let (ctx2, b_fe) = make_fe::<2, 4>(ctx1, 7u128, 0u128);
    let params = ECParams::<2, 4> { a: a_fe, b: b_fe };

    let q = AffinePoint::<2, 4> { x: fe_zero::<2, 4>(), y: fe_one::<2, 4>(), is_inf: 1u128 };

    let k0 = BigNum::<2>::from_limbs([0u128, 0u128]);
    let (ctx3, r0) = scalar_mul_trivial_point::<2, 4>(ctx2, params, q, k0);
    println("scalar_mul: fast path (Q=inf or k=0)");
    let zero_bn = BigNum::<2>::from_limbs([0u128, 0u128]);
    assert(bignum_ct_eq::<2>(r0.Z.mont, zero_bn) == 1u128);

    let k1 = BigNum::<2>::from_limbs([1u128, 0u128]);
    let (ctx4, r1) = scalar_mul_trivial_point::<2, 4>(ctx3, params, q, k1);
    println("scalar_mul: fast path (k=1)");
    assert(bignum_ct_eq::<2>(r1.Z.mont, zero_bn) == 1u128);

    let k2 = BigNum::<2>::from_limbs([2u128, 0u128]);
    let (_ctx5, r2) = scalar_mul_trivial_point::<2, 4>(ctx4, params, q, k2);
    println("scalar_mul: fast path (Q=inf covers all k)");
    assert(bignum_ct_eq::<2>(r2.Z.mont, zero_bn) == 1u128);
}


#[test]
fn ec_scalar_mul_matches_repeated_add() {
    let ctx0 = MontCtx::<2, 4> {
        p: test_p_N2(),
        r2: test_r2_N2(),
        ninv: test_ninv(),
    };

    let (ctx1, a_fe) = make_fe::<2, 4>(ctx0, 0u128, 0u128); // a = 0
    let (ctx2, b_fe) = make_fe::<2, 4>(ctx1, 7u128, 0u128); // b = 7
    let params = ECParams::<2, 4> { a: a_fe, b: b_fe };

    let (ctx3, qx) = make_fe::<2, 4>(ctx2, 5u128, 0u128);
    let (ctx4, qy) = make_fe::<2, 4>(ctx3, 9u128, 0u128);
    let q = AffinePoint::<2, 4> {
        x: qx,
        y: qy,
        is_inf: 0u128,
    };

    let k3 = BigNum::<2>::from_limbs([3u128, 0u128]);

    let (ctx5, r_mul) = scalar_mul_small::<2, 4>(ctx4, params, q, k3);

    let (ctx6, qj) = affine_to_jacobian::<2, 4>(ctx5, q);
    let (ctx7, two_q) = point_double::<2, 4>(ctx6, params, qj);
    let (ctx8, three_q) = point_add_mixed::<2, 4>(ctx7, params, two_q, q);

    let (_ctx9, eq) = jacobian_equal_proj::<2, 4>(ctx8, r_mul, three_q);
    assert(eq == 1u128);
}

