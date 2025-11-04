use crate::montgomery::MontCtx;
use crate::field::{FieldElt, field_add, field_sub, field_mul, field_square};
use crate::cmp::bignum_ct_eq;
use crate::bignum::BigNum;
use crate::ec_types::{ECParams, AffinePoint, JacobianPoint};

fn fe_zero<let N: u32, let P: u32>() -> FieldElt<N, P> {
    FieldElt::<N, P> { mont: BigNum::<N>::from_limbs([0u128; N]) }
}

fn fe_one<let N: u32, let P: u32>() -> FieldElt<N, P> {
    let mut limbs: [u128; N] = [0u128; N];
    limbs[0] = 1u128;
    FieldElt::<N, P> { mont: BigNum::<N>::from_limbs(limbs) }
}

pub fn affine_to_jacobian<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: AffinePoint<N, P>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    let x_m = BigNum::<N>::ct_select(fe_zero::<N, P>().mont, a.x.mont, ((1u128 - a.is_inf) & 1u128));
    let y_m = BigNum::<N>::ct_select(fe_one::<N, P>().mont, a.y.mont, ((1u128 - a.is_inf) & 1u128));
    let z_m = BigNum::<N>::ct_select(fe_zero::<N, P>().mont, fe_one::<N, P>().mont, ((1u128 - a.is_inf) & 1u128));
    (ctx, JacobianPoint::<N, P> {
        X: FieldElt::<N, P> { mont: x_m },
        Y: FieldElt::<N, P> { mont: y_m },
        Z: FieldElt::<N, P> { mont: z_m },
    })
}

pub fn jacobian_is_inf<let N: u32, let P: u32>(p: JacobianPoint<N, P>) -> u128 {
    bignum_ct_eq(p.Z.mont, BigNum::<N>::from_limbs([0u128; N]))
}

pub fn jacobian_ct_select<let N: u32, let P: u32>(
    a: JacobianPoint<N, P>, b: JacobianPoint<N, P>, flag: u128
) -> JacobianPoint<N, P> {
    let X = FieldElt::<N, P> { mont: BigNum::<N>::ct_select(a.X.mont, b.X.mont, flag) };
    let Y = FieldElt::<N, P> { mont: BigNum::<N>::ct_select(a.Y.mont, b.Y.mont, flag) };
    let Z = FieldElt::<N, P> { mont: BigNum::<N>::ct_select(a.Z.mont, b.Z.mont, flag) };
    JacobianPoint::<N, P> { X, Y, Z }
}

pub fn jacobian_equal_proj<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, p: JacobianPoint<N, P>, q: JacobianPoint<N, P>
) -> (MontCtx<N, P>, u128) {
    let (ctx1, z1z1) = field_square::<N, P>(ctx, p.Z);
    let (ctx2, z2z2) = field_square::<N, P>(ctx1, q.Z);
    let (ctx3, z1z1z1) = field_mul::<N, P>(ctx2, z1z1, p.Z);
    let (ctx4, z2z2z2) = field_mul::<N, P>(ctx3, z2z2, q.Z);

    let (ctx5, x1n) = field_mul::<N, P>(ctx4, p.X, z2z2);
    let (ctx6, x2n) = field_mul::<N, P>(ctx5, q.X, z1z1);
    let (ctx7, y1n) = field_mul::<N, P>(ctx6, p.Y, z2z2z2);
    let (ctx8, y2n) = field_mul::<N, P>(ctx7, q.Y, z1z1z1);

    let eqx: u128 = bignum_ct_eq(x1n.mont, x2n.mont);
    let eqy: u128 = bignum_ct_eq(y1n.mont, y2n.mont);
    (ctx8, (eqx & eqy) & 1u128)
}

pub fn point_double<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, params: ECParams<N, P>, p: JacobianPoint<N, P>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    let is_inf = jacobian_is_inf::<N, P>(p);

    let (ctx1, A) = field_square::<N, P>(ctx, p.X);
    let (ctx2, B) = field_square::<N, P>(ctx1, p.Y);
    let (ctx3, C) = field_square::<N, P>(ctx2, B);

    let (ctx4, Z2) = field_square::<N, P>(ctx3, p.Z);
    let (ctx5, Z4) = field_square::<N, P>(ctx4, Z2);

    let (ctx6, threeA) = {
        let (c1, t1) = field_add::<N, P>(ctx5, A, A);
        let (c2, t2) = field_add::<N, P>(c1, t1, A);
        (c2, t2)
    };
    let (ctx7, aZ4) = field_mul::<N, P>(ctx6, params.a, Z4);
    let (ctx8, E) = field_add::<N, P>(ctx7, threeA, aZ4);

    let (ctx9, XplusB) = field_add::<N, P>(ctx8, p.X, B);
    let (ctx10, XplusB2) = field_square::<N, P>(ctx9, XplusB);
    let (ctx11, tD1) = field_sub::<N, P>(ctx10, XplusB2, A);
    let (ctx12, tD2) = field_sub::<N, P>(ctx11, tD1, C);
    let (ctx13, D) = field_add::<N, P>(ctx12, tD2, tD2);

    let (ctx14, F) = field_square::<N, P>(ctx13, E);
    let (ctx15, twoD) = field_add::<N, P>(ctx14, D, D);
    let (ctx16, X3) = field_sub::<N, P>(ctx15, F, twoD);

    let (ctx17, DmX3) = field_sub::<N, P>(ctx16, D, X3);
    let (ctx18, E_DmX3) = field_mul::<N, P>(ctx17, E, DmX3);
    let (ctx19, eightC) = {
        let (c1, t1) = field_add::<N, P>(ctx18, C, C);
        let (c2, t2) = field_add::<N, P>(c1, t1, t1);
        field_add::<N, P>(c2, t2, t2)
    };
    let (ctx20, Y3) = field_sub::<N, P>(ctx19, E_DmX3, eightC);

    let (ctx21, YZ) = field_mul::<N, P>(ctx20, p.Y, p.Z);
    let (ctx22, Z3) = field_add::<N, P>(ctx21, YZ, YZ);

    let infP = JacobianPoint::<N, P> { X: fe_zero::<N, P>(), Y: fe_one::<N, P>(), Z: fe_zero::<N, P>() };
    let R = JacobianPoint::<N, P> { X: X3, Y: Y3, Z: Z3 };
    let out = jacobian_ct_select::<N, P>(infP, R, ((1u128 - is_inf) & 1u128));
    (ctx22, out)
}

pub fn point_add_mixed<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, _params: ECParams<N, P>,
    p: JacobianPoint<N, P>, q: AffinePoint<N, P>
) -> (MontCtx<N, P>, JacobianPoint<N, P>) {
    let q_is_inf = q.is_inf & 1u128;
    let p_is_inf = jacobian_is_inf::<N, P>(p);

    let (ctx0, qj) = affine_to_jacobian::<N, P>(ctx, q);

    let (ctx1, Z1Z1) = field_square::<N, P>(ctx0, p.Z);
    let (ctx2, U2) = field_mul::<N, P>(ctx1, qj.X, Z1Z1);
    let (ctx3, Z1_cubed) = field_mul::<N, P>(ctx2, Z1Z1, p.Z);
    let (ctx4, S2) = field_mul::<N, P>(ctx3, qj.Y, Z1_cubed);

    let U1 = p.X;
    let S1 = p.Y;

    let (ctx5, H) = field_sub::<N, P>(ctx4, U2, U1);
    let (ctx6, _HH) = field_square::<N, P>(ctx5, H);
    let (ctx7, I) = {
        let (c1, twoH) = field_add::<N, P>(ctx6, H, H);
        field_square::<N, P>(c1, twoH)
    };
    let (ctx8, J) = field_mul::<N, P>(ctx7, H, I);
    let (ctx9, S2mS1) = field_sub::<N, P>(ctx8, S2, S1);
    let (ctx10, r) = field_add::<N, P>(ctx9, S2mS1, S2mS1);
    let (ctx11, V) = field_mul::<N, P>(ctx10, U1, I);

    let (ctx12, r2) = field_square::<N, P>(ctx11, r);
    let (ctx13, t) = field_sub::<N, P>(ctx12, r2, J);
    let (ctx14, twoV) = field_add::<N, P>(ctx13, V, V);
    let (ctx15, X3) = field_sub::<N, P>(ctx14, t, twoV);

    let (ctx16, VmX3) = field_sub::<N, P>(ctx15, V, X3);
    let (ctx17, rVmX3) = field_mul::<N, P>(ctx16, r, VmX3);
    let (ctx18, S1J) = field_mul::<N, P>(ctx17, S1, J);
    let (ctx19, twoS1J) = field_add::<N, P>(ctx18, S1J, S1J);
    let (ctx20, Y3) = field_sub::<N, P>(ctx19, rVmX3, twoS1J);

    let (ctx21, Z1plus1) = field_add::<N, P>(ctx20, p.Z, fe_one::<N, P>());
    let (ctx22, Z1plus1_sq) = field_square::<N, P>(ctx21, Z1plus1);
    let (ctx23, t2) = field_sub::<N, P>(ctx22, Z1plus1_sq, Z1Z1);
    let (ctx24, t3) = field_sub::<N, P>(ctx23, t2, fe_one::<N, P>());
    let (ctx25, Z3) = field_mul::<N, P>(ctx24, t3, H);

    let R = JacobianPoint::<N, P> { X: X3, Y: Y3, Z: Z3 };

    let infP = JacobianPoint::<N, P> { X: fe_zero::<N, P>(), Y: fe_one::<N, P>(), Z: fe_zero::<N, P>() };
    let useP = jacobian_ct_select::<N, P>(infP, R, ((1u128 - p_is_inf) & 1u128));
    let out1 = jacobian_ct_select::<N, P>(qj, useP, ((1u128 - p_is_inf) & 1u128));
    let out2 = jacobian_ct_select::<N, P>(R, out1, q_is_inf);

    (ctx25, out2)
}

