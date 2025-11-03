use crate::bignum::BigNum;
use crate::limbs::{add_no_mod, sub_no_mod};

pub fn bignum_add_no_mod<let N: u32>(a: BigNum<N>, b: BigNum<N>) -> (BigNum<N>, u128) {
    let (sum, carry) = add_no_mod::<N>(a.limbs, b.limbs);
    (BigNum::from_limbs(sum), carry)
}

pub fn bignum_sub_no_mod<let N: u32>(a: BigNum<N>, b: BigNum<N>) -> (BigNum<N>, u128) {
    let (diff, borrow) = sub_no_mod::<N>(a.limbs, b.limbs);
    (BigNum::from_limbs(diff), borrow)
}

