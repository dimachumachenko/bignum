use crate::bignum::BigNum;
use crate::limbs::MASK_120;

pub fn test_p_N2() -> BigNum<2> {
    let mut limbs: [u128; 2] = [0u128; 2];
    limbs[0] = MASK_120;
    limbs[1] = MASK_120;
    BigNum::<2>::from_limbs(limbs)
}

pub fn test_r2_N2() -> BigNum<2> {
    BigNum::<2>::from_limbs([1u128, 0u128])
}

pub fn test_ninv() -> u128 {
    1u128
}
pub fn test_p_minus_2_N2() -> BigNum<2> {
    let a0 = (MASK_120 - 2u128) & MASK_120;
    let a1 = MASK_120;
    BigNum::<2>::from_limbs([a0, a1])
}
