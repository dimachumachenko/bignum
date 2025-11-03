use crate::bignum::BigNum;
use crate::limbs::{MASK_120, normalize_limb};

fn bn<let N: u32>(xs: [u128; N]) -> BigNum<N> {
    BigNum::from_limbs(xs)
}

#[test]
fn zero_and_from_limbs_normalizes() {
    let z = BigNum::<2>::zero();
    assert(z.limbs[0] == 0u128);
    assert(z.limbs[1] == 0u128);

    let hi_bit: u128 = 1u128 << 127;                // above 120
    let low_bits: u128 = 0x1234_5678_9ABC_DEF0u128; // within low 64
    let x = bn::<3>([low_bits | hi_bit, MASK_120, hi_bit]);

    assert(x.limbs[0] == normalize_limb(low_bits | hi_bit));
    assert(x.limbs[1] == MASK_120);
    assert(x.limbs[2] == 0u128);
}

#[test]
fn endian_le_single_limb_to_from_bytes_15() {
    let v: u128 =
          0x00u128
        | (0x11u128 << 0)
        | (0x22u128 << 8)
        | (0x33u128 << 16)
        | (0x44u128 << 24)
        | (0x55u128 << 32)
        | (0x66u128 << 40)
        | (0x77u128 << 48)
        | (0x88u128 << 56)
        | (0x99u128 << 64)
        | (0xAAu128 << 72)
        | (0xBBu128 << 80)
        | (0xCCu128 << 88)
        | (0xDDu128 << 96)
        | (0xEEu128 << 104);

    let a = bn::<1>([v]);

    let bytes = a.to_bytes_le::<15>();
    assert(bytes[0]  == 0x11u8);
    assert(bytes[1]  == 0x22u8);
    assert(bytes[2]  == 0x33u8);
    assert(bytes[3]  == 0x44u8);
    assert(bytes[4]  == 0x55u8);
    assert(bytes[5]  == 0x66u8);
    assert(bytes[6]  == 0x77u8);
    assert(bytes[7]  == 0x88u8);
    assert(bytes[8]  == 0x99u8);
    assert(bytes[9]  == 0xAAu8);
    assert(bytes[10] == 0xBBu8);
    assert(bytes[11] == 0xCCu8);
    assert(bytes[12] == 0xDDu8);
    assert(bytes[13] == 0xEEu8);
    assert(bytes[14] == 0x00u8); 

    let b = BigNum::<1>::from_bytes_le(bytes);
    assert(b.limbs[0] == normalize_limb(v));
}

#[test]
fn roundtrip_bytes_le_two_limbs_30() {
    let bytes: [u8; 30] = [
        0,1,2,3,4,5,6,7,8,9,
        10,11,12,13,14,15,16,17,18,19,
        20,21,22,23,24,25,26,27,28,29
    ];
    let a = BigNum::<2>::from_bytes_le(bytes);
    let bb = a.to_bytes_le::<30>();
    assert(bb == bytes);
}

#[test]
fn roundtrip_bytes_le_three_limbs_45() {
    let bytes: [u8; 45] = [
        0,  3,  6,  9, 12, 15, 18, 21, 24, 27,
        30, 33, 36, 39, 42, 45, 48, 51, 54, 57,
        60, 63, 66, 69, 72, 75, 78, 81, 84, 87,
        90, 93, 96, 99, 102,105,108,111,114,117,
        120,123,126,129,132
    ];
    let a = BigNum::<3>::from_bytes_le(bytes);
    let bb = a.to_bytes_le::<45>();
    assert(bb == bytes);
}

#[test]
fn ct_select_works() {
    let a = bn::<2>([1u128, 2u128]);
    let b = bn::<2>([3u128, 4u128]);

    let pick_a = BigNum::<2>::ct_select(a, b, 0u128);
    let pick_b = BigNum::<2>::ct_select(a, b, 1u128);

    assert(pick_a.limbs[0] == 1u128);
    assert(pick_a.limbs[1] == 2u128);
    assert(pick_b.limbs[0] == 3u128);
    assert(pick_b.limbs[1] == 4u128);
}

#[test]
fn from_bytes_does_not_overflow_limbs() {
    let src: [u8; 32] = [
        0xAA,0xBB,0xCC,0xDD,0xEE,0xFF,0x10,0x20,
        0x30,0x40,0x50,0x60,0x70,0x80,0x90,0xA0,
        0xB0,0xC0,0xD0,0xE0,0xF0,0x01,0x02,0x03,
        0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B
    ];
    let a = BigNum::<2>::from_bytes_le(src);
    let back = a.to_bytes_le::<30>();
    let expect: [u8; 30] = [
        0xAA,0xBB,0xCC,0xDD,0xEE,0xFF,0x10,0x20,
        0x30,0x40,0x50,0x60,0x70,0x80,0x90,0xA0,
        0xB0,0xC0,0xD0,0xE0,0xF0,0x01,0x02,0x03,
        0x04,0x05,0x06,0x07,0x08,0x09
    ];
    assert(back == expect);
}


